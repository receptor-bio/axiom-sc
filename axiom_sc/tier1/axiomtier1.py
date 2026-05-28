"""
AXIOMTier1: MLP ensemble for cell type classification.

Architecture based on: Fischer et al. (2024) Nature Communications 15:6611
  "scTab: Scaling cross-tissue supervised cell type annotation models"
  doi: 10.1038/s41467-024-46499-y

Training data: CELLxGENE Census (CC BY 4.0)
  Tabula Sapiens Consortium (2022) Science 376:eabl4896
  doi: 10.1126/science.abl4896

Our trained weights: Apache 2.0
  Released at: huggingface.co/receptor-bio/axiomtier1

This module implements the inference side only. Training code is in
axiom_sc/tier1/training/train.py (run once on GPU cluster).

License: Apache 2.0
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    import pandas as pd

logger = logging.getLogger(__name__)

# ── constants (Phase 1 validated thresholds — do not change without re-validation) ──

CONFIDENCE_THRESHOLD_ACCEPT = 0.85   # route: Tier 1 terminal accept
CONFIDENCE_THRESHOLD_VERIFY = 0.50   # 0.50–0.85: Tier 2 verification of Tier 1 call
# below 0.50: Tier 2 full search, Tier 1 call not trusted

N_MODELS = 10          # ensemble size
N_GENES = 2000         # top HVGs used as input
HIDDEN_DIM = 512       # MLP hidden layer width
HUGGINGFACE_REPO = "receptor-bio/axiomtier1"


# ── network architecture ──────────────────────────────────────────────────────

def _build_net(n_genes: int, n_classes: int, hidden: int = HIDDEN_DIM):
    """
    Build a single MLP in the ensemble.
    Lazy-imports torch so the package installs without it.
    """
    try:
        import torch.nn as nn
    except ImportError:
        raise ImportError(
            "torch required for AXIOMTier1Net. "
            "Install with: pip install 'axiom-sc[science]'"
        )

    return nn.Sequential(
        nn.Linear(n_genes, hidden), nn.BatchNorm1d(hidden),
        nn.ReLU(), nn.Dropout(0.1),
        nn.Linear(hidden, hidden), nn.BatchNorm1d(hidden),
        nn.ReLU(), nn.Dropout(0.1),
        nn.Linear(hidden, 256), nn.BatchNorm1d(256),
        nn.ReLU(), nn.Dropout(0.1),
        nn.Linear(256, n_classes),
    )


class AXIOMTier1Net:
    """
    Thin wrapper around the sequential MLP so callers get a consistent API
    regardless of whether they use the nn.Sequential directly or this class.

    Parameters
    ----------
    n_genes : int   — number of input HVGs (must match training)
    n_classes : int — number of cell types (must match training label encoder)
    hidden : int    — hidden layer width (default 512)
    """

    def __init__(self, n_genes: int, n_classes: int, hidden: int = HIDDEN_DIM):
        self.n_genes = n_genes
        self.n_classes = n_classes
        self.hidden = hidden
        self._net = _build_net(n_genes, n_classes, hidden)

    def forward(self, x):
        return self._net(x)

    def __call__(self, x):
        return self._net(x)

    def eval(self):
        self._net.eval()
        return self

    def load_state_dict(self, state_dict, strict: bool = True):
        self._net.load_state_dict(state_dict, strict=strict)

    def state_dict(self):
        return self._net.state_dict()

    def parameters(self):
        return self._net.parameters()

    @property
    def net(self):
        return self._net


# ── confidence calculation ────────────────────────────────────────────────────

def _compute_confidence(mean_probs: "np.ndarray") -> "np.ndarray":
    """
    Confidence = 1 - normalised Shannon entropy of the mean softmax distribution.

    H(p) = -sum(p * log(p))
    H_norm = H(p) / log(n_classes)  ∈ [0, 1]
    confidence = 1 - H_norm         ∈ [0, 1]

    High confidence → peaked distribution → low entropy → confidence near 1.
    Uniform distribution → maximum entropy → confidence = 0.
    """
    # Clip to avoid log(0)
    p = np.clip(mean_probs, 1e-10, 1.0)
    n_classes = p.shape[-1]
    h = -np.sum(p * np.log(p), axis=-1)
    h_norm = h / np.log(n_classes)
    return 1.0 - h_norm


def _route(confidence: "np.ndarray") -> list[str]:
    """Map confidence scores to routing decisions."""
    routes = []
    for c in confidence:
        if c >= CONFIDENCE_THRESHOLD_ACCEPT:
            routes.append("accept")
        elif c >= CONFIDENCE_THRESHOLD_VERIFY:
            routes.append("tier2_verify")
        else:
            routes.append("tier2_full")
    return routes


# ── ensemble ──────────────────────────────────────────────────────────────────

class AXIOMTier1Ensemble:
    """
    Ensemble of N=10 independently trained AXIOMTier1Net models.

    Inference pipeline
    ------------------
    1. Each model produces a softmax distribution over cell types.
    2. Mean softmax is taken across the 10 models.
    3. Confidence = 1 - normalised_entropy(mean_softmax).
    4. Routing:
       confidence ≥ 0.85  → 'accept'        (Tier 1 terminal)
       0.50 ≤ conf < 0.85 → 'tier2_verify'  (Tier 2 verifies the Tier 1 call)
       conf < 0.50        → 'tier2_full'     (Tier 2 full KG search)

    Parameters
    ----------
    weights_path : str | None
        Local directory containing 10 model weight files
        (axiomtier1_model_0.pt … axiomtier1_model_9.pt) and a
        label_encoder.json file mapping class indices to cell type names.
        If None: downloads receptor-bio/axiomtier1 from HuggingFace on first use.
    gene_list_path : str | None
        Path to the HVG gene list (axiomtier1_genes.txt).
        Required for aligning input adata to model's feature space.
        If None: loaded from the weights directory.
    _models : list | None
        For testing only — pre-loaded model objects bypass weight loading.
    _label_encoder : dict | None
        For testing only — {int → cell_type_str} mapping.
    _gene_list : list[str] | None
        For testing only — ordered gene list.
    """

    N_MODELS = N_MODELS
    CONFIDENCE_THRESHOLD_ACCEPT = CONFIDENCE_THRESHOLD_ACCEPT
    CONFIDENCE_THRESHOLD_VERIFY = CONFIDENCE_THRESHOLD_VERIFY

    def __init__(
        self,
        weights_path: str | None = None,
        gene_list_path: str | None = None,
        *,
        _models: list | None = None,
        _label_encoder: dict | None = None,
        _gene_list: list[str] | None = None,
    ):
        self._weights_path = weights_path
        self._gene_list_path = gene_list_path

        # Injected in tests / after weight download
        self._models: list | None = _models
        self._label_encoder: dict | None = _label_encoder   # {idx: cell_type_str}
        self._gene_list: list[str] | None = _gene_list
        self._loaded = _models is not None

    # ── weight loading ────────────────────────────────────────────────────────

    def load_weights(self, weights_path: str | None = None) -> None:
        """
        Load all 10 model weights from a local directory.

        Expected layout::
            <weights_path>/
                axiomtier1_model_0.pt … axiomtier1_model_9.pt
                label_encoder.json    — {str(idx): cell_type_name}
                axiomtier1_genes.txt  — one HVG per line

        Parameters
        ----------
        weights_path : str
            Directory containing model files. Uses self._weights_path if None.
        """
        try:
            import torch
        except ImportError:
            raise ImportError(
                "torch required. Install with: pip install 'axiom-sc[science]'"
            )
        import json

        path = Path(weights_path or self._weights_path)
        if not path.is_dir():
            raise FileNotFoundError(f"Weights directory not found: {path}")

        # Load gene list
        gene_file = Path(self._gene_list_path) if self._gene_list_path else path / "axiomtier1_genes.txt"
        if not gene_file.exists():
            raise FileNotFoundError(f"Gene list not found: {gene_file}")
        self._gene_list = gene_file.read_text().strip().splitlines()

        # Load label encoder
        encoder_file = path / "label_encoder.json"
        if not encoder_file.exists():
            raise FileNotFoundError(f"label_encoder.json not found in {path}")
        raw = json.loads(encoder_file.read_text())
        self._label_encoder = {int(k): v for k, v in raw.items()}

        n_genes = len(self._gene_list)
        n_classes = len(self._label_encoder)

        # Read architecture config (hidden dim may differ from default in tests)
        config_file = path / "model_config.json"
        hidden = HIDDEN_DIM  # default
        if config_file.exists():
            import json as _json
            cfg = _json.loads(config_file.read_text())
            hidden = cfg.get("hidden", HIDDEN_DIM)
            n_models = cfg.get("n_models", self.N_MODELS)
        else:
            n_models = self.N_MODELS

        # Load each model
        self._models = []
        for i in range(n_models):
            model_file = path / f"axiomtier1_model_{i}.pt"
            if not model_file.exists():
                raise FileNotFoundError(f"Model file not found: {model_file}")
            net = AXIOMTier1Net(n_genes=n_genes, n_classes=n_classes, hidden=hidden)
            sd = torch.load(model_file, map_location="cpu", weights_only=True)
            net.load_state_dict(sd)
            net.eval()
            self._models.append(net)

        self._loaded = True
        logger.info(
            "AXIOMTier1Ensemble: loaded %d models | %d genes | %d cell types",
            len(self._models), n_genes, n_classes,
        )

    def download_from_huggingface(self) -> Path:
        """
        Download weights from HuggingFace receptor-bio/axiomtier1.
        Caches to ~/.cache/axiom_sc/axiomtier1/.
        Returns the local cache directory path.
        """
        try:
            from huggingface_hub import snapshot_download
        except ImportError:
            raise ImportError(
                "huggingface_hub required: pip install huggingface-hub"
            )
        cache_dir = Path.home() / ".cache" / "axiom_sc" / "axiomtier1"
        cache_dir.mkdir(parents=True, exist_ok=True)
        local_path = snapshot_download(
            repo_id=HUGGINGFACE_REPO,
            cache_dir=str(cache_dir),
        )
        logger.info("AXIOMTier1: downloaded weights to %s", local_path)
        return Path(local_path)

    def _ensure_loaded(self) -> None:
        """Load weights on first predict() call if not already loaded."""
        if self._loaded:
            return
        if self._weights_path and Path(self._weights_path).is_dir():
            self.load_weights(self._weights_path)
        else:
            logger.info("Downloading AXIOMTier1 weights from HuggingFace …")
            local = self.download_from_huggingface()
            self.load_weights(str(local))

    # ── inference ─────────────────────────────────────────────────────────────

    def predict(self, adata) -> dict:
        """
        Predict cell types for all cells in an AnnData object.

        Parameters
        ----------
        adata : anndata.AnnData
            Must be:
            - Normalised to 10k counts (sc.pp.normalize_total, target_sum=1e4)
            - log1p transformed (sc.pp.log1p)
            - .var_names must overlap with the model's gene list

        Returns
        -------
        dict with keys:
            'label'      : pd.Series[str], index = adata.obs_names
            'confidence' : pd.Series[float], index = adata.obs_names
            'route'      : pd.Series[str], values ∈ {'accept','tier2_verify','tier2_full'}
            'probs'      : np.ndarray shape (n_cells, n_classes) — mean softmax
        """
        import pandas as pd

        self._ensure_loaded()
        assert self._models is not None
        assert self._label_encoder is not None
        assert self._gene_list is not None

        X = self._align_features(adata)              # (n_cells, n_genes) float32
        all_probs = self._run_ensemble(X)             # (n_models, n_cells, n_classes)
        mean_probs = all_probs.mean(axis=0)           # (n_cells, n_classes)
        confidence = _compute_confidence(mean_probs)  # (n_cells,)
        pred_idx = mean_probs.argmax(axis=1)          # (n_cells,)

        labels = [self._label_encoder[i] for i in pred_idx]
        routes = _route(confidence)

        obs_names = adata.obs_names
        return {
            "label":      pd.Series(labels, index=obs_names, name="tier1_label"),
            "confidence": pd.Series(confidence.tolist(), index=obs_names, name="tier1_confidence"),
            "route":      pd.Series(routes, index=obs_names, name="tier1_route"),
            "probs":      mean_probs,
        }

    def _align_features(self, adata) -> "np.ndarray":
        """
        Align adata.X to the model's expected gene order.
        Genes missing from adata are zeroed; extra genes are dropped.

        Returns float32 numpy array shape (n_cells, n_genes).
        """
        import scipy.sparse as sp

        n_cells = adata.n_obs
        n_genes = len(self._gene_list)
        X_aligned = np.zeros((n_cells, n_genes), dtype=np.float32)

        var_names = list(adata.var_names)
        gene_to_idx_adata = {g: i for i, g in enumerate(var_names)}

        for j, gene in enumerate(self._gene_list):
            if gene in gene_to_idx_adata:
                adata_j = gene_to_idx_adata[gene]
                if sp.issparse(adata.X):
                    X_aligned[:, j] = np.asarray(
                        adata.X[:, adata_j].todense()
                    ).ravel()
                else:
                    X_aligned[:, j] = adata.X[:, adata_j]

        missing = sum(1 for g in self._gene_list if g not in gene_to_idx_adata)
        if missing > 0:
            pct = 100 * missing / n_genes
            logger.warning(
                "AXIOMTier1: %d/%d model genes (%.1f%%) not in adata.var_names — "
                "set to zero. Accuracy may be reduced if >20%% missing.",
                missing, n_genes, pct,
            )

        return X_aligned

    def _run_ensemble(self, X: "np.ndarray") -> "np.ndarray":
        """
        Run all N models on input X.
        Returns softmax probabilities shape (N_MODELS, n_cells, n_classes).
        """
        import torch
        import torch.nn.functional as F

        tensor = torch.from_numpy(X)
        all_probs = []

        with torch.no_grad():
            for model in self._models:
                logits = model(tensor)
                probs = F.softmax(logits, dim=-1).numpy()
                all_probs.append(probs)

        return np.stack(all_probs, axis=0)  # (N_MODELS, n_cells, n_classes)

    def __repr__(self) -> str:
        status = "loaded" if self._loaded else "not loaded"
        n_types = len(self._label_encoder) if self._label_encoder else "?"
        return (
            f"AXIOMTier1Ensemble(n_models={self.N_MODELS}, "
            f"n_cell_types={n_types}, status={status!r})"
        )
