"""
Tier 1 tests — Day 5 requirement.

Tests run without GPU or real weights:
  - AXIOMTier1Net architecture (CPU forward pass)
  - AXIOMTier1Ensemble with injected mock models (no HuggingFace required)
  - Confidence calculation and routing logic
  - CensusAnnotator instantiation and predict() with injected kNN index
  - Tier1Ensemble fusion logic
  - Feature alignment (gene list mismatch handling)

All tests must pass in CI (ubuntu, CPU-only, no network calls).
"""
from __future__ import annotations

import json
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest


# ── fixtures ──────────────────────────────────────────────────────────────────

N_GENES   = 64    # small for CPU tests (real model uses 2000)
N_CLASSES = 10
N_CELLS   = 20
HIDDEN    = 32


@pytest.fixture(scope="module")
def mock_torch_models():
    """10 tiny AXIOMTier1Net models with random weights — no GPU needed."""
    torch = pytest.importorskip("torch")
    from axiom_sc.tier1.axiomtier1 import AXIOMTier1Net
    models = []
    for i in range(10):
        torch.manual_seed(i)
        net = AXIOMTier1Net(n_genes=N_GENES, n_classes=N_CLASSES, hidden=HIDDEN)
        net.eval()
        models.append(net)
    return models


@pytest.fixture(scope="module")
def label_encoder():
    return {i: f"CellType_{i}" for i in range(N_CLASSES)}


@pytest.fixture(scope="module")
def gene_list():
    return [f"GENE_{i:04d}" for i in range(N_GENES)]


@pytest.fixture(scope="module")
def mock_ensemble(mock_torch_models, label_encoder, gene_list):
    """AXIOMTier1Ensemble pre-loaded with mock models — no weight files needed."""
    from axiom_sc.tier1.axiomtier1 import AXIOMTier1Ensemble
    return AXIOMTier1Ensemble(
        _models=mock_torch_models,
        _label_encoder=label_encoder,
        _gene_list=gene_list,
    )


@pytest.fixture(scope="module")
def mock_adata(gene_list):
    """Tiny synthetic AnnData (N_CELLS × N_GENES) with log1p-normalised values."""
    anndata = pytest.importorskip("anndata")
    rng = np.random.default_rng(0)
    X = rng.exponential(scale=0.5, size=(N_CELLS, N_GENES)).astype(np.float32)
    import pandas as pd
    obs = pd.DataFrame(index=[f"cell_{i}" for i in range(N_CELLS)])
    var = pd.DataFrame(index=gene_list)
    return anndata.AnnData(X=X, obs=obs, var=var)


# ── AXIOMTier1Net architecture tests ─────────────────────────────────────────

def test_axiomtier1net_builds():
    """AXIOMTier1Net builds without error and has correct layer structure."""
    torch = pytest.importorskip("torch")
    from axiom_sc.tier1.axiomtier1 import AXIOMTier1Net

    net = AXIOMTier1Net(n_genes=N_GENES, n_classes=N_CLASSES, hidden=HIDDEN)
    assert net.n_genes == N_GENES
    assert net.n_classes == N_CLASSES


def test_axiomtier1net_forward_pass(mock_torch_models):
    """Single model forward pass produces correct output shape."""
    torch = pytest.importorskip("torch")

    net = mock_torch_models[0]
    X = torch.randn(N_CELLS, N_GENES)
    with torch.no_grad():
        logits = net(X)
    assert logits.shape == (N_CELLS, N_CLASSES), (
        f"Expected logits shape ({N_CELLS}, {N_CLASSES}), got {logits.shape}"
    )


def test_axiomtier1net_output_unnormalised(mock_torch_models):
    """Model outputs raw logits (not softmax). Softmax applied in ensemble."""
    torch = pytest.importorskip("torch")
    import torch.nn.functional as F

    net = mock_torch_models[0]
    X = torch.randn(5, N_GENES)
    with torch.no_grad():
        logits = net(X)
    # Logits should NOT sum to 1
    row_sums = logits.sum(dim=-1).numpy()
    assert not np.allclose(row_sums, 1.0, atol=0.01), (
        "Model output looks like softmax — should be raw logits"
    )
    # But softmax of logits should
    probs = F.softmax(logits, dim=-1).numpy()
    assert np.allclose(probs.sum(axis=-1), 1.0, atol=1e-5)


# ── confidence calculation tests ──────────────────────────────────────────────

def test_compute_confidence_uniform():
    """Uniform distribution → minimum confidence (maximum entropy)."""
    from axiom_sc.tier1.axiomtier1 import _compute_confidence

    # Uniform: maximum entropy → confidence should be near 0
    probs = np.ones((5, N_CLASSES)) / N_CLASSES
    conf = _compute_confidence(probs)
    assert conf.shape == (5,)
    assert np.all(conf < 0.05), f"Uniform dist should have confidence ≈ 0, got {conf}"


def test_compute_confidence_peaked():
    """Peaked distribution → high confidence (near 1)."""
    from axiom_sc.tier1.axiomtier1 import _compute_confidence

    # Near-deterministic: put 0.999 on one class
    probs = np.full((5, N_CLASSES), 0.001 / (N_CLASSES - 1))
    probs[:, 0] = 0.999
    conf = _compute_confidence(probs)
    assert np.all(conf > 0.90), f"Peaked dist should have high confidence, got {conf}"


def test_compute_confidence_range():
    """Confidence must always be in [0, 1]."""
    from axiom_sc.tier1.axiomtier1 import _compute_confidence

    rng = np.random.default_rng(42)
    # Random valid probability distributions
    raw = rng.dirichlet(np.ones(N_CLASSES), size=100)
    conf = _compute_confidence(raw)
    assert np.all(conf >= 0.0) and np.all(conf <= 1.0), (
        f"Confidence out of [0,1] range: min={conf.min():.4f}, max={conf.max():.4f}"
    )


# ── routing tests ─────────────────────────────────────────────────────────────

def test_route_accept():
    from axiom_sc.tier1.axiomtier1 import _route
    assert _route(np.array([0.90, 0.85])) == ["accept", "accept"]


def test_route_tier2_verify():
    from axiom_sc.tier1.axiomtier1 import _route
    assert _route(np.array([0.75, 0.50])) == ["tier2_verify", "tier2_verify"]


def test_route_tier2_full():
    from axiom_sc.tier1.axiomtier1 import _route
    assert _route(np.array([0.49, 0.10, 0.0])) == [
        "tier2_full", "tier2_full", "tier2_full"
    ]


def test_route_boundary_conditions():
    """Exact threshold values route correctly (accept at 0.85, verify at 0.50)."""
    from axiom_sc.tier1.axiomtier1 import _route
    routes = _route(np.array([0.85, 0.50, 0.84999, 0.49999]))
    assert routes[0] == "accept"
    assert routes[1] == "tier2_verify"
    assert routes[2] == "tier2_verify"
    assert routes[3] == "tier2_full"


# ── AXIOMTier1Ensemble predict tests ──────────────────────────────────────────

def test_ensemble_predict_shape(mock_ensemble, mock_adata):
    """Ensemble predict() returns correct shapes for all output keys."""
    result = mock_ensemble.predict(mock_adata)

    assert "label" in result
    assert "confidence" in result
    assert "route" in result
    assert "probs" in result

    assert len(result["label"]) == N_CELLS
    assert len(result["confidence"]) == N_CELLS
    assert len(result["route"]) == N_CELLS
    assert result["probs"].shape == (N_CELLS, N_CLASSES)


def test_ensemble_predict_labels_valid(mock_ensemble, mock_adata, label_encoder):
    """All predicted labels must be known cell type names."""
    result = mock_ensemble.predict(mock_adata)
    valid_labels = set(label_encoder.values())
    for label in result["label"]:
        assert label in valid_labels, f"Unknown label: {label}"


def test_ensemble_predict_confidence_in_range(mock_ensemble, mock_adata):
    """All confidence values must be in [0, 1]."""
    result = mock_ensemble.predict(mock_adata)
    confs = result["confidence"].values
    assert np.all(confs >= 0.0) and np.all(confs <= 1.0), (
        f"Confidence out of range: min={confs.min():.4f}, max={confs.max():.4f}"
    )


def test_ensemble_predict_routes_valid(mock_ensemble, mock_adata):
    """All routes must be one of the three valid values."""
    result = mock_ensemble.predict(mock_adata)
    valid_routes = {"accept", "tier2_verify", "tier2_full"}
    for route in result["route"]:
        assert route in valid_routes, f"Invalid route: {route}"


def test_ensemble_predict_obs_names_preserved(mock_ensemble, mock_adata):
    """Output Series index must match adata.obs_names."""
    result = mock_ensemble.predict(mock_adata)
    assert list(result["label"].index) == list(mock_adata.obs_names)
    assert list(result["confidence"].index) == list(mock_adata.obs_names)


def test_ensemble_repr(mock_ensemble):
    """Ensemble repr includes expected fields."""
    r = repr(mock_ensemble)
    assert "AXIOMTier1Ensemble" in r
    assert "loaded" in r


# ── feature alignment tests ───────────────────────────────────────────────────

def test_feature_alignment_missing_genes(mock_ensemble):
    """Genes not in adata.var_names are zero-filled without error."""
    anndata = pytest.importorskip("anndata")
    import pandas as pd

    # adata with only HALF the genes the model expects
    half_genes = mock_ensemble._gene_list[:N_GENES // 2]
    rng = np.random.default_rng(1)
    X = rng.exponential(size=(5, len(half_genes))).astype(np.float32)
    obs = pd.DataFrame(index=[f"c_{i}" for i in range(5)])
    var = pd.DataFrame(index=half_genes)
    adata_partial = anndata.AnnData(X=X, obs=obs, var=var)

    # Should not raise — missing genes zero-filled
    result = mock_ensemble.predict(adata_partial)
    assert len(result["label"]) == 5


def test_feature_alignment_extra_genes(mock_ensemble):
    """Extra genes in adata (not in model's list) are ignored gracefully."""
    anndata = pytest.importorskip("anndata")
    import pandas as pd

    extra_genes = mock_ensemble._gene_list + [f"EXTRA_{i}" for i in range(20)]
    rng = np.random.default_rng(2)
    X = rng.exponential(size=(5, len(extra_genes))).astype(np.float32)
    obs = pd.DataFrame(index=[f"c_{i}" for i in range(5)])
    var = pd.DataFrame(index=extra_genes)
    adata_extra = anndata.AnnData(X=X, obs=obs, var=var)

    result = mock_ensemble.predict(adata_extra)
    assert len(result["label"]) == 5


def test_feature_alignment_sparse_input(mock_ensemble, gene_list):
    """Sparse AnnData input (scipy.sparse) aligns correctly."""
    anndata = pytest.importorskip("anndata")
    import scipy.sparse as sp
    import pandas as pd

    rng = np.random.default_rng(3)
    X_dense = rng.exponential(size=(5, N_GENES)).astype(np.float32)
    X_sparse = sp.csr_matrix(X_dense)
    obs = pd.DataFrame(index=[f"c_{i}" for i in range(5)])
    var = pd.DataFrame(index=gene_list)
    adata_sparse = anndata.AnnData(X=X_sparse, obs=obs, var=var)

    result = mock_ensemble.predict(adata_sparse)
    assert len(result["label"]) == 5


# ── weight loading tests (from temp files) ───────────────────────────────────

def test_load_weights_from_disk(mock_torch_models, label_encoder, gene_list):
    """load_weights() correctly restores models, gene list, and label encoder."""
    torch = pytest.importorskip("torch")
    from axiom_sc.tier1.axiomtier1 import AXIOMTier1Ensemble

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)

        # Save model weights
        for i, net in enumerate(mock_torch_models):
            torch.save(net.state_dict(), tmp / f"axiomtier1_model_{i}.pt")

        # Save gene list, label encoder, and architecture config
        (tmp / "axiomtier1_genes.txt").write_text("\n".join(gene_list))
        (tmp / "label_encoder.json").write_text(
            json.dumps({str(k): v for k, v in label_encoder.items()})
        )
        # model_config.json tells load_weights() the hidden dim used in tests
        (tmp / "model_config.json").write_text(
            json.dumps({"hidden": HIDDEN, "n_models": 10, "n_genes": N_GENES})
        )

        # Load fresh ensemble from disk
        ensemble = AXIOMTier1Ensemble()
        ensemble.load_weights(str(tmp))

        assert ensemble._loaded
        assert len(ensemble._models) == 10
        assert ensemble._gene_list == gene_list
        assert ensemble._label_encoder == label_encoder


def test_load_weights_missing_dir_raises():
    """load_weights() raises FileNotFoundError on missing directory."""
    from axiom_sc.tier1.axiomtier1 import AXIOMTier1Ensemble
    pytest.importorskip("torch")

    ensemble = AXIOMTier1Ensemble(weights_path="/nonexistent/path")
    with pytest.raises(FileNotFoundError):
        ensemble.load_weights("/nonexistent/path")


# ── CensusAnnotator tests ─────────────────────────────────────────────────────

def test_census_annotator_instantiation():
    """CensusAnnotator instantiates with default params."""
    from axiom_sc.tier1.census_annotator import CensusAnnotator

    ann = CensusAnnotator(tissue_filter=["blood", "lung"], n_neighbors=20)
    assert ann.tissue_filter == ["blood", "lung"]
    assert ann.n_neighbors == 20
    r = repr(ann)
    assert "CensusAnnotator" in r


def test_census_annotator_predict_without_reference_raises():
    """predict() before build_reference() raises RuntimeError."""
    from axiom_sc.tier1.census_annotator import CensusAnnotator
    anndata = pytest.importorskip("anndata")
    import pandas as pd

    ann = CensusAnnotator(use_existing_obsm="X_scvi")
    # No reference built
    obs = pd.DataFrame(index=["c0"])
    adata = anndata.AnnData(X=np.zeros((1, 10)), obs=obs)

    with pytest.raises(RuntimeError, match="reference not built"):
        ann.predict(adata)


def test_census_annotator_predict_with_injected_index():
    """CensusAnnotator.predict() with a pre-built sklearn index (no Census needed)."""
    sklearn = pytest.importorskip("sklearn")
    anndata = pytest.importorskip("anndata")
    import pandas as pd
    from sklearn.neighbors import NearestNeighbors
    from axiom_sc.tier1.census_annotator import CensusAnnotator

    # Build a tiny reference
    rng = np.random.default_rng(0)
    n_ref, emb_dim = 100, 16
    ref_emb = rng.normal(size=(n_ref, emb_dim)).astype(np.float32)
    ref_labels = [f"CellType_{i % 5}" for i in range(n_ref)]

    # Inject reference + index
    ann = CensusAnnotator(n_neighbors=10, use_existing_obsm="X_emb")
    ann._ref_embeddings = ref_emb
    ann._ref_labels = ref_labels
    ann._build_knn_index(ref_emb)

    # Query
    n_query = 8
    query_emb = rng.normal(size=(n_query, emb_dim)).astype(np.float32)
    obs = pd.DataFrame(index=[f"q_{i}" for i in range(n_query)])
    adata = anndata.AnnData(
        X=np.zeros((n_query, 1), dtype=np.float32),
        obs=obs,
        obsm={"X_emb": query_emb},
    )

    result = ann.predict(adata)

    assert len(result["label"]) == n_query
    assert len(result["confidence"]) == n_query
    assert len(result["route"]) == n_query
    assert list(result["label"].index) == list(adata.obs_names)

    # All confidence values in [0, 1]
    confs = result["confidence"].values
    assert np.all(confs >= 0.0) and np.all(confs <= 1.0)


# ── Tier1Ensemble fusion tests ─────────────────────────────────────────────────

def test_tier1_ensemble_mlp_only(mock_ensemble, mock_adata):
    """Tier1Ensemble with mlp-only falls through to MLP predict."""
    from axiom_sc.tier1.ensemble import Tier1Ensemble

    t1 = Tier1Ensemble(mlp=mock_ensemble)
    result = t1.predict(mock_adata)
    assert len(result["label"]) == N_CELLS


def test_tier1_ensemble_no_annotators_raises():
    """Tier1Ensemble with neither mlp nor census raises ValueError."""
    from axiom_sc.tier1.ensemble import Tier1Ensemble

    with pytest.raises(ValueError, match="At least one"):
        Tier1Ensemble(mlp=None, census=None)


def test_tier1_ensemble_acceptance_rate(mock_ensemble, mock_adata):
    """acceptance_rate() returns float in [0, 1]."""
    from axiom_sc.tier1.ensemble import Tier1Ensemble

    t1 = Tier1Ensemble(mlp=mock_ensemble)
    result = t1.predict(mock_adata)
    rate = t1.acceptance_rate(result)
    assert 0.0 <= rate <= 1.0, f"Acceptance rate {rate} out of [0, 1]"


def test_tier1_ensemble_repr(mock_ensemble):
    from axiom_sc.tier1.ensemble import Tier1Ensemble
    t1 = Tier1Ensemble(mlp=mock_ensemble)
    assert "Tier1Ensemble" in repr(t1)


# ── CensusDataset tests (training data loader) ────────────────────────────────

def test_census_dataset_loads_npz():
    """CensusDataset reads .npz and returns correct (X, y) tensors."""
    torch = pytest.importorskip("torch")
    from axiom_sc.tier1.training.data_loader import CensusDataset

    with tempfile.TemporaryDirectory() as tmpdir:
        npz_path = Path(tmpdir) / "test.npz"
        rng = np.random.default_rng(0)
        X = rng.normal(size=(50, N_GENES)).astype(np.float32)
        y = rng.integers(0, N_CLASSES, size=50).astype(np.int32)
        np.savez_compressed(str(npz_path), X=X, y=y)

        ds = CensusDataset(str(npz_path))
        assert len(ds) == 50

        x0, y0 = ds[0]
        assert x0.shape == (N_GENES,)
        assert y0.dtype == torch.long


def test_stratified_split_min_per_type():
    """_stratified_split respects min_per_type guarantee."""
    from axiom_sc.tier1.training.data_loader import _stratified_split

    rng = np.random.default_rng(0)
    # 5 classes, 3000 cells each
    labels = np.repeat(np.arange(5), 3000)
    rng.shuffle(labels)

    train_idx, val_idx, test_idx = _stratified_split(
        labels, train_frac=0.80, val_frac=0.10,
        min_per_type=2000, random_seed=42,
    )

    # Each class must have ≥ 2000 in training
    for cls in range(5):
        n_cls_train = (labels[train_idx] == cls).sum()
        assert n_cls_train >= 2000, (
            f"Class {cls} has only {n_cls_train} training cells, expected ≥ 2000"
        )

    # No overlap between splits
    train_set, val_set, test_set = set(train_idx), set(val_idx), set(test_idx)
    assert not (train_set & val_set), "Train/val overlap"
    assert not (train_set & test_set), "Train/test overlap"
    assert not (val_set & test_set), "Val/test overlap"
