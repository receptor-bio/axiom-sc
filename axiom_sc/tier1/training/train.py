"""
One-time training pipeline for AXIOMTier1 MLP ensemble.

Run on GPU cluster::
    axiom-train-tier1 --output weights/ --n-models 10

Compute estimate:
  Single A100 (80GB): ~5h per model, ~50h total for 10 models
  Cost: ~$75 at Lambda Labs / Vast.ai
  Run once; weights pushed to HuggingFace receptor-bio/axiomtier1

Architecture reference:
  Fischer et al. (2024) Nature Communications 15:6611
  doi: 10.1038/s41467-024-46499-y

Training data:
  CELLxGENE Census (MIT API, CC BY 4.0 data)
  doi: 10.1126/science.abl4896

Our trained weights: Apache 2.0
  huggingface.co/receptor-bio/axiomtier1

License: Apache 2.0
"""
from __future__ import annotations

import argparse
import json
import logging
import time
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)

# ── training hyper-parameters ─────────────────────────────────────────────────

DEFAULTS = {
    "n_models":      10,
    "epochs":        50,
    "batch_size":    512,
    "lr":            1e-3,
    "weight_decay":  1e-4,
    "hidden":        512,
    "dropout":       0.1,
    "patience":      5,     # early stopping patience (epochs without val improvement)
}

TARGET_MACRO_F1 = 0.80   # minimum acceptable macro-F1 on test set


# ── main entry point ──────────────────────────────────────────────────────────

def main(argv: list[str] | None = None) -> None:
    """
    CLI entry point: axiom-train-tier1

    Steps
    -----
    1. Build (or load) training dataset from CELLxGENE Census
    2. Train N_MODELS independently with different random seeds
    3. Evaluate ensemble on held-out test set (target: macro-F1 ≥ 0.80)
    4. Save all weights + gene list + label encoder to output directory
    5. (Optional) Push to HuggingFace receptor-bio/axiomtier1
    """
    args = _parse_args(argv)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    # ── step 1: dataset ────────────────────────────────────────────────────
    data_dir = out / "data"
    if not (data_dir / "train.npz").exists() or args.rebuild_data:
        logger.info("Step 1: downloading Census training data …")
        from axiom_sc.tier1.training.data_loader import build_training_dataset
        stats = build_training_dataset(
            output_dir=str(data_dir),
            census_version=args.census_version,
            tissue_filter=args.tissue_filter or None,
            max_cells=args.max_cells or None,
            n_hvg=args.n_hvg,
            random_seed=args.seed,
        )
        logger.info(
            "Dataset: %d train / %d val / %d test | %d types | %d genes",
            stats["n_cells_train"], stats["n_cells_val"], stats["n_cells_test"],
            stats["n_cell_types"], stats["n_genes"],
        )
    else:
        logger.info("Step 1: reusing existing dataset in %s", data_dir)

    # Load metadata
    with open(data_dir / "label_encoder.json") as f:
        label_encoder = {int(k): v for k, v in json.load(f).items()}
    gene_list = (data_dir / "axiomtier1_genes.txt").read_text().strip().splitlines()
    n_genes = len(gene_list)
    n_classes = len(label_encoder)

    logger.info("n_genes=%d  n_classes=%d", n_genes, n_classes)

    # ── step 2: train N_MODELS ─────────────────────────────────────────────
    val_f1_scores = []
    for model_idx in range(args.n_models):
        seed_i = args.seed + model_idx * 137   # deterministic per-model seed
        logger.info("─" * 60)
        logger.info("Training model %d/%d (seed=%d) …", model_idx + 1, args.n_models, seed_i)
        t0 = time.time()

        val_f1 = _train_single_model(
            model_idx=model_idx,
            seed=seed_i,
            data_dir=data_dir,
            output_dir=out,
            n_genes=n_genes,
            n_classes=n_classes,
            hidden=args.hidden,
            epochs=args.epochs,
            batch_size=args.batch_size,
            lr=args.lr,
            weight_decay=args.weight_decay,
            patience=args.patience,
        )
        elapsed = (time.time() - t0) / 3600
        logger.info(
            "Model %d done: val macro-F1=%.4f | %.1f h",
            model_idx, val_f1, elapsed,
        )
        val_f1_scores.append(val_f1)

    # ── step 3: evaluate ensemble on test set ──────────────────────────────
    logger.info("─" * 60)
    logger.info("Step 3: evaluating ensemble on held-out test set …")
    test_f1 = _evaluate_ensemble(
        data_dir=data_dir,
        weights_dir=out,
        n_models=args.n_models,
        n_genes=n_genes,
        n_classes=n_classes,
        hidden=args.hidden,
        label_encoder=label_encoder,
    )
    logger.info("Ensemble test macro-F1: %.4f (target ≥ %.2f)", test_f1, TARGET_MACRO_F1)

    if test_f1 < TARGET_MACRO_F1:
        logger.warning(
            "⚠ Ensemble F1=%.4f below target %.2f. "
            "Consider: more epochs, more data, or architecture changes.",
            test_f1, TARGET_MACRO_F1,
        )
    else:
        logger.info("✓ Target F1 met.")

    # ── step 4: copy gene list + label encoder to output root ─────────────
    import shutil
    shutil.copy(data_dir / "axiomtier1_genes.txt", out / "axiomtier1_genes.txt")
    shutil.copy(data_dir / "label_encoder.json", out / "label_encoder.json")

    # Save training summary
    summary = {
        "n_models": args.n_models,
        "n_genes": n_genes,
        "n_cell_types": n_classes,
        "val_f1_per_model": val_f1_scores,
        "mean_val_f1": float(np.mean(val_f1_scores)),
        "test_macro_f1": float(test_f1),
        "target_f1": TARGET_MACRO_F1,
        "hyperparams": vars(args),
    }
    (out / "training_summary.json").write_text(json.dumps(summary, indent=2))
    logger.info("Training complete. Weights saved to %s", out)

    # ── step 5: optional HuggingFace push ─────────────────────────────────
    if args.push_to_hub:
        _push_to_huggingface(out, repo_id=args.hub_repo)


# ── single model training ─────────────────────────────────────────────────────

def _train_single_model(
    model_idx: int,
    seed: int,
    data_dir: Path,
    output_dir: Path,
    n_genes: int,
    n_classes: int,
    hidden: int,
    epochs: int,
    batch_size: int,
    lr: float,
    weight_decay: float,
    patience: int,
) -> float:
    """Train one MLP and save weights. Returns validation macro-F1."""
    try:
        import torch
        import torch.nn as nn
        from torch.utils.data import DataLoader
        from sklearn.metrics import f1_score
    except ImportError:
        raise ImportError("torch + scikit-learn required: pip install 'axiom-sc[science]'")

    from axiom_sc.tier1.training.data_loader import CensusDataset
    from axiom_sc.tier1.axiomtier1 import AXIOMTier1Net

    torch.manual_seed(seed)
    np.random.seed(seed)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    logger.info("Device: %s", device)

    # DataLoaders
    train_ds = CensusDataset(str(data_dir / "train.npz"))
    val_ds   = CensusDataset(str(data_dir / "val.npz"))
    train_loader = DataLoader(train_ds, batch_size=batch_size, shuffle=True,
                              num_workers=4, pin_memory=True)
    val_loader   = DataLoader(val_ds, batch_size=batch_size * 4, shuffle=False,
                              num_workers=4, pin_memory=True)

    # Model + optimiser
    net = AXIOMTier1Net(n_genes=n_genes, n_classes=n_classes, hidden=hidden)
    net.net.to(device)
    optimizer = torch.optim.AdamW(net.parameters(), lr=lr, weight_decay=weight_decay)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=epochs)
    criterion = nn.CrossEntropyLoss()

    best_val_f1 = 0.0
    best_state = None
    patience_count = 0

    for epoch in range(epochs):
        # Train
        net.net.train()
        for X_batch, y_batch in train_loader:
            X_batch, y_batch = X_batch.to(device), y_batch.to(device)
            optimizer.zero_grad()
            logits = net(X_batch)
            loss = criterion(logits, y_batch)
            loss.backward()
            optimizer.step()
        scheduler.step()

        # Validate
        net.eval()
        all_preds, all_true = [], []
        with torch.no_grad():
            for X_batch, y_batch in val_loader:
                X_batch = X_batch.to(device)
                preds = net(X_batch).argmax(dim=1).cpu().numpy()
                all_preds.extend(preds)
                all_true.extend(y_batch.numpy())

        val_f1 = f1_score(all_true, all_preds, average="macro", zero_division=0)

        if (epoch + 1) % 5 == 0:
            logger.info("  epoch %3d/%d  val_macro_f1=%.4f", epoch + 1, epochs, val_f1)

        if val_f1 > best_val_f1:
            best_val_f1 = val_f1
            best_state = {k: v.cpu().clone() for k, v in net.state_dict().items()}
            patience_count = 0
        else:
            patience_count += 1
            if patience_count >= patience:
                logger.info("  Early stopping at epoch %d", epoch + 1)
                break

    # Save best weights
    weight_file = output_dir / f"axiomtier1_model_{model_idx}.pt"
    torch.save(best_state, weight_file)
    logger.info("  Saved weights → %s  (best val F1=%.4f)", weight_file, best_val_f1)

    return best_val_f1


# ── ensemble evaluation ───────────────────────────────────────────────────────

def _evaluate_ensemble(
    data_dir: Path,
    weights_dir: Path,
    n_models: int,
    n_genes: int,
    n_classes: int,
    hidden: int,
    label_encoder: dict,
) -> float:
    """Load all N models and evaluate ensemble on the test set. Returns macro-F1."""
    try:
        import torch
        import torch.nn.functional as F
        from torch.utils.data import DataLoader
        from sklearn.metrics import f1_score, classification_report
    except ImportError:
        raise ImportError("torch + scikit-learn required.")

    from axiom_sc.tier1.training.data_loader import CensusDataset
    from axiom_sc.tier1.axiomtier1 import AXIOMTier1Net

    test_ds = CensusDataset(str(data_dir / "test.npz"))
    test_loader = DataLoader(test_ds, batch_size=1024, shuffle=False, num_workers=4)

    # Load all models
    models = []
    for i in range(n_models):
        net = AXIOMTier1Net(n_genes=n_genes, n_classes=n_classes, hidden=hidden)
        sd = torch.load(weights_dir / f"axiomtier1_model_{i}.pt",
                        map_location="cpu", weights_only=True)
        net.load_state_dict(sd)
        net.eval()
        models.append(net)

    # Ensemble predict
    all_mean_probs = []
    all_true = []

    with torch.no_grad():
        for X_batch, y_batch in test_loader:
            batch_probs = []
            for net in models:
                logits = net(X_batch)
                probs = F.softmax(logits, dim=-1).numpy()
                batch_probs.append(probs)
            mean_probs = np.stack(batch_probs, axis=0).mean(axis=0)  # (batch, n_classes)
            all_mean_probs.append(mean_probs)
            all_true.extend(y_batch.numpy())

    mean_probs = np.vstack(all_mean_probs)
    preds = mean_probs.argmax(axis=1)
    macro_f1 = f1_score(all_true, preds, average="macro", zero_division=0)

    # Per-class report
    target_names = [label_encoder[i] for i in range(n_classes)]
    report = classification_report(
        all_true, preds, target_names=target_names, zero_division=0
    )
    logger.info("Classification report:\n%s", report)

    return float(macro_f1)


# ── HuggingFace push ──────────────────────────────────────────────────────────

def _push_to_huggingface(weights_dir: Path, repo_id: str) -> None:
    """Upload weights directory to HuggingFace Hub under Apache 2.0 license."""
    try:
        from huggingface_hub import HfApi
    except ImportError:
        raise ImportError("huggingface-hub required: pip install huggingface-hub")

    api = HfApi()
    logger.info("Pushing to HuggingFace %s …", repo_id)
    api.upload_folder(
        folder_path=str(weights_dir),
        repo_id=repo_id,
        repo_type="model",
        commit_message=f"AXIOMTier1 weights — Apache 2.0",
    )
    logger.info("✓ Pushed to https://huggingface.co/%s", repo_id)


# ── CLI argument parser ───────────────────────────────────────────────────────

def _parse_args(argv):
    p = argparse.ArgumentParser(
        description="Train AXIOMTier1 MLP ensemble on CELLxGENE Census."
    )
    p.add_argument("--output", default="weights/",
                   help="Output directory for weights and metadata.")
    p.add_argument("--n-models", type=int, default=DEFAULTS["n_models"],
                   help="Number of independent models in ensemble.")
    p.add_argument("--epochs", type=int, default=DEFAULTS["epochs"])
    p.add_argument("--batch-size", type=int, default=DEFAULTS["batch_size"])
    p.add_argument("--lr", type=float, default=DEFAULTS["lr"])
    p.add_argument("--weight-decay", type=float, default=DEFAULTS["weight_decay"])
    p.add_argument("--hidden", type=int, default=DEFAULTS["hidden"])
    p.add_argument("--patience", type=int, default=DEFAULTS["patience"])
    p.add_argument("--n-hvg", type=int, default=2000,
                   help="Number of HVGs for feature selection.")
    p.add_argument("--census-version", default="stable",
                   help="CELLxGENE Census version.")
    p.add_argument("--tissue-filter", nargs="+", default=None,
                   help="Limit to specific tissues (space-separated).")
    p.add_argument("--max-cells", type=int, default=None,
                   help="Hard cap on total cells (for dev runs).")
    p.add_argument("--rebuild-data", action="store_true",
                   help="Re-download Census data even if cached.")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--push-to-hub", action="store_true",
                   help="Push weights to HuggingFace after training.")
    p.add_argument("--hub-repo", default="receptor-bio/axiomtier1",
                   help="HuggingFace repo ID for push.")
    return p.parse_args(argv)


if __name__ == "__main__":
    main()


# ── streaming Census training (Colab / low-RAM workflow) ──────────────────────
#
# Public API used by the Google Colab notebook:
#
#   from axiom_sc.tier1.training.train import train_single_model
#   val_acc = train_single_model(model_idx=0, output_path="weights/model_0.pt",
#                                 output_dir="weights/")
#
# Design: streams raw counts from Census via tiledbsoma-ml (ExperimentDataset),
# normalises on-the-fly (10k + log1p), avoids storing multi-GB matrices on disk.
# Gene vocabulary (HVGs) and label encoder are computed once from a ~500k-cell
# sample and saved to output_dir; all subsequent model runs reuse them.

_VOCAB_SAMPLE_CELLS = 500_000   # cells used to compute HVGs (fits in ~6 GB RAM)
_EXCLUDED_LABELS = {            # same as data_loader.EXCLUDED_LABELS
    "unknown", "native cell", "", "cell", "eukaryotic cell",
}
_MIN_CELLS_INCLUDE = 100        # drop types with fewer cells in Census


def build_gene_vocabulary(
    output_dir: str | Path,
    census_version: str = "stable",
    n_hvg: int = 2000,
    sample_cells: int = _VOCAB_SAMPLE_CELLS,
    random_seed: int = 42,
) -> dict:
    """
    Build HVG gene list and label encoder from a Census sample.

    Saves four files to output_dir (idempotent — skipped if they exist):
      axiomtier1_genes.txt        — one HVG symbol per line
      axiomtier1_var_ids.json     — Census var soma_joinid for each HVG
      label_encoder.json          — {str(idx): cell_type_name}
      model_config.json           — architecture config for load_weights()

    Parameters
    ----------
    output_dir : str | Path
    census_version : str
    n_hvg : int
    sample_cells : int     — cells drawn for HVG computation
    random_seed : int

    Returns
    -------
    dict with keys: n_genes, n_classes, gene_list, label_encoder, var_ids
    """
    try:
        import cellxgene_census
        import scanpy as sc
    except ImportError:
        raise ImportError("pip install cellxgene-census scanpy")

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    gene_file    = out / "axiomtier1_genes.txt"
    var_id_file  = out / "axiomtier1_var_ids.json"
    encoder_file = out / "label_encoder.json"
    config_file  = out / "model_config.json"

    # Fast path: reuse if all four files exist
    if all(f.exists() for f in [gene_file, var_id_file, encoder_file, config_file]):
        logger.info("Gene vocabulary already built — reusing files in %s", out)
        gene_list       = gene_file.read_text().strip().splitlines()
        var_ids         = json.loads(var_id_file.read_text())
        label_encoder   = {int(k): v for k, v in json.loads(encoder_file.read_text()).items()}
        return dict(
            n_genes=len(gene_list), n_classes=len(label_encoder),
            gene_list=gene_list, label_encoder=label_encoder, var_ids=var_ids,
        )

    logger.info("Building gene vocabulary from Census (sample=%d cells) …", sample_cells)

    with cellxgene_census.open_soma(census_version=census_version) as census:
        # ── 1. Get obs metadata to build label encoder ─────────────────────
        obs_df = census["census_data"]["homo_sapiens"].obs.read(
            column_names=["soma_joinid", "cell_type", "is_primary_data"],
            value_filter="is_primary_data == True",
        ).concat().to_pandas()

        ct_counts = (
            obs_df.loc[~obs_df["cell_type"].isin(_EXCLUDED_LABELS), "cell_type"]
            .value_counts()
        )
        keep_types = set(ct_counts[ct_counts >= _MIN_CELLS_INCLUDE].index)
        obs_df = obs_df[obs_df["cell_type"].isin(keep_types)]

        sorted_types = sorted(keep_types)
        label_encoder = {i: ct for i, ct in enumerate(sorted_types)}

        logger.info("Label encoder: %d cell types", len(label_encoder))

        # ── 2. Sample cells for HVG computation ────────────────────────────
        rng = np.random.default_rng(random_seed)
        sample_ids = rng.choice(
            obs_df["soma_joinid"].values,
            min(sample_cells, len(obs_df)),
            replace=False,
        )

        logger.info("Fetching %d cells for HVG computation …", len(sample_ids))
        adata = cellxgene_census.get_anndata(
            census=census,
            organism="Homo sapiens",
            measurement_name="RNA",
            obs_coords=sample_ids.tolist(),
            obs_column_names=["cell_type"],
        )

    # ── 3. Preprocess + HVG selection ─────────────────────────────────────
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, flavor="seurat_v3")
    adata = adata[:, adata.var["highly_variable"]].copy()
    gene_list = list(adata.var_names)
    logger.info("HVG selection done: %d genes", len(gene_list))

    # ── 4. Get var soma_joinids for the HVGs ──────────────────────────────
    # Needed for the tiledbsoma-ml ExperimentDataset var_query coords.
    # adata.var index is feature_id (Ensembl); we need soma_joinid.
    with cellxgene_census.open_soma(census_version=census_version) as census:
        var_df = census["census_data"]["homo_sapiens"]["ms"]["RNA"]["var"].read(
            column_names=["soma_joinid", "feature_id", "feature_name"],
        ).concat().to_pandas()

    gene_set = set(gene_list)
    # gene_list may be feature_id or feature_name depending on Census version
    hvg_var = var_df[
        var_df["feature_id"].isin(gene_set) | var_df["feature_name"].isin(gene_set)
    ].copy()
    # Build ordered list matching gene_list order
    id_to_joinid = dict(zip(hvg_var["feature_id"], hvg_var["soma_joinid"]))
    name_to_joinid = dict(zip(hvg_var["feature_name"], hvg_var["soma_joinid"]))
    var_ids = []
    final_gene_list = []
    for g in gene_list:
        jid = id_to_joinid.get(g) or name_to_joinid.get(g)
        if jid is not None:
            var_ids.append(int(jid))
            final_gene_list.append(g)
    gene_list = final_gene_list
    logger.info("Mapped %d/%d HVGs to Census var soma_joinids", len(var_ids), n_hvg)

    # ── 5. Save files ──────────────────────────────────────────────────────
    gene_file.write_text("\n".join(gene_list))
    var_id_file.write_text(json.dumps(var_ids))
    encoder_file.write_text(json.dumps({str(k): v for k, v in label_encoder.items()}, indent=2))
    save_config(out, n_genes=len(gene_list), n_classes=len(label_encoder))
    logger.info("Vocabulary saved to %s", out)

    return dict(
        n_genes=len(gene_list), n_classes=len(label_encoder),
        gene_list=gene_list, label_encoder=label_encoder, var_ids=var_ids,
    )


def save_config(
    output_dir: str | Path,
    n_genes: int,
    n_classes: int,
    hidden: int = 512,
    n_models: int = 10,
) -> None:
    """Save model_config.json so AXIOMTier1Ensemble.load_weights() can reconstruct architecture."""
    cfg = {"hidden": hidden, "n_models": n_models, "n_genes": n_genes, "n_classes": n_classes}
    (Path(output_dir) / "model_config.json").write_text(json.dumps(cfg, indent=2))


def build_streaming_datapipe(
    census_version: str,
    var_soma_joinids: list[int],
    label_to_idx: dict,
    batch_size: int = 512,
    seed: int = 42,
    shuffle: bool = True,
):
    """
    Stream Census batches via tiledbsoma-ml (official ExperimentDataPipe successor).

    Requires: pip install tiledbsoma-ml

    Returns (loader, census) — caller MUST call census.close() when done
    (or use the context-manager helper _open_streaming_loader instead).

    Each batch from loader is a dict:
      batch["X"]         : torch.Tensor float32 (batch, n_genes) — raw counts
      batch["cell_type"] : list[str]

    Normalisation (10k + log1p) is applied inside train_single_model after
    receiving each batch, keeping this function simple and composable.

    Parameters
    ----------
    census_version : str
    var_soma_joinids : list[int]  — Census var soma_joinids for the HVG genes
    label_to_idx : dict           — {cell_type_str: int class index}
    batch_size : int
    seed : int
    shuffle : bool
    """
    try:
        from tiledbsoma_ml import ExperimentDataset, experiment_dataloader
        import cellxgene_census
        from tiledbsoma import AxisQuery
    except ImportError:
        raise ImportError(
            "pip install tiledbsoma-ml cellxgene-census tiledbsoma"
        )

    census = cellxgene_census.open_soma(census_version=census_version)
    experiment = census["census_data"]["homo_sapiens"]

    query = experiment.axis_query(
        measurement_name="RNA",
        obs_query=AxisQuery(value_filter="is_primary_data == True"),
        var_query=AxisQuery(coords=(var_soma_joinids,)),
    )

    dataset = ExperimentDataset(
        query,
        layer_name="raw",
        obs_column_names=["cell_type"],
        batch_size=batch_size,
        shuffle=shuffle,
        seed=seed,
    )

    loader = experiment_dataloader(dataset, num_workers=2)
    return loader, census   # caller must census.close() when done


def train_single_model(
    model_idx: int,
    output_path: str,
    output_dir: str,
    census_version: str = "stable",
    n_hvg: int = 2000,
    epochs: int = 5,
    batch_size: int = 512,
    lr: float = 1e-3,
    weight_decay: float = 1e-4,
    hidden: int = 512,
    patience: int = 3,
    seed: int = 42,
) -> float:
    """
    Train one model in the AXIOMTier1 ensemble using streaming Census data.

    This is the Colab-facing entry point. Call from Cell 3:

        val_acc = train_single_model(
            model_idx=MODEL_IDX,
            output_path=f"{WEIGHTS_DIR}/model_{MODEL_IDX}.pt",
            output_dir=WEIGHTS_DIR,
        )

    On the first run it builds the HVG vocabulary (once, ~15 min).
    Subsequent models reuse the saved vocabulary files.
    Progress is logged every 500 batches.

    Parameters
    ----------
    model_idx : int
        Index of this model in the ensemble (0–9).
    output_path : str
        Path to save the weight file (.pt).
    output_dir : str
        Directory for vocabulary files (axiomtier1_genes.txt etc.).
        Also used as a fallback to load pre-built npz splits if present.
    census_version : str
    n_hvg : int
    epochs : int
    batch_size : int
    lr : float
    weight_decay : float
    hidden : int
    patience : int   — epochs without val F1 improvement before early stopping
    seed : int

    Returns
    -------
    float — validation macro-F1 of the best checkpoint.
    """
    try:
        import torch
        import torch.nn as nn
        from sklearn.metrics import f1_score
    except ImportError:
        raise ImportError("pip install torch scikit-learn")

    from axiom_sc.tier1.axiomtier1 import AXIOMTier1Net

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    model_seed = seed + model_idx * 137
    torch.manual_seed(model_seed)
    np.random.seed(model_seed)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    logger.info("train_single_model idx=%d  device=%s  seed=%d", model_idx, device, model_seed)

    # ── 1. Vocabulary (build once, reuse across models) ───────────────────
    vocab = build_gene_vocabulary(
        output_dir=out_dir,
        census_version=census_version,
        n_hvg=n_hvg,
        random_seed=seed,
    )
    n_genes   = vocab["n_genes"]
    n_classes = vocab["n_classes"]
    label_encoder  = vocab["label_encoder"]          # {int: cell_type_str}
    label_to_idx   = {v: k for k, v in label_encoder.items()}
    var_ids        = vocab["var_ids"]

    logger.info("n_genes=%d  n_classes=%d", n_genes, n_classes)
    save_config(out_dir, n_genes=n_genes, n_classes=n_classes, hidden=hidden)

    # ── 2. Model + optimiser ──────────────────────────────────────────────
    net = AXIOMTier1Net(n_genes=n_genes, n_classes=n_classes, hidden=hidden)
    net.net.to(device)
    optimizer = torch.optim.AdamW(net.parameters(), lr=lr, weight_decay=weight_decay)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=epochs)
    criterion = nn.CrossEntropyLoss()

    # ── 3. Training loop ──────────────────────────────────────────────────
    best_val_f1    = 0.0
    best_state     = None
    patience_count = 0

    import torch.nn.functional as F

    def _normalise(X_raw: "torch.Tensor") -> "torch.Tensor":
        """10k normalise + log1p on a raw-count batch tensor."""
        X = X_raw.float()
        row_sums = X.sum(dim=1, keepdim=True).clamp(min=1.0)
        return torch.log1p(X * (1e4 / row_sums))

    def _labels_from_batch(batch: dict) -> "torch.Tensor | None":
        """Map cell_type strings to integer indices; return None if all unknown."""
        cell_types = batch["cell_type"]
        idxs = [label_to_idx.get(ct, -1) for ct in cell_types]
        y = torch.tensor(idxs, dtype=torch.long)
        if (y < 0).all():
            return None
        return y

    for epoch in range(epochs):
        # ── train ──────────────────────────────────────────────────────────
        net.net.train()
        train_loader, train_census = build_streaming_datapipe(
            census_version=census_version,
            var_soma_joinids=var_ids,
            label_to_idx=label_to_idx,
            batch_size=batch_size,
            seed=model_seed + epoch,   # different shuffle each epoch
            shuffle=True,
        )
        try:
            train_loss_sum = 0.0
            n_batches = 0
            for batch_idx, batch in enumerate(train_loader):
                y_batch = _labels_from_batch(batch)
                if y_batch is None:
                    continue
                X_batch = _normalise(batch["X"])
                valid = y_batch >= 0
                X_batch, y_batch = X_batch[valid].to(device), y_batch[valid].to(device)

                optimizer.zero_grad()
                loss = criterion(net(X_batch), y_batch)
                loss.backward()
                optimizer.step()

                train_loss_sum += loss.item()
                n_batches += 1
                if (batch_idx + 1) % 500 == 0:
                    logger.info(
                        "  epoch %d  batch %d  avg_loss=%.4f",
                        epoch + 1, batch_idx + 1, train_loss_sum / n_batches,
                    )
        finally:
            train_census.close()   # always release the Census connection

        scheduler.step()

        # ── validate (no shuffle, first _MAX_VAL_BATCHES batches) ─────────
        net.eval()
        val_loader, val_census = build_streaming_datapipe(
            census_version=census_version,
            var_soma_joinids=var_ids,
            label_to_idx=label_to_idx,
            batch_size=batch_size * 4,
            seed=model_seed,
            shuffle=False,
        )
        all_preds, all_true = [], []
        _MAX_VAL_BATCHES = 100   # ~50k cells — enough for reliable macro-F1
        try:
            with torch.no_grad():
                for _vi, batch in enumerate(val_loader):
                    if _vi >= _MAX_VAL_BATCHES:
                        break
                    y_batch = _labels_from_batch(batch)
                    if y_batch is None:
                        continue
                    X_batch = _normalise(batch["X"])
                    valid = y_batch >= 0
                    preds = net(X_batch[valid].to(device)).argmax(dim=1).cpu().numpy()
                    all_preds.extend(preds)
                    all_true.extend(y_batch[valid].numpy())
        finally:
            val_census.close()

        val_f1 = (
            f1_score(all_true, all_preds, average="macro", zero_division=0)
            if all_true else 0.0
        )
        logger.info(
            "epoch %d/%d  val_macro_f1=%.4f  train_loss=%.4f  n_val_cells=%d",
            epoch + 1, epochs, val_f1,
            train_loss_sum / max(n_batches, 1), len(all_true),
        )

        if val_f1 > best_val_f1:
            best_val_f1    = val_f1
            best_state     = {k: v.cpu().clone() for k, v in net.state_dict().items()}
            patience_count = 0
        else:
            patience_count += 1
            if patience_count >= patience:
                logger.info("Early stopping at epoch %d", epoch + 1)
                break

    # ── 4. Save best weights ──────────────────────────────────────────────
    weight_path = Path(output_path)
    weight_path.parent.mkdir(parents=True, exist_ok=True)
    torch.save(best_state, weight_path)
    logger.info("Saved model_%d → %s  (best val_f1=%.4f)", model_idx, weight_path, best_val_f1)

    return best_val_f1
