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
