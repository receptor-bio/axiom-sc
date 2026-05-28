"""
Post-training evaluation and calibration checks for AXIOMTier1.

Generates:
  - Macro-F1 and per-class F1 breakdown
  - Confidence calibration plot (confidence vs. accuracy)
  - Acceptance rate at each confidence threshold
  - Comparison baseline: Census kNN annotator

Usage::
    python -m axiom_sc.tier1.training.evaluate \\
        --weights weights/ \\
        --test-data weights/data/test.npz \\
        --output evaluation/

License: Apache 2.0
"""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)


def evaluate(
    weights_dir: str,
    test_npz: str,
    output_dir: str,
    n_bins: int = 10,
) -> dict:
    """
    Full evaluation of the AXIOMTier1Ensemble on a held-out test set.

    Parameters
    ----------
    weights_dir : str
        Directory with model weights, gene list, and label encoder.
    test_npz : str
        Path to test.npz split file.
    output_dir : str
        Where to save evaluation outputs.
    n_bins : int
        Number of confidence bins for calibration plot.

    Returns
    -------
    dict with macro_f1, calibration, acceptance_rates, per_class_f1
    """
    try:
        import torch
        import torch.nn.functional as F
        from sklearn.metrics import f1_score
    except ImportError:
        raise ImportError("torch + scikit-learn required.")

    from axiom_sc.tier1.axiomtier1 import AXIOMTier1Ensemble
    from axiom_sc.tier1.training.data_loader import CensusDataset
    from torch.utils.data import DataLoader

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # Load ensemble
    ensemble = AXIOMTier1Ensemble(weights_path=weights_dir)
    ensemble.load_weights(weights_dir)

    # Load label encoder
    with open(Path(weights_dir) / "label_encoder.json") as f:
        label_encoder = {int(k): v for k, v in json.load(f).items()}
    gene_list = (Path(weights_dir) / "axiomtier1_genes.txt").read_text().strip().splitlines()

    # Load test data
    test_ds = CensusDataset(test_npz)
    test_loader = DataLoader(test_ds, batch_size=1024, shuffle=False)

    # Collect predictions
    all_mean_probs, all_true = [], []
    with torch.no_grad():
        for X_batch, y_batch in test_loader:
            batch_probs = ensemble._run_ensemble(X_batch.numpy())  # (N, batch, C)
            mean_probs = batch_probs.mean(axis=0)
            all_mean_probs.append(mean_probs)
            all_true.extend(y_batch.numpy())

    mean_probs = np.vstack(all_mean_probs)   # (n_test, n_classes)
    y_true = np.array(all_true)
    y_pred = mean_probs.argmax(axis=1)

    from axiom_sc.tier1.axiomtier1 import _compute_confidence
    confidence = _compute_confidence(mean_probs)

    # ── overall metrics ────────────────────────────────────────────────────
    macro_f1 = float(f1_score(y_true, y_pred, average="macro", zero_division=0))
    per_class_f1 = f1_score(y_true, y_pred, average=None, zero_division=0)
    per_class = {
        label_encoder[i]: float(per_class_f1[i])
        for i in range(len(label_encoder))
    }

    logger.info("Macro-F1: %.4f", macro_f1)

    # ── confidence calibration ─────────────────────────────────────────────
    bin_edges = np.linspace(0, 1, n_bins + 1)
    calibration = []
    for lo, hi in zip(bin_edges[:-1], bin_edges[1:]):
        mask = (confidence >= lo) & (confidence < hi)
        if mask.sum() == 0:
            continue
        acc = (y_pred[mask] == y_true[mask]).mean()
        calibration.append({
            "bin_lower": float(lo),
            "bin_upper": float(hi),
            "mean_confidence": float(confidence[mask].mean()),
            "accuracy": float(acc),
            "n_cells": int(mask.sum()),
        })

    # ── acceptance rate by threshold ───────────────────────────────────────
    thresholds = [0.5, 0.6, 0.7, 0.75, 0.80, 0.85, 0.90, 0.95]
    acceptance_rates = {}
    for thresh in thresholds:
        accepted = (confidence >= thresh).mean()
        accepted_correct = (
            (confidence >= thresh) & (y_pred == y_true)
        ).sum() / max((confidence >= thresh).sum(), 1)
        acceptance_rates[str(thresh)] = {
            "acceptance_rate": float(accepted),
            "precision_at_threshold": float(accepted_correct),
        }

    # Log the two key operating points
    for t in [0.50, 0.85]:
        ar = acceptance_rates[str(t)]
        logger.info(
            "  thresh=%.2f → accept_rate=%.1f%%  precision=%.1f%%",
            t, 100 * ar["acceptance_rate"], 100 * ar["precision_at_threshold"],
        )

    # ── save results ───────────────────────────────────────────────────────
    results = {
        "macro_f1": macro_f1,
        "target_f1": 0.80,
        "meets_target": macro_f1 >= 0.80,
        "n_cell_types": len(label_encoder),
        "n_test_cells": len(y_true),
        "calibration": calibration,
        "acceptance_rates": acceptance_rates,
        "per_class_f1": per_class,
    }

    (out / "evaluation_results.json").write_text(json.dumps(results, indent=2))
    logger.info("Evaluation saved → %s", out / "evaluation_results.json")

    return results


def _parse_args():
    p = argparse.ArgumentParser(description="Evaluate AXIOMTier1 ensemble.")
    p.add_argument("--weights", required=True, help="Weights directory.")
    p.add_argument("--test-data", required=True, help="Path to test.npz.")
    p.add_argument("--output", default="evaluation/", help="Output directory.")
    p.add_argument("--n-bins", type=int, default=10, help="Calibration bins.")
    return p.parse_args()


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    args = _parse_args()
    evaluate(
        weights_dir=args.weights,
        test_npz=args.test_data,
        output_dir=args.output,
        n_bins=args.n_bins,
    )
