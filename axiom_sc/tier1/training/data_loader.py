"""
CELLxGENE Census data loader for AXIOMTier1 training.

Downloads human single-cell RNA-seq data from CELLxGENE Census,
applies preprocessing, and produces stratified train/val/test splits.

Data source:
  CELLxGENE Census: Chanzuckerberg Initiative (2023)
  MIT license API; CC BY 4.0 underlying data
  doi: 10.1126/science.abl4896

Preprocessing pipeline:
  1. Normalize to 10k counts (scanpy normalize_total)
  2. log1p transform
  3. Select top 2000 highly variable genes (HVGs)
  4. Stratified train/val/test split (≥2000 cells per type in training)

License: Apache 2.0
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)

# ── constants ─────────────────────────────────────────────────────────────────

N_HVG = 2000                    # highly variable genes selected
TARGET_SUM = 1e4                # normalize_total target
MIN_CELLS_PER_TYPE_TRAIN = 2000 # minimum training cells per cell type
MIN_CELLS_PER_TYPE_INCLUDE = 100 # exclude types with fewer cells in Census
TRAIN_FRAC = 0.80
VAL_FRAC   = 0.10
TEST_FRAC  = 0.10

CENSUS_ORGANISM = "Homo sapiens"
CELL_TYPE_COLUMN = "cell_type"  # harmonised CL label in Census obs

# Exclude low-quality or ambiguous labels
EXCLUDED_LABELS = {
    "unknown",
    "native cell",
    "",
    "cell",
    "eukaryotic cell",
}


# ── main download + preprocess function ───────────────────────────────────────

def build_training_dataset(
    output_dir: str,
    census_version: str = "stable",
    tissue_filter: list[str] | None = None,
    max_cells: int | None = None,
    n_hvg: int = N_HVG,
    random_seed: int = 42,
) -> dict:
    """
    Download Census data, preprocess, and save train/val/test splits.

    Parameters
    ----------
    output_dir : str
        Directory to save split files and metadata.
    census_version : str
        Census version ("stable" or date string).
    tissue_filter : list[str] | None
        Limit to specific tissues. None = all tissues.
    max_cells : int | None
        Hard cap on total cells (for development runs). None = all.
    n_hvg : int
        Number of HVGs to select.
    random_seed : int

    Returns
    -------
    dict with keys:
        'n_cells_train', 'n_cells_val', 'n_cells_test',
        'n_cell_types', 'n_genes', 'gene_list', 'label_encoder',
        'output_dir'
    """
    try:
        import cellxgene_census
        import scanpy as sc
        import anndata as ad
        import pandas as pd
        import scipy.sparse as sp
    except ImportError:
        raise ImportError(
            "Science stack required: pip install 'axiom-sc[science]'"
        )

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    logger.info("Opening Census version=%s …", census_version)

    with cellxgene_census.open_soma(census_version=census_version) as census:
        # Build obs filter
        obs_value_filter = "is_primary_data == True"
        if tissue_filter:
            tissue_str = " or ".join(f'tissue_general == "{t}"' for t in tissue_filter)
            obs_value_filter += f" and ({tissue_str})"

        logger.info("Fetching AnnData from Census …")
        adata = cellxgene_census.get_anndata(
            census=census,
            organism=CENSUS_ORGANISM,
            measurement_name="RNA",
            obs_value_filter=obs_value_filter,
            obs_column_names=[CELL_TYPE_COLUMN, "tissue_general", "dataset_id"],
            var_value_filter=None,  # all genes initially
        )

    logger.info("Downloaded %d cells × %d genes", adata.n_obs, adata.n_vars)

    # ── filter cell types ──────────────────────────────────────────────────
    adata = _filter_cell_types(adata, CELL_TYPE_COLUMN, MIN_CELLS_PER_TYPE_INCLUDE)
    logger.info("After cell-type filter: %d cells, %d types",
                adata.n_obs, adata.obs[CELL_TYPE_COLUMN].nunique())

    # ── optional hard cap ──────────────────────────────────────────────────
    if max_cells and adata.n_obs > max_cells:
        idx = np.random.default_rng(random_seed).choice(
            adata.n_obs, max_cells, replace=False
        )
        adata = adata[idx].copy()
        logger.info("Subsampled to %d cells (max_cells=%d)", adata.n_obs, max_cells)

    # ── preprocessing ──────────────────────────────────────────────────────
    logger.info("Normalising and selecting HVGs …")
    sc.pp.normalize_total(adata, target_sum=TARGET_SUM)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, flavor="seurat_v3")
    adata = adata[:, adata.var["highly_variable"]].copy()
    gene_list = list(adata.var_names)
    logger.info("HVG selection: %d genes", len(gene_list))

    # ── label encoder ──────────────────────────────────────────────────────
    cell_types = sorted(adata.obs[CELL_TYPE_COLUMN].unique())
    label_encoder = {i: ct for i, ct in enumerate(cell_types)}
    label_to_idx = {ct: i for i, ct in label_encoder.items()}
    adata.obs["label_idx"] = adata.obs[CELL_TYPE_COLUMN].map(label_to_idx).astype(int)

    # ── stratified split ───────────────────────────────────────────────────
    logger.info("Building stratified train/val/test split …")
    train_idx, val_idx, test_idx = _stratified_split(
        adata.obs["label_idx"].values,
        train_frac=TRAIN_FRAC, val_frac=VAL_FRAC,
        min_per_type=MIN_CELLS_PER_TYPE_TRAIN,
        random_seed=random_seed,
    )

    # ── save splits ────────────────────────────────────────────────────────
    _save_split(adata, train_idx, out / "train.npz", "train")
    _save_split(adata, val_idx,   out / "val.npz",   "val")
    _save_split(adata, test_idx,  out / "test.npz",  "test")

    # Save gene list and label encoder
    import json
    (out / "axiomtier1_genes.txt").write_text("\n".join(gene_list))
    (out / "label_encoder.json").write_text(
        json.dumps({str(k): v for k, v in label_encoder.items()}, indent=2)
    )

    stats = {
        "n_cells_train": len(train_idx),
        "n_cells_val":   len(val_idx),
        "n_cells_test":  len(test_idx),
        "n_cell_types":  len(cell_types),
        "n_genes":       len(gene_list),
        "gene_list":     gene_list,
        "label_encoder": label_encoder,
        "output_dir":    str(out),
    }
    logger.info(
        "Dataset ready: %d train / %d val / %d test | %d types | %d genes",
        stats["n_cells_train"], stats["n_cells_val"], stats["n_cells_test"],
        stats["n_cell_types"], stats["n_genes"],
    )
    return stats


# ── helpers ───────────────────────────────────────────────────────────────────

def _filter_cell_types(adata, col: str, min_cells: int):
    """Remove excluded labels and types with fewer than min_cells cells."""
    mask = ~adata.obs[col].isin(EXCLUDED_LABELS)
    adata = adata[mask].copy()
    counts = adata.obs[col].value_counts()
    keep = counts[counts >= min_cells].index
    mask2 = adata.obs[col].isin(keep)
    return adata[mask2].copy()


def _stratified_split(
    labels: np.ndarray,
    train_frac: float,
    val_frac: float,
    min_per_type: int,
    random_seed: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Stratified split ensuring at least min_per_type cells per type in training.

    Returns arrays of integer indices into the original array.
    """
    rng = np.random.default_rng(random_seed)
    unique_labels = np.unique(labels)

    train_idx, val_idx, test_idx = [], [], []

    for label in unique_labels:
        idx = np.where(labels == label)[0]
        rng.shuffle(idx)

        n = len(idx)
        # Guarantee at least min_per_type in training, then split remaining
        n_train = max(min_per_type, int(n * train_frac))
        n_train = min(n_train, n)   # can't exceed total
        n_remaining = n - n_train
        n_val = int(n_remaining * (val_frac / (val_frac + (1 - train_frac - val_frac))))
        n_val = max(1, min(n_val, n_remaining))
        n_test = n_remaining - n_val

        train_idx.extend(idx[:n_train])
        val_idx.extend(idx[n_train:n_train + n_val])
        test_idx.extend(idx[n_train + n_val:n_train + n_val + n_test])

    return (
        np.array(train_idx),
        np.array(val_idx),
        np.array(test_idx),
    )


def _save_split(adata, indices: np.ndarray, path: Path, split_name: str) -> None:
    """Save one split as a compressed .npz file."""
    import scipy.sparse as sp

    subset = adata[indices]
    X = subset.X
    if sp.issparse(X):
        X = X.toarray()
    y = subset.obs["label_idx"].values.astype(np.int32)

    np.savez_compressed(str(path), X=X.astype(np.float32), y=y)
    logger.info("Saved %s split: %d cells → %s", split_name, len(indices), path)


# ── dataset class for PyTorch DataLoader ─────────────────────────────────────

class CensusDataset:
    """
    PyTorch Dataset wrapping a .npz split file.

    Usage::
        ds = CensusDataset("data/train.npz")
        loader = DataLoader(ds, batch_size=256, shuffle=True)
    """

    def __init__(self, npz_path: str):
        data = np.load(npz_path)
        self._X = data["X"]   # (n_cells, n_genes) float32
        self._y = data["y"]   # (n_cells,) int32

    def __len__(self) -> int:
        return len(self._y)

    def __getitem__(self, idx):
        try:
            import torch
        except ImportError:
            raise ImportError("torch required for DataLoader usage.")
        return (
            torch.from_numpy(self._X[idx]),
            torch.tensor(int(self._y[idx]), dtype=torch.long),
        )
