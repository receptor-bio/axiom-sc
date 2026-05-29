# CHANGELOG

All notable changes to `axiom-sc` are documented here.
Versions are published to [PyPI](https://pypi.org/project/axiom-sc/).

---

## 0.2.6 (2026-05-29)

### Fixed
- **`_labels_from_batch()` TypeError** — `tiledbsoma-ml` yields `MiniBatch`
  tuples `(X, obs_df)`, not dicts. Fixed `batch["cell_type"]` → `batch[1]["cell_type"]`
  and `batch["X"]` → `batch[0]` throughout the training loop.
- **`_normalise()` TypeError** — `batch[0]` arrives as `np.ndarray` or
  `scipy.csr_matrix`, not a `torch.Tensor`. Fixed to handle all three types
  (sparse → dense → tensor) before applying 10k normalisation + log1p.

### Upgrade notes
```python
# Colab Cell 1
!pip install "axiom-sc[science]==0.2.6" --upgrade --quiet
```

---

## 0.2.5 (2026-05-29)

### Fixed
- **`tiledbsoma-ml` install failure** — version constraint was `>=1.0` but the
  package only exists at `0.1.0` on PyPI. Fixed to `>=0.1`.
- **Redundant dependency** — removed `tiledbsoma>=1.12` from `[science]` extra;
  it is already a transitive dependency of `cellxgene-census` and does not need
  to be listed separately.

---

## 0.2.4 (2026-05-29)

### Fixed
- **Stale comments** — removed two references to `ExperimentDataPipe` in
  `train.py` comments left over from the 0.2.2 implementation.
- **Missing `[science]` deps** — `tiledbsoma-ml` was added to the `[science]`
  optional extra in `pyproject.toml`. Without this, `pip install "axiom-sc[science]"`
  would succeed but `build_streaming_datapipe()` would immediately raise
  `ImportError: No module named 'tiledbsoma_ml'`.

---

## 0.2.3 (2026-05-29)

### Fixed
- **`ModuleNotFoundError: No module named 'torchdata.datapipes'`** —
  `cellxgene_census.experimental.ml.ExperimentDataPipe` (used in 0.2.2) depends
  on `torchdata` which was removed from PyTorch and is not available on
  Python 3.12. Replaced with `tiledbsoma-ml` (`ExperimentDataset` +
  `experiment_dataloader`), the official successor maintained by the
  TileDB-SOMA team.
- **Census connection leak** — `train_single_model()` now wraps each epoch's
  train and validation loops in `try/finally` blocks that call `census.close()`,
  ensuring the Census connection is always released even if training fails
  mid-epoch.

---

## 0.2.2 (2026-05-29)

### Added
- **Streaming Census training pipeline** for Google Colab and low-RAM GPU
  environments. No full-dataset download required — raw counts stream directly
  from CELLxGENE Census and are normalised on-the-fly (10k + log1p).

- **`train_single_model(model_idx, output_path, output_dir, ...)`** — public
  Colab entry point. Trains one MLP in the AXIOMTier1 ensemble. Call from
  Colab Cell 3:
  ```python
  from axiom_sc.tier1.training.train import train_single_model
  val_acc = train_single_model(
      model_idx=MODEL_IDX,
      output_path=f"{WEIGHTS_DIR}/model_{MODEL_IDX}.pt",
      output_dir=WEIGHTS_DIR,
  )
  ```

- **`build_gene_vocabulary(output_dir, ...)`** — computes HVG gene list and
  label encoder from a ~500k-cell Census sample. Saves four files to
  `output_dir` (`axiomtier1_genes.txt`, `axiomtier1_var_ids.json`,
  `label_encoder.json`, `model_config.json`). Idempotent — reuses saved files
  on subsequent model runs so HVG computation only runs once.

- **`save_config(output_dir, n_genes, n_classes, hidden, n_models)`** — writes
  `model_config.json` required by `AXIOMTier1Ensemble.load_weights()` to
  reconstruct the correct architecture.

- **`build_streaming_datapipe(...)`** — streams Census data via
  `tiledbsoma-ml`. Returns `(loader, census)` — caller closes census when done.

### Known issues in 0.2.2 (fixed in 0.2.3)
- `build_streaming_datapipe()` used `ExperimentDataPipe` which requires
  `torchdata` — not available on Python 3.12 / modern PyTorch.

---

## 0.2.1 (2026-05-28)

### Fixed
- **Missing KG data in wheel** — `oracle_kg_v0.2.0.json` was not bundled inside
  the PyPI wheel in 0.2.0, causing `FileNotFoundError` on import. Fixed by
  adding the `kg_data/` directory to the hatchling build include list.

---

## 0.2.0 (2026-05-27)

### Added
- AXIOM-SC v0.2.0 — full 5-tier mechanistic cell type annotation SDK.
- **Knowledge Graph**: 640 ACTIVE rules / 198 cell types across 28 tissues.
  All PMIDs verified via NCBI Entrez. Negative and circuit rules are the
  core AXIOM innovation (proof-by-contradiction).
- **Profile system**: `oss-mit`, `oss-apache` (default), `commercial` built-in
  profiles + custom `axiom_profile.json` support.
- **Tier 1**: AXIOMTier1 MLP ensemble inference code + CensusAnnotator kNN
  fallback (weights pending GPU training job).
- **Tier 2**: AxiomAnnotator with all 6 Phase 1 corrections (circuit
  completeness, support counting reform, PROVEN co-confirmation, PDC/ILC3
  fixes). pySCENIC via GPL-isolated subprocess (`scenic-env` conda env).
- **Tier 3**: scVelo velocity, Signac chromatin, COMMOT communication,
  CELLama spatial niche, OrthoFinder cross-species — multi-stream convergence.
- **Tier 4**: LLM ensemble (Anthropic + OpenAI) + CellMarker 2.0 RAG.
- **Tier 5**: Novel attractor discovery + Tier 5 → KG feedback loop.
- CI via GitHub Actions; OIDC trusted publishing to PyPI (no stored secrets).

---

## 0.1.0 (Phase 1 baseline)

- AXIOM KG engine with proof-by-contradiction: 41 rules, 18 cell types.
- Validated on Human Thymus Cell Atlas, Lung Cell Atlas, Tabula Sapiens.
- KG in-scope accuracy: 22.2% vs 5.6% (CASSIA) and 0% (mLLMCelltype).
- PROVEN precision: 67% (2/3 independently validated).
