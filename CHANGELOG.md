# CHANGELOG

## 0.2.5 (2026-05-29)

### Fixed
- Drop redundant `tiledbsoma>=1.12` from `[science]` extra — already a
  transitive dependency of `cellxgene-census`. Keep only `tiledbsoma-ml>=0.1`
  which is the actual new dependency for `ExperimentDataset`.

## 0.2.4 (2026-05-29)

### Fixed
- Removed stale `ExperimentDataPipe` references from comments in `train.py`.
- Added `tiledbsoma>=1.12` and `tiledbsoma-ml>=1.0` to `[science]` optional
  extra in `pyproject.toml` — previously missing, causing `ImportError` when
  running `pip install "axiom-sc[science]"` and then calling
  `build_streaming_datapipe()`.

## 0.2.3 (2026-05-29)

### Fixed
- `build_streaming_datapipe()` rewritten to use `tiledbsoma-ml`
  (`ExperimentDataset` + `experiment_dataloader`) instead of the removed
  `cellxgene_census.experimental.ml.ExperimentDataPipe` which required
  `torchdata.datapipes` (not available in Python 3.12 / modern PyTorch).
- `train_single_model()` updated: applies 10k + log1p normalisation inline
  per batch, filters unknown cell types via `_labels_from_batch()`, and
  always closes the Census connection in a `finally` block after each epoch.

## 0.2.2 (2026-05-29)

### Added
- `train_single_model()` — public Colab-facing entry point for training one
  model in the AXIOMTier1 ensemble using streaming CELLxGENE Census data.
  No full-dataset download required; raw counts are normalised on-the-fly
  (10k + log1p) via `ExperimentDataPipe`.
- `build_gene_vocabulary()` — builds HVG gene list and label encoder from a
  configurable Census sample (~500k cells default). Files are saved to
  `output_dir` and reused across all 10 model runs.
- `save_config()` — writes `model_config.json` so `AXIOMTier1Ensemble.load_weights()`
  can reconstruct the architecture from disk.
- `build_streaming_datapipe()` — wraps `cellxgene_census.experimental.ml.ExperimentDataPipe`
  with per-batch 10k normalisation, log1p, and cell-type label filtering.
  Falls back gracefully when the experimental API is unavailable.

## 0.2.1 (2026-05-28)

### Fixed
- Bundle `oracle_kg_v0.2.0.json` inside the PyPI wheel (was missing from 0.2.0).

## 0.2.0 (2026-05-27)

### Added
- AXIOM-SC v0.2.0 — 5-tier mechanistic cell type annotation SDK.
- KG: 640 ACTIVE rules / 198 cell types (all PMIDs verified).
- Profile system: oss-mit, oss-apache (default), commercial + custom JSON profiles.
- AXIOMTier1 ensemble inference code (weights pending GPU training).
- CensusAnnotator kNN fallback (Tier 1 primary until weights arrive).
- pySCENIC subprocess isolation (GPL-clean, scenic-env conda env).
- Tier 3: velocity, chromatin, communication, spatial_niche, cross_species convergence.
- Tier 4: LLM ensemble (Anthropic + OpenAI) + CellMarker 2.0 RAG.
- Tier 5: novel attractor discovery + Tier 5 → KG feedback loop.
- CI via GitHub Actions; OIDC trusted publishing to PyPI.

## 0.1.0 (Phase 1 baseline)

- Phase 1: AXIOM KG engine with proof-by-contradiction (41 rules, 18 cell types).
- Validated on Human Thymus Cell Atlas, Lung Cell Atlas, Tabula Sapiens.
- 22.2% in-scope accuracy vs 5.6% CASSIA, 0% mLLMCelltype.
