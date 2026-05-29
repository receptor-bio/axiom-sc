# CHANGELOG

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
