# AXIOM-SC

**Mechanistic cell type annotation for single-cell RNA-seq via proof-by-contradiction**

[![PyPI](https://img.shields.io/pypi/v/axiom-sc)](https://pypi.org/project/axiom-sc/)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue)](LICENSE)
[![CI](https://github.com/receptor-bio/axiom-sc/actions/workflows/ci.yml/badge.svg)](https://github.com/receptor-bio/axiom-sc/actions/workflows/ci.yml)
[![Coverage](https://img.shields.io/badge/coverage-84%25-green)]()

---

## What is AXIOM-SC?

AXIOM-SC is a **5-tier routing system** for single-cell RNA-seq cell type annotation. It is not a classifier — each cell is routed to the *minimum tier required* to annotate it confidently.

**The core innovation:** proof-by-contradiction using mechanistic biological rules. A single hard rule violation eliminates a candidate cell type, regardless of how many positive markers support it. This resolves a fundamental failure shared by all LLM-based annotation systems (CASSIA, mLLMCelltype): they rely exclusively on positive marker matching and cannot rule out candidates mechanistically.

> *One violated rule eliminates a candidate. No amount of supporting evidence rescues it.*

---

## The 5-Tier Routing Pipeline

```
┌─────────────────────────────────────────────────────────────────────────┐
│                         Input: scRNA-seq AnnData                        │
│              (clustered h5ad · raw counts · CellRanger .h5 · marker CSV)│
└─────────────────────────────────┬───────────────────────────────────────┘
                                  │
                    ┌─────────────▼──────────────┐
                    │         TIER 1             │  ~75–78% of cells
                    │      AXIOMTier1            │  ──────────────────
                    │   MLP Ensemble (×10)       │  Trained on 22M+ cells
                    │   CELLxGENE Census kNN     │  (CELLxGENE Census)
                    │                            │
                    │  confidence ≥ 0.85 → ACCEPT│
                    │  0.50–0.85 → Tier 2 verify │
                    │  < 0.50   → Tier 2 full    │
                    └──────────┬─────────────────┘
                               │ uncertain / low-confidence
                    ┌──────────▼──────────────────┐
                    │         TIER 2              │  ~15% of cells
                    │   AXIOM KG Engine           │  ────────────────
                    │   + pySCENIC (subprocess)   │  Mechanistic rules
                    │                             │  640 rules · 198 types
                    │  Proof-by-contradiction:    │  Orthogonal to Tier 1 —
                    │  PROVEN · UNCERTAIN ·       │  one violation = ruled out
                    │  CONTRADICTED               │  regardless of markers
                    └──────────┬──────────────────┘
                               │ UNCERTAIN
                    ┌──────────▼──────────────────┐
                    │         TIER 3              │  ~5–7% of cells
                    │   Multi-Stream Convergence  │  ────────────────
                    │                             │  6 orthogonal streams:
                    │  velocity  · chromatin      │  4+ streams agree → PROVEN
                    │  L-R comm  · spatial niche  │  resolves exhausted vs
                    │  cross-species · scType     │  progenitor T cells,
                    │                             │  FOXP3/AIRE locus access
                    └──────────┬──────────────────┘
                               │ still uncertain
                    ┌──────────▼──────────────────┐
                    │         TIER 4              │  ~2–3% of cells
                    │   LLM Elite Ensemble        │  ────────────────
                    │                             │  Full evidence bundle:
                    │  Claude · GPT · Kimi · Grok │  Tier 1–3 results +
                    │  + CellMarker 2.0 RAG       │  velocity + chromatin +
                    │                             │  spatial context + rules
                    └──────────┬──────────────────┘
                               │ still unresolved
                    ┌──────────▼──────────────────┐
                    │         TIER 5              │  < 1% of cells
                    │   Novel Attractor Discovery │  ────────────────
                    │                             │  Characterized, not
                    │  GRN attractors · velocity  │  discarded. Generates
                    │  sinks · Cell Ontology dist │  falsifiable predictions
                    │  → KG rule candidates       │  + feeds back to Tier 2
                    └─────────────────────────────┘
```

---

## Tier Strengths at a Glance

| Tier | Mechanism | Key Strength | Cells handled |
|------|-----------|-------------|---------------|
| **1 · AXIOMTier1** | MLP ensemble trained on 22M+ cells | Fast, broad coverage of common types | ~75–78% |
| **2 · KG Engine** | Proof-by-contradiction, 640 mechanistic rules | Eliminates false positives that markers alone cannot catch | ~15% |
| **3 · Convergence** | 6 orthogonal evidence streams | Resolves ambiguous subtypes: exhausted vs progenitor T, Treg vs activated T | ~5–7% |
| **4 · LLM Ensemble** | Reasoning over full Tier 1–3 evidence bundle | Open-vocabulary; handles rare/novel types not in KG | ~2–3% |
| **5 · Discovery** | GRN attractor characterization | Turns unknown states into characterized discoveries with falsifiable predictions | <1% |

### Why proof-by-contradiction changes everything

Traditional LLM annotators send marker gene lists to an LLM and pick the most-mentioned cell type. AXIOM-SC Tier 2 instead asks: *can this candidate be **ruled out**?*

```
Candidate: ILC3
  ✓ PASS    ILC3_CIRCUIT_001  RORC regulon active + NCR2/NCR3 present  →  circuit satisfied
  ✗ FAIL    ILC3_NEG_001      TRAC present  →  CONTRADICTED
            Mechanistic basis: ILCs lack VDJ recombination; TRAC = T cell receptor constant

Candidate: Myofibroblast (same dataset, RORC active via circadian regulation)
  — N/T    ILC3_CIRCUIT_001  NCR2 absent  →  circuit NOT SATISFIED
            → False positive eliminated. Myofibroblast correctly retained as candidate.
```

Tier 2 recovers signal that SCENIC misses at low cell counts: FOXP3, AIRE, TBX21 regulons recovered at 50k cells with NES threshold 2.0 (vs published default 3.0 which gives zero recovery on immune master TFs).

---

## Benchmark Results (Phase 1)

Evaluated on Human Thymus Cell Atlas, Lung Cell Atlas, Tabula Sapiens (immune).

| Metric | CASSIA | mLLMCelltype | **AXIOM-SC** |
|--------|--------|-------------|-------------|
| All 100 low-QS clusters | 3% | 0% | **4%** |
| KG in-scope (18 clusters) | 5.6% | 0% | **22.2%** |
| PROVEN precision | — | — | **67%** |

> AXIOM-SC's advantage is largest precisely where other methods fail: ambiguous clusters with low quality scores that LLMs cannot confidently resolve from markers alone.

---

## Installation

```bash
pip install axiom-sc
```

**Python 3.10+** required. For GPU-accelerated pySCENIC (Tier 2), see [environment setup](docs/setup.md).

---

## Quick Start

```python
import axiom_sc

# Load your clustered AnnData (normalized, log1p, Leiden clusters in obs)
import anndata
adata = anndata.read_h5ad("my_dataset.h5ad")

# Annotate with Tiers 1 + 2 (default — ~5 min for 50k cells)
annotator = axiom_sc.AXIOMAnnotator(profile="oss-apache", tiers=[1, 2])
result = annotator.annotate(adata)

# Results per cluster
print(result.summary())
# cluster  label       verdict      confidence  tier
# 0        CD4_Tcm     PROVEN       0.94        1
# 1        pDC         PROVEN       0.81        2
# 2        ILC3        PROVEN       0.76        2
# 3        Unknown_3   UNCERTAIN    0.41        2   ← routes to Tier 3
```

### Run all 5 tiers

```python
annotator = axiom_sc.AXIOMAnnotator(tiers=[1, 2, 3, 4, 5])
result = annotator.annotate(adata)
```

### Inspect which rules fired per cluster

```python
# See exactly why a cluster was PROVEN or CONTRADICTED
for firing in result.clusters["thy-22"].rule_firings:
    print(f"{firing.verdict:12s}  {firing.rule_id}  —  {firing.mechanistic_basis}")

# PASS          PDC_CIRCUIT_001  —  IRF7 regulon active (z=3.7): master pDC TF
# PASS          PDC_POS_001      —  SIGLEC1 high: pDC surface marker
# NOT_TESTABLE  PDC_NEG_001      —  PAX5 regulon: not in SCENIC output
# Verdict: PROVEN (confidence 0.81)
```

### Use a custom profile (academic, with PanglaoDB)

```python
# axiom_profile.json
{
  "name": "academic-full",
  "base_profile": "oss-apache",
  "additional_components": ["kg_panglao"],
  "license_acknowledgements": {
    "kg_panglao": "Used for academic non-commercial research only per CC BY-NC 4.0"
  }
}
```

```python
annotator = axiom_sc.AXIOMAnnotator(profile_path="axiom_profile.json")
```

### Query the knowledge graph

```python
# Browse rules for a cell type
treg_rules = axiom_sc.list_kg_rules(cell_type="Treg", rule_type="negative")
for r in treg_rules:
    print(f"{r['rule_id']}  —  {r['mechanistic_basis']}  [PMID:{r['pmid']}]")

# TREG_NEG_001  —  IL2 high contradicts Treg: FOXP3 suppresses IL2 transcription  [PMID:7584460]
# TREG_NEG_002  —  IFNG regulon active contradicts Treg identity  [PMID:15790681]

# Add a new pending rule (queued for expert review before activation)
axiom_sc.add_kg_rule({
    "cell_type": "My_Cell_Type",
    "rule_id": "MYC_NEG_001",
    "rule_type": "negative",
    "evidence_source": "marker_genes",
    "gene_or_regulon": ["CD3D"],
    "direction": "high",
    "mechanistic_basis": "CD3D high contradicts non-T cell: TCR complex subunit",
    "pmid": "1698053",
    "confidence": "high",
    "tissue_context": ["blood"],
    "source_db": "manual",
    "status": "ACTIVE",
    "added_in_version": "0.2.0"
})
```

---

## Profile System

Three built-in profiles control which components are enabled:

| Profile | Use case | What's included |
|---------|----------|-----------------|
| `oss-mit` | Strictest open source — MIT/Apache only | Tier 1 + KG engine |
| `oss-apache` | **Default** — all permissive licenses | All tiers; pySCENIC subprocess-isolated |
| `commercial` | Production deployments | All tiers; excludes CC BY-NC data sources |

```python
# Explicitly set profile
annotator = axiom_sc.AXIOMAnnotator(profile="commercial")

# Or via environment variable (overrides code)
# export AXIOM_PROFILE=oss-apache
```

---

## Knowledge Graph

The bundled KG (`kg_data/oracle_kg_v0.2.0.json`) contains **640 ACTIVE rules across 198 cell types**, derived from primary literature and seeded from CellMarker 2.0 (CC BY 4.0).

Every rule has:
- A **verified PubMed ID** — no rule is accepted without a primary source
- A **mechanistic basis** — human-readable explanation of the biological logic
- A **rule type**: `positive` · `negative` · `circuit` · `spatial`
- **Tissue context** and **confidence** level

The KG is released separately under **CC BY 4.0** — scientific knowledge should be maximally open regardless of code licensing.

Full citations: [REFERENCES.md](REFERENCES.md)

---

## Supported Input Formats

AXIOM-SC auto-detects input type and applies only the preprocessing steps needed:

| Input | Auto-detected as | Preprocessing |
|-------|-----------------|---------------|
| Clustered `.h5ad` (normalized, Leiden in `obs`) | `clustered_h5ad` | DE markers only |
| Unclustered `.h5ad` (normalized, no clusters) | `unclustered_h5ad` | Leiden clustering → DE |
| Raw counts `.h5ad` | `raw_count_h5ad` | Normalize → log1p → HVG → PCA → Leiden |
| CellRanger `.h5` | `cellranger_h5` | Full preprocessing pipeline |
| Marker gene `.csv` (scanpy / Seurat / CASSIA format) | `marker_csv` | Skip to Tier 2 directly |

---

## License

Apache 2.0 — see [LICENSE](LICENSE).

`pySCENIC` (GPL v3) is called via subprocess isolation only and is **never imported** into `axiom_sc`. `PanglaoDB` (CC BY-NC 4.0) is not bundled and only available via user-provided custom profiles for academic use.

---

## Citation

If you use AXIOM-SC in your research, please cite:

```
[Author et al.] AXIOM-SC: mechanistic cell type annotation via proof-by-contradiction.
bioRxiv (2026). doi: [doi]
```

Key tools used by AXIOM-SC (please also cite):

| Tool | Reference |
|------|-----------|
| CellMarker 2.0 (KG seed) | Hu et al. (2023) *Nucleic Acids Research* 51:D870 |
| pySCENIC | Van de Sande et al. (2020) *Nature Protocols* 15:2247 |
| scVelo | Bergen et al. (2020) *Nature Biotechnology* 38:1408 |
| Signac | Stuart et al. (2021) *Nature Methods* 18:1272 |
| COMMOT | Cang et al. (2023) *Nature Communications* 14:7706 |
| CELLxGENE Census | Tabula Sapiens Consortium (2022) *Science* 376:eabl4896 |

Full reference list: [REFERENCES.md](REFERENCES.md)
