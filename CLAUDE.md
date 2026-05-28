# AXIOM-SC — Phase 2
## CLAUDE.md — Definitive implementation spec for Claude Code

**Package:** axiom-sc (PyPI) | **Import:** axiom_sc | **License:** Apache 2.0
**GitHub:** receptor.bio/axiom-sc (SDK) · receptor.bio/axiom-playground (UI)
**Version:** 0.2.0

---

## 1. What AXIOM-SC is

AXIOM-SC is a 5-tier routing system for single-cell RNA-seq cell type annotation.
It is NOT a classifier. Each cell is routed to the minimum tier required to annotate
it confidently. Tiers use orthogonal evidence so failures cannot cascade.

**Core innovation (Phase 1 proven):** proof-by-contradiction using mechanistic
biological rules. One hard rule violation eliminates a candidate cell type regardless
of how many positive markers support it. This resolves failures all LLM-based systems
(CASSIA, mLLMCelltype) share because they only use positive marker matching.

**Five tiers:**

| Tier | Tool | Cells handled | Mechanism |
|---|---|---|---|
| 1 | AXIOMTier1 (our trained MLP) | ~75–78% | Statistical deep learning on 22M+ cells |
| 2 | AXIOM KG engine + SCENIC+ | ~15% | Mechanistic constraint proof-by-contradiction |
| 3 | scVelo + Signac + COMMOT + spatial | ~5–7% | Multi-stream orthogonal convergence |
| 4 | LLM elite ensemble (API) | ~2–3% | Open-vocabulary reasoning over full evidence |
| 5 | Novel attractor discovery | <1% | Unknown cell states characterized as discoveries |

---

## 2. Phase 1 results (carry as baseline — Phase 2 must beat these)

| Metric | CASSIA | mLLMCelltype | AXIOM Phase 1 |
|---|---|---|---|
| All 100 low-QS clusters | 3% | 0% | 4% |
| KG in-scope (18 clusters) | 5.6% | 0% | **22.2%** |
| PROVEN precision | — | — | **67% (2/3)** |

**Phase 1 datasets:** Human Thymus Cell Atlas · Lung Cell Atlas · Tabula Sapiens (immune)
**Phase 1 KG:** 41 rules · 18 cell types · 8 failure categories

---

## 3. License and profile system

### 3a. Apache 2.0 compliance rules (NON-NEGOTIABLE)

| Tool | License | Action |
|---|---|---|
| pySCENIC | GPL v3 | Subprocess isolation only. Never import into axiom_sc. |
| ArchR | Non-commercial | Replace with Signac (MIT) everywhere. |
| scType (R code) | GPL v3 | Reimplement scoring algorithm in Python (~200 lines). Do not import. |
| scTypeDB (data) | GPL v3 repo | Use CellMarker 2.0 (CC BY 4.0) for KG seeding instead. |
| PanglaoDB | CC BY-NC | Not bundled. Available only in academic profile as user-provided data. |
| scTab weights | No license filed | Not used. Build AXIOMTier1 with own trained weights (see Section 8). |

### 3b. Profile system — full specification

Profiles control which components are loaded and run. Three built-in profiles.
Profile is set at SDK instantiation and propagated to the playground UI.

```python
# axiom_sc/profiles/registry.py
# License: Apache 2.0

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional

class ProfileType(str, Enum):
    OSS_MIT     = "oss-mit"
    OSS_APACHE  = "oss-apache"   # DEFAULT — covers all permissive licenses
    COMMERCIAL  = "commercial"   # adds licensed restricted tools

@dataclass
class ComponentMeta:
    component_id: str
    display_name: str
    license: str               # SPDX identifier
    commercial_ok: bool        # safe for commercial use without extra licensing
    requires_permission: bool  # needs explicit author permission
    bundled: bool              # included in pip install (vs user-provided)
    tier: int                  # which AXIOM-SC tier uses this
    paper_doi: str             # primary paper DOI for citation
    notes: str = ""

COMPONENT_REGISTRY: dict[str, ComponentMeta] = {

    # ── Tier 1 ───────────────────────────────────────────────────────────────
    "tier1_axiomtier1": ComponentMeta(
        "tier1_axiomtier1", "AXIOMTier1 (own trained MLP)", "Apache-2.0",
        True, False, True, 1,
        "10.1038/s41467-024-46499-y",   # scTab paper (architecture reference)
        "Our own trained weights on CELLxGENE Census. Apache 2.0 weights."
    ),
    "tier1_census": ComponentMeta(
        "tier1_census", "CELLxGENE Census API", "MIT",
        True, False, True, 1,
        "10.1126/science.abl4896",
        "Used for nearest-neighbor annotation fallback alongside AXIOMTier1."
    ),

    # ── Tier 2 ───────────────────────────────────────────────────────────────
    "tier2_axiom_kg": ComponentMeta(
        "tier2_axiom_kg", "AXIOM KG Engine", "Apache-2.0",
        True, False, True, 2,
        "axiom-sc-phase1",   # our own paper
        "Proof-by-contradiction mechanistic KG. Phase 1 validated."
    ),
    "tier2_scenic": ComponentMeta(
        "tier2_scenic", "pySCENIC (subprocess)", "GPL-3.0",
        False, False, False, 2,
        "10.1038/s41596-020-0336-2",
        "GPL isolated via subprocess. Install in separate scenic-env conda env."
    ),
    "tier2_sctype_scoring": ComponentMeta(
        "tier2_sctype_scoring", "scType scoring (reimplemented)", "Apache-2.0",
        True, False, True, 2,
        "10.1038/s41467-022-28803-w",
        "Algorithm reimplemented in Python. Original code is GPL v3 — not imported."
    ),
    "tier2_gnn_propagation": ComponentMeta(
        "tier2_gnn_propagation", "kNN verdict propagation (ScInfeR-style)", "Apache-2.0",
        True, False, True, 2,
        "10.1093/bioinformatics/btad680",
        "Algorithm reimplemented. Original ScInfeR code not imported."
    ),

    # ── Tier 3 ───────────────────────────────────────────────────────────────
    "tier3_scvelo": ComponentMeta(
        "tier3_scvelo", "scVelo", "BSD-3-Clause",
        True, False, False, 3,
        "10.1038/s41587-020-0591-3", ""
    ),
    "tier3_signac": ComponentMeta(
        "tier3_signac", "Signac (ATAC)", "MIT",
        True, False, False, 3,
        "10.1038/s41592-021-01282-5",
        "R package accessed via rpy2. Replaces ArchR (non-commercial)."
    ),
    "tier3_commot": ComponentMeta(
        "tier3_commot", "COMMOT (cell communication)", "MIT",
        True, False, False, 3,
        "10.1038/s41467-023-43600-9", ""
    ),
    "tier3_ortho": ComponentMeta(
        "tier3_ortho", "OrthoFinder (cross-species)", "MIT",
        True, False, False, 3,
        "10.1186/s13059-019-1832-y", ""
    ),
    "tier3_cellama_niche": ComponentMeta(
        "tier3_cellama_niche", "CELLama niche encoder", "MIT",
        True, False, False, 3,
        "10.1002/advs.202413514",
        "Gene-rank sentence encoding for spatial niche context. Cheap, no API."
    ),

    # ── Tier 4 ───────────────────────────────────────────────────────────────
    "tier4_anthropic": ComponentMeta(
        "tier4_anthropic", "Anthropic API (Claude)", "Commercial API",
        True, False, False, 4,
        "N/A", "User provides ANTHROPIC_API_KEY."
    ),
    "tier4_openai": ComponentMeta(
        "tier4_openai", "OpenAI API (GPT)", "Commercial API",
        True, False, False, 4,
        "N/A", "User provides OPENAI_API_KEY."
    ),
    "tier4_cellmarker_rag": ComponentMeta(
        "tier4_cellmarker_rag", "CellMarker 2.0 RAG", "CC-BY-4.0",
        True, False, True, 4,
        "10.1093/nar/gkac947",
        "Attribution required. 2578 cell types, 26915 markers."
    ),

    # ── Tier 5 ───────────────────────────────────────────────────────────────
    "tier5_discovery": ComponentMeta(
        "tier5_discovery", "Novel attractor discovery", "Apache-2.0",
        True, False, True, 5,
        "axiom-sc-phase1", "Own code. GRN attractor characterization."
    ),

    # ── KG data sources ──────────────────────────────────────────────────────
    "kg_cellmarker2": ComponentMeta(
        "kg_cellmarker2", "CellMarker 2.0 (KG seed)", "CC-BY-4.0",
        True, False, False, 2,
        "10.1093/nar/gkac947",
        "Primary KG seeding source. Download and process offline."
    ),
    "kg_panglao": ComponentMeta(
        "kg_panglao", "PanglaoDB (optional, academic)", "CC-BY-NC-4.0",
        False, False, False, 2,
        "10.1093/database/baz046",
        "NOT bundled. Academic use only. User provides path if desired."
    ),
    "kg_dorothea": ComponentMeta(
        "kg_dorothea", "DoRothEA / OmniPath TF regulons", "MIT",
        True, False, False, 2,
        "10.1186/s13073-019-0668-6", ""
    ),
}

# Profile → allowed component set
OSS_MIT_COMPONENTS: set[str] = {
    k for k, v in COMPONENT_REGISTRY.items()
    if v.license in ("MIT", "Apache-2.0") and v.commercial_ok
}

OSS_APACHE_COMPONENTS: set[str] = {
    k for k, v in COMPONENT_REGISTRY.items()
    if v.commercial_ok and not v.requires_permission
    and v.license not in ("GPL-3.0",)   # GPL excluded even from OSS-Apache
}

COMMERCIAL_COMPONENTS: set[str] = {
    k for k, v in COMPONENT_REGISTRY.items()
    if k not in ("kg_panglao",)   # PanglaoDB CC BY-NC excluded from commercial
}

PROFILE_COMPONENT_MAP: dict[ProfileType, set[str]] = {
    ProfileType.OSS_MIT:    OSS_MIT_COMPONENTS,
    ProfileType.OSS_APACHE: OSS_APACHE_COMPONENTS,
    ProfileType.COMMERCIAL: COMMERCIAL_COMPONENTS,
}
```

**Profile loading at SDK instantiation:**

```python
# axiom_sc/__init__.py
from axiom_sc.profiles.registry import ProfileType, PROFILE_COMPONENT_MAP
from axiom_sc.profiles.loader import load_profile_from_env

def AXIOMAnnotator(
    profile: str | ProfileType = "oss-apache",
    profile_path: str | None = None,  # path to custom axiom_profile.json
    tiers: list[int] | None = None,
    **kwargs
):
    """
    Main entry point.

    profile: "oss-mit" | "oss-apache" (default) | "commercial"
             OR path to a custom axiom_profile.json
    tiers:   subset of [1,2,3,4,5] to run. None = all tiers in profile.

    Environment override: AXIOM_PROFILE=oss-apache
    """
```

**Custom profile JSON schema** (`axiom_profile.json`):

```json
{
  "name": "my-custom-profile",
  "description": "Academic project, all tools available",
  "base_profile": "oss-apache",
  "additional_components": ["kg_panglao"],
  "disabled_components": ["tier4_openai"],
  "tiers_enabled": [1, 2, 3, 4, 5],
  "license_acknowledgements": {
    "kg_panglao": "Used for academic non-commercial research only per CC BY-NC 4.0",
    "tier2_scenic": "pySCENIC GPL v3 run via subprocess isolation"
  },
  "created_at": "2026-05-23",
  "created_by": "researcher@institution.edu"
}
```

---

## 4. Repository structure

### axiom-sc (SDK + PyPI)

```
axiom-sc/
├── LICENSE                              # Apache 2.0
├── README.md                            # see Section 14 for citation spec
├── REFERENCES.md                        # auto-generated from KG PMIDs
├── CHANGELOG.md                         # semver: 0.1.0 (Phase 1) → 0.2.0 (Phase 2)
├── pyproject.toml
├── axiom_sc/
│   ├── __init__.py                      # AXIOMAnnotator, __version__, kg_version
│   ├── version.py                       # __version__ = "0.2.0"
│   ├── profiles/
│   │   ├── registry.py                  # ComponentMeta + profile maps (above)
│   │   ├── loader.py                    # load profile from env/file/string
│   │   └── validator.py                 # check component compatibility
│   ├── tier1/
│   │   ├── axiomtier1.py               # our trained MLP ensemble
│   │   ├── census_annotator.py          # CELLxGENE Census kNN fallback
│   │   ├── ensemble.py                  # vote + confidence fusion
│   │   └── training/
│   │       ├── train.py                 # training pipeline (run once)
│   │       ├── data_loader.py           # Census download + stratified split
│   │       └── evaluate.py             # calibration + F1 benchmarks
│   ├── tier2/
│   │   ├── axiom_annotator.py          # AxiomAnnotator (Phase 1 extended)
│   │   ├── evidence.py                  # EvidenceBundle dataclass
│   │   ├── kg_loader.py                 # loads + validates KG JSON
│   │   ├── scenic_runner.py             # GPL-isolated pySCENIC subprocess
│   │   ├── scenic_worker.py             # runs IN scenic-env (GPL side)
│   │   ├── sctype_scorer.py             # scType algorithm reimplemented
│   │   └── gnn_propagator.py            # kNN verdict propagation
│   ├── tier3/
│   │   ├── velocity.py                  # scVelo RNA velocity stream
│   │   ├── chromatin.py                 # Signac ATAC locus accessibility
│   │   ├── communication.py             # COMMOT L-R stream
│   │   ├── spatial_niche.py             # CELLama sentence encoder
│   │   ├── cross_species.py             # OrthoFinder conservation oracle
│   │   └── convergence.py              # multi-stream vote + confidence fusion
│   ├── tier4/
│   │   ├── llm_ensemble.py             # LLM routing + consensus
│   │   ├── evidence_bundler.py          # formats Tier 1-3 evidence for LLM
│   │   ├── rag.py                       # CellMarker 2.0 RAG retriever
│   │   └── prompts.py                   # reproducible prompt templates
│   ├── tier5/
│   │   ├── attractor_discovery.py       # novel state characterization
│   │   ├── candidate_state.py           # CandidateState dataclass
│   │   └── feedback.py                  # Tier 5 → Tier 2 KG loop
│   ├── kg/
│   │   ├── schema.json                  # KG rule JSON schema
│   │   ├── loader.py                    # validates + indexes rules
│   │   ├── seeder.py                    # CellMarker 2.0 → candidate rules
│   │   ├── review_cli.py                # interactive rule review tool
│   │   └── references.py               # PMID lookup + REFERENCES.md generator
│   ├── pipelines/
│   │   ├── preprocess.py               # scanpy pipeline (Phase 1 validated)
│   │   ├── scenic_pipeline.py           # full SCENIC+ orchestrator
│   │   └── evaluate.py                 # comparison metrics + plots
│   └── utils/
│       ├── anndata_utils.py
│       ├── stratified_sample.py         # ≥300 cells/type
│       └── zscore_auc.py               # Phase 1 validated normalization
├── kg_data/
│   └── oracle_kg_v0.2.0.json           # KG snapshot CC BY 4.0 (see note)
├── model_weights/
│   └── README.md                        # points to HuggingFace: receptor-bio/axiomtier1
├── tests/
│   ├── fixtures/                        # small synthetic AnnData (50 cells)
│   ├── test_profiles.py
│   ├── test_tier1.py
│   ├── test_tier2.py
│   ├── test_tier3.py
│   ├── test_tier4.py
│   ├── test_tier5.py
│   └── test_kg.py
└── .github/
    └── workflows/
        ├── ci.yml                       # pytest on every PR
        └── publish.yml                  # tag v* → PyPI
```

**Note on kg_data license:** The KG JSON file contains rules derived from primary
literature (our writing) seeded from CellMarker 2.0 (CC BY 4.0). The compiled KG
is released under CC BY 4.0 separately from the Apache 2.0 code. This is intentional
— the scientific knowledge should be maximally open regardless of code licensing.

### axiom-playground (FastAPI + React)

```
axiom-playground/
├── LICENSE                              # Apache 2.0
├── README.md
├── docker-compose.yml                   # LOCAL DEV ONLY — production uses ECS Fargate (Section 27)
├── config.yaml                          # sdk_version: "0.2.0"
├── .env.example
├── .aws/                                # ECS task definition templates (Section 27)
├── .github/workflows/deploy.yml         # OIDC CI/CD → ECR → ECS (Section 27)
├── api/
│   ├── main.py
│   ├── config.py                        # SDK version validation
│   ├── Dockerfile
│   ├── Dockerfile.worker                # Dramatiq worker image
│   ├── requirements.txt                 # axiom-sc=={AXIOM_SDK_VERSION}
│   ├── tasks.py                         # Dramatiq jobs (NOT Celery — see Section 22/26)
│   └── routes/
│       ├── annotate.py
│       ├── jobs.py
│       ├── kg.py
│       ├── profiles.py                  # GET/POST /profiles
│       └── health.py
└── app/
    ├── package.json
    ├── vite.config.ts
    ├── Dockerfile
    └── src/
        ├── App.tsx
        ├── api/                         # typed API client
        ├── store/                       # Zustand
        └── components/
            ├── layout/
            ├── upload/
            ├── pipeline/
            ├── results/
            ├── kg/
            ├── profiles/                # ProfileSelector + ProfileCreator
            └── shared/
```

---

## 5. Environment setup

### 5a. Main conda environment (Apache 2.0 clean — axiom-env)

```bash
conda create -n axiom-env python=3.10 -y
conda activate axiom-env

# Core science
pip install scanpy anndata leidenalg python-igraph
pip install scvelo commot
pip install cellxgene-census
pip install pyomnipath
pip install torch torchvision --index-url https://download.pytorch.org/whl/cu118

# API + async
pip install fastapi "uvicorn[standard]" redis "dramatiq[redis]>=1.16" python-multipart aiofiles
pip install httpx pytest-asyncio boto3   # boto3 for S3 (RunPod/Modal loom exchange) + Secrets Manager

# KG tooling
pip install jsonschema biopython requests openpyxl

# Dev
pip install pytest pytest-cov black ruff mypy pre-commit hatch

# R dependencies for Signac (MIT)
conda install -c conda-forge r-base r-signac r-seurat bioconductor-bsgenome.hsapiens.ucsc.hg38
pip install rpy2

# Install axiom-sc in editable mode
cd axiom-sc && pip install -e ".[dev]"

conda env export > environment.yml
```

### 5b. SCENIC subprocess environment (GPL isolated — scenic-env)

```bash
conda create -n scenic-env python=3.10 -y
conda activate scenic-env
pip install pyscenic loompy
# This env is ONLY called via subprocess from axiom_sc/tier2/scenic_runner.py
# Never activate scenic-env manually during development
# Never import pyscenic inside axiom-env
```

### 5c. Node environment (React playground)

```bash
nvm install 20 && nvm use 20
cd axiom-playground/app
npm install
# Key packages: vite, react, typescript, @tanstack/react-query,
#               zustand, plotly.js-dist-min, react-dropzone, tailwindcss,
#               @visx/visx, react-hot-toast, lucide-react
```

### 5d. .env file (copy from .env.example — NEVER COMMIT)

```
AXIOM_SDK_VERSION=0.2.0
AXIOM_PROFILE=oss-apache
AXIOM_KG_PATH=./kg_data/oracle_kg_v0.2.0.json
SCENIC_CONDA_ENV=scenic-env
REDIS_URL=redis://localhost:6379
ANTHROPIC_API_KEY=           # required for Tier 4
OPENAI_API_KEY=              # required for Tier 4
HUGGINGFACE_TOKEN=           # required for AXIOMTier1 weight download
```

---

## 6. pyproject.toml

```toml
[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[project]
name = "axiom-sc"
version = "0.2.0"
description = "AXIOM-SC: mechanistic cell type annotation via knowledge graph constraints"
readme = "README.md"
license = {text = "Apache-2.0"}
requires-python = ">=3.10"
authors = [{name = "receptor.bio", email = "oss@receptor.bio"}]
keywords = ["single-cell", "rna-seq", "cell-type", "annotation", "bioinformatics",
            "knowledge-graph", "mechanistic", "attractor"]
classifiers = [
    "License :: OSI Approved :: Apache Software License",
    "Programming Language :: Python :: 3.10",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]
dependencies = [
    "scanpy>=1.9", "anndata>=0.9",
    "scvelo>=0.3", "commot>=0.0.3",
    "cellxgene-census>=1.9",
    "pyomnipath>=0.1",
    "jsonschema>=4.17",
    "pandas>=2.0", "numpy>=1.24",
    "scikit-learn>=1.3", "scipy>=1.11",
    "torch>=2.0",
    "biopython>=1.81",
    "requests>=2.31",
    "openpyxl>=3.1",
]
[project.optional-dependencies]
dev = ["pytest", "pytest-cov", "black", "ruff", "mypy", "pre-commit", "hatch"]
gpu = ["cupy-cuda12x"]

[project.scripts]
axiom-annotate    = "axiom_sc.cli:annotate"
axiom-kg          = "axiom_sc.cli:kg"
axiom-benchmark   = "axiom_sc.cli:benchmark"
axiom-train-tier1 = "axiom_sc.tier1.training.train:main"
```

---

## 7. SDK versioning and release pipeline

### .github/workflows/ci.yml

```yaml
name: CI
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with: {python-version: "3.10"}
      - run: pip install -e ".[dev]"
      - run: pytest tests/ -v --cov=axiom_sc --cov-report=xml --cov-fail-under=80
      - uses: codecov/codecov-action@v4
```

### .github/workflows/publish.yml

```yaml
name: Publish to PyPI
on:
  push:
    tags: ["v*"]
jobs:
  publish:
    runs-on: ubuntu-latest
    environment: pypi
    permissions:
      id-token: write   # OIDC trusted publishing — no stored token needed
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with: {python-version: "3.10"}
      - run: pip install hatch && hatch build
      - uses: pypa/gh-action-pypi-publish@release/v1
```

### SDK version management in playground

```python
# axiom-playground/api/config.py
import yaml, importlib.metadata, os

def load_config() -> dict:
    with open("config.yaml") as f:
        cfg = yaml.safe_load(f)
    configured = os.getenv("AXIOM_SDK_VERSION") or cfg.get("sdk_version", "0.2.0")
    installed  = importlib.metadata.version("axiom-sc")
    if installed != configured:
        raise RuntimeError(
            f"SDK mismatch: config wants axiom-sc=={configured}, "
            f"but {installed} is installed. Run:\n"
            f"  pip install axiom-sc=={configured}"
        )
    return cfg
```

### Release workflow (researcher runs this)

```bash
# 1. Bump version
sed -i 's/__version__ = "0.2.0"/__version__ = "0.2.1"/' axiom_sc/version.py
sed -i 's/version = "0.2.0"/version = "0.2.1"/' pyproject.toml

# 2. Update CHANGELOG.md

# 3. Tag and push → GitHub Action publishes to PyPI automatically
git add -A && git commit -m "release: v0.2.1"
git tag v0.2.1 && git push origin main --tags

# 4. Update playground to use new version
echo "sdk_version: 0.2.1" > axiom-playground/config.yaml
git add axiom-playground/config.yaml && git commit -m "chore: bump sdk to 0.2.1"
git push origin main
# → GitHub Actions CI/CD builds new images → ECS rolling update → /health returns sdk_version: 0.2.1
# Local dev only: docker run --rm -p 6379:6379 redis:7-alpine & uvicorn main:app --reload
```

---

## 8. AXIOMTier1 — our own trained cell type model

**Why we train our own:** scTab has no license filed (all-rights-reserved by default).
We replicate the core approach — MLP ensemble on CELLxGENE data — using our own
training pipeline. Released as Apache 2.0 weights on HuggingFace under receptor-bio org.

**Architecture (MLP ensemble — paper-described, no patent):**

```python
# axiom_sc/tier1/axiomtier1.py
"""
AXIOMTier1: MLP ensemble for cell type classification.

Architecture based on: Fischer et al. (2024) Nature Communications 15:6611
doi: 10.1038/s41467-024-46499-y

Training data: CELLxGENE Census (CC BY 4.0)
Our trained weights: Apache 2.0, released at HuggingFace receptor-bio/axiomtier1
"""
import torch
import torch.nn as nn

class AXIOMTier1Net(nn.Module):
    """
    Single MLP in the 10-model ensemble.
    Input: log1p normalized expression (10k), ~2000 HVGs
    Output: softmax over n_cell_types
    """
    def __init__(self, n_genes: int, n_classes: int, hidden: int = 512):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(n_genes, hidden), nn.BatchNorm1d(hidden),
            nn.ReLU(), nn.Dropout(0.1),
            nn.Linear(hidden, hidden), nn.BatchNorm1d(hidden),
            nn.ReLU(), nn.Dropout(0.1),
            nn.Linear(hidden, 256), nn.BatchNorm1d(256),
            nn.ReLU(), nn.Dropout(0.1),
            nn.Linear(256, n_classes)
        )
    def forward(self, x): return self.net(x)


class AXIOMTier1Ensemble:
    """
    Ensemble of N=10 independently trained AXIOMTier1Net models.
    Confidence = 1 - entropy(mean_softmax_across_10_models)
    """
    N_MODELS = 10
    CONFIDENCE_THRESHOLD_ACCEPT  = 0.85  # route: Tier 1 accept
    CONFIDENCE_THRESHOLD_VERIFY  = 0.50  # below: Tier 2 full search

    def __init__(self, weights_path: str | None = None):
        """
        weights_path: local path or HuggingFace repo ID
        If None: downloads receptor-bio/axiomtier1 from HuggingFace
        """

    def predict(self, adata) -> dict:
        """
        Returns {
            'label': pd.Series per cell,
            'confidence': pd.Series per cell,
            'route': pd.Series — 'accept' | 'tier2_verify' | 'tier2_full'
        }
        Requires: adata normalized to 10k counts, log1p transformed.
        """
```

**Training pipeline:**

```python
# axiom_sc/tier1/training/train.py
"""
One-time training of AXIOMTier1 on CELLxGENE Census.
Run: axiom-train-tier1 --n-models 10 --output weights/

TRAINING DATA:
  Source: CELLxGENE Census (cellxgene_census Python API, MIT license)
  Data license: CC BY 4.0 per dataset
  n_cells: ~22M+ (full Census) or configurable subset
  n_cell_types: 164+ from Census harmonized labels
  Preprocessing: 10k normalization, log1p, top 2000 HVG, stratified split

COMPUTE ESTIMATE:
  Single A100 (80GB): ~5h per model, 50h for 10 models
  Cost: ~$75 total at $1.50/hr on Lambda Labs
  Run once. Weights stored in HuggingFace receptor-bio/axiomtier1.

LICENSE:
  Our trained weights: Apache 2.0
  Released at: huggingface.co/receptor-bio/axiomtier1
"""
import cellxgene_census
import torch
from torch.utils.data import DataLoader

def main():
    # 1. Download Census training data
    # 2. Stratified split (≥2000 cells per cell type for rare types)
    # 3. Train 10 independent models with different random seeds
    # 4. Evaluate ensemble: target macro-F1 ≥ 0.80 on held-out test set
    # 5. Push weights to HuggingFace receptor-bio/axiomtier1
    pass
```

---

## 9. Tier 2 — AXIOM KG engine

**Carries forward all Phase 1 validated code and corrections.**

### Phase 1 lessons — HARDCODED as non-negotiable defaults

```python
# axiom_sc/tier2/scenic_pipeline.py
"""
pySCENIC pipeline orchestrator.
Calls pySCENIC via subprocess isolation (GPL v3 — not imported directly).

References:
  Van de Sande et al. (2020) Nature Protocols 15:2247
  doi: 10.1038/s41596-020-0336-2
"""

SCENIC_DEFAULTS = {
    # Phase 1 lesson: NES 3.0 (default) recovers ZERO immune master TFs at 20k cells
    # NES 2.0 recovers IRF7, GATA3, RORC
    "nes_threshold":    2.0,    # NOT the published default of 3.0

    # Phase 1 lesson: default min_genes=30 too strict for rare TFs
    "min_genes":        5,

    # Phase 1 lesson: HVG selection excludes master TFs (low variance in mixed datasets)
    # These MUST be added to loom input regardless of HVG rank
    "forced_genes": [
        "FOXP3","TBX21","GATA3","RORC","IRF7","PAX5","EBF1",
        "AIRE","DLL4","PSMB11","TOX","TCF7","NR4A1","RUNX3",
        "IKZF2","IKZF3","BCL6","PRDM1","XBP1","CEBPA","SPI1",
        "IRF4","IRF8","BATF3","ID2","EOMES",
    ],

    # Phase 1 lesson: raw AUC values not comparable across runs/datasets
    # (thymus: 0-12, lung: 0.01-0.35)
    "z_score_auc":      True,   # ALWAYS z-score per regulon per dataset
    "z_active_thresh":  1.5,    # z ≥ 1.5 = regulon active
    "z_inactive_thresh":-0.5,   # z ≤ -0.5 = regulon inactive

    # Phase 1 lesson: 20k cells insufficient to recover FOXP3/TBX21/AIRE
    # (expressed in <2% of cells → GRN signal too weak)
    "min_cells_per_type": 300,  # stratified sampling requirement
    "recommended_total": 50000, # phase 2 uses 50k (not 20k)

    # SCENIC database URLs
    "tf_list_url":    "https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt",
    "rankings_db_url":"https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather",
    "motif_db_url":   "https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl",
}
```

### Phase 1 annotator corrections — ALL must be preserved

```python
# axiom_sc/tier2/axiom_annotator.py
"""
AxiomAnnotator: mechanistic cell type annotation via proof-by-contradiction.

Phase 1 validated corrections (DO NOT revert):

1. FULL KG TEST MODE: always test all KG cell types, not just CASSIA suggestions.
   CASSIA free-text labels fail to map to KG keys. Exhaustive testing
   (n_clusters × n_kg_types × n_rules ≈ 73,800 evaluations) completes in <3 min.

2. CIRCUIT COMPLETENESS: if any component gene is not in evidence bundle,
   circuit rule returns NOT_TESTABLE (not partial PASS).
   if rule_type == "circuit" and len(testable_genes) < len(genes):
       return NOT_TESTABLE

3. SUPPORT COUNTING REFORM: only positive/circuit PASS counts as support.
   Negative rule PASS = "not contradicted", not affirmative evidence.
   supports = [f for f in firings
               if f.verdict == "PASS"
               and f.rule_type in ("positive", "circuit")]

4. PROVEN CO-CONFIRMATION: circuit PASS alone is not sufficient for PROVEN.
   Requires at least one positive rule PASS (if positive rules exist for type).
   if circuit_passes and (positive_passes or not has_positive_rules):
       verdict = "PROVEN"
   This eliminates false PROVEN calls from non-immune RORC expression.

5. PDC_NEG_001 CORRECTION: check PAX5 REGULON z-score < -0.5 (not IGKC marker absence).
   pDCs naturally express IGKC (developmental B-lineage relationship).
   Checking IGKC absence fires false contradictions on true pDCs.
   Reference: Rodrigues et al. 2018 doi:10.1016/j.celrep.2018.11.094

6. ILC3_CIRCUIT_001 PHASE 2 FIX: add NCR2/NCR3 co-requirement.
   Phase 1 false positive (myofibroblast→ILC3) caused by non-immune RORC expression
   (circadian regulation). NCR2/NCR3 innate receptors absent in fibroblasts.
   Updated rule: RORC active + NCR2/NCR3 present + NO TRAC + NO GATA3.
"""
```

---

## 10. Tier 3 — multi-stream convergence

```python
# axiom_sc/tier3/convergence.py
"""
MultiStreamConvergence: 6 orthogonal evidence streams.

References:
  scVelo: Bergen et al. (2020) Nature Biotechnology 38:1408
          doi: 10.1038/s41587-020-0591-3
  Signac: Stuart et al. (2021) Nature Methods 18:1272
          doi: 10.1038/s41592-021-01282-5
  COMMOT: Cang et al. (2023) Nature Communications 14:7706
          doi: 10.1038/s41467-023-43600-9
  OrthoFinder: Emms & Kelly (2019) Genome Biology 20:238
               doi: 10.1186/s13059-019-1832-y
  CELLama: Lv et al. (2025) Advanced Science
           doi: 10.1002/advs.202413514

CONVERGENCE CRITERION:
  4+ of 6 streams agree → PROVEN (Tier 3)
  3 streams agree → HIGH_CONFIDENCE
  <3 agree → UNCERTAIN (Tier 4 handles in Phase 2+)

MINIMUM VIABLE (non-spatial data): velocity + chromatin (2 streams)
FULL (Visium/Multiome): all 6 streams
"""

STREAMS = {
    "velocity":     "scVelo trajectory direction",
    "chromatin":    "Signac locus accessibility",
    "communication":"COMMOT L-R signaling",
    "spatial_niche":"CELLama neighbor composition",
    "cross_species":"OrthoFinder conservation",
    "spatial_type": "scType spatial (reimplemented)",
}

# Key Tier 3 discriminations:
# - Exhausted vs progenitor-exhausted T: OPPOSITE velocity vectors
# - FOXP3 locus accessibility → Treg identity (SCENIC miss in Phase 1 fixed here)
# - AIRE locus accessibility → mature mTEC (SCENIC miss in Phase 1 fixed here)
# - Monocyte vs tissue-resident Mφ: spatial position relative to vessels
# - CAF vs SMC: COMMOT TGFβ signaling context
```

---

## 11. Tier 4 — LLM elite ensemble

```python
# axiom_sc/tier4/llm_ensemble.py
"""
LLM ensemble for hard cases: rare/novel types not in KG.

Model selection based on DeepCellSeek benchmark (2025):
  Elite ensemble: Kimi-k2, GPT-5, Claude-4, Grok-4
  Key finding: ensemble outperforms any single model on subtype annotation.

CRITICAL DIFFERENCE FROM CASSIA:
  CASSIA: sends marker genes to LLM (single evidence stream)
  AXIOM-SC Tier 4: sends FULL evidence bundle (Tier 1-3 results +
  velocity direction + chromatin state + spatial context + AXIOM rules fired)
  LLM reasons over all evidence simultaneously, not just marker genes.

API costs: user provides own API keys.
No data is sent to LLM APIs except the evidence bundle for the specific cluster.
"""
```

---

## 12. Tier 5 — novel attractor discovery

```python
# axiom_sc/tier5/attractor_discovery.py
"""
Novel attractor characterization.
Cells failing all 4 tiers are characterized, not discarded.

OUTPUT per novel candidate:
  candidate_id:    str (hash of dataset + cluster_id)
  active_circuits: list[str]  (which TF programs running from SCENIC)
  velocity_sink:   str | None (nearest known attractor in velocity space)
  spatial_context: str        (tissue niche description)
  cross_species:   dict       (OrthoFinder conservation scores)
  nearest_ontology_term: str  (closest Cell Ontology term by distance)
  perturbation_predictions: list[dict]  (falsifiable experimental predictions)
  status: "CANDIDATE" | "VALIDATED" | "PUBLISHED"

FEEDBACK LOOP (Tier 5 → Tier 2):
  Validated discoveries auto-generate KG rule candidates
  → expert review queue (axiom_sc/kg/review_cli.py)
  → merged to axiom_kg.json with new pmid
"""
```

---

## 13. KG expansion — 350 rules, 80+ cell types

### 13a. CellMarker 2.0 seeding pipeline

```python
# axiom_sc/kg/seeder.py
"""
KG seeder from CellMarker 2.0.

License: CellMarker 2.0 paper is CC BY 4.0 (unrestricted including commercial).
Source: Hu et al. (2023) Nucleic Acids Research 51:D870
        doi: 10.1093/nar/gkac947
WARNING: CellMarker 1.0 was CC BY-NC 4.0 — ensure using v2.0 data specifically.

DO NOT USE for seeding:
  - PanglaoDB: CC BY-NC 4.0 (non-commercial only) — excluded from all bundled profiles
  - scTypeDB: in GPL v3 repository (use CellMarker 2.0 instead)
"""

CELLMARKER2_URL = (
    "http://bio-bigdata.hrbmu.edu.cn/CellMarker/"
    "CellMarker_download_files/file/Cell_marker_Human.xlsx"
)

# After seeding: all rules have status="PENDING_REVIEW"
# Researcher must set status="ACTIVE" after verifying mechanistic_basis and pmid
# Rules with status != "ACTIVE" are not used by AxiomAnnotator
```

### 13b. KG JSON schema (updated from Phase 1)

```json
{
  "$schema": "http://json-schema.org/draft-07/schema",
  "required": ["cell_type","rule_id","rule_type","evidence_source",
               "gene_or_regulon","direction","mechanistic_basis",
               "pmid","confidence","status","tissue_context","source_db"],
  "properties": {
    "cell_type":          {"type": "string"},
    "rule_id":            {"type": "string", "pattern": "^[A-Z0-9_]+$"},
    "rule_type":          {"enum": ["positive","negative","circuit","spatial"]},
    "evidence_source":    {"enum": ["marker_genes","regulon","trajectory",
                                    "spatial","communication","chromatin"]},
    "gene_or_regulon":    {"type": "array", "items": {"type": "string"}, "minItems": 1},
    "direction":          {"enum": ["high","low","absent","present","active","inactive"]},
    "paired_with":        {"type": "array", "items": {"type": "string"}},
    "incompatible_with":  {"type": "array", "items": {"type": "string"}},
    "mechanistic_basis":  {"type": "string", "minLength": 20},
    "pmid":               {"type": "string", "minLength": 1},
    "confidence":         {"enum": ["high","medium","low"]},
    "tissue_context":     {"type": "array", "items": {"type": "string"}},
    "source_db":          {"type": "string"},
    "status":             {"enum": ["ACTIVE","PENDING_REVIEW","DEPRECATED"]},
    "added_in_version":   {"type": "string"}
  }
}
```

### 13c. Priority order (Phase 2 target: 350 ACTIVE rules)

| Priority | Cell types | Rules | Source |
|---|---|---|---|
| 1 (Week 1) | All immune (T/B/NK/ILC/myeloid) | ~170 | CellMarker 2.0 + literature |
| 2 (Week 2) | Lung/liver/gut epithelial | ~80 | CellMarker 2.0 + literature |
| 3 (Week 3) | Stromal/endothelial | ~50 | CellMarker 2.0 + literature |
| 4 (Week 4) | CNS/endocrine/rare | ~50 | CellMarker 2.0 + literature |

---

## 14. Citation system

### 14a. KG rule citations (enforced by schema)

Every rule in `axiom_kg.json` MUST have a non-empty `pmid` field.
If no single paper covers the rule, cite the most relevant mechanistic paper.
Claude Code: use NCBI Entrez API to look up PMIDs when writing rules.

```python
# axiom_sc/kg/references.py
"""
Generates REFERENCES.md from KG rule PMIDs.
Run: axiom-kg references --output REFERENCES.md
"""
import Bio.Entrez as Entrez

def fetch_citation(pmid: str) -> str:
    """Fetches formatted citation from PubMed."""
    Entrez.email = "oss@receptor.bio"
    handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="text")
    # parse → APA format
```

### 14b. Module-level docstrings (required in every tier module)

Every file in `axiom_sc/tier*/` MUST have a module docstring citing:
1. The primary paper for the tool/algorithm used
2. The doi
3. What the module does

Example format (already shown in Tier 2 scenic_pipeline.py above).

### 14c. README citation section

```markdown
## Citation

If you use AXIOM-SC in your research, please cite:

**AXIOM-SC:**
> [Author et al.] AXIOM-SC: mechanistic cell type annotation via proof-by-contradiction.
> bioRxiv (2026). doi: [doi]

**Key tools used by AXIOM-SC** (please also cite):
> AXIOMTier1 training data: Tabula Sapiens Consortium (2022) Science 376:eabl4896
> CellMarker 2.0 (KG source): Hu et al. (2023) Nucleic Acids Research 51:D870
> pySCENIC: Van de Sande et al. (2020) Nature Protocols 15:2247
> scVelo: Bergen et al. (2020) Nature Biotechnology 38:1408
> Signac: Stuart et al. (2021) Nature Methods 18:1272
> COMMOT: Cang et al. (2023) Nature Communications 14:7706

Full reference list: [REFERENCES.md](REFERENCES.md) (auto-generated from KG PMIDs)
```

---

## 15. FastAPI backend

```python
# axiom-playground/api/main.py
from fastapi import FastAPI
from contextlib import asynccontextmanager
import axiom_sc
from config import load_config

@asynccontextmanager
async def lifespan(app: FastAPI):
    cfg = load_config()              # hard-fails if SDK version mismatch
    app.state.sdk_version    = axiom_sc.__version__
    app.state.kg_version     = axiom_sc.kg_version
    app.state.kg_rule_count  = axiom_sc.kg_rule_count()
    app.state.active_profile = cfg.get("default_profile", "oss-apache")
    yield

app = FastAPI(title="AXIOM-SC Playground", version="0.2.0", lifespan=lifespan)
```

### Routes spec

```
GET  /health
     → {sdk_version, kg_version, kg_rule_count, active_profile, tiers_available}

POST /annotate
     body: {file: .h5ad upload, tiers: [1,2,3], modalities: ["rna"], profile: "oss-apache"}
     → {job_id: str}

GET  /jobs/{job_id}
     → {status, progress_pct, current_step, results?, figures?}

GET  /jobs/{job_id}/stream
     → Server-Sent Events stream:
       {"event": "tier1_complete", "data": {"n_accepted": 750, "n_routed_tier2": 200}}
       {"event": "scenic_progress", "data": {"pct": 45, "regulons_so_far": 12}}
       {"event": "tier2_complete",  "data": {"n_proven": 12, "n_uncertain": 188}}
       {"event": "complete",        "data": {"accuracy": 0.27, "job_id": "..."}}

GET  /kg
     query: ?cell_type=Treg&tissue=thymus&rule_type=negative&status=ACTIVE
     → [array of KG rules]

POST /kg/rules
     body: KGRule (validated against schema)
     → {rule_id, status: "PENDING_REVIEW"}

DELETE /kg/rules/{rule_id}     → sets status="DEPRECATED" (soft delete only)

GET  /profiles
     → [list of built-in + custom profiles with component details]

POST /profiles
     body: custom profile JSON
     → {profile_id, validation_warnings}

GET  /profiles/{name}/components
     → [{component_id, display_name, license, commercial_ok, enabled, notes}, ...]

GET  /datasets
     → pre-loaded CELLxGENE demo datasets available for playground
```

### Task queue: Dramatiq (replaces Celery — see Sections 22 and 26 for full implementation)

```python
# axiom-playground/api/tasks.py — DRAMATIQ (not Celery)
# Full implementation in Section 22. Key points:
import dramatiq
from dramatiq.brokers.redis import RedisBroker
# @dramatiq.actor(max_retries=0, time_limit=3_600_000)
# Worker launch: dramatiq tasks --processes 1 --threads 2
# Result storage: Redis key axiom:result:{job_id} (not Celery result backend)
# Progress events: redis_client.publish(f"axiom:job:{job_id}", ...) → SSE stream
```

### Local development only (docker-compose.yml)

For local dev, a minimal docker-compose.yml is useful to start Redis and run the Dramatiq
worker alongside the API. **Production deployment uses AWS ECS Fargate — see Section 27.**

```yaml
# docker-compose.yml — LOCAL DEV ONLY
version: "3.9"
services:
  api:
    build: ./api
    ports: ["8000:8000"]
    env_file: .env
    environment:
      PYTHONPATH: /axiom-sc
    volumes:
      - ../axiom-sc:/axiom-sc:ro   # mount SDK in dev mode
    depends_on: [redis]
    command: uvicorn main:app --host 0.0.0.0 --port 8000 --reload

  worker:
    build: ./api
    env_file: .env
    environment:
      PYTHONPATH: /axiom-sc
    volumes:
      - ../axiom-sc:/axiom-sc:ro
    depends_on: [redis]
    command: dramatiq tasks --processes 1 --threads 2   # Dramatiq, not Celery

  app:
    build: ./app
    ports: ["3000:3000"]
    environment:
      - VITE_API_URL=http://localhost:8000
    depends_on: [api]

  redis:
    image: redis:7-alpine
    ports: ["6379:6379"]
```

---

## 16. React playground — component spec

### Components to build (in this exact order)

**1. `SDKVersionBadge`** — always visible in header
```tsx
// GET /health on mount
// Renders: axiom-sc v0.2.0  |  KG: 312 rules  |  Profile: oss-apache
// Click → ChangelogModal showing CHANGELOG.md content
```

**2. `ProfileSelector`** — choose active profile
```tsx
// Dropdown: oss-mit | oss-apache | commercial | custom
// Shows component count per profile: "42 components enabled"
// Link → ProfileCreator
```

**3. `ProfileCreator`** — the key UX innovation
```tsx
// Wizard: "What's your use case?"
//   → Academic research / Startup / Enterprise
//
// Component table (fetched from GET /profiles/{name}/components):
//   | Component           | License       | Commercial OK | Toggle |
//   | AXIOMTier1          | Apache-2.0    | ✓ Yes         | [ON]   |
//   | pySCENIC (subprocess)| GPL-3.0      | ⚠ No          | [ON]   |
//   | CELLama niche       | MIT           | ✓ Yes         | [ON]   |
//   | PanglaoDB (optional) | CC BY-NC 4.0 | ✗ No          | [OFF]  |
//
// License badge colors: MIT/Apache=green | BSD=blue | GPL=orange | NC=red
//
// Validation panel: "Your profile enables 3 non-commercial components.
//                    Commercial use would require: [list]"
//
// Export button: downloads axiom_profile.json
// Save to playground: POST /profiles
```

**4. `FileUpload`**
```tsx
// Two modes:
// A) Upload .h5ad (drag-drop, max 2GB)
// B) Demo dataset: Human Thymus / Lung Cell Atlas / Tabula Sapiens
//    (streamed from CELLxGENE, cached in api/cache/)
```

**5. `TierSelector`**
```tsx
// Tier 1 only       → fastest, common types, ~30s
// Tier 1 + 2        → adds mechanistic verification, ~5 min (DEFAULT)
// Tier 1 + 2 + 3    → multi-stream, requires velocity/ATAC data, ~20 min
// Tier 1 + 2 + 3 + 4 → adds LLM, requires API keys, ~30 min
// All 5 tiers        → full AXIOM-SC, ~45 min
// Shows: estimated runtime + API cost per option
// Disables Tier 3 options if no ATAC/velocity in uploaded adata.obsm
```

**6. `ProgressStream`** — SSE consumer
```tsx
// Real-time pipeline events:
// ● Tier 1 running...
// ✓ Tier 1 complete — 750/1000 cells accepted (75%)
// ● SCENIC+ running... ████░░░░ 45%
// ✓ SCENIC+ complete — 29 regulons recovered (IRF7, GATA3, RORC ✓)
// ● Tier 2 running...
// ✓ Tier 2 complete — 12 PROVEN · 188 UNCERTAIN · 0 CONTRADICTED
```

**7. `UMAPViewer`** — interactive UMAP
```tsx
// Color by: original_labels | tier1_label | tier2_label | confidence | verdict
// Click cluster → opens VerdictPanel
// Hover cell → shows label + confidence tooltip
// Downloadable as SVG
```

**8. `VerdictCards`** — per-cluster annotation explanation
```tsx
// Badge: PROVEN (green) | UNCERTAIN (amber) | CONTRADICTED (red)
// Confidence bar
// CASSIA annotation (if run) vs AXIOM-SC verdict
// Rules fired list:
//   ✓ PASS   ILC1_POS_001 (positive): TBX21 regulon active (z=2.8)
//            Mechanistic basis: TBX21 drives IFN-gamma in ILC1
//            [PubMed: 29958055]
//   ✗ FAIL   ILC1_NEG_001 (negative): TRAC absent — but TRAC found in top markers!
//            Mechanistic basis: ILCs lack VDJ recombination
//            [PubMed: 30033366]
//   — N/T    ILC1_NEG_002 (negative): GATA3 regulon — not in SCENIC output
```

**9. `RuleFireHeatmap`** — Figure 2 from Phase 1, interactive
```tsx
// Rows: clusters. Columns: KG rules.
// Color: PASS=green | FAIL=red | NOT_TESTABLE=grey
// Click cell → shows rule details + PMID link
// Filter: "show only FAIL rules" | "show only PROVEN clusters"
// Export as PNG/SVG
```

**10. `KGBrowser`** — explore and edit the knowledge graph
```tsx
// Table: cell_type | rule_id | rule_type | confidence | tissue | status | actions
// Filter by: cell_type, tissue, rule_type, confidence, status
// Click rule → inline editor (mechanistic_basis, pmid, genes, direction)
// Add rule button → rule creation form with PMID lookup
// Export KG as JSON (versioned filename)
// Stats: "312 rules · 82 cell types · 28 tissues"
```

**11. `ComparisonTable`** — accuracy comparison
```tsx
// CASSIA | mLLMCelltype | AXIOM-SC Tier 1+2 | AXIOM-SC All Tiers
// Rows: per failure category
// Shows improvement over Phase 1 baseline
// Download as CSV
```

**12. `ConfidenceCalibration`** — Figure 3 from Phase 1
```tsx
// Scatter: AXIOM confidence vs correctness
// Shows: "PROVEN at 0.8+ confidence = 67% precision"
// CASSIA QS for comparison
```

---

## 17. Day-by-day execution plan (20 working days)

### WEEK 1: SDK foundation + KG + AXIOMTier1 training launch

**Day 1 — Scaffold both repos**

```bash
# Claude Code: execute these steps exactly

# 1. Create axiom-sc repo structure from Section 4
mkdir -p axiom-sc/axiom_sc/{profiles,tier1/training,tier2,tier3,tier4,tier5,kg,pipelines,utils}
mkdir -p axiom-sc/{kg_data,model_weights,tests/fixtures,.github/workflows}

# 2. Write pyproject.toml from Section 6
# 3. Write axiom_sc/__init__.py exporting AXIOMAnnotator, __version__, kg_version, kg_rule_count
# 4. Write axiom_sc/version.py: __version__ = "0.2.0"
# 5. Write axiom_sc/profiles/registry.py (full COMPONENT_REGISTRY from Section 3b)
# 6. Write axiom_sc/profiles/loader.py and validator.py
# 7. Write .github/workflows/ci.yml and publish.yml from Section 7
# 8. pip install -e ".[dev]"
# 9. Write 3 smoke tests in tests/test_profiles.py:
#    - load oss-apache profile → assert tier2_axiom_kg in components
#    - load oss-mit profile → assert tier3_scvelo in components (BSD is OK in oss-mit? No)
#    - create custom profile → assert validation works
# 10. pytest tests/test_profiles.py — all pass
```

**Day 2 — Copy Phase 1 code into axiom_sc/tier2/**

```
Copy and refactor:
  Phase1: axiom/annotator.py      → axiom_sc/tier2/axiom_annotator.py
  Phase1: axiom/evidence.py       → axiom_sc/tier2/evidence.py
  Phase1: axiom/kg_loader.py      → axiom_sc/tier2/kg_loader.py (add pmid/status fields)
  Phase1: axiom/utils.py          → axiom_sc/utils/anndata_utils.py
  Phase1: kg/axiom_kg.json        → kg_data/oracle_kg_v0.1.0.json
  Phase1: kg/kg_schema.json       → axiom_sc/kg/schema.json (add pmid/status/tissue_context)

Apply all Phase 1 corrections from Section 9 (circuit completeness, support counting, etc.)
Update ILC3_CIRCUIT_001 to require NCR2/NCR3 co-expression (Phase 2 fix from Section 9)
Update PDC_NEG_001 to use PAX5 regulon z-score (Phase 1 correction)

Write tests/test_tier2.py with these 6 test cases:
  1. Treg with IL2 high → CONTRADICTED (TREG_NEG_001)
  2. ILC1 with TRAC present → CONTRADICTED (ILC1_NEG_001)
  3. pDC with IRF7 active + SIGLEC1 high + PAX5 regulon inactive → PROVEN
  4. ILC3 with RORC active + NCR2 present → PROVEN (new Phase 2 fix)
  5. Myofibroblast with RORC active + no NCR2 → NOT PROVEN (Phase 2 fix works)
  6. Full-KG mode: exhaust all 18 types, return highest-support verdict

pytest tests/test_tier2.py — all 6 pass
```

**Day 3 — KG seeding pipeline + start immune KG**

```
1. Implement axiom_sc/kg/seeder.py from Section 13a
2. Implement axiom_sc/kg/review_cli.py interactive review tool
3. Implement axiom_sc/kg/references.py PMID lookup via biopython Entrez
4. Run seeder for 5 tissue types → kg_data/candidates/cellmarker2_candidates.json
5. Start review queue: review + accept 60 immune T-cell rules (researcher task)
   Claude Code assists: look up PMIDs, generate mechanistic_basis text,
   cross-check against Phase 1 validated rules for conflicts
```

**Day 4 — Continue immune KG (researcher + Claude Code)**

```
Review and accept:
  - B cell maturation rules (30 target)
  - NK/ILC/NKT rules (40 target)
  - Myeloid (40 target)

End of Day 4 target: ≥170 ACTIVE rules in kg_data/oracle_kg_v0.2.0.json
Run axiom-kg references → REFERENCES.md generated
```

**Day 5 — AXIOMTier1 training (GPU job + Tier 1 inference code)**

```
Training (run on Vast.ai A100 — takes ~50 GPU hours, start now, runs overnight/weekend):
1. Write axiom_sc/tier1/training/train.py from Section 8
2. Write axiom_sc/tier1/training/data_loader.py (Census download, stratified split)
3. Submit training job to Vast.ai: 10 models × ~5h each
4. Cost: ~$75

In parallel (local, CPU):
5. Write axiom_sc/tier1/axiomtier1.py inference class (loads weights)
6. Write axiom_sc/tier1/census_annotator.py (Census kNN fallback while weights train)
7. Write tests/test_tier1.py with mock weight loading
8. Tier 1 works end-to-end with Census kNN while AXIOMTier1 weights training
```

---

### WEEK 2: SCENIC+ + Tier 3 + epithelial KG

**Day 6 — pySCENIC subprocess isolation**

```
1. Create scenic-env conda env: conda create -n scenic-env python=3.10 && pip install pyscenic
2. Write axiom_sc/tier2/scenic_runner.py (subprocess wrapper, Apache 2.0 side)
3. Write axiom_sc/tier2/scenic_worker.py (runs inside scenic-env, GPL side)
4. Write axiom_sc/tier2/scenic_pipeline.py with SCENIC_DEFAULTS from Section 9
5. Write axiom_sc/utils/stratified_sample.py (≥300 cells/type)
6. Write axiom_sc/utils/zscore_auc.py (Phase 1 validated normalization)

Test: run scenic_runner on Phase 1 thymus h5ad at 50k cells
Expected: FOXP3 recovered (was zero at 20k — this is the Phase 2 validation test)
Print: regulon recovery report for all 26 forced genes
```

**Days 7–8 — Tier 3 streams**

```
Day 7:
1. Write axiom_sc/tier3/velocity.py (scVelo, BSD)
   Key: exhausted vs progenitor-exhausted T velocity discrimination
2. Write axiom_sc/tier3/chromatin.py (Signac via rpy2, MIT)
   Key: FOXP3/AIRE locus accessibility (fixes Phase 1 SCENIC misses)
3. Write axiom_sc/tier3/convergence.py framework
4. Write tests/test_tier3.py

Day 8:
5. Write axiom_sc/tier3/communication.py (COMMOT, MIT)
6. Write axiom_sc/tier3/spatial_niche.py (CELLama, MIT)
7. Write axiom_sc/tier3/cross_species.py (OrthoFinder, MIT)
8. Hook all streams into convergence.py
9. Run Tier 3 on Phase 1 datasets — measure stream agreement rates
```

**Days 9–10 — Epithelial + stromal KG + Tiers 4 & 5 scaffold**

```
Day 9: KG expansion
  - Review + accept ~80 lung epithelial rules (AT1, AT2, Club, Ciliated, etc.)
  - Review + accept ~30 liver rules (hepatocyte zones 1-3, LSEC, Kupffer)
  - Total KG target end of Day 9: ≥250 ACTIVE rules

Day 10: Tier 4 + 5 scaffold
1. Write axiom_sc/tier4/llm_ensemble.py (multi-model API calls)
2. Write axiom_sc/tier4/evidence_bundler.py (formats Tier 1-3 evidence for LLM)
3. Write axiom_sc/tier4/rag.py (CellMarker 2.0 RAG for rare types)
4. Write axiom_sc/tier4/prompts.py (reproducible prompt templates)
5. Write axiom_sc/tier5/attractor_discovery.py
6. Write axiom_sc/tier5/candidate_state.py
7. Write axiom_sc/tier5/feedback.py (Tier 5 → KG loop)
```

---

### WEEK 3: FastAPI + React playground + AWS infrastructure

**Day 11 — FastAPI backend (Dramatiq)**

```
Full implementation spec in Section 22. Summary:
1. Scaffold axiom-playground/api/ with Dramatiq tasks.py (NOT Celery)
2. Write all 5 route files: annotate, jobs, kg, profiles, health
3. Write config.py with SDK version validation
4. Add .aws/ task definition templates and .github/workflows/deploy.yml
5. Local verification:
   docker run -d -p 6379:6379 redis:7-alpine
   PYTHONPATH=../axiom-sc uvicorn main:app --reload
   dramatiq tasks --processes 1 --threads 2   ← Dramatiq worker
   curl http://localhost:8000/health
   → {sdk_version, kg_rule_count: ≥495, tier1_backend: "census_knn"}
```

**Day 12 — API tests + AWS infrastructure**

```
Morning: write tests/test_api.py
   GET /health → sdk_version present, kg_rule_count ≥495
   GET /kg?cell_type=Treg → TREG_NEG_001 and TREG_CIRCUIT_001 present
   POST /annotate (upload fixture h5ad) → job_id returned
   GET /jobs/{id}/stream → SSE events flow
   Verify Verdict serializes as "PROVEN" not "Verdict.PROVEN" (Python 3.10 gotcha)

Afternoon: one-time AWS setup (Section 27)
   aws ecr create-repository (×3: api, worker, app)
   aws ecs create-cluster
   aws elasticache create-serverless-cache
   aws efs create-file-system
   aws s3 mb s3://axiom-sc-jobs + lifecycle rule (1-day TTL)
   aws secretsmanager create-secret (anthropic-api-key, openai-api-key)
   Create IAM roles: OIDC deploy role + ECS task/execution roles
   Register task definitions from .aws/*.json
   Create ECS services (desired_count=0 initially)
   Set GitHub secrets: AWS_ACCOUNT_ID, ALB_DNS_NAME
```

**Days 13–15 — React frontend + GitHub Actions wiring**

```
Day 13: GitHub Actions CI/CD + React foundation
  Write .github/workflows/deploy.yml (Section 27 — OIDC, no stored keys)
  Push to main → verify 3 ECR images build → ECS services start → ALB reachable
  React: SDKVersionBadge, ProfileSelector, ProfileCreator (license table + validation panel)
  Zustand store + src/api/client.ts

Day 14: Upload + pipeline + results
  FileUpload (drag-drop .h5ad/.h5/.csv + demo dataset selector)
  InspectionPanel — shown immediately after upload: detected input type, cell count,
    active streams, estimated time, preprocessing steps (from POST /annotate response)
  TierSelector — disables Tier 3 if no velocity/ATAC; disables Tier 4 if no API keys
  ProgressStream — SSE consumer showing "SCENIC+ via RunPod: 69 regulons, FOXP3 ✓ z=5.74"
  UMAPViewer + VerdictCards (rules fired with PMID links)
  Each push → GitHub Actions deploys to ECS staging automatically

Day 15: KG browser + comparison + polish
  RuleFireHeatmap, KGBrowser (inline editor, add rule form with PMID validation),
  ComparisonTable, ConfidenceCalibration
  Full e2e via ALB DNS: upload thymus h5ad → PROVEN pDC in VerdictCards
  Fix bugs
```

---

### WEEK 4: AXIOMTier1 weights + integration + release

**Day 16 — AXIOMTier1 weights arrive (training completes)**

```
1. Download trained weights from Vast.ai job output
2. Evaluate ensemble on held-out test set
   Target: macro-F1 ≥ 0.80 (vs scTab published 0.83)
   If F1 < 0.75: extend training, check data quality
3. Upload weights to HuggingFace: receptor-bio/axiomtier1 (Apache 2.0)
4. Update axiom_sc/tier1/axiomtier1.py to download from HuggingFace on first use
5. Replace Census kNN fallback with AXIOMTier1 as primary Tier 1
6. Re-run end-to-end: verify Tier 1 acceptance rate ≥70% on Phase 1 datasets
```

**Days 17–18 — Full integration test + RunPod GPU wiring**

```
Day 17:
1. Full pipeline: thymus h5ad → Tiers 1+2+3 → verify results
2. Expected:
   - Tier 1: ≥70% cells accepted at ≥0.85 confidence
   - Tier 2: pDC (thy-22) PROVEN (IRF7 z=3.7, SIGLEC1 high, PMID link visible in UI)
   - Tier 2: ILC3 (tab-20) PROVEN (RORC z=3.37)
   - Tier 2: false PROVEN myofibroblast eliminated (NCR2/NCR3 fix works)
   - FOXP3 regulon recovered at 50k cells (Phase 2 target)
3. Compare all accuracy metrics to Phase 1 baseline — must show improvement

Day 18: RunPod/Modal GPU dispatch setup (Section 28)
4. Build and push RunPod Docker image for pySCENIC handler
   (SCENIC resources — rankings.feather 1.3GB + motif DB — baked in at build time)
5. Create RunPod serverless endpoint → get RUNPOD_SCENIC_ENDPOINT_ID
6. Add to AWS Secrets Manager: runpod-api-key + runpod-scenic-endpoint-id
7. Add Secrets Manager ARNs to ECS worker task definition
8. Test: upload thymus h5ad → verify pySCENIC dispatches to RunPod
   ProgressStream shows "SCENIC+ via RunPod RTX 4090: 69 regulons, FOXP3 ✓ z=5.74"
9. Tier 4: test LLM ensemble on 5 hard-case clusters
10. Tier 5: test novel attractor on 2 UNCERTAIN clusters
11. ProfileCreator e2e: create academic profile → verify license validation works
```

**Days 19–20 — KG extension + release**

```
Day 19:
1. Begin KG universe extension (Section 25): retina, bone, gonads priority
2. Final KG target: ≥500 ACTIVE rules, ≥155 cell types
   (stretch: ≥640 rules / ≥198 types per Section 25 full target)
3. Run pytest: must achieve ≥80% coverage, 0 failures
4. Run axiom-kg references → regenerate REFERENCES.md
5. Update CHANGELOG.md

Day 20 — SDK release 0.2.0:
6. In axiom-sc repo:
   git tag v0.2.0 && git push origin v0.2.0
   → GitHub Action publishes axiom-sc==0.2.0 to PyPI automatically

7. In axiom-playground repo:
   echo "sdk_version: 0.2.0" > config.yaml
   git add config.yaml && git commit -m "chore: bump to axiom-sc 0.2.0"
   git push origin main
   → GitHub Actions: requires production environment approval
   → Approve → ECS rolling update, zero downtime
   → /health returns sdk_version: 0.2.0 from PyPI package
   → React SDKVersionBadge shows v0.2.0 (green background, not amber "dev")
```

---

## 18. Definition of Phase 2 success

| Criterion | Target | Strong claim |
|---|---|---|
| axiom-sc published to PyPI | 0.2.0 ✓ | — |
| Profile system working | 3 built-in + custom creator ✓ | — |
| AXIOMTier1 trained weights | HuggingFace receptor-bio/axiomtier1, Apache 2.0 | F1 ≥ 0.80 |
| Tier 1 acceptance rate | ≥70% cells at ≥0.85 confidence | ≥78% |
| Tier 2 in-scope accuracy | ≥25% (Phase 1 baseline: 22.2%) | ≥40% |
| FOXP3 regulon recovered | In ≥2 datasets at 50k cells | All 3 |
| KG rules | ≥500 ACTIVE (patch: 495) | ≥640 (Section 25 full extension) |
| Cell types | ≥155 (patch: 153) | ≥198 (Section 25) |
| Playground running | ECS Fargate via ALB — upload → results < 5 min (Tier 1+2) | — |
| Task queue | Dramatiq (not Celery) — zero task-duplication bugs | — |
| GPU dispatch | pySCENIC via RunPod/Modal — cost ~$0.20–0.60/job | — |
| Secrets | AWS Secrets Manager (not SSM) for all confidential config | — |
| CI/CD | OIDC GitHub Actions → ECR → ECS — no stored AWS keys | — |
| Tier 4 LLM | Works on 5 test hard cases | — |
| Tier 5 discovery | Returns characterized output on UNCERTAIN cells | — |
| Apache 2.0 clean | Zero GPL imports in axiom_sc/ (subprocess only) | — |
| Citation system | Every KG rule has pmid, REFERENCES.md generated | — |
| SDK version loop | Tag → PyPI → playground config → /health verified | — |

---

## 19. What NOT to implement in Phase 2

Defer to Phase 3:
- DISCO API integration (no commercial license confirmed for DISCO platform)
- scTab model weights (no license — use AXIOMTier1 instead)
- ArchR integration (non-commercial — use Signac instead, already specified)
- Human-in-the-loop annotation review queue (Tier 4 extension)
- Perturb-seq perturbation prediction module (Tier 5 extension)
- Multi-species training of AXIOMTier1 (Phase 2 trains human only)
- CellMaster / DeepCellSeek API integration (Phase 3 Tier 4 option)

---

## 20. Cost estimate

| Item | Cost |
|---|---|
| AXIOMTier1 training (10× A100, ~50h total) | ~$75 |
| pySCENIC GPU Phase 2 (50k × 3 datasets, RTX 4090) | ~$12 |
| RunPod Docker image build + SCENIC resources download | ~$5 |
| CASSIA comparison baseline runs | ~$5 |
| LLM API testing (Tier 4, 20 test clusters) | ~$8 |
| AWS ECS Fargate (api + worker + app, 1 month dev) | ~$30 |
| AWS ElastiCache Redis serverless (1 month dev) | ~$15 |
| AWS ECR storage (3 images) | ~$3 |
| Playground annotation cost per job (production) | ~$0.20–0.60/job |
| Claude.ai Pro + Cursor Pro (1 month) | $40 |
| **Total direct costs (Phase 2 build)** | **~$193** |
| **Calendar time (solo researcher + Claude Code)** | **~20 working days** |

---

## 21. Implementation notes — Days 1–5 (coding record)

This section records actual implementation decisions made during coding.
These **override or clarify** the spec above where they differ.

### 21a. Environment / tooling

- **Python interpreter**: always use `/Users/sivamoturi/miniconda/bin/python3.10`.
  System python3 on this machine is 3.9 and lacks project dependencies.
- **pyomnipath does not exist on PyPI** — correct package name is `omnipath`.
  It was moved to `[science]` optional extra (not core deps) to avoid install failures.
- **biopython** is in `[science]` optional extra; tests mock network calls.
- **anndata** must be pip-installed separately (`pip install anndata`) — not auto-installed with base deps.

### 21b. pyproject.toml additions (not in original spec)

```toml
[tool.pytest.ini_options]
testpaths = ["tests"]
addopts = "--tb=short"

[tool.coverage.run]
source = ["axiom_sc"]
omit = [
    "axiom_sc/annotator.py",           # Day 10: full pipeline stub
    "axiom_sc/tier1/*",                 # Day 5+: remove once training complete
    "axiom_sc/tier2/scenic_runner.py",  # Day 6
    "axiom_sc/tier2/scenic_worker.py",  # Day 6
    "axiom_sc/tier2/sctype_scorer.py",  # Day 6
    "axiom_sc/tier2/gnn_propagator.py", # Day 7
    "axiom_sc/tier3/*",                 # Days 7-8
    "axiom_sc/tier4/*",                 # Day 10
    "axiom_sc/tier5/*",                 # Day 10
    "axiom_sc/utils/anndata_utils.py",  # Day 6
    "axiom_sc/pipelines/*",             # Day 6
    "axiom_sc/kg/review_cli.py",        # interactive terminal CLI
    "axiom_sc/cli.py",
]

[tool.coverage.report]
fail_under = 80
show_missing = true
exclude_lines = ["pragma: no cover", "raise NotImplementedError", ...]
```

CI workflow updated to use `--cov-config=pyproject.toml` (no `--cov-fail-under` on the command line).

### 21c. KG biological decisions (non-negotiable)

**Pre_B — intentionally has NO circuit rule.**
VPREB1 and IGLL1 (surrogate light chain subunits) are co-expressed in BOTH Pro_B AND Pre_B.
A VPREB1+IGLL1 circuit would fire PROVEN on Pro_B cells → false positive.
The biologically correct Pre_B circuit is `["IGHM", "VPREB1"]` (μ-heavy chain expression confirms
completed VDJ recombination). This requires adding IGHM as a positive rule first.
Current correct behaviour: Pre_B returns UNCERTAIN from Tier 2; Tier 3 scVelo handles the boundary.

**ILC3_CIRCUIT_001 requires NCR2** — `gene_or_regulon: ["RORC", "NCR2"]`.
Without NCR2, myofibroblasts with circadian RORC expression get false PROVEN as ILC3.

**PDC_NEG_001 checks PAX5 REGULON inactive** (z ≤ -0.5), not IGKC absence.
pDCs naturally express IGKC due to developmental B-lineage relationship.

**M1 ↔ M2 cross-contradiction is intentional.**
M1MAC_NEG_001: MRC1 high → M1 CONTRADICTED. M2MAC_NEG_001: IDO1 high → M2 CONTRADICTED.

### 21d. Tier 2 annotator: tested cell type coverage (end of Day 4)

220 ACTIVE rules, 43 cell types. All have PROVEN paths tested except Pre_B (intentional — see above).
PROVEN requires: circuit_pass AND (positive_pass OR no positive rules for that type).

### 21e. Tier 1 architecture decisions

**AXIOMTier1Net** outputs raw logits. Softmax is applied in `AXIOMTier1Ensemble._run_ensemble()`.
Do not add softmax inside the model's forward() — training uses CrossEntropyLoss which expects logits.

**Confidence formula**: `1 - H(p)/log(n_classes)` where H = Shannon entropy.
Uniform distribution → confidence ≈ 0. Peaked distribution → confidence ≈ 1.

**model_config.json** must be saved alongside model weight files:
```json
{"hidden": 512, "n_models": 10, "n_genes": 2000}
```
`load_weights()` reads this to reconstruct the correct architecture. Without it, defaults to hidden=512
which will fail for any model trained with a different hidden dim.

**Training not yet run.** `axiom-train-tier1` CLI is written and ready to submit to GPU cluster.
While waiting for weights: `CensusAnnotator` (kNN over Census scVI embeddings) is Tier 1 primary.

**Routing thresholds** (hardcoded, do not change without re-validation):
- ≥ 0.85 → "accept" (Tier 1 terminal, ~75-78% of cells)
- 0.50–0.85 → "tier2_verify" (Tier 2 checks the Tier 1 call)
- < 0.50 → "tier2_full" (Tier 2 full KG search, Tier 1 call not trusted)

### 21f. Test infrastructure

Tests always use `/Users/sivamoturi/miniconda/bin/python3.10 -m pytest`.
`pytest.importorskip("anndata")` / `pytest.importorskip("torch")` used for optional deps.
Tier 1 tests inject mock models (N_GENES=64, HIDDEN=32, N_CLASSES=10) — run on CPU, no weights.

**Current totals (end of Day 5):** 102 passing, 0 failing, 85.84% coverage.

### 21g. Day 6 priorities (next session)

Per CLAUDE.md Section 17 (Day 6):
1. Create `scenic-env` conda env: `conda create -n scenic-env python=3.10 && pip install pyscenic`
2. Implement `axiom_sc/tier2/scenic_runner.py` (subprocess wrapper — Apache 2.0 side)
3. Implement `axiom_sc/tier2/scenic_worker.py` (runs inside scenic-env — GPL side)
4. Implement `axiom_sc/tier2/scenic_pipeline.py` with all `SCENIC_DEFAULTS` from Section 9
5. Implement `axiom_sc/utils/stratified_sample.py` (≥300 cells/type)
6. Implement `axiom_sc/utils/zscore_auc.py` (Phase 1 validated normalization)
7. Validation test: run scenic_runner on Phase 1 thymus h5ad at 50k cells → FOXP3 recovered

Remove `axiom_sc/tier2/scenic_runner.py`, `scenic_worker.py`, `sctype_scorer.py` from
`[tool.coverage.run] omit` in pyproject.toml after implementing them.

---

## 22. Day 11 — FastAPI backend (Dramatiq, not Celery)

### Task queue decision: Dramatiq replaces Celery

Celery has documented task-duplication bugs in advanced Canvas features, its original
maintainer is no longer active, and issues are routinely left unresolved. Dramatiq
gives identical SSE pub/sub semantics, ~10× better throughput, and a simpler API.
All Dramatiq code uses the same Redis broker and identical Redis pub/sub for SSE.

### Playground repo setup

```bash
cd /Users/sivamoturi/Documents/siva-wrk/
mkdir axiom-playground && cd axiom-playground
git init && cp ../axiom-sc/LICENSE .
mkdir api app .aws
```

### api/requirements.txt

```
# axiom-sc installed from local mount in dev; PyPI in production
axiom-sc @ file:///Users/sivamoturi/Documents/siva-wrk/axiom-sc
fastapi>=0.110
uvicorn[standard]>=0.27
redis>=5.0
dramatiq[redis]>=1.16          # NOT celery
python-multipart>=0.0.9
aiofiles>=23.2
httpx>=0.26
pyyaml>=6.0
boto3>=1.34                     # for S3 (loom ↔ RunPod) and Secrets Manager
```

### config.yaml (project root)

```yaml
sdk_version: dev               # change to "0.2.0" once published to PyPI
default_profile: oss-apache
demo_datasets:
  - name: "Human Thymus Cell Atlas (Phase 1)"
    path: ../axiom-sc/tests/fixtures/thymus_mini.h5ad
  - name: "Lung Cell Atlas subset"
    path: ../axiom-sc/tests/fixtures/lung_mini.h5ad
max_upload_bytes: 2147483648   # 2 GB
```

### .env.example (copy to .env — never commit)

```
AXIOM_SDK_VERSION=dev
AXIOM_PROFILE=oss-apache
REDIS_URL=redis://localhost:6379
# Secrets — populated from AWS Secrets Manager in production; set locally for dev:
ANTHROPIC_API_KEY=
OPENAI_API_KEY=
# GPU dispatch (leave blank for local subprocess mode):
RUNPOD_API_KEY=
RUNPOD_SCENIC_ENDPOINT_ID=
MODAL_TOKEN_ID=
MODAL_TOKEN_SECRET=
# S3 for loom/AUC exchange with GPU workers:
AXIOM_S3_BUCKET=axiom-sc-jobs
```

### api/config.py

```python
import yaml, os, importlib.metadata
from pathlib import Path

def load_config() -> dict:
    with open(Path(__file__).parent.parent / "config.yaml") as f:
        cfg = yaml.safe_load(f)
    configured = os.getenv("AXIOM_SDK_VERSION") or cfg.get("sdk_version", "dev")
    if configured != "dev":
        installed = importlib.metadata.version("axiom-sc")
        if installed != configured:
            raise RuntimeError(
                f"SDK mismatch: config={configured}, installed={installed}\n"
                f"Run: pip install axiom-sc=={configured}"
            )
    return cfg
```

### api/tasks.py (Dramatiq)

```python
import dramatiq
from dramatiq.brokers.redis import RedisBroker
import json, redis as redis_lib, traceback, os

REDIS_URL = os.getenv("REDIS_URL", "redis://localhost:6379/0")
broker = RedisBroker(url=REDIS_URL)
dramatiq.set_broker(broker)
redis_client = redis_lib.Redis.from_url(REDIS_URL)

def emit(job_id: str, event: str, data: dict):
    redis_client.publish(f"axiom:job:{job_id}", json.dumps({"event": event, **data}))

@dramatiq.actor(max_retries=0, time_limit=3_600_000)  # 60-min hard limit; no retries
def run_annotation_job(job_id, h5ad_path, tiers, inspection_dict, profile):
    import axiom_sc, anndata
    from axiom_sc.input.detector import InputInspection, InputType
    from axiom_sc.input.preprocessor import prepare_for_annotation

    emit(job_id, "loading", {"message": "Loading dataset...", "progress_pct": 5})
    try:
        adata = anndata.read_h5ad(h5ad_path)
        inspection = InputInspection(**inspection_dict)

        emit(job_id, "preprocessing",
             {"message": f"Preprocessing {adata.n_obs} cells...", "progress_pct": 10})
        adata = prepare_for_annotation(adata, inspection)

        annotator = axiom_sc.AXIOMAnnotator(profile=profile, tiers=tiers)
        result = annotator.annotate(
            adata, progress_callback=lambda ev, d: emit(job_id, ev, d)
        )
        output = result.to_serializable_dict()
        emit(job_id, "complete", {"progress_pct": 100, "summary": output.get("summary", {})})
        redis_client.setex(f"axiom:result:{job_id}", 3600, json.dumps(output))
        return output
    except Exception as e:
        emit(job_id, "failed", {"error": str(e), "traceback": traceback.format_exc()})
        redis_client.setex(f"axiom:result:{job_id}", 3600,
                           json.dumps({"status": "failed", "error": str(e)}))
        raise
```

### api/main.py

```python
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from contextlib import asynccontextmanager
from config import load_config
import axiom_sc

@asynccontextmanager
async def lifespan(app: FastAPI):
    cfg = load_config()
    app.state.cfg            = cfg
    app.state.sdk_version    = axiom_sc.__version__
    app.state.kg_version     = axiom_sc.kg_version
    app.state.kg_rule_count  = axiom_sc.kg_rule_count()
    app.state.active_profile = cfg.get("default_profile", "oss-apache")
    yield

app = FastAPI(title="AXIOM-SC Playground", version="0.2.0", lifespan=lifespan)
app.add_middleware(CORSMiddleware,
    allow_origins=["http://localhost:3000", "http://localhost:5173"],
    allow_credentials=True, allow_methods=["*"], allow_headers=["*"])

from routes import health, annotate, jobs, kg, profiles
for r in [health.router, annotate.router, jobs.router, kg.router, profiles.router]:
    app.include_router(r)
```

### api/routes/health.py

```python
from fastapi import APIRouter, Request
router = APIRouter()

@router.get("/health")
async def health(request: Request):
    return {
        "status": "ok",
        "sdk_version":    request.app.state.sdk_version,
        "kg_version":     request.app.state.kg_version,
        "kg_rule_count":  request.app.state.kg_rule_count,
        "active_profile": request.app.state.active_profile,
        "tiers_available": [1, 2, 3, 4, 5],
        "tier1_backend":  "census_knn",   # → "axiomtier1" once GPU weights arrive
    }
```

### api/routes/annotate.py

```python
from fastapi import APIRouter, UploadFile, Form, Request, HTTPException
from pathlib import Path
import uuid, os, aiofiles
from tasks import run_annotation_job
from axiom_sc.input.detector import inspect_input

router = APIRouter(prefix="/annotate")
UPLOAD_DIR = "/tmp/axiom_uploads"
os.makedirs(UPLOAD_DIR, exist_ok=True)

@router.post("")
async def annotate(
    request: Request, file: UploadFile,
    tiers:   str = Form(default="1,2"),
    profile: str = Form(default="oss-apache"),
):
    content = await file.read()
    if not any(file.filename.endswith(e) for e in [".h5ad", ".h5", ".csv"]):
        raise HTTPException(400, "Accepted: .h5ad, .h5 (CellRanger), .csv (markers)")
    max_bytes = request.app.state.cfg.get("max_upload_bytes", 2_147_483_648)
    if len(content) > max_bytes:
        raise HTTPException(413, f"File too large (max {max_bytes/1e9:.1f} GB)")

    job_id = str(uuid.uuid4())[:8]
    file_path = f"{UPLOAD_DIR}/{job_id}{Path(file.filename).suffix}"
    async with aiofiles.open(file_path, "wb") as f:
        await f.write(content)

    inspection = inspect_input(file_path)   # <2s, runs before job is queued
    run_annotation_job.send(
        job_id, file_path, [int(t) for t in tiers.split(",")],
        inspection.__dict__, profile
    )
    return {
        "job_id": job_id, "status": "queued",
        "inspection": {
            "input_type":          inspection.input_type,
            "n_cells":             inspection.n_cells,
            "n_clusters":          inspection.n_clusters,
            "preprocessing_needed":inspection.preprocessing_needed,
            "estimated_minutes":   inspection.estimated_minutes,
            "warnings":            inspection.warnings,
        }
    }
```

### api/routes/jobs.py

```python
from fastapi import APIRouter
from fastapi.responses import StreamingResponse
import json, asyncio, redis as redis_lib, os

router = APIRouter(prefix="/jobs")
redis_client = redis_lib.Redis.from_url(os.getenv("REDIS_URL","redis://localhost:6379"))

@router.get("/{job_id}")
async def get_job(job_id: str):
    raw = redis_client.get(f"axiom:result:{job_id}")
    if raw:
        return {"job_id": job_id, "status": "complete", **json.loads(raw)}
    return {"job_id": job_id, "status": "running", "progress_pct": 0}

@router.get("/{job_id}/stream")
async def stream_job(job_id: str):
    async def gen():
        ps = redis_client.pubsub()
        ps.subscribe(f"axiom:job:{job_id}")
        try:
            while True:
                msg = ps.get_message(ignore_subscribe_messages=True)
                if msg:
                    data = json.loads(msg["data"])
                    yield f"data: {json.dumps(data)}\n\n"
                    if data.get("event") in ("complete", "failed"):
                        break
                await asyncio.sleep(0.2)
        finally:
            ps.unsubscribe()
    return StreamingResponse(gen(), media_type="text/event-stream",
        headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"})
```

### api/routes/kg.py

```python
from fastapi import APIRouter, HTTPException
from typing import Optional
import axiom_sc

router = APIRouter(prefix="/kg")

@router.get("")
async def get_kg(cell_type: Optional[str]=None, tissue: Optional[str]=None,
                 rule_type: Optional[str]=None, status: str="ACTIVE"):
    rules = axiom_sc.list_kg_rules(cell_type=cell_type, tissue=tissue,
                                    rule_type=rule_type, status=status)
    return {"rules": rules, "total": len(rules)}

@router.post("/rules")
async def add_rule(rule: dict):
    try:
        r = axiom_sc.add_kg_rule(rule)
        return {"rule_id": r["rule_id"], "status": "PENDING_REVIEW"}
    except ValueError as e:
        raise HTTPException(422, str(e))

@router.delete("/rules/{rule_id}")
async def deprecate_rule(rule_id: str):
    axiom_sc.deprecate_kg_rule(rule_id)
    return {"rule_id": rule_id, "status": "DEPRECATED"}
```

### api/routes/profiles.py

```python
from fastapi import APIRouter
import axiom_sc

router = APIRouter(prefix="/profiles")

@router.get("")
async def list_profiles():
    return axiom_sc.list_profiles()

@router.get("/{name}/components")
async def get_profile_components(name: str):
    return axiom_sc.get_profile_components(name)

@router.post("")
async def create_profile(profile: dict):
    return axiom_sc.create_custom_profile(profile)
```

### Day 11 local verification

```bash
cd axiom-playground/api
cp .env.example .env
docker run -d -p 6379:6379 redis:7-alpine   # local Redis
pip install -r requirements.txt
PYTHONPATH=/Users/sivamoturi/Documents/siva-wrk/axiom-sc \
  uvicorn main:app --reload --port 8000

# Verify:
curl http://localhost:8000/health
# → {"status":"ok","sdk_version":"0.2.0","kg_rule_count":495,...}

curl "http://localhost:8000/kg?cell_type=Treg&status=ACTIVE"
# → {"rules":[...TREG_NEG_001, TREG_CIRCUIT_001...],"total":5}

# Start Dramatiq worker:
PYTHONPATH=/Users/sivamoturi/Documents/siva-wrk/axiom-sc \
  dramatiq tasks --processes 1 --threads 2
```

### Day 12 — API tests

```python
# tests/test_api.py
import pytest
from fastapi.testclient import TestClient
from main import app

client = TestClient(app)

def test_health():
    r = client.get("/health")
    assert r.status_code == 200
    d = r.json()
    assert "sdk_version" in d
    assert "kg_rule_count" in d
    assert d["kg_rule_count"] >= 495

def test_kg_returns_treg_rules():
    r = client.get("/kg?cell_type=Treg&status=ACTIVE")
    assert r.status_code == 200
    rules = r.json()["rules"]
    rule_ids = [r["rule_id"] for r in rules]
    assert "TREG_NEG_001" in rule_ids     # FOXP3 + IL2 contradiction rule
    assert "TREG_CIRCUIT_001" in rule_ids

def test_verdict_serialization():
    # Verdict enum must serialize as "PROVEN" not "Verdict.PROVEN" (Python 3.10)
    from axiom_sc.tier2.axiom_annotator import AxiomAnnotator, EvidenceBundle
    annotator = AxiomAnnotator()
    # ...build minimal evidence... verify .to_serializable_dict()["verdict"] == "PROVEN"

def test_annotate_upload_and_stream(tmp_path):
    # Upload thymus mini h5ad, verify job queued, events stream
    import anndata, numpy as np
    adata = anndata.AnnData(X=np.random.rand(50, 100).astype("float32"))
    h5ad = tmp_path / "test.h5ad"
    adata.write_h5ad(h5ad)
    with open(h5ad, "rb") as f:
        r = client.post("/annotate", files={"file": ("test.h5ad", f, "application/octet-stream")})
    assert r.status_code == 200
    assert "job_id" in r.json()
```

---

## 23. Days 13–17 — React playground

### Setup

```bash
cd axiom-playground/app
npm create vite@latest . -- --template react-ts
npm install @tanstack/react-query zustand plotly.js-dist-min react-dropzone \
  tailwindcss lucide-react react-hot-toast @visx/visx
npm install -D typescript @types/react @types/node
```

### src/api/client.ts

```typescript
const API = import.meta.env.VITE_API_URL ?? "http://localhost:8000";

export const getHealth = () => fetch(`${API}/health`).then(r => r.json());
export const getKG = (p?: Record<string,string>) =>
  fetch(`${API}/kg${p ? "?"+new URLSearchParams(p) : ""}`).then(r => r.json());
export const postAnnotate = (file: File, tiers: number[], profile: string) => {
  const form = new FormData();
  form.append("file", file);
  form.append("tiers", tiers.join(","));
  form.append("profile", profile);
  return fetch(`${API}/annotate`, { method:"POST", body: form }).then(r => r.json());
};
export const streamJob = (id: string, onEvent: (e:any) => void) => {
  const es = new EventSource(`${API}/jobs/${id}/stream`);
  es.onmessage = e => onEvent(JSON.parse(e.data));
  es.onerror = () => es.close();
  return () => es.close();
};
```

### src/store/index.ts

```typescript
import { create } from "zustand";
interface AppStore {
  sdkVersion: string; kgRuleCount: number; activeProfile: string;
  currentJobId: string|null; jobStatus: string; jobProgress: number;
  jobEvents: any[]; results: any|null; selectedCluster: string|null;
  inspection: any|null;
  setHealth: (h:any) => void;
  setJob: (id:string) => void;
  pushEvent: (e:any) => void;
  selectCluster: (id:string|null) => void;
}
export const useStore = create<AppStore>(set => ({
  sdkVersion:"", kgRuleCount:0, activeProfile:"",
  currentJobId:null, jobStatus:"idle", jobProgress:0,
  jobEvents:[], results:null, selectedCluster:null, inspection:null,
  setHealth: h => set({ sdkVersion:h.sdk_version, kgRuleCount:h.kg_rule_count,
                         activeProfile:h.active_profile }),
  setJob: id => set({ currentJobId:id, jobStatus:"queued", jobProgress:0,
                       jobEvents:[], results:null }),
  pushEvent: e => set(s => ({
    jobEvents: [...s.jobEvents, e],
    jobProgress: e.progress_pct ?? s.jobProgress,
    jobStatus: e.event === "complete" ? "complete"
             : e.event === "failed"   ? "failed" : "running",
    results: e.event === "complete" ? e.results ?? s.results : s.results,
  })),
  selectCluster: id => set({ selectedCluster: id }),
}));
```

### Component build order (Day 13–15)

Build strictly in this sequence — each component depends on the previous:

**Day 13:**
1. `SDKVersionBadge` — GET /health, shows `axiom-sc v0.2.0 | KG: 495 rules | oss-apache`
2. `ProfileSelector` — dropdown of built-in + custom profiles
3. `ProfileCreator` — component table with license badges (Apache=green, GPL=orange, NC=red),
   commercial-ok toggles, validation panel showing non-commercial components count.
   Exports `axiom_profile.json`. This is the UX centrepiece.
4. `FileUpload` — drag-drop .h5ad/.h5/.csv + demo dataset selector.
   Shows `InspectionPanel` immediately after upload (input type, cell count, cluster count,
   estimated minutes, active streams, warnings).
5. `TierSelector` — Tier 1 / 1+2 (default) / 1+2+3 / all 5 with runtime estimates.
   Disables Tier 3 if no velocity/ATAC in inspection; disables Tier 4 if no API keys.

**Day 14:**
6. `ProgressStream` — SSE consumer showing timeline:
   ```
   ✓ Tier 1: 33,421/45,230 cells accepted (73.9%)
   ● SCENIC+ running... ████░░ 62% (via RunPod RTX 4090)
   ✓ SCENIC+: 69 regulons — FOXP3 ✓ z=5.74
   ✓ Tier 2: 8 PROVEN · 12,809 UNCERTAIN
   ```
7. `UMAPViewer` — plotly.js scatter, color by tier1/tier2/verdict/confidence.
   Click cluster → `VerdictPanel`.
8. `VerdictCards` — per-cluster: PROVEN/UNCERTAIN/CONTRADICTED badge,
   rules fired list with mechanistic_basis + PubMed link, CASSIA comparison.
9. `InspectionPanel` — shown after upload before Run button. Displays detected
   input type, modalities, preprocessing steps that will run, stream count.

**Day 15:**
10. `RuleFireHeatmap` — Phase 1 Figure 2 reproduced interactively.
    Rows: clusters. Cols: KG rules (paginated 30 at a time).
    Colors: PASS=green, FAIL=red, NOT_TESTABLE=grey. Click cell → tooltip with PMID link.
11. `KGBrowser` — searchable rule table with inline editor.
    Add rule form enforces: pmid must be all-digits ≥7 chars, mechanistic_basis ≥20 chars.
    Shows "PENDING_REVIEW — verify PMID via PubMed before activating."
12. `ComparisonTable` — CASSIA / mLLMCelltype / AXIOM accuracy table with Phase 1 baseline.

### App.tsx layout

```tsx
// Two views: pre-results (upload + pipeline) and post-results (tabs)
// Header: <Logo /> AXIOM-SC | <SDKVersionBadge /> | <ProfileSelector />
// Pre-results: <FileUpload /> → <InspectionPanel /> → <TierSelector /> → <ProgressStream />
// Post-results tabs: UMAP+Verdicts | Rule Firing | KG Browser | Comparison
```

---

## 24. Pre-Day-11 patch phase — KG + Tier 3 resources

### Status (completed — do not re-run patch scripts)

- **KG: 495 ACTIVE rules / 153 cell types** (scripts: patch_a_kg_expansion.py,
  patch_a_supplement.py — do NOT re-run on existing KG, will duplicate)
- **Tier 3 resources**: 8 external databases auto-load from `~/.axiom_sc/resources/`
  when present. CLI: `axiom-sc download-resources`
- **Tests: 379 passing / 84.51% coverage**
- New tissue classes: CNS, skin, kidney, heart, stromal/EC, endocrine,
  hematopoietic, cancer states (TEX_progenitor/TEX_terminal), trophoblast

### Resource auto-loading (Patch B — already implemented)

| Module | When resource present | Fallback (resource absent) |
|---|---|---|
| `communication.py` | CellPhoneDB 2,300+ L-R genes | 40 hardcoded L-R pairs |
| `cross_species.py` | ENSEMBL genome-wide conservation | ~300 hardcoded conserved genes |
| `spatial_niche.py` | CELLama reference embeddings (~50MB) | kNN composition mode |
| `scenic_pipeline.py` | Rankings feather (1.3GB) + motif DB | Adjacency-only fallback |

### RNA-only auto-detection (already in convergence.py)

`run_and_converge(auto_detect=True)` → 4/6 streams on RNA-only data
(communication + cross_species + spatial_niche + SCENIC).
velocity and chromatin activate automatically when layers present in adata.

### Gate check (run before Day 11)

```bash
/Users/sivamoturi/miniconda/bin/python3.10 -m pytest tests/ -v \
  --cov=axiom_sc --cov-config=pyproject.toml 2>&1 | tail -5
# Expected: ≥379 passing, 0 failing, ≥84% coverage

/Users/sivamoturi/miniconda/bin/python3.10 -c "
from axiom_sc.kg.loader import KGLoader
kg = KGLoader().load('kg_data/oracle_kg_v0.2.0.json')
rules = [r for r in kg if r.get('status') == 'ACTIVE']
types = {r['cell_type'] for r in rules}
missing = [r['rule_id'] for r in rules if not r.get('pmid') or r['pmid']=='NEEDS_REVIEW']
print(f'ACTIVE rules: {len(rules)}  (need ≥450)')
print(f'Cell types:   {len(types)}  (need ≥150)')
print(f'Unverified PMIDs: {len(missing)}  (need 0)')
"
```

---

## 25. KG universe extension — day-by-day Claude Code execution

**When:** Run in parallel with playground Days 14–19 (KG work is CPU-only,
independent of the API/React build). Can be interleaved: KG session in morning,
playground session in afternoon. Or as dedicated focused days after Day 20.

**Duration:** ~4 focused days to reach ≥640 rules / ≥198 cell types.

**PMID discipline (mandatory — do not skip):**
Every PMID must be found and verified via 3-step Entrez before writing any rule:
```python
from Bio import Entrez
Entrez.email = "oss@receptor.bio"
# Step 1: search
handle = Entrez.esearch(db="pubmed", term="YOUR SPECIFIC QUERY")
record = Entrez.read(handle)
# Step 2: fetch and read title
handle2 = Entrez.efetch(db="pubmed", id=record["IdList"][0], rettype="medline")
print(handle2.read()[:500])
# Step 3: confirm title directly supports the gene/cell-type claim → then write rule
```
One anchor paper can cover multiple rules for the same cell type.
After every batch: run `axiom-kg check-pmids` → must return empty list.

---

### KG Extension Day A — Retina (10 types, ~30 rules)

**Claude Code instructions:**

Step 1: Seed candidates from CellMarker 2.0
```python
from axiom_sc.kg.seeder import seed_from_cellmarker2
seed_from_cellmarker2(
    output_path="kg_data/candidates/ext_retina.json",
    tissues=["Retina"]
)
```

Step 2: Find and verify anchor PMIDs FIRST before writing any rule
```python
# Search these queries via Entrez.esearch, verify titles, record PMIDs:
queries = [
    "human retina single cell transcriptomics atlas photoreceptor",  # Menon 2019 or Lukowski 2019
    "retinal ganglion cell ATOH7 BRN3 single cell",
    "Muller glia RLBP1 VIM single cell retina",
    "NRL rod photoreceptor transcription factor development",
    "cone photoreceptor RXRG ARR3 single cell",
    "RPE RPE65 BEST1 single cell retinal pigment epithelium",
]
```

Step 3: Write rules directly into `kg_data/oracle_kg_v0.2.0.json`
using `axiom-kg add` CLI or direct JSON append. Rules for this batch:

| rule_id | cell_type | rule_type | gene_or_regulon | direction | critical_note |
|---|---|---|---|---|---|
| ROD_POS_001 | Rod_photoreceptor | positive | ["NRL","RHO","CNGB1"] | high | NRL regulon must be active |
| ROD_NEG_001 | Rod_photoreceptor | **negative** | ["ARR3"] | absent | ARR3 = cone arrestin, defines mutual exclusion |
| CONE_POS_001 | Cone_photoreceptor | positive | ["RXRG","ARR3","OPN1SW"] | high | |
| CONE_NEG_001 | Cone_photoreceptor | **negative** | ["NRL"] | absent | NRL = rod master TF — absent in all cone types |
| CONE_CIRCUIT_001 | Cone_photoreceptor | circuit | ["RXRG","ARR3"] | active+high | co-required |
| MUL_POS_001 | Muller_glia | positive | ["RLBP1","VIM","SLC1A3","GLUL"] | high | GLUL = glutamine synthetase |
| RGC_CIRCUIT_001 | Retinal_ganglion_cell | circuit | ["ATOH7","SNCG","RBPMS"] | active+high | ATOH7 as regulon |
| RPE_POS_001 | RPE | positive | ["RPE65","BEST1","RLBP1"] | high | |
| RPE_NEG_001 | RPE | **negative** | ["RHO"] | absent | RPE shares RLBP1 with Muller glia; RHO distinguishes |
| HORIZ_CIRCUIT_001 | Horizontal_cell | circuit | ["LHX1","ONECUT1","CALB1"] | active+high | LHX1 as regulon |
| AMACRINE_POS_001 | Amacrine_cell | positive | ["GAD1","TFAP2A"] | high | inhibitory subtype |
| BIPOLAR_CIRCUIT_001 | Bipolar_cell | circuit | ["VSX2"] | active | VSX2 regulon = bipolar identity |
| CIL_EPI_POS_001 | Ciliary_epithelium | positive | ["AQP1","PLCL1"] | high | |
| ENDO_RET_POS_001 | Endothelium_retinal | positive | ["CDH5","NRP1"] | high | |

Step 4: Write tests for this batch
```python
# Add to tests/test_kg_extension_a.py
def test_rod_cone_mutual_exclusion():
    # ARR3 high → Rod CONTRADICTED; NRL high → Cone CONTRADICTED
    ev_arr3 = make_evidence({"ARR3": 2.5})
    assert annotator.annotate_candidate("Rod_photoreceptor", ev_arr3).verdict == "CONTRADICTED"
    ev_nrl  = make_evidence({"NRL": 2.5})
    assert annotator.annotate_candidate("Cone_photoreceptor", ev_nrl).verdict == "CONTRADICTED"

def test_rpe_muller_distinction():
    # Both express RLBP1 — RHO distinguishes: present=rod/not RPE
    ev_rpe = make_evidence({"RPE65": 2.0, "RLBP1": 2.0, "RHO": 0.05})
    assert annotator.annotate_candidate("RPE", ev_rpe).verdict in ("PROVEN","UNCERTAIN")
    ev_rod_contamination = make_evidence({"RPE65": 2.0, "RLBP1": 2.0, "RHO": 2.5})
    assert annotator.annotate_candidate("RPE", ev_rod_contamination).verdict == "CONTRADICTED"

def test_retina_batch_pmids():
    rules = kg_loader.get_rules_by_prefix("ROD_,CONE_,MUL_,RGC_,RPE_,HORIZ_,AMACRINE_,BIPOLAR_")
    missing = [r["rule_id"] for r in rules if not r.get("pmid") or r["pmid"] == "NEEDS_REVIEW"]
    assert missing == [], f"Unverified PMIDs: {missing}"
```

Step 5: Run tests + KG integrity check
```bash
/Users/sivamoturi/miniconda/bin/python3.10 -m pytest tests/test_kg_extension_a.py -v
axiom-kg check-pmids   # must return empty
```

**End of Day A:** ≥14 new rules, 10 new cell types, all PMIDs verified.

---

### KG Extension Day B — Bone/Cartilage + Gonads (13 types, ~45 rules)

**Claude Code instructions:**

Step 1: Verify anchor PMIDs for both tissue classes
```python
queries = [
    "osteoblast osteoclast single cell scRNA-seq bone marrow",       # Baryawno 2019 or Zhong 2020
    "RUNX2 SP7 osterix osteoblast transcription factor",
    "cathepsin K CTSK osteoclast marker",
    "chondrocyte SOX9 COL2A1 single cell cartilage",
    "human testis single cell spermatogenesis atlas",                # Wang 2018 Cell Research
    "Sertoli cell SOX9 AMH single cell testis",
    "granulosa cell FSHR CYP19A1 aromatase ovary single cell",
    "theca cell CYP17A1 steroidogenesis ovary",
]
```

Step 2: Write bone/cartilage rules

| rule_id | cell_type | rule_type | gene_or_regulon | direction |
|---|---|---|---|---|
| OSTBL_CIRCUIT_001 | Osteoblast | circuit | ["RUNX2","SP7","COL1A1"] | active+high |
| OSTBL_NEG_001 | Osteoblast | **negative** | ["CTSK"] | absent |
| OSTCL_POS_001 | Osteoclast | positive | ["ACP5","CTSK","MMP9"] | high |
| OSTCL_NEG_001 | Osteoclast | **negative** | ["COL1A1","RUNX2"] | absent |
| OSTCYTE_POS_001 | Osteocyte | positive | ["PHEX","DMP1","SOST"] | high |
| CHOND_CIRCUIT_001 | Chondrocyte | circuit | ["SOX9","COL2A1","ACAN"] | active+high |
| CHOND_NEG_001 | Chondrocyte | **negative** | ["COL1A1"] | absent |
| HYPCHOND_POS_001 | Hypertrophic_chondrocyte | positive | ["COL10A1","MMP13"] | high |
| HYPCHOND_NEG_001 | Hypertrophic_chondrocyte | **negative** | ["SOX9"] | inactive |

Step 3: Write gonad rules

| rule_id | cell_type | rule_type | gene_or_regulon | direction |
|---|---|---|---|---|
| SSC_POS_001 | Spermatogonial_SC | positive | ["ID4","GFRA1","ZBTB16"] | high |
| SPCYTE_POS_001 | Spermatocyte | positive | ["SYCP3","SYCP1","MLH1"] | high |
| SPTID_POS_001 | Spermatid | positive | ["PRM1","TNP1","ACR"] | high |
| SERTOLI_CIRCUIT_001 | Sertoli_cell | circuit | ["SOX9","AMH","CLDN11","WT1"] | active+high |
| LEYDIG_POS_001 | Leydig_cell | positive | ["CYP17A1","STAR","HSD3B2"] | high |
| GRANUL_CIRCUIT_001 | Granulosa_cell | circuit | ["FOXL2","FSHR","CYP19A1"] | active+high |
| GRANUL_NEG_001 | Granulosa_cell | **negative** | ["CYP17A1"] | absent |
| THECA_POS_001 | Theca_cell | positive | ["CYP17A1","LHR","STAR"] | high |
| THECA_NEG_001 | Theca_cell | **negative** | ["CYP19A1"] | absent |
| OOCYTE_POS_001 | Oocyte | positive | ["ZP1","ZP2","GDF9","BMP15"] | high |
| SERTOLI_NEG_001 | Sertoli_cell | **negative** | ["CYP17A1"] | absent |

Step 4: Mandatory mutual-exclusion tests
```python
# tests/test_kg_extension_b.py
def test_osteoblast_osteoclast_mutual_exclusion():
    ev_ctsk = make_evidence({"CTSK": 2.8, "RUNX2": 2.0, "COL1A1": 1.5})
    assert annotator.annotate_candidate("Osteoblast", ev_ctsk).verdict == "CONTRADICTED"
    ev_col1a1 = make_evidence({"ACP5": 2.0, "CTSK": 2.5, "COL1A1": 2.5})
    assert annotator.annotate_candidate("Osteoclast", ev_col1a1).verdict == "CONTRADICTED"

def test_granulosa_theca_mutual_exclusion():
    # Granulosa: CYP17A1 must be absent
    ev_cyp17 = make_evidence({"FSHR": 2.0, "CYP19A1": 1.8, "CYP17A1": 2.5})
    assert annotator.annotate_candidate("Granulosa_cell", ev_cyp17).verdict == "CONTRADICTED"
    # Theca: CYP19A1 must be absent
    ev_cyp19 = make_evidence({"CYP17A1": 2.5, "STAR": 2.0, "CYP19A1": 2.5})
    assert annotator.annotate_candidate("Theca_cell", ev_cyp19).verdict == "CONTRADICTED"

def test_chondrocyte_not_fibroblast():
    # COL1A1 high → Chondrocyte CONTRADICTED (fibrocartilage contamination)
    ev_col1 = make_evidence({"SOX9_regulon": 2.0, "COL2A1": 2.0, "COL1A1": 2.5})
    assert annotator.annotate_candidate("Chondrocyte", ev_col1).verdict == "CONTRADICTED"
```

**End of Day B:** ≥20 new rules, 13 new cell types.

---

### KG Extension Day C — Inner Ear + Cornea + Salivary + Synovium (16 types, ~49 rules)

**Claude Code instructions:**

Step 1: Verify anchor PMIDs for all four tissue classes
```python
queries = [
    "inner ear hair cell single cell transcriptomics cochlea",       # Kolla 2020 or Li 2018
    "ATOH1 hair cell regeneration transcription factor",
    "prestin SLC26A5 outer hair cell electromotility",
    "corneal epithelium KRT12 KRT3 single cell scRNA-seq",
    "corneal endothelium SLC4A11 ZEB1 single cell",
    "salivary gland acinar ductal serous mucous single cell",
    "synovial fibroblast PRG4 lubricin single cell scRNA-seq",
]
```

Step 2: Write inner ear rules

| rule_id | cell_type | rule_type | gene_or_regulon | direction |
|---|---|---|---|---|
| IHC_CIRCUIT_001 | Inner_Hair_cell | circuit | ["ATOH1","MYO7A","ESPN"] | active+high |
| OHC_POS_001 | Outer_Hair_cell | positive | ["SLC26A5","OCM","ONCOMODULIN"] | high |
| OHC_NEG_001 | Outer_Hair_cell | **negative** | ["ATOH1"] | inactive |
| SUPCELL_POS_001 | Supporting_cell | positive | ["SOX2","PROX1"] | high |
| SUPCELL_NEG_001 | Supporting_cell | **negative** | ["MYO7A"] | absent |
| SGN_POS_001 | Spiral_ganglion_neuron | positive | ["NEFM","RBFOX3","SLC17A6"] | high |
| STRIA_POS_001 | Stria_vascularis | positive | ["KCNJ10","SLC12A2"] | high |

Step 3: Write cornea rules

| rule_id | cell_type | rule_type | gene_or_regulon | direction |
|---|---|---|---|---|
| CORNEP_POS_001 | Corneal_epithelium | positive | ["KRT12","KRT3"] | high |
| CORNSTRO_POS_001 | Corneal_keratocyte | positive | ["KERA","ALDH3A1"] | high |
| CORNSTRO_NEG_001 | Corneal_keratocyte | **negative** | ["ACTA2"] | absent |
| CORNENDO_POS_001 | Corneal_endothelium | positive | ["SLC4A11","ZEB1","COL8A2"] | high |
| TM_POS_001 | Trabecular_meshwork | positive | ["AQP1","MGP","MYOC"] | high |

Step 4: Write salivary gland + synovium rules

| rule_id | cell_type | rule_type | gene_or_regulon | direction |
|---|---|---|---|---|
| ACINAR_SER_POS_001 | Acinar_serous | positive | ["BPIFB2","ZG16B","PRSS2"] | high |
| ACINAR_MUC_POS_001 | Acinar_mucous | positive | ["MUC5B","CLCA1"] | high |
| ACINAR_MUC_NEG_001 | Acinar_mucous | **negative** | ["PRSS2"] | absent |
| DUCT_STRI_POS_001 | Striated_duct | positive | ["AQP5","KRT14"] | high |
| MYOEPI_POS_001 | Myoepithelial | positive | ["ACTA2","KRT14","MYLK"] | high |
| SYNFIB_POS_001 | Synovial_fibroblast | positive | ["PRG4","THY1","PDPN"] | high |
| SYNMAC_POS_001 | Synovial_macrophage | positive | ["CD68","FOLR2","LYVE1"] | high |
| CHONDPROG_POS_001 | Chondroprogenitor | positive | ["PDGFRA","NES"] | high |
| CHONDPROG_NEG_001 | Chondroprogenitor | **negative** | ["MBP","COL2A1"] | absent |

Step 5: Batch tests
```python
# tests/test_kg_extension_c.py
def test_inner_outer_hair_cell_mutual_exclusion():
    # Outer hair cell: ATOH1 regulon must be INACTIVE (ATOH1 drives IHC not OHC fate)
    ev_atoh1 = make_evidence({"ATOH1_regulon": 2.5, "SLC26A5": 2.0})
    assert annotator.annotate_candidate("Outer_Hair_cell", ev_atoh1).verdict == "CONTRADICTED"

def test_supporting_not_hair_cell():
    # MYO7A absent in supporting cells
    ev_myo7a = make_evidence({"SOX2": 2.0, "PROX1": 1.5, "MYO7A": 2.5})
    assert annotator.annotate_candidate("Supporting_cell", ev_myo7a).verdict == "CONTRADICTED"

def test_acinar_serous_mucous_distinction():
    # Mucous: PRSS2 absent (serous digestive enzyme)
    ev_prss2 = make_evidence({"MUC5B": 2.0, "CLCA1": 1.8, "PRSS2": 2.5})
    assert annotator.annotate_candidate("Acinar_mucous", ev_prss2).verdict == "CONTRADICTED"
```

**End of Day C:** ≥33 new rules, 16 new cell types.

---

### KG Extension Day D — Additional vascular subtypes (8 types, ~25 rules) + validation sprint

**Claude Code instructions:**

Step 1: Verify anchor PMIDs
```python
queries = [
    "arterial venous endothelial single cell scRNA-seq NOTCH4 APLNR",  # Trimm 2021 or similar
    "high endothelial venule ACKR1 GLYCAM1 single cell lymphocyte trafficking",
    "lymphatic endothelial PROX1 LYVE1 single cell",
    "tip cell ESM1 ANGPT2 DLL4 single cell angiogenesis",
    "pericyte brain blood-brain-barrier RGS5 KCNJ8 single cell",
    "synthetic smooth muscle cell COL1A1 vascular remodeling",
    "sinusoidal endothelial liver STAB1 STAB2 CLEC4G single cell",
    "lung capillary aerocyte CA4 EDNRB single cell",
]
```

Step 2: Write vascular rules

| rule_id | cell_type | rule_type | gene_or_regulon | direction |
|---|---|---|---|---|
| ART_EC_POS_001 | Arterial_Endothelial | positive | ["NOTCH4","GJA4","CXCL12"] | high |
| ART_EC_NEG_001 | Arterial_Endothelial | **negative** | ["APLNR"] | absent |
| VEN_EC_POS_001 | Venous_Endothelial | positive | ["APLNR","NR2F2"] | high |
| VEN_EC_NEG_001 | Venous_Endothelial | **negative** | ["GJA4"] | absent |
| HEV_POS_001 | High_endothelial_venule | positive | ["ACKR1","GLYCAM1","CCL21"] | high |
| LYMPH_EC_CIRCUIT_001 | Lymphatic_Endothelial | circuit | ["PROX1","LYVE1","PDPN"] | active+high |
| TIP_CIRCUIT_001 | Tip_cell | circuit | ["ESM1","DLL4","KDR"] | high |
| TIP_NEG_001 | Tip_cell | **negative** | ["NR2F2"] | absent |
| PERI_BRAIN_POS_001 | Pericyte_brain | positive | ["PDGFRB","RGS5","KCNJ8","SLC7A5"] | high |
| SYNTH_SMC_POS_001 | Synthetic_SMC | positive | ["COL1A1","FN1"] | high |
| SYNTH_SMC_NEG_001 | Synthetic_SMC | **negative** | ["MYH11"] | absent |
| LSEC_POS_001 | Sinusoidal_EC | positive | ["STAB1","STAB2","CLEC4G","LYVE1"] | high |
| LUNG_CAP_POS_001 | Lung_capillary | positive | ["CA4","EDNRB"] | high |

Step 3: Final KG validation sprint
```bash
# 1. Full PMID audit
/Users/sivamoturi/miniconda/bin/python3.10 -c "
from axiom_sc.kg.loader import KGLoader
kg = KGLoader().load('kg_data/oracle_kg_v0.2.0.json')
rules = [r for r in kg if r.get('status') == 'ACTIVE']
types = {r['cell_type'] for r in rules}
missing_pmid = [r['rule_id'] for r in rules
                if not r.get('pmid') or r['pmid'] in ('NEEDS_REVIEW','')]
print(f'ACTIVE rules : {len(rules)}')
print(f'Cell types   : {len(types)}')
print(f'Missing PMID : {len(missing_pmid)}')
assert len(missing_pmid) == 0, f'Fix these: {missing_pmid}'
print('✓ All PMIDs verified')
"

# 2. Full test suite
/Users/sivamoturi/miniconda/bin/python3.10 -m pytest tests/ -v \
  --cov=axiom_sc --cov-config=pyproject.toml 2>&1 | tail -5
# Expected: ≥415 passing, 0 failing, ≥84% coverage

# 3. Regenerate REFERENCES.md
axiom-kg references --output REFERENCES.md
echo "REFERENCES.md now has $(grep -c 'doi\|PMID\|pubmed' REFERENCES.md) citations"

# 4. Final KG stats
/Users/sivamoturi/miniconda/bin/python3.10 -c "
from axiom_sc.kg.loader import KGLoader
kg = KGLoader().load('kg_data/oracle_kg_v0.2.0.json')
rules = [r for r in kg if r.get('status') == 'ACTIVE']
types = {r['cell_type'] for r in rules}
neg_rules = [r for r in rules if r['rule_type'] == 'negative']
circuit_rules = [r for r in rules if r['rule_type'] == 'circuit']
print(f'ACTIVE rules: {len(rules)} target ≥640')
print(f'Cell types:   {len(types)} target ≥198')
print(f'Negative:     {len(neg_rules)} (these are the core AXIOM innovation)')
print(f'Circuit:      {len(circuit_rules)}')
print(f'Positive:     {len(rules)-len(neg_rules)-len(circuit_rules)}')
"
```

Step 4: Tests for vascular extension
```python
# tests/test_kg_extension_d.py
def test_arterial_venous_mutual_exclusion():
    # APLNR high → Arterial CONTRADICTED; GJA4 high → Venous CONTRADICTED
    ev_aplnr = make_evidence({"NOTCH4": 2.0, "APLNR": 2.5})
    assert annotator.annotate_candidate("Arterial_Endothelial", ev_aplnr).verdict == "CONTRADICTED"
    ev_gja4 = make_evidence({"APLNR": 2.0, "NR2F2": 1.8, "GJA4": 2.5})
    assert annotator.annotate_candidate("Venous_Endothelial", ev_gja4).verdict == "CONTRADICTED"

def test_synthetic_smc_not_contractile():
    # MYH11 high → Synthetic_SMC CONTRADICTED (has lost contractile program)
    ev_myh11 = make_evidence({"COL1A1": 2.0, "FN1": 1.8, "MYH11": 2.5})
    assert annotator.annotate_candidate("Synthetic_SMC", ev_myh11).verdict == "CONTRADICTED"

def test_final_kg_scale():
    rules = kg_loader.list_rules(status="ACTIVE")
    types = {r["cell_type"] for r in rules}
    assert len(rules) >= 640, f"Only {len(rules)} rules — need ≥640"
    assert len(types) >= 198, f"Only {len(types)} types — need ≥198"
    neg = [r for r in rules if r["rule_type"] == "negative"]
    assert len(neg) >= 120, "Need ≥120 negative rules (core AXIOM innovation)"
```

**End of Day D:** ≥640 rules, ≥198 cell types, all PMIDs verified, full test suite green.

---

### Scheduling within the Day 11–20 plan

These KG extension days can be slotted as follows (choose based on energy/pace):

| Option | When | Notes |
|---|---|---|
| **Parallel** (recommended) | Morning KG / Afternoon playground | KG is CPU-only, independent of playground build |
| **Dedicated block** | Days 21–24 after playground release | After Day 20 SDK release, KG extension as v0.2.1 prep |
| **Interleaved** | KG Day A on Day 14, B on Day 16, C on Day 18, D on Day 19 | Lighter playground days (React polish) free up focus |

Each KG day produces a commit to `kg_data/oracle_kg_v0.2.0.json` and
`tests/test_kg_extension_{a,b,c,d}.py`. REFERENCES.md regenerated after Day D.

### KG ceiling

~200 types is the mechanistic ceiling. Beyond this, cell types lack sufficient
primary literature for non-trivial negative/circuit rules. Tier 4 LLM handles
novel/rare types beyond this ceiling dynamically via CellMarker 2.0 RAG retrieval.

---

## 26. Task queue: Dramatiq (full spec)

All code in Section 22 uses Dramatiq. Reference spec only:

```bash
# Worker launch (Dramatiq — replaces Celery)
dramatiq tasks --processes 1 --threads 2
# --processes 1: annotation jobs are memory-heavy (h5ad can be 2GB)
# --threads 2:   one active job + one waiting for GPU result

# requirements.txt: dramatiq[redis]>=1.16
# NOT: celery>=5.3
```

No retry on annotation jobs (`max_retries=0`) — bioinformatics pipelines are not
idempotent. If a job fails, the user re-submits with corrected data.

---

## 27. AWS ECS Fargate deployment + GitHub Actions OIDC CI/CD

### Secrets: AWS Secrets Manager (not SSM Parameter Store)

Only confidential values go in Secrets Manager. Non-sensitive config uses ECS
environment variables directly. AXIOM_SDK_VERSION is a Docker build arg baked
into the image — it does NOT go in Secrets Manager or SSM.

```bash
# Create secrets once:
aws secretsmanager create-secret \
  --name axiom-sc/prod/anthropic-api-key \
  --secret-string "sk-ant-YOUR_KEY"

aws secretsmanager create-secret \
  --name axiom-sc/prod/openai-api-key \
  --secret-string "sk-proj-YOUR_KEY"

# GPU dispatch secrets (when RunPod/Modal in use):
aws secretsmanager create-secret \
  --name axiom-sc/prod/runpod-api-key \
  --secret-string "rp_YOUR_KEY"

aws secretsmanager create-secret \
  --name axiom-sc/prod/runpod-scenic-endpoint-id \
  --secret-string "your_endpoint_id"
```

### AWS infrastructure (one-time setup)

```bash
# ECR — one repo per service
aws ecr create-repository --repository-name axiom-sc/api
aws ecr create-repository --repository-name axiom-sc/worker
aws ecr create-repository --repository-name axiom-sc/app

# ECS cluster
aws ecs create-cluster --cluster-name axiom-sc-cluster

# ElastiCache Redis (serverless)
aws elasticache create-serverless-cache \
  --serverless-cache-name axiom-redis --engine redis

# EFS for shared uploads (api + worker mount the same volume)
aws efs create-file-system --tags Key=Name,Value=axiom-uploads

# S3 for loom/AUC exchange with RunPod/Modal GPU workers
aws s3 mb s3://axiom-sc-jobs
aws s3api put-bucket-lifecycle-configuration \
  --bucket axiom-sc-jobs \
  --lifecycle-configuration '{"Rules":[{"ID":"ttl","Status":"Enabled",
    "Filter":{"Prefix":"jobs/"},"Expiration":{"Days":1}}]}'
```

### IAM OIDC role for GitHub Actions (no stored secrets)

```json
{
  "Statement": [{
    "Effect": "Allow",
    "Principal": {
      "Federated": "arn:aws:iam::ACCOUNT_ID:oidc-provider/token.actions.githubusercontent.com"
    },
    "Action": "sts:AssumeRoleWithWebIdentity",
    "Condition": {
      "StringEquals": {"token.actions.githubusercontent.com:aud": "sts.amazonaws.com"},
      "StringLike":   {"token.actions.githubusercontent.com:sub":
                        "repo:receptor-bio/axiom-playground:*"}
    }
  }]
}
```

Attach policy: `ecr:*`, `ecs:RegisterTaskDefinition`, `ecs:UpdateService`,
`ecs:DescribeServices`, `iam:PassRole` for the ECS task role.

### ECS task definition — corrected secrets block

```json
"secrets": [
  {"name": "ANTHROPIC_API_KEY",
   "valueFrom": "arn:aws:secretsmanager:REGION:ACCOUNT_ID:secret:axiom-sc/prod/anthropic-api-key"},
  {"name": "OPENAI_API_KEY",
   "valueFrom": "arn:aws:secretsmanager:REGION:ACCOUNT_ID:secret:axiom-sc/prod/openai-api-key"},
  {"name": "RUNPOD_API_KEY",
   "valueFrom": "arn:aws:secretsmanager:REGION:ACCOUNT_ID:secret:axiom-sc/prod/runpod-api-key"},
  {"name": "RUNPOD_SCENIC_ENDPOINT_ID",
   "valueFrom": "arn:aws:secretsmanager:REGION:ACCOUNT_ID:secret:axiom-sc/prod/runpod-scenic-endpoint-id"}
],
"environment": [
  {"name": "REDIS_URL",     "value": "redis://ELASTICACHE_ENDPOINT:6379"},
  {"name": "AXIOM_PROFILE", "value": "oss-apache"},
  {"name": "AXIOM_S3_BUCKET","value": "axiom-sc-jobs"}
]
```

Note: NO AXIOM_SDK_VERSION at runtime — it is baked into the Docker image
at build time via `--build-arg AXIOM_SDK_VERSION=$(cat config.yaml | ...)`.

### ECS task role IAM policy

```json
{"Statement": [
  {"Effect":"Allow","Action":["secretsmanager:GetSecretValue","secretsmanager:DescribeSecret"],
   "Resource":"arn:aws:iam::ACCOUNT_ID:secret:axiom-sc/prod/*"},
  {"Effect":"Allow","Action":["s3:PutObject","s3:GetObject","s3:DeleteObject"],
   "Resource":"arn:aws:s3:::axiom-sc-jobs/jobs/*"},
  {"Effect":"Allow","Action":["elasticfilesystem:ClientMount","elasticfilesystem:ClientWrite"],
   "Resource":"arn:aws:elasticfilesystem:REGION:ACCOUNT_ID:file-system/EFS_ID"},
  {"Effect":"Allow","Action":["logs:CreateLogStream","logs:PutLogEvents"],
   "Resource":"arn:aws:logs:REGION:ACCOUNT_ID:log-group:/axiom-sc/*"}
]}
```

### Dockerfiles

**api/Dockerfile:**
```dockerfile
FROM python:3.10-slim
WORKDIR /app
RUN apt-get update && apt-get install -y libhdf5-dev curl && rm -rf /var/lib/apt/lists/*
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
ARG AXIOM_SDK_VERSION=0.2.0
RUN pip install --no-cache-dir axiom-sc==${AXIOM_SDK_VERSION}
COPY . .
EXPOSE 8000
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000", "--timeout-keep-alive", "120"]
```

**api/Dockerfile.worker:**
```dockerfile
FROM python:3.10-slim
WORKDIR /app
RUN apt-get update && apt-get install -y libhdf5-dev && rm -rf /var/lib/apt/lists/*
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
ARG AXIOM_SDK_VERSION=0.2.0
RUN pip install --no-cache-dir axiom-sc==${AXIOM_SDK_VERSION}
COPY . .
CMD ["dramatiq", "tasks", "--processes", "1", "--threads", "2"]
```

**app/Dockerfile (multi-stage):**
```dockerfile
FROM node:20-alpine AS builder
WORKDIR /app
COPY package*.json .
RUN npm ci
COPY . .
ARG VITE_API_URL
ENV VITE_API_URL=${VITE_API_URL}
RUN npm run build

FROM nginx:alpine
COPY --from=builder /app/dist /usr/share/nginx/html
COPY nginx.conf /etc/nginx/conf.d/default.conf
EXPOSE 80
CMD ["nginx", "-g", "daemon off;"]
```

**app/nginx.conf:**
```nginx
server {
    listen 80;
    root /usr/share/nginx/html;
    index index.html;
    location / { try_files $uri $uri/ /index.html; }
    location ~* \.(js|css|png|svg|woff2)$ { expires 1y; add_header Cache-Control "public, immutable"; }
}
```

### GitHub Actions CI/CD (.github/workflows/deploy.yml)

```yaml
name: Deploy to AWS ECS
on:
  push:
    branches: [main, develop]

permissions:
  id-token: write   # OIDC — no stored AWS keys needed
  contents: read

env:
  AWS_REGION: us-east-1
  ECR_REGISTRY: ${{ secrets.AWS_ACCOUNT_ID }}.dkr.ecr.us-east-1.amazonaws.com
  ECS_CLUSTER: axiom-sc-cluster

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with: {python-version: "3.10"}
      - run: |
          SDK=$(grep sdk_version axiom-playground/config.yaml | awk '{print $2}')
          pip install axiom-sc==${SDK} 2>/dev/null || pip install -e axiom-sc/
          cd axiom-playground/api && pip install -r requirements.txt pytest httpx
          pytest tests/ -v --tb=short
      - uses: actions/setup-node@v4
        with: {node-version: "20"}
      - run: cd axiom-playground/app && npm ci && npm run type-check

  build-and-deploy:
    needs: test
    runs-on: ubuntu-latest
    strategy:
      matrix:
        include:
          - service: api
            ecr_repo: axiom-sc/api
            dockerfile: api/Dockerfile
            context: axiom-playground/api
            ecs_service: axiom-api-service
            task_def: .aws/task-def-api.json
            container: api
          - service: worker
            ecr_repo: axiom-sc/worker
            dockerfile: api/Dockerfile.worker
            context: axiom-playground/api
            ecs_service: axiom-worker-service
            task_def: .aws/task-def-worker.json
            container: worker
          - service: app
            ecr_repo: axiom-sc/app
            dockerfile: app/Dockerfile
            context: axiom-playground/app
            ecs_service: axiom-app-service
            task_def: .aws/task-def-app.json
            container: app
    steps:
      - uses: actions/checkout@v4

      - name: Configure AWS via OIDC (no secrets in GitHub)
        uses: aws-actions/configure-aws-credentials@v4
        with:
          role-to-assume: arn:aws:iam::${{ secrets.AWS_ACCOUNT_ID }}:role/axiom-github-deploy-role
          aws-region: ${{ env.AWS_REGION }}

      - uses: aws-actions/amazon-ecr-login@v2

      - name: Get SDK version
        id: sdk
        run: echo "version=$(grep sdk_version axiom-playground/config.yaml | awk '{print $2}')" >> $GITHUB_OUTPUT

      - name: Build and push image
        uses: docker/build-push-action@v5
        with:
          context: ${{ matrix.context }}
          file: axiom-playground/${{ matrix.dockerfile }}
          push: true
          build-args: |
            AXIOM_SDK_VERSION=${{ steps.sdk.outputs.version }}
            VITE_API_URL=https://${{ secrets.ALB_DNS_NAME }}
          tags: |
            ${{ env.ECR_REGISTRY }}/${{ matrix.ecr_repo }}:${{ github.sha }}
            ${{ env.ECR_REGISTRY }}/${{ matrix.ecr_repo }}:latest
          cache-from: type=registry,ref=${{ env.ECR_REGISTRY }}/${{ matrix.ecr_repo }}:latest
          cache-to: type=inline

      - name: Update ECS task definition
        id: task-def
        uses: aws-actions/amazon-ecs-render-task-definition@v1
        with:
          task-definition: ${{ matrix.task_def }}
          container-name: ${{ matrix.container }}
          image: ${{ env.ECR_REGISTRY }}/${{ matrix.ecr_repo }}:${{ github.sha }}

      - name: Deploy to ECS (rolling, zero-downtime)
        uses: aws-actions/amazon-ecs-deploy-task-definition@v1
        with:
          task-definition: ${{ steps.task-def.outputs.task-definition }}
          service: ${{ matrix.ecs_service }}
          cluster: ${{ env.ECS_CLUSTER }}
          wait-for-service-stability: true

environments:
  staging:    # develop branch → no approval gate
  production: # main branch → requires 1 reviewer approval
```

### GitHub Secrets required (set once, no AWS credentials)

```
AWS_ACCOUNT_ID   → 12-digit account ID (non-sensitive — used for ECR URL construction)
ALB_DNS_NAME     → your ALB DNS name for VITE_API_URL
```

### SDK version promotion (one command)

```bash
# Bump axiom-sc in config.yaml → triggers full rebuild + deploy automatically
echo "sdk_version: 0.2.1" > axiom-playground/config.yaml
git add config.yaml && git commit -m "chore: bump sdk to 0.2.1"
git push origin main
# → CI/CD rebuilds all 3 images with axiom-sc==0.2.1 from PyPI
# → ECS rolling update, zero downtime
# → /health returns sdk_version: 0.2.1
```

---

## 28. GPU task isolation — RunPod / Modal serverless

### Which tasks need GPU

| Task | GPU needed | Where runs |
|---|---|---|
| All preprocessing, KG checks, Tier 3-5 | No | ECS Fargate CPU |
| **pySCENIC (GRNBoost2 + AUCell)** | **YES — major** | RunPod/Modal RTX 4090 |
| **scVelo dynamical mode** | GPU-beneficial | RunPod/Modal (optional; falls back to CPU) |

### Auto-dispatch in ScenicRunner

```python
# axiom_sc/tier2/scenic_runner.py  — dispatch priority
def run(self, loom_path, config, output_dir):
    if os.getenv("MODAL_TOKEN_ID"):
        return self._run_modal(loom_path, config, output_dir)
    elif os.getenv("RUNPOD_API_KEY"):
        return self._run_runpod(loom_path, config, output_dir)
    else:
        return self._run_subprocess(loom_path, config, output_dir)  # dev/test fallback
```

### RunPod handler (container-based)

Key decisions:
- SCENIC resources (rankings DB 1.3GB + motif DB) baked into Docker image at build time
  → no cold-start download latency
- loom files exchanged via S3 (Dramatiq worker uploads → RunPod downloads → uploads AUC)
- S3 bucket `axiom-sc-jobs` has 1-day lifecycle rule (privacy + cost)

```python
# runpod_handlers/scenic_handler.py
import runpod, boto3, json, subprocess, tempfile, os

s3 = boto3.client("s3")

def handler(job):
    inp = job["input"]
    # Download loom from S3
    with tempfile.TemporaryDirectory() as d:
        loom = f"{d}/input.loom"
        bucket, key = inp["loom_s3_uri"].replace("s3://","").split("/",1)
        s3.download_file(bucket, key, loom)

        config = {**inp["config"], "loom_path": loom, "output_dir": d,
                  "rankings_db_path": "/scenic_resources/hg38_rankings.feather",
                  "motif_db_path":    "/scenic_resources/motifs-v9.tbl",
                  "tf_list_path":     "/scenic_resources/allTFs_hg38.txt"}

        result = json.loads(subprocess.check_output(
            ["python", "/app/scenic_worker.py", json.dumps(config)], text=True))

        auc_key = inp["output_s3_prefix"] + "auc_matrix.csv"
        s3.upload_file(f"{d}/auc_matrix.csv", bucket, auc_key)
        return {"auc_s3_uri": f"s3://{bucket}/{auc_key}", **result}

runpod.serverless.start({"handler": handler})
```

### Modal alternative (Python-native, recommended for greenfield)

```python
# axiom_sc/tier2/scenic_modal.py
import modal

scenic_image = (
    modal.Image.debian_slim(python_version="3.10")
    .pip_install("setuptools<81", "dask[dataframe]==2024.2.0", "pyscenic", "loompy", "h5py")
    .run_commands(
        "python /app/patch_pyscenic_numpy.py",  # apply numpy 2.x patches
        "wget -q -O /scenic_resources/hg38_rankings.feather <URL>",  # bake in at build
        "wget -q -O /scenic_resources/motifs-v9.tbl <URL>",
        "wget -q -O /scenic_resources/allTFs_hg38.txt <URL>",
    )
)
app = modal.App("axiom-scenic")

@app.function(gpu="A10G", timeout=7200, retries=0, image=scenic_image)
def run_scenic(loom_bytes: bytes, config: dict) -> dict:
    import subprocess, json, tempfile
    with tempfile.TemporaryDirectory() as d:
        open(f"{d}/input.loom","wb").write(loom_bytes)
        config.update({"loom_path":f"{d}/input.loom","output_dir":d,
                        "rankings_db_path":"/scenic_resources/hg38_rankings.feather"})
        result = json.loads(subprocess.check_output(
            ["python","/app/scenic_worker.py",json.dumps(config)],text=True))
        return {**result, "auc_matrix_csv": open(f"{d}/auc_matrix.csv","rb").read()}
```

### Cost per annotation job (all tiers)

| Component | Service | Cost |
|---|---|---|
| Preprocessing + KG + Tier 3-5 | ECS Fargate (1vCPU) | ~$0.002 |
| pySCENIC | RunPod RTX 4090 (~30 min) | ~$0.15–0.40 |
| scVelo (optional) | RunPod RTX 4090 (~15 min) | ~$0.05–0.15 |
| LLM (Tier 4, optional) | External API | ~$0.01–0.05 |
| **Total** | | **~$0.20–0.60** |

---

## 29. Input format flexibility — real compbio workflows

### Typical researcher workflow reaching annotation

```
CellRanger → filter/QC → normalize (10k) → HVG → PCA → batch correct
→ kNN graph → Leiden cluster → DE markers → ANNOTATION ← AXIOM-SC
```

Researchers arrive with a **clustered h5ad ~65% of the time**.
They want each cluster labeled, not their preprocessing re-done.

### Five input scenarios

| Input | Freq | What it has | AXIOM action |
|---|---|---|---|
| Clustered h5ad | 65% | X normalized, obs['leiden'], X_umap | Compute DE markers → annotate |
| Unclustered processed h5ad | 15% | X normalized, X_pca | Auto-Leiden (res=0.5) → annotate |
| Raw count h5ad | 5% | Integer X, no preprocessing | Full preprocessing → annotate |
| Marker gene CSV | 10% | cluster_id, gene, logFC | Skip to Tier 2 directly |
| CellRanger .h5 | 5% | barcodes + features + matrix | Read 10x → full preprocessing |

### Input detection (axiom_sc/input/detector.py)

```python
def inspect_input(path: str) -> InputInspection:
    """Auto-detects input type in <2 seconds. Result shown in UI before job starts."""
    p = Path(path)
    if p.suffix == ".csv":
        return InputInspection(input_type=InputType.MARKER_CSV, ...)
    if p.suffix == ".h5":
        return InputInspection(input_type=InputType.CELLRANGER_H5, ...)
    # h5ad: check if raw (max X > 100), check for cluster key, check modality layers
    adata = anndata.read_h5ad(path, backed="r")
    is_raw = np.max(adata.X[:100]) > 100
    cluster_key = next((k for k in ["leiden","louvain","seurat_clusters","cell_type"]
                        if k in adata.obs.columns), None)
    has_velocity = "velocity" in adata.layers or "spliced" in adata.layers
    has_atac     = any(k in adata.obsm for k in ["X_atac","gene_activities"])
    has_spatial  = "spatial" in adata.obsm
    ...
    return InputInspection(input_type=..., n_cells=adata.n_obs, ...)
```

### InspectionPanel (shown in UI immediately after upload, before Run)

```
Detected: clustered h5ad
45,230 cells · 24 Leiden clusters · velocity detected · no ATAC
Estimated: ~8 min · Streams: 4/6 active (velocity ✓, ATAC ✗)
⚙ No preprocessing needed
⚠ Velocity layers found → scVelo stream will run (Tier 3)
```

### Mixed cluster detection

When Tier 1 cell-level predictions show <60% agreement within a cluster:

```python
# In annotator.py — tier1_cell_to_cluster():
if agreement_rate < 0.60:
    results[cluster_id] = {
        "label": f"Mixed: {top1} ({pct1:.0%}) / {top2} ({pct2:.0%})",
        "confidence": 0.2,
        "mixed": True,
        # Possible doublets or genuine cell type boundary
    }
```

### Marker CSV fast path

Accepts scanpy (`logfoldchanges` column), Seurat (`avg_log2FC`), or CASSIA-compatible
(`markers` column as comma-separated string). Skips Tier 1 and SCENIC entirely.
Completes in ~2 minutes at near-zero cost.

---

## 30. Updated Day 11–20 execution plan

### AXIOMTier1 GPU training — start today (parallel with playground development)

```bash
conda activate axiom-env
axiom-train-tier1 --n-models 10 --output model_weights/ \
  --census-version "2023-12-15" --n-cells 2000000 --epochs 5 --gpu-device cuda:0
# Submit to Vast.ai A100 — takes ~50 GPU hours total (~$75)
# Weights arrive ~Day 16 as planned
```

### Day 11 — FastAPI + Dramatiq (see Section 22)
Complete all steps in Section 22. End-of-day check:
```bash
curl http://localhost:8000/health
# → {sdk_version, kg_rule_count: ≥495, tier1_backend: "census_knn"}
```

### Day 12 — API tests + AWS infrastructure setup
Morning: write tests/test_api.py (see Section 22, Day 12 tests).
Afternoon: run all AWS CLI commands from Section 27 (ECR, ECS cluster, ElastiCache,
EFS, S3, Secrets Manager, IAM roles). Set GitHub secrets (AWS_ACCOUNT_ID, ALB_DNS_NAME).

### Day 13 — GitHub Actions CI/CD wiring
Write .github/workflows/deploy.yml (Section 27). Push to main.
Expected: all 3 images build → push to ECR → ECS services start → /health reachable via ALB.

### Days 14–16 — React components (Section 23)
Day 14: Components 1–5 (SDKVersionBadge through InspectionPanel)
Day 15: Components 6–9 (ProgressStream through VerdictCards)
Day 16: Components 10–12 (RuleFireHeatmap, KGBrowser, ComparisonTable)
Push each day → GitHub Actions deploys to ECS staging automatically.

### Day 17 — AXIOMTier1 weights + integration
Weights from GPU training job should be ready.
Upload to HuggingFace receptor-bio/axiomtier1 (Apache 2.0).
Update axiomtier1.py → /health returns "tier1_backend": "axiomtier1".

### Day 18 — RunPod/Modal setup + GPU dispatch wiring
Write runpod_handlers/scenic_handler.py (Section 28).
Build and push RunPod Docker image (with SCENIC resources baked in).
Add RunPod secrets to Secrets Manager.
Test: upload thymus h5ad → verify pySCENIC dispatches to RunPod → AUC returns.

### Day 19 — Full integration test + KG extension
End-to-end: upload Phase 1 thymus h5ad → all tiers → verify:
- pDC (thy-22) PROVEN with IRF7 circuit (visible in VerdictCards + PMID link)
- ILC3 (tab-20) PROVEN with RORC circuit
- ProgressStream shows "SCENIC+ via RunPod RTX 4090: 69 regulons, FOXP3 ✓ z=5.74"
Begin retina/bone/gonad KG extension (Section 25) in parallel.

### Day 20 — SDK release 0.2.0 + production deploy
```bash
# axiom-sc repo:
git tag v0.2.0 && git push origin v0.2.0
# → PyPI publish workflow fires → axiom-sc==0.2.0 on PyPI

# axiom-playground repo:
echo "sdk_version: 0.2.0" > config.yaml
git add config.yaml && git commit -m "chore: bump to axiom-sc 0.2.0"
git push origin main
# → GitHub Actions: requires production environment approval
# → Approve → ECS rolling update with axiom-sc==0.2.0 from PyPI
# → /health returns sdk_version: 0.2.0
```

---

## 31. Architectural summary — all decisions

| Decision | Choice | Reason |
|---|---|---|
| Task queue | **Dramatiq** | Celery has maintenance issues and task-duplication bugs |
| GPU compute | **RunPod serverless / Modal** | pySCENIC only; billed per-second, zero idle cost |
| Secrets | **AWS Secrets Manager** | Confidential values only (API keys). Not SSM. |
| Non-sensitive config | ECS environment variables | REDIS_URL, AXIOM_PROFILE |
| SDK version | Docker build arg from config.yaml | Not a runtime secret; baked at build time |
| Container registry | **AWS ECR** | 3 repos: api / worker / app |
| CI/CD auth | **GitHub Actions OIDC** | No stored AWS keys in GitHub |
| Shared storage | **AWS EFS** | h5ad upload shared between api + worker Fargate tasks |
| Redis | **AWS ElastiCache serverless** | Broker + result store + SSE pub/sub |
| Cluster-level vs cell-level | **Hybrid**: Tier 1 cell → aggregate → Tier 2 cluster | Best of both |
| Input formats | 5 accepted (h5ad, clustered/unclustered/raw, .h5, .csv) | Real compbio workflow coverage |
| Mixed clusters | Flagged distinctly (not forced label) | Honest output for doublets/boundaries |
| KG ceiling | ~200 types / ~650 rules | Beyond this, Tier 4 LLM handles dynamically |
| PMID verification | 3-step Entrez mandatory | Prevents Day-4-style PMID incident recurrence |
