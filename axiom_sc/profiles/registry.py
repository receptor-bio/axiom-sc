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
        "10.1038/s41467-024-46499-y",
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
        "axiom-sc-phase1",
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
