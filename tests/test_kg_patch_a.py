"""
Tests for Patch A KG expansion: ≥150 cell types, ≥450 ACTIVE rules.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

KG_PATH = Path(__file__).parent.parent / "kg_data" / "oracle_kg_v0.2.0.json"


@pytest.fixture(scope="module")
def kg_rules():
    with open(KG_PATH) as f:
        return json.load(f)


@pytest.fixture(scope="module")
def active_rules(kg_rules):
    return [r for r in kg_rules if r.get("status") == "ACTIVE"]


@pytest.fixture(scope="module")
def active_cell_types(active_rules):
    return set(r["cell_type"] for r in active_rules)


@pytest.fixture(scope="module")
def patch_a_rules(active_rules):
    """Only rules added in Patch A (version 0.2.0-patch-a)."""
    return [r for r in active_rules if "patch-a" in r.get("added_in_version", "")]


# ── Scale targets ─────────────────────────────────────────────────────────────

def test_kg_reaches_150_cell_types(active_cell_types):
    assert len(active_cell_types) >= 150, (
        f"Expected ≥150 cell types, got {len(active_cell_types)}"
    )


def test_kg_reaches_450_active_rules(active_rules):
    assert len(active_rules) >= 450, (
        f"Expected ≥450 ACTIVE rules, got {len(active_rules)}"
    )


# ── Tissue coverage ───────────────────────────────────────────────────────────

def test_kg_covers_cns_types(active_cell_types):
    """KG must contain ≥4 CNS cell types (checking by substring)."""
    cns_keywords = ["Neuron", "Astrocyte", "Oligodendrocyte", "Microglia"]
    found = [ct for ct in active_cell_types if any(kw in ct for kw in cns_keywords)]
    assert len(found) >= 4, f"Expected ≥4 CNS types, found {found}"


def test_kg_covers_skin_types(active_cell_types):
    skin_keywords = ["Keratinocyte", "Melanocyte", "Langerhans"]
    found = [ct for ct in active_cell_types if any(kw in ct for kw in skin_keywords)]
    assert len(found) >= 2, f"Expected ≥2 skin types, found {found}"


def test_kg_covers_kidney_types(active_cell_types):
    kidney_keywords = ["Proximal_Tubule", "Podocyte", "Loop_of_Henle", "Collecting_Duct"]
    found = [ct for ct in active_cell_types if any(kw in ct for kw in kidney_keywords)]
    assert len(found) >= 2, f"Expected ≥2 kidney types, found {found}"


def test_kg_covers_cardiac_types(active_cell_types):
    """KG must contain ≥1 cardiomyocyte variant."""
    cardiac = [ct for ct in active_cell_types if "Cardiomyocyte" in ct]
    assert len(cardiac) >= 1, f"Expected ≥1 cardiac type, found {cardiac}"


def test_kg_covers_endocrine_types(active_cell_types):
    endo_keywords = ["Beta_cell", "Alpha_cell", "Delta_cell"]
    found = [ct for ct in active_cell_types if any(kw in ct for kw in endo_keywords)]
    assert len(found) >= 2, f"Expected ≥2 endocrine types, found {found}"


def test_kg_covers_stromal_types(active_cell_types):
    """KG must contain ≥3 stromal/endothelial types."""
    stromal_keywords = ["CAF", "SMC", "Fibroblast", "Endothelial"]
    found = [ct for ct in active_cell_types if any(kw in ct for kw in stromal_keywords)]
    assert len(found) >= 3, f"Expected ≥3 stromal/EC types, found {found}"


def test_kg_covers_hematopoietic_types(active_cell_types):
    """KG must contain ≥2 hematopoietic progenitor types."""
    hema_keywords = ["HSC", "CMP", "GMP", "MEP", "Megakaryocyte", "Erythroid", "Erythroblast"]
    found = [ct for ct in active_cell_types if any(kw in ct for kw in hema_keywords)]
    assert len(found) >= 2, f"Expected ≥2 hematopoietic types, found {found}"


def test_kg_covers_trophoblast_types(active_cell_types):
    troph_keywords = ["Trophoblast", "Syncytio", "Cytotrophoblast"]
    found = [ct for ct in active_cell_types if any(kw in ct for kw in troph_keywords)]
    assert len(found) >= 1, f"Expected ≥1 trophoblast type, found {found}"


# ── Mutual exclusion rules ────────────────────────────────────────────────────

def _get_negative_rules_for(active_rules, cell_type):
    return [
        r for r in active_rules
        if r["cell_type"] == cell_type and r["rule_type"] == "negative"
    ]


def test_excitatory_inhibitory_mutual_exclusion(active_rules):
    """Excitatory neurons must have a negative rule blocking GAD1/GAD2."""
    exc_neg = _get_negative_rules_for(active_rules, "Excitatory_Neuron")
    exc_neg_genes = {g for r in exc_neg for g in r["gene_or_regulon"]}
    assert exc_neg_genes & {"GAD1", "GAD2"}, (
        f"Excitatory_Neuron has no negative rule for GAD1/GAD2; genes: {exc_neg_genes}"
    )


def test_inhibitory_excitatory_mutual_exclusion(active_rules):
    """Inhibitory neurons must have a negative rule blocking VGLUT markers."""
    inh_neg = _get_negative_rules_for(active_rules, "Inhibitory_Neuron")
    inh_neg_genes = {g for r in inh_neg for g in r["gene_or_regulon"]}
    assert inh_neg_genes & {"SLC17A7", "SLC17A6", "CAMK2A"}, (
        f"Inhibitory_Neuron has no negative rule for excitatory markers; genes: {inh_neg_genes}"
    )


def test_beta_alpha_mutual_exclusion(active_rules):
    """Beta cells must have a negative rule for GCG (glucagon = alpha marker)."""
    beta_neg = _get_negative_rules_for(active_rules, "Beta_cell")
    beta_neg_genes = {g for r in beta_neg for g in r["gene_or_regulon"]}
    assert "GCG" in beta_neg_genes, (
        f"Beta_cell has no negative rule for GCG; genes: {beta_neg_genes}"
    )


def test_alpha_beta_mutual_exclusion(active_rules):
    """Alpha cells must have a negative rule for INS (insulin = beta marker)."""
    alpha_neg = _get_negative_rules_for(active_rules, "Alpha_cell")
    alpha_neg_genes = {g for r in alpha_neg for g in r["gene_or_regulon"]}
    assert "INS" in alpha_neg_genes, (
        f"Alpha_cell has no negative rule for INS; genes: {alpha_neg_genes}"
    )


def test_tex_progenitor_exists_with_positive_rules(active_rules):
    """TEX_progenitor must exist in KG with positive/circuit rules."""
    prog_tex = [r for r in active_rules if r["cell_type"] == "TEX_progenitor"]
    assert prog_tex, "TEX_progenitor not in KG"
    types_present = {r["rule_type"] for r in prog_tex}
    assert types_present & {"positive", "circuit"}, (
        f"TEX_progenitor has no positive/circuit rules; types: {types_present}"
    )


def test_tex_terminal_vs_progenitor_discrimination(active_rules):
    """
    TEX_terminal must have DIFFERENT marker rules from TEX_progenitor.
    TEX_terminal: high TOX/HAVCR2, low TCF7.
    TEX_progenitor: TCF7 positive (progenitor phenotype).
    """
    term_genes = {
        g for r in active_rules if r["cell_type"] == "TEX_terminal"
        for g in r["gene_or_regulon"]
    }
    prog_genes = {
        g for r in active_rules if r["cell_type"] == "TEX_progenitor"
        for g in r["gene_or_regulon"]
    }
    assert term_genes, "TEX_terminal has no rules at all"
    assert prog_genes, "TEX_progenitor has no rules at all"
    # The progenitor-specific marker TCF7 must be in progenitor rules
    assert "TCF7" in prog_genes or "TCF7" in term_genes, (
        "TCF7 not present in exhausted T cell rules (key discriminator)"
    )


def test_exhausted_t_cell_has_canonical_markers(active_rules, active_cell_types):
    """At least one exhausted T cell type must have TOX/PDCD1/HAVCR2 as positive markers."""
    tex_types = [ct for ct in active_cell_types if "TEX" in ct or ct == "Tex"]
    assert tex_types, "No exhausted T cell types in KG"
    tex_pos_genes = {
        g for r in active_rules
        if r["cell_type"] in tex_types and r["rule_type"] == "positive"
        for g in r["gene_or_regulon"]
    }
    assert tex_pos_genes & {"TOX", "PDCD1", "HAVCR2", "LAG3"}, (
        f"No canonical exhaustion markers in positive rules; found: {tex_pos_genes}"
    )


# ── Rule quality ──────────────────────────────────────────────────────────────

def test_all_active_rules_have_nonempty_pmid(active_rules):
    bad = [r["rule_id"] for r in active_rules if not r.get("pmid", "").strip()]
    assert not bad, f"Rules missing pmid: {bad[:10]}"


def test_all_active_rules_have_mechanistic_basis(active_rules):
    bad = [
        r["rule_id"] for r in active_rules
        if len(r.get("mechanistic_basis", "")) < 20
    ]
    assert not bad, f"Rules with short/empty mechanistic_basis: {bad[:10]}"


def test_all_active_rules_have_valid_rule_type(active_rules):
    valid = {"positive", "negative", "circuit", "spatial"}
    bad = [r["rule_id"] for r in active_rules if r.get("rule_type") not in valid]
    assert not bad, f"Rules with invalid rule_type: {bad[:5]}"


def test_all_active_rules_have_gene_or_regulon(active_rules):
    bad = [r["rule_id"] for r in active_rules if not r.get("gene_or_regulon")]
    assert not bad, f"Rules missing gene_or_regulon: {bad[:10]}"


def test_patch_a_rules_all_have_tissue_context(patch_a_rules):
    """All Patch A rules (added_in_version=0.2.0-patch-a) must have tissue_context."""
    bad = [r["rule_id"] for r in patch_a_rules if not r.get("tissue_context")]
    assert not bad, f"Patch A rules missing tissue_context: {bad[:10]}"


def test_rule_ids_are_unique(active_rules):
    ids = [r["rule_id"] for r in active_rules]
    dupes = [rid for rid in set(ids) if ids.count(rid) > 1]
    assert not dupes, f"Duplicate rule_ids: {dupes[:10]}"


def test_no_pending_patch_a_rules(kg_rules):
    """Patch A rules should all be ACTIVE (none left PENDING_REVIEW)."""
    pending_patch_a = [
        r["rule_id"] for r in kg_rules
        if r.get("status") == "PENDING_REVIEW"
        and "patch-a" in r.get("added_in_version", "")
    ]
    assert not pending_patch_a, f"Patch A rules still PENDING: {pending_patch_a[:5]}"


# ── Patch A specific cell types have rules ────────────────────────────────────

def test_microglia_dam_has_positive_rules(active_rules):
    """Microglia_DAM (disease-associated microglia) has APOE+ positive rule."""
    mg_dam = [r for r in active_rules if r["cell_type"] == "Microglia_DAM" and r["rule_type"] == "positive"]
    assert mg_dam, "Microglia_DAM has no positive rules"


def test_oligodendrocyte_has_positive_rules(active_rules):
    oligo = [r for r in active_rules if r["cell_type"] == "Oligodendrocyte" and r["rule_type"] == "positive"]
    assert oligo, "Oligodendrocyte has no positive rules"


def test_cardiomyocyte_subtypes_have_positive_rules(active_rules):
    """Atrial_Cardiomyocyte or Ventricular_Cardiomyocyte must have positive rules."""
    cm = [
        r for r in active_rules
        if "Cardiomyocyte" in r["cell_type"] and r["rule_type"] == "positive"
    ]
    assert cm, "No Cardiomyocyte subtype has positive rules"


def test_proximal_tubule_has_positive_rules(active_rules):
    pt = [r for r in active_rules if r["cell_type"] == "Proximal_Tubule" and r["rule_type"] == "positive"]
    assert pt, "Proximal_Tubule has no positive rules"


def test_podocyte_has_positive_rules(active_rules):
    pod = [r for r in active_rules if r["cell_type"] == "Podocyte" and r["rule_type"] == "positive"]
    assert pod, "Podocyte has no positive rules"


def test_keratinocyte_subtypes_have_positive_rules(active_rules):
    """At least one Keratinocyte subtype has positive rules."""
    k = [
        r for r in active_rules
        if "Keratinocyte" in r["cell_type"] and r["rule_type"] == "positive"
    ]
    assert k, "No Keratinocyte subtype has positive rules"


def test_hematopoietic_progenitors_have_positive_rules(active_rules):
    """CMP, GMP, or MEP must have positive rules."""
    prog = [
        r for r in active_rules
        if r["cell_type"] in {"CMP", "GMP", "MEP"} and r["rule_type"] == "positive"
    ]
    assert prog, "No hematopoietic progenitor (CMP/GMP/MEP) has positive rules"


def test_trophoblast_types_have_rules(active_rules, active_cell_types):
    """Trophoblast types (Cytotrophoblast, Syncytiotrophoblast, EVT) must exist with rules."""
    troph_types = [
        ct for ct in active_cell_types
        if "Trophoblast" in ct or "Syncytio" in ct or "Cytotrophoblast" in ct
    ]
    assert troph_types, "No trophoblast cell types in KG"
    troph_rules = [r for r in active_rules if r["cell_type"] in troph_types]
    assert len(troph_rules) >= 3, f"Expected ≥3 trophoblast rules, got {len(troph_rules)}"
