"""
Tests for KG universe extension Days A-D.

Covers:
  - Scale targets: ≥640 ACTIVE rules, ≥198 cell types
  - All 4-day tissue coverage (retina, bone/gonads, inner ear, cornea,
    salivary gland, synovium, vascular specializations, misc)
  - Mutual exclusion rules (negative rules fire correctly)
  - PMID non-empty for all new rules
  - No duplicate rule IDs
"""
from __future__ import annotations

import pytest

from axiom_sc.tier2.kg_loader import KGLoader, KGRule

from pathlib import Path
import axiom_sc as _axiom_sc
KG_PATH = str(Path(_axiom_sc.__file__).parent / "kg_data" / "oracle_kg_v0.2.0.json")

DAYS_A_D_RULE_PREFIXES = [
    # Day A — retina
    "ROD_", "CONE_", "MUL_", "RGC_", "RPE_", "HORIZ_", "AMACRINE_",
    "BIPOLAR_", "CILRET_", "RETENDO_",
    # Day B — bone/gonads
    "OSTCYTE_", "OSTBL_", "OSTCL_", "CHOND_", "HYPCHOND_",
    "SSC_", "SPCYTE_", "SPTID_", "SERTOLI_", "LEYDIG2_", "GRANUL_", "THECA_", "OOCYTE_",
    # Day C — inner ear, cornea, salivary, synovium
    "IHC_", "OHC_", "COCH_SUPP_", "SGN_", "STRIAV_",
    "CORNEP_", "CORNKER_", "TM_", "CORNEND_NEG",
    "ACINAR_SER_", "ACINAR_MUC_", "STRI_DUCT_", "MYOEPI_",
    "SYNFIB_", "SYNMAC_",
    # Day D — vascular/misc
    "HEV_", "PERICB_", "SYNTH_SMC_", "AEROCYTE_", "PORTALENDO_",
    "PRBASSAL_", "PRLUM_", "UROBASE_", "UROSURF_",
    "RPMAC_", "MZB_", "FAP_", "ICC_", "CPE_",
]

DAY_A_TYPES = [
    "Rod_photoreceptor", "Cone_photoreceptor", "Muller_glia",
    "Retinal_Ganglion_cell", "Retinal_Pigment_Epithelium", "Horizontal_cell",
    "Amacrine_cell", "Bipolar_cell",
]
DAY_B_TYPES = [
    "Osteocyte", "Hypertrophic_chondrocyte",
    "Spermatogonial_SC", "Spermatocyte", "Spermatid",
    "Granulosa_cell", "Theca_cell", "Oocyte",
]
DAY_C_TYPES = [
    "Inner_Hair_cell", "Outer_Hair_cell", "Supporting_cell_cochlear",
    "Spiral_ganglion_neuron", "Stria_vascularis",
    "Corneal_epithelium", "Corneal_keratocyte", "Trabecular_meshwork",
    "Acinar_serous", "Acinar_mucous", "Striated_duct", "Myoepithelial",
    "Synovial_fibroblast", "Synovial_macrophage", "Chondroprogenitor",
]
DAY_D_TYPES = [
    "High_endothelial_venule", "Pericyte_brain", "Synthetic_SMC",
    "Lung_aerocyte", "Portal_endothelial",
    "Prostate_basal", "Prostate_luminal",
    "Urothelial_basal", "Urothelial_superficial",
    "Red_pulp_macrophage", "Marginal_zone_B",
    "Fibro_adipogenic_progenitor", "Interstitial_cell_Cajal",
    "Choroid_plexus_epithelium",
]


@pytest.fixture(scope="module")
def kg() -> KGLoader:
    return KGLoader.from_file(KG_PATH)


@pytest.fixture(scope="module")
def active_rules(kg: KGLoader) -> list[KGRule]:
    return kg.active_rules()


@pytest.fixture(scope="module")
def active_types(kg: KGLoader) -> set[str]:
    return set(kg.active_cell_types())


# ── Scale targets ─────────────────────────────────────────────────────────────

def test_active_rule_count_hits_640(active_rules):
    assert len(active_rules) >= 640, (
        f"Only {len(active_rules)} ACTIVE rules — need ≥640"
    )


def test_cell_type_count_hits_198(active_types):
    assert len(active_types) >= 198, (
        f"Only {len(active_types)} cell types — need ≥198"
    )


def test_negative_rules_at_least_120(active_rules):
    neg = [r for r in active_rules if r.rule_type == "negative"]
    assert len(neg) >= 120, (
        f"Only {len(neg)} negative rules — these are the core AXIOM innovation"
    )


def test_circuit_rules_at_least_90(active_rules):
    circ = [r for r in active_rules if r.rule_type == "circuit"]
    assert len(circ) >= 90, f"Only {len(circ)} circuit rules"


# ── No PMID missing ───────────────────────────────────────────────────────────

def test_no_unverified_pmids(active_rules):
    missing = [
        r.rule_id for r in active_rules
        if not r.pmid or r.pmid in ("", "NEEDS_REVIEW")
    ]
    assert missing == [], f"Rules with missing/unverified PMIDs: {missing}"


# ── No duplicate rule IDs ─────────────────────────────────────────────────────

def test_no_duplicate_rule_ids(kg: KGLoader):
    all_ids = [r.rule_id for r in kg]
    seen: set[str] = set()
    dupes = [rid for rid in all_ids if rid in seen or seen.add(rid)]
    assert dupes == [], f"Duplicate rule IDs: {dupes}"


# ── Day A: Retina coverage ────────────────────────────────────────────────────

def test_day_a_cell_types_present(active_types):
    missing = [ct for ct in DAY_A_TYPES if ct not in active_types]
    assert missing == [], f"Missing retina types: {missing}"


def test_rod_photoreceptor_has_positive_and_negative(kg: KGLoader):
    rules = kg.active_rules("Rod_photoreceptor")
    pos = [r for r in rules if r.rule_type == "positive"]
    neg = [r for r in rules if r.rule_type == "negative"]
    assert pos, "Rod_photoreceptor needs a positive rule"
    assert neg, "Rod_photoreceptor needs a negative rule (ARR3 absent)"
    neg_genes = {g for r in neg for g in r.gene_or_regulon}
    assert "ARR3" in neg_genes, "ARR3 should be the Rod→Cone exclusion gene"


def test_cone_photoreceptor_has_positive_and_negative(kg: KGLoader):
    rules = kg.active_rules("Cone_photoreceptor")
    pos = [r for r in rules if r.rule_type == "positive"]
    neg = [r for r in rules if r.rule_type == "negative"]
    assert pos, "Cone_photoreceptor needs a positive rule"
    assert neg, "Cone_photoreceptor needs a negative rule (NRL absent)"
    neg_genes = {g for r in neg for g in r.gene_or_regulon}
    assert "NRL" in neg_genes, "NRL should be the Cone→Rod exclusion gene"


def test_rpe_distinguished_from_muller_glia(kg: KGLoader):
    rpe_rules = kg.active_rules("Retinal_Pigment_Epithelium")
    neg = [r for r in rpe_rules if r.rule_type == "negative"]
    assert neg, "RPE must have a negative rule to distinguish from Muller_glia (both RLBP1+)"
    neg_genes = {g for r in neg for g in r.gene_or_regulon}
    # RPE should be RHO-negative (rod marker) — prevents false PROVEN on rod-contaminated RPE
    assert neg_genes.intersection({"RHO", "NRL", "PRKCZ"}), (
        f"RPE negative rule should exclude rod/Muller markers, got: {neg_genes}"
    )


def test_rgc_has_circuit_rule(kg: KGLoader):
    rules = kg.active_rules("Retinal_Ganglion_cell")
    circ = [r for r in rules if r.rule_type == "circuit"]
    assert circ, "RGC must have a circuit rule (RBPMS + SNCG co-requirement)"


def test_horizontal_cell_has_circuit_or_positive(kg: KGLoader):
    rules = kg.active_rules("Horizontal_cell")
    assert rules, "Horizontal_cell must have at least one rule"
    all_genes = {g for r in rules for g in r.gene_or_regulon}
    assert all_genes.intersection({"LHX1", "ONECUT1", "CALB1"}), (
        f"Horizontal_cell should express LHX1/ONECUT1/CALB1, got: {all_genes}"
    )


# ── Day B: Bone/Gonads coverage ───────────────────────────────────────────────

def test_day_b_cell_types_present(active_types):
    missing = [ct for ct in DAY_B_TYPES if ct not in active_types]
    assert missing == [], f"Missing bone/gonad types: {missing}"


def test_osteocyte_has_specific_markers(kg: KGLoader):
    rules = kg.active_rules("Osteocyte")
    all_genes = {g for r in rules for g in r.gene_or_regulon}
    assert all_genes.intersection({"SOST", "DMP1", "PHEX"}), (
        f"Osteocyte should have SOST/DMP1/PHEX markers, got: {all_genes}"
    )


def test_hypertrophic_chondrocyte_has_negative_sox9(kg: KGLoader):
    rules = kg.active_rules("Hypertrophic_chondrocyte")
    neg = [r for r in rules if r.rule_type == "negative"]
    assert neg, "Hypertrophic_chondrocyte must have a negative rule"
    neg_genes = {g for r in neg for g in r.gene_or_regulon}
    assert "SOX9" in neg_genes or "SOX9_regulon" in neg_genes, (
        "SOX9 must be inactive in hypertrophic chondrocytes (lost during hypertrophy)"
    )


def test_granulosa_theca_mutual_exclusion(kg: KGLoader):
    granul = kg.active_rules("Granulosa_cell")
    granul_neg_genes = {g for r in granul if r.rule_type == "negative" for g in r.gene_or_regulon}
    theca = kg.active_rules("Theca_cell")
    theca_neg_genes = {g for r in theca if r.rule_type == "negative" for g in r.gene_or_regulon}

    assert "CYP17A1" in granul_neg_genes, (
        "Granulosa_cell must exclude CYP17A1 (theca steroidogenesis enzyme)"
    )
    assert "CYP19A1" in theca_neg_genes, (
        "Theca_cell must exclude CYP19A1 (aromatase — granulosa specific)"
    )


def test_sertoli_leydig_rules_present(kg: KGLoader):
    for ct in ("Sertoli_cell", "Leydig_cell"):
        rules = kg.active_rules(ct)
        assert len(rules) >= 1, f"{ct} must have at least one rule"


def test_spermatocyte_has_meiosis_markers(kg: KGLoader):
    rules = kg.active_rules("Spermatocyte")
    all_genes = {g for r in rules for g in r.gene_or_regulon}
    assert all_genes.intersection({"SYCP3", "SYCP1", "MLH1"}), (
        f"Spermatocyte should express meiotic markers, got: {all_genes}"
    )


def test_oocyte_has_zona_pellucida_markers(kg: KGLoader):
    rules = kg.active_rules("Oocyte")
    all_genes = {g for r in rules for g in r.gene_or_regulon}
    assert all_genes.intersection({"ZP1", "ZP2", "GDF9", "BMP15"}), (
        f"Oocyte should have zona pellucida/oocyte-specific genes, got: {all_genes}"
    )


# ── Day C: Inner ear, cornea, salivary, synovium coverage ────────────────────

def test_day_c_cell_types_present(active_types):
    missing = [ct for ct in DAY_C_TYPES if ct not in active_types]
    assert missing == [], f"Missing Day C types: {missing}"


def test_inner_outer_hair_cell_mutual_exclusion(kg: KGLoader):
    ihc = kg.active_rules("Inner_Hair_cell")
    ohc = kg.active_rules("Outer_Hair_cell")

    ihc_genes = {g for r in ihc for g in r.gene_or_regulon}
    ohc_neg_genes = {g for r in ohc if r.rule_type == "negative" for g in r.gene_or_regulon}

    assert ihc_genes.intersection({"MYO7A", "ESPN", "ATOH1"}), (
        f"IHC should express MYO7A/ESPN/ATOH1, got: {ihc_genes}"
    )
    assert ohc_neg_genes.intersection({"ATOH1", "ESPN"}), (
        f"OHC must exclude IHC markers, got neg: {ohc_neg_genes}"
    )


def test_outer_hair_cell_has_slc26a5(kg: KGLoader):
    rules = kg.active_rules("Outer_Hair_cell")
    all_genes = {g for r in rules for g in r.gene_or_regulon}
    assert "SLC26A5" in all_genes, (
        "Outer_Hair_cell must include SLC26A5/prestin (electromotility gene)"
    )


def test_stria_vascularis_has_kcnj10(kg: KGLoader):
    rules = kg.active_rules("Stria_vascularis")
    all_genes = {g for r in rules for g in r.gene_or_regulon}
    assert "KCNJ10" in all_genes, "Stria_vascularis KCNJ10 required for endocochlear potential"


def test_corneal_epithelium_has_krt12_and_negative(kg: KGLoader):
    rules = kg.active_rules("Corneal_epithelium")
    pos_genes = {g for r in rules if r.rule_type == "positive" for g in r.gene_or_regulon}
    neg = [r for r in rules if r.rule_type == "negative"]
    assert "KRT12" in pos_genes or "KRT3" in pos_genes, (
        f"Corneal_epithelium should be KRT12+/KRT3+, got pos: {pos_genes}"
    )
    assert neg, "Corneal_epithelium must have a negative rule"


def test_corneal_keratocyte_negative_acta2(kg: KGLoader):
    rules = kg.active_rules("Corneal_keratocyte")
    neg = [r for r in rules if r.rule_type == "negative"]
    assert neg, "Corneal_keratocyte must have a negative rule"
    neg_genes = {g for r in neg for g in r.gene_or_regulon}
    assert "ACTA2" in neg_genes, (
        "Corneal_keratocyte must exclude ACTA2 (myofibroblast activation marker)"
    )


def test_acinar_serous_mucous_mutual_exclusion(kg: KGLoader):
    serous = kg.active_rules("Acinar_serous")
    mucous = kg.active_rules("Acinar_mucous")

    serous_genes = {g for r in serous if r.rule_type == "positive" for g in r.gene_or_regulon}
    mucous_neg_genes = {g for r in mucous if r.rule_type == "negative" for g in r.gene_or_regulon}

    assert serous_genes.intersection({"PRSS2", "BPIFB2", "ZG16B"}), (
        f"Acinar_serous should express serous enzymes, got: {serous_genes}"
    )
    assert mucous_neg_genes.intersection({"PRSS2", "BPIFB2"}), (
        f"Acinar_mucous must exclude serous enzymes, got neg: {mucous_neg_genes}"
    )


def test_synovial_fibroblast_has_prg4(kg: KGLoader):
    rules = kg.active_rules("Synovial_fibroblast")
    all_genes = {g for r in rules for g in r.gene_or_regulon}
    assert "PRG4" in all_genes, "Synovial_fibroblast must include PRG4/lubricin"


# ── Day D: Vascular and misc coverage ────────────────────────────────────────

def test_day_d_cell_types_present(active_types):
    missing = [ct for ct in DAY_D_TYPES if ct not in active_types]
    assert missing == [], f"Missing Day D types: {missing}"


def test_high_endothelial_venule_has_ackr1(kg: KGLoader):
    rules = kg.active_rules("High_endothelial_venule")
    all_genes = {g for r in rules for g in r.gene_or_regulon}
    assert "ACKR1" in all_genes, "HEV must express ACKR1/DARC (chemokine scavenger)"


def test_synthetic_smc_negative_myh11(kg: KGLoader):
    rules = kg.active_rules("Synthetic_SMC")
    neg = [r for r in rules if r.rule_type == "negative"]
    assert neg, "Synthetic_SMC must have a negative rule"
    neg_genes = {g for r in neg for g in r.gene_or_regulon}
    assert "MYH11" in neg_genes, (
        "Synthetic_SMC must exclude MYH11 (lost contractile program)"
    )


def test_pericyte_brain_has_kcnj8(kg: KGLoader):
    rules = kg.active_rules("Pericyte_brain")
    all_genes = {g for r in rules for g in r.gene_or_regulon}
    assert "KCNJ8" in all_genes or "PDGFRB" in all_genes, (
        f"Pericyte_brain should have KCNJ8/PDGFRB, got: {all_genes}"
    )


def test_lung_aerocyte_has_ca4(kg: KGLoader):
    rules = kg.active_rules("Lung_aerocyte")
    all_genes = {g for r in rules for g in r.gene_or_regulon}
    assert "CA4" in all_genes, "Lung_aerocyte must include CA4 (carbonic anhydrase IV)"


def test_lung_aerocyte_has_negative_rule(kg: KGLoader):
    rules = kg.active_rules("Lung_aerocyte")
    neg = [r for r in rules if r.rule_type == "negative"]
    assert neg, "Lung_aerocyte must have a negative rule to distinguish from other capillary ECs"


def test_portal_endothelial_negative_stab2(kg: KGLoader):
    rules = kg.active_rules("Portal_endothelial")
    neg = [r for r in rules if r.rule_type == "negative"]
    assert neg, "Portal_endothelial must have a negative rule"
    neg_genes = {g for r in neg for g in r.gene_or_regulon}
    assert "STAB2" in neg_genes, (
        "Portal_endothelial must exclude STAB2 (sinusoidal EC marker)"
    )


def test_red_pulp_macrophage_has_negative(kg: KGLoader):
    rules = kg.active_rules("Red_pulp_macrophage")
    neg = [r for r in rules if r.rule_type == "negative"]
    assert neg, "Red_pulp_macrophage must have a negative rule"
    neg_genes = {g for r in neg for g in r.gene_or_regulon}
    assert "SIGLEC1" in neg_genes, (
        "Red_pulp_macrophage must exclude SIGLEC1 (metallophilic marginal-zone macrophage)"
    )


def test_marginal_zone_b_rules_present(kg: KGLoader):
    rules = kg.active_rules("Marginal_zone_B")
    assert len(rules) >= 2, "Marginal_zone_B should have positive and negative rules"


def test_fibro_adipogenic_progenitor_negative_myod1(kg: KGLoader):
    rules = kg.active_rules("Fibro_adipogenic_progenitor")
    neg = [r for r in rules if r.rule_type == "negative"]
    assert neg, "FAP must have a negative rule"
    neg_genes = {g for r in neg for g in r.gene_or_regulon}
    assert "MYOD1" in neg_genes, "FAP must exclude MYOD1 (myoblast master TF)"


def test_corneal_endothelial_negative_krt12(kg: KGLoader):
    rules = kg.active_rules("Corneal_Endothelial")
    neg = [r for r in rules if r.rule_type == "negative"]
    assert neg, "Corneal_Endothelial must have a negative rule"
    neg_genes = {g for r in neg for g in r.gene_or_regulon}
    assert "KRT12" in neg_genes, (
        "Corneal_Endothelial must exclude KRT12 (corneal epithelium keratin)"
    )


# ── Full KG integrity ─────────────────────────────────────────────────────────

def test_all_new_types_have_at_least_two_rules(active_rules):
    # Note: Retinal_Ganglion_cell and Retinal_Pigment_Epithelium were enriched from
    # pre-existing single-rule stubs; exclude from the "new types" check since they
    # were in the KG before Days A-D.
    all_new_types = DAY_A_TYPES + DAY_B_TYPES + DAY_C_TYPES + DAY_D_TYPES
    by_type: dict[str, list] = {}
    for r in active_rules:
        by_type.setdefault(r.cell_type, []).append(r)
    weak = [ct for ct in all_new_types
            if ct in by_type and len(by_type[ct]) < 2]
    assert weak == [], (
        f"These new cell types have fewer than 2 rules (need pos+neg): {weak}"
    )


def test_pmid_format_all_new_rules(active_rules):
    day_a_d_rules = [
        r for r in active_rules
        if r.cell_type in set(DAY_A_TYPES + DAY_B_TYPES + DAY_C_TYPES + DAY_D_TYPES)
    ]
    bad = [r.rule_id for r in day_a_d_rules if not r.pmid.isdigit()]
    assert bad == [], f"Non-numeric PMIDs in Days A-D rules: {bad}"


def test_final_stats(active_rules):
    types = {r.cell_type for r in active_rules}
    neg = [r for r in active_rules if r.rule_type == "negative"]
    circ = [r for r in active_rules if r.rule_type == "circuit"]
    assert len(active_rules) >= 640
    assert len(types) >= 198
    assert len(neg) >= 120
    assert len(circ) >= 90
