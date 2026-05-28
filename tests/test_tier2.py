"""
Tier 2 annotator tests — Day 2 requirement.

Six mandatory test cases (all must pass):
  1. Treg with IL2 high → CONTRADICTED (TREG_NEG_001)
  2. ILC1 with TRAC present → CONTRADICTED (ILC1_NEG_001)
  3. pDC with IRF7 active + SIGLEC1 high + PAX5 regulon inactive → PROVEN
  4. ILC3 with RORC active + NCR2 present → PROVEN (Phase 2 fix)
  5. Myofibroblast with RORC active + no NCR2 → NOT PROVEN (circuit NOT_TESTABLE)
  6. Full-KG mode: exhaust all 18 types, return highest-support verdict (pDC)

Additional tests verify each Phase 1 correction individually.
"""
from pathlib import Path

import pytest

from axiom_sc.tier2.axiom_annotator import AxiomAnnotator, Verdict, RuleVerdict
from axiom_sc.tier2.evidence import EvidenceBundle
from axiom_sc.tier2.kg_loader import KGLoader

# ── fixtures ──────────────────────────────────────────────────────────────────

KG_PATH = Path(__file__).parent.parent / "kg_data" / "oracle_kg_v0.2.0.json"


@pytest.fixture(scope="session")
def kg() -> KGLoader:
    return KGLoader.from_file(KG_PATH)


@pytest.fixture(scope="session")
def annotator(kg) -> AxiomAnnotator:
    return AxiomAnnotator(kg=kg)


# ── required test 1 ───────────────────────────────────────────────────────────

def test_treg_il2_high_contradicted(annotator):
    """
    TREG_NEG_001: Treg should NOT produce IL-2.
    If IL2 is high in a cluster → that cluster is CONTRADICTED as Treg.
    """
    evidence = EvidenceBundle(
        cluster_id="test-treg-contradicted",
        marker_genes={
            "IL2": 3.0,     # high IL2 → TREG_NEG_001 fires → CONTRADICTED
            "IL2RA": 2.0,   # high CD25 (would normally support Treg)
            "FOXP3": 2.5,   # high FOXP3 (would normally support Treg)
        },
        regulons={
            "FOXP3": 2.2,   # FOXP3 regulon active (positive support, irrelevant once contradicted)
        },
    )

    verdicts = annotator.annotate_cluster(evidence)

    assert "Treg" in verdicts, "Treg must be tested in full-KG mode"
    v = verdicts["Treg"]
    assert v.verdict == Verdict.CONTRADICTED, (
        f"Expected CONTRADICTED for Treg with IL2=3.0, got {v.verdict}. "
        f"Rule firings: {[(f.rule_id, f.verdict) for f in v.rule_firings]}"
    )
    # Verify the specific rule that fired
    contradicting_rules = [f.rule_id for f in v.failed_negative_firings()]
    assert "TREG_NEG_001" in contradicting_rules, (
        f"Expected TREG_NEG_001 in contradictions, got {contradicting_rules}"
    )


# ── required test 2 ───────────────────────────────────────────────────────────

def test_ilc1_trac_present_contradicted(annotator):
    """
    ILC1_NEG_001: ILCs lack TCR. TRAC present → NOT an ILC1.
    """
    evidence = EvidenceBundle(
        cluster_id="test-ilc1-contradicted",
        marker_genes={
            "TRAC": 2.5,   # TRAC present → ILC1_NEG_001 fires → CONTRADICTED
            "NKG7": 3.0,   # NK/ILC effector marker (would support ILC1)
        },
        regulons={
            "TBX21": 2.5,  # TBX21 active (would support ILC1)
        },
    )

    verdicts = annotator.annotate_cluster(evidence)

    assert "ILC1" in verdicts
    v = verdicts["ILC1"]
    assert v.verdict == Verdict.CONTRADICTED, (
        f"Expected CONTRADICTED for ILC1 with TRAC=2.5, got {v.verdict}"
    )
    contradicting_rules = [f.rule_id for f in v.failed_negative_firings()]
    assert "ILC1_NEG_001" in contradicting_rules


# ── required test 3 ───────────────────────────────────────────────────────────

def test_pdc_irf7_siglec1_pax5_inactive_proven(annotator):
    """
    pDC PROVEN: IRF7 regulon active + SIGLEC1 high + PAX5 regulon INACTIVE.

    Phase 1 PDC_NEG_001 correction: check PAX5 REGULON z-score < -0.5
    (not IGKC absence — pDCs naturally express IGKC).
    PAX5 inactive → pDC not contradicted → PROVEN.
    """
    evidence = EvidenceBundle(
        cluster_id="thy-22",
        marker_genes={
            "SIGLEC1": 3.2,   # high → PDC_POS_001 PASS
        },
        regulons={
            "IRF7": 3.7,      # active → PDC_CIRCUIT_001 PASS (circuit)
            "PAX5": -1.0,     # inactive (z < -0.5) → PDC_NEG_001 PASS (not contradicted)
        },
    )

    verdicts = annotator.annotate_cluster(evidence)

    assert "pDC" in verdicts
    v = verdicts["pDC"]
    assert v.verdict == Verdict.PROVEN, (
        f"Expected PROVEN for pDC (IRF7 active + SIGLEC1 high + PAX5 inactive), "
        f"got {v.verdict}. "
        f"Firings: {[(f.rule_id, f.verdict.value) for f in v.rule_firings]}"
    )
    assert v.support_count >= 2, (
        f"Expected support_count ≥ 2 (circuit + positive), got {v.support_count}"
    )


def test_pdc_pax5_active_contradicted(annotator):
    """
    PDC_NEG_001 (corrected): active PAX5 regulon → pDC CONTRADICTED.
    This verifies the correction fires correctly in the positive direction.
    """
    evidence = EvidenceBundle(
        cluster_id="test-pdc-b-cell",
        marker_genes={"SIGLEC1": 3.2},
        regulons={
            "IRF7": 3.7,
            "PAX5": 2.0,   # ACTIVE → PDC_NEG_001 fires → CONTRADICTED
        },
    )

    verdicts = annotator.annotate_cluster(evidence)
    v = verdicts["pDC"]
    assert v.verdict == Verdict.CONTRADICTED, (
        f"Expected pDC CONTRADICTED when PAX5 regulon active, got {v.verdict}"
    )


# ── required test 4 ───────────────────────────────────────────────────────────

def test_ilc3_rorc_ncr2_present_proven(annotator):
    """
    ILC3 PROVEN: RORC regulon active + NCR2 present (Phase 2 circuit fix).
    Both required for ILC3_CIRCUIT_001 PASS.
    """
    evidence = EvidenceBundle(
        cluster_id="tab-20-ilc3",
        marker_genes={
            "NCR2": 2.0,   # present → satisfies 'active' direction (z ≥ 1.5)
        },
        regulons={
            "RORC": 3.37,  # active → ILC3_POS_001 PASS + part of circuit
        },
    )

    verdicts = annotator.annotate_cluster(evidence)

    assert "ILC3" in verdicts
    v = verdicts["ILC3"]
    assert v.verdict == Verdict.PROVEN, (
        f"Expected ILC3 PROVEN (RORC active + NCR2 present), got {v.verdict}. "
        f"Firings: {[(f.rule_id, f.verdict.value, f.notes) for f in v.rule_firings]}"
    )

    # Verify the circuit rule passed
    circuit_fires = [f for f in v.rule_firings if f.rule_type == "circuit"]
    assert any(f.verdict == RuleVerdict.PASS for f in circuit_fires), (
        "ILC3_CIRCUIT_001 must PASS when both RORC and NCR2 are present"
    )


# ── required test 5 ───────────────────────────────────────────────────────────

def test_myofibroblast_rorc_no_ncr2_not_proven(annotator):
    """
    Phase 2 fix verified: myofibroblast with RORC active but NO NCR2.
    ILC3_CIRCUIT_001 is NOT_TESTABLE (NCR2 absent from evidence) → ILC3 NOT PROVEN.

    This was the Phase 1 false-positive: non-immune RORC expression (circadian
    regulation) caused myofibroblasts to spuriously call as ILC3 in Phase 1.
    """
    evidence = EvidenceBundle(
        cluster_id="test-myofibroblast",
        marker_genes={
            # NCR2 deliberately absent — myofibroblasts don't express innate receptors
            "ACTA2": 3.5,   # alpha-SMA — myofibroblast marker
        },
        regulons={
            "RORC": 3.0,   # RORC active (circadian expression, not ILC3-specific)
        },
    )

    verdicts = annotator.annotate_cluster(evidence)

    assert "ILC3" in verdicts
    v = verdicts["ILC3"]
    assert v.verdict != Verdict.PROVEN, (
        f"ILC3 must NOT be PROVEN when NCR2 is absent "
        f"(Phase 2 circuit fix). Got {v.verdict}. "
        f"Firings: {[(f.rule_id, f.verdict.value, f.notes) for f in v.rule_firings]}"
    )

    # Verify the circuit is NOT_TESTABLE (NCR2 missing)
    circuit_fires = [f for f in v.rule_firings
                     if f.rule_type == "circuit" and f.rule_id == "ILC3_CIRCUIT_001"]
    assert circuit_fires, "ILC3_CIRCUIT_001 must be in firings"
    assert circuit_fires[0].verdict == RuleVerdict.NOT_TESTABLE, (
        f"ILC3_CIRCUIT_001 must be NOT_TESTABLE when NCR2 absent, "
        f"got {circuit_fires[0].verdict}: {circuit_fires[0].notes}"
    )


# ── required test 6 ───────────────────────────────────────────────────────────

def test_full_kg_exhaustive_pdc_best(annotator, kg):
    """
    FULL KG TEST MODE (Phase 1 correction 1):
    annotate_cluster tests ALL cell types in KG.
    With pDC evidence, all 18 types are tested and pDC wins.
    """
    evidence = EvidenceBundle(
        cluster_id="thy-22-full",
        marker_genes={
            "SIGLEC1": 3.2,
        },
        regulons={
            "IRF7": 3.7,
            "PAX5": -1.0,
        },
    )

    verdicts = annotator.annotate_cluster(evidence)

    # All active cell types in KG must be tested (count grows as KG expands)
    expected_count = len(kg.active_cell_types())
    assert expected_count >= 18, f"KG must have ≥ 18 active cell types, has {expected_count}"
    assert len(verdicts) == expected_count, (
        f"Full-KG mode must test all {expected_count} cell types, "
        f"tested {len(verdicts)}: {sorted(verdicts.keys())}"
    )

    # pDC must be PROVEN
    assert verdicts["pDC"].verdict == Verdict.PROVEN, (
        f"pDC must be PROVEN in full-KG mode, got {verdicts['pDC'].verdict}"
    )

    # get_best_verdict must return pDC
    best = annotator.get_best_verdict(evidence)
    assert best.cell_type == "pDC", (
        f"Best verdict must be pDC, got {best.cell_type} ({best.verdict})"
    )
    assert best.verdict == Verdict.PROVEN


# ── Phase 1 correction tests ──────────────────────────────────────────────────

def test_correction_2_circuit_completeness(annotator):
    """
    Correction 2 (CIRCUIT COMPLETENESS): circuit rule with any missing gene
    returns NOT_TESTABLE, not partial PASS.
    """
    # Treg circuit requires FOXP3 AND IL2RA
    # If only FOXP3 present but IL2RA absent → NOT_TESTABLE (not partial PASS)
    evidence = EvidenceBundle(
        cluster_id="test-circuit-incomplete",
        marker_genes={
            "FOXP3": 3.0,   # present
            # IL2RA deliberately absent
        },
        regulons={"FOXP3": 2.5},
    )

    verdicts = annotator.annotate_cluster(evidence)
    treg_firings = verdicts["Treg"].rule_firings

    circuit_firings = [f for f in treg_firings if f.rule_id == "TREG_CIRCUIT_001"]
    assert circuit_firings, "TREG_CIRCUIT_001 must be in firings"
    assert circuit_firings[0].verdict == RuleVerdict.NOT_TESTABLE, (
        f"Incomplete circuit must be NOT_TESTABLE, got {circuit_firings[0].verdict}: "
        f"{circuit_firings[0].notes}"
    )


def test_correction_3_support_counting(annotator):
    """
    Correction 3: negative rule PASS must NOT add to support_count.
    """
    # ILC1 evidence: TBX21 active (positive), TRAC absent (negative PASS)
    evidence = EvidenceBundle(
        cluster_id="test-support-count",
        marker_genes={
            "TRAC": -1.5,  # absent → ILC1_NEG_001 PASS (not contradicted)
        },
        regulons={
            "TBX21": 2.5,  # active → ILC1_POS_001 PASS
        },
    )

    verdicts = annotator.annotate_cluster(evidence)
    v = verdicts["ILC1"]

    # support_count should only include positive + circuit passes, NOT negative PASS
    positive_circuit_passes = sum(
        1 for f in v.rule_firings
        if f.rule_type in ("positive", "circuit") and f.verdict == RuleVerdict.PASS
    )
    assert v.support_count == positive_circuit_passes, (
        f"support_count={v.support_count} must equal positive/circuit PASSes={positive_circuit_passes}"
    )

    # The negative rule PASS should NOT be counted
    negative_passes = [
        f for f in v.rule_firings
        if f.rule_type == "negative" and f.verdict == RuleVerdict.PASS
    ]
    # Verify negative PASS is not included in support_count
    if negative_passes:
        total_all_passes = sum(
            1 for f in v.rule_firings if f.verdict == RuleVerdict.PASS
        )
        assert v.support_count < total_all_passes, (
            "Negative rule PASS must not be counted in support_count"
        )


def test_correction_4_proven_requires_positive_pass(annotator):
    """
    Correction 4: circuit PASS alone is NOT sufficient for PROVEN.
    Need at least one positive PASS too (if positive rules exist).
    """
    # ILC1 has: ILC1_CIRCUIT_001 (TBX21 circuit) + ILC1_POS_001 (TBX21 positive)
    # If we make TBX21 NOT active, circuit FAILS and positive FAILS → UNCERTAIN
    # This test checks: circuit PASS + no positive PASS → NOT PROVEN (UNCERTAIN)

    # ILC1: ILC1_NEG_001 (TRAC) + ILC1_POS_001 (TBX21 regulon) + ILC1_CIRCUIT_001 (TBX21 circuit)
    # ILC3 has: circuit (RORC+NCR2) + positive (RORC) + negative rules
    # To test correction 4: need a cell type whose circuit fires but positive doesn't.
    # NK has NK_CIRCUIT_001 (EOMES circuit) and NK_POS_001 (NCAM1 marker).
    # Provide EOMES active but NO NCAM1 → circuit PASS but no positive PASS → UNCERTAIN (not PROVEN)

    evidence = EvidenceBundle(
        cluster_id="test-circuit-no-positive",
        marker_genes={
            "TRAC": -1.5,   # absent → NK_NEG_001 PASS (not contradicted)
        },
        regulons={
            "EOMES": 2.5,  # active → NK_CIRCUIT_001 PASS
            # NCAM1 not in regulons; check marker:
        },
    )
    # Note: NK_POS_001 checks NCAM1 in marker_genes, which is absent → FAIL
    # So: circuit PASS (EOMES) but no positive PASS (NCAM1 absent) → correction 4 → NOT PROVEN

    verdicts = annotator.annotate_cluster(evidence)
    v = verdicts["NK"]

    circuit_passes = [f for f in v.rule_firings if f.rule_type == "circuit" and f.verdict == RuleVerdict.PASS]
    positive_passes = [f for f in v.rule_firings if f.rule_type == "positive" and f.verdict == RuleVerdict.PASS]

    if circuit_passes and not positive_passes:
        # Correction 4: circuit PASS but no positive PASS → UNCERTAIN (not PROVEN)
        assert v.verdict != Verdict.PROVEN, (
            f"NK must not be PROVEN when circuit PASS but no positive PASS. "
            f"Got {v.verdict}. Correction 4 violated."
        )


# ── KGLoader tests ────────────────────────────────────────────────────────────

def test_kg_loads_correctly(kg):
    """KG loads, validates against schema, and has the correct minimum rule count."""
    assert kg.rule_count("ACTIVE") > 0
    # KG grows across days; assert the Phase 1 baseline of 18 types is present
    assert len(kg.active_cell_types()) >= 18


def test_kg_active_cell_types(kg):
    """All 18 Phase 1 cell types remain present in the KG after expansion."""
    # Day 3 adds 9 new T-cell subtypes; original 18 must still be present
    phase1_types = {
        "Treg", "ILC1", "ILC2", "ILC3", "NK", "CD4_T", "CD8_T",
        "pDC", "cDC1", "cDC2", "Monocyte", "Macrophage",
        "B_cell", "Plasma_cell", "mTEC", "cTEC", "NKT", "Neutrophil",
    }
    actual = set(kg.active_cell_types())
    missing = phase1_types - actual
    assert not missing, (
        f"Phase 1 cell types missing after Day 3 expansion: {missing}"
    )


def test_kg_all_rules_have_pmid(kg):
    """Every ACTIVE rule must have a non-empty pmid (citation required)."""
    for rule in kg:
        if rule.status == "ACTIVE":
            assert rule.pmid, f"Rule {rule.rule_id} missing pmid"


def test_kg_tissue_filtering(kg, annotator):
    """Tissue filter restricts rules to matching + universal tissue rules."""
    # thymus-specific rules: mTEC, cTEC have tissue_context=["thymus"]
    thymus_rules = kg.active_rules(cell_type="mTEC", tissue="thymus")
    assert len(thymus_rules) > 0, "mTEC must have thymus-context rules"

    # lung tissue: thymus-specific rules should be excluded
    lung_rules = kg.active_rules(cell_type="mTEC", tissue="lung")
    assert len(lung_rules) == 0, "mTEC thymus rules must not appear when tissue=lung"


# ── EvidenceBundle direct method coverage ────────────────────────────────────

def test_evidence_get_marker_and_get_regulon():
    """get_marker / get_regulon return z-scores or None; repr works."""
    eb = EvidenceBundle(
        cluster_id="test-direct",
        marker_genes={"CD3E": 2.5, "CD4": -1.0},
        regulons={"TBX21": 3.2},
    )
    # get_marker
    assert eb.get_marker("CD3E") == 2.5
    assert eb.get_marker("CD4") == -1.0
    assert eb.get_marker("NOTHERE") is None

    # get_regulon
    assert eb.get_regulon("TBX21") == 3.2
    assert eb.get_regulon("FOXP3") is None

    # __repr__
    r = repr(eb)
    assert "test-direct" in r
    assert "markers=2" in r


def test_evidence_satisfies_direction_all_directions():
    """satisfies_direction covers high/present/active and low/absent/inactive."""
    eb = EvidenceBundle(cluster_id="test-direction", marker_genes={}, regulons={})
    # high / present / active — all mean z >= Z_ACTIVE_THRESH (1.5)
    assert eb.satisfies_direction(2.0, "high") is True
    assert eb.satisfies_direction(0.5, "high") is False
    assert eb.satisfies_direction(1.5, "present") is True
    assert eb.satisfies_direction(2.0, "active") is True

    # low / absent / inactive — all mean z <= Z_INACTIVE_THRESH (-0.5)
    assert eb.satisfies_direction(-1.0, "low") is True
    assert eb.satisfies_direction(0.0, "low") is False
    assert eb.satisfies_direction(-0.6, "absent") is True
    assert eb.satisfies_direction(-1.5, "inactive") is True

    # Unknown direction raises
    with pytest.raises(ValueError, match="Unknown direction"):
        eb.satisfies_direction(1.0, "sideways")


def test_evidence_lookup_any_prefers_regulon():
    """lookup_any returns regulon z-score when gene appears in both dicts."""
    eb = EvidenceBundle(
        cluster_id="test-lookup-any",
        marker_genes={"FOXP3": 1.0},   # lower z in marker genes
        regulons={"FOXP3": 3.5},        # higher z in regulons
    )
    # lookup_any should prefer regulon
    val = eb.lookup_any("FOXP3")
    assert val == 3.5, f"lookup_any must prefer regulon z-score, got {val}"
