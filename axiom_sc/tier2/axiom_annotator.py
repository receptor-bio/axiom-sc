"""
AxiomAnnotator: proof-by-contradiction cell type annotation via mechanistic KG.

Core principle: one hard rule violation eliminates a candidate cell type regardless
of how many positive markers support it. This differs from all LLM-based systems
(CASSIA, mLLMCelltype) which only use positive marker matching.

Phase 1 corrections — ALL preserved here (do not revert):

  1. FULL KG TEST MODE
     Always test all KG cell types, not just CASSIA suggestions.
     CASSIA free-text labels fail to map to KG keys. Exhaustive testing completes
     in <3 min (n_clusters × n_kg_types × n_rules ≈ 73,800 evaluations).

  2. CIRCUIT COMPLETENESS
     If any component gene is not in evidence bundle, circuit rule returns
     NOT_TESTABLE (not partial PASS).

  3. SUPPORT COUNTING REFORM
     Only positive/circuit PASS counts as support.
     Negative rule PASS = "not contradicted" — not affirmative evidence.

  4. PROVEN CO-CONFIRMATION
     Circuit PASS alone is not sufficient for PROVEN.
     Requires at least one positive rule PASS (if positive rules exist for type).

  5. PDC_NEG_001 CORRECTION
     Check PAX5 REGULON z-score < -0.5 (not IGKC marker absence).
     pDCs naturally express IGKC (developmental B-lineage relationship).
     This is encoded in the KG JSON, not in code.

  6. ILC3_CIRCUIT_001 PHASE 2 FIX
     Requires NCR2/NCR3 co-expression with RORC.
     Non-immune RORC (circadian) is excluded because NCR2/NCR3 absent in fibroblasts.
     This is encoded in the KG JSON, not in code.

License: Apache 2.0
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional

from axiom_sc.tier2.evidence import EvidenceBundle
from axiom_sc.tier2.kg_loader import KGLoader, KGRule


# ── enums and result types ────────────────────────────────────────────────────

class RuleVerdict(str, Enum):
    PASS = "PASS"
    FAIL = "FAIL"
    NOT_TESTABLE = "NOT_TESTABLE"


class Verdict(str, Enum):
    PROVEN = "PROVEN"
    CONTRADICTED = "CONTRADICTED"
    UNCERTAIN = "UNCERTAIN"
    NOT_TESTABLE = "NOT_TESTABLE"   # no rules were testable at all


@dataclass
class RuleFiring:
    """Result of firing one KG rule against one evidence bundle."""

    rule_id: str
    cell_type: str
    rule_type: str       # positive | negative | circuit | spatial
    verdict: RuleVerdict
    genes_tested: list[str]
    gene_values: dict[str, float | None]
    mechanistic_basis: str
    pmid: str
    notes: str = ""      # human-readable explanation of the verdict


@dataclass
class CellTypeVerdict:
    """
    Aggregated verdict for one candidate cell type against one evidence bundle.
    """

    cell_type: str
    verdict: Verdict
    support_count: int          # positive + circuit PASSes (correction 3)
    contradiction_count: int    # negative rule FAILes
    rule_firings: list[RuleFiring]
    confidence: float = 0.0     # 0–1 score derived from support_count

    def best_positive_firings(self) -> list[RuleFiring]:
        return [
            f for f in self.rule_firings
            if f.rule_type in ("positive", "circuit") and f.verdict == RuleVerdict.PASS
        ]

    def failed_negative_firings(self) -> list[RuleFiring]:
        return [
            f for f in self.rule_firings
            if f.rule_type == "negative" and f.verdict == RuleVerdict.FAIL
        ]


# ── main annotator ────────────────────────────────────────────────────────────

class AxiomAnnotator:
    """
    Phase 1-validated proof-by-contradiction annotator.

    Parameters
    ----------
    kg : KGLoader
        Loaded knowledge graph.
    z_active_thresh : float
        z-score threshold above which a gene/regulon is considered active/high/present.
        Phase 1 default: 1.5
    z_inactive_thresh : float
        z-score threshold below which a gene/regulon is considered inactive/low/absent.
        Phase 1 default: -0.5

    Usage::

        kg = KGLoader.from_file("kg_data/oracle_kg_v0.2.0.json")
        annotator = AxiomAnnotator(kg=kg)
        verdicts = annotator.annotate_cluster(evidence_bundle)
        best = annotator.get_best_verdict(evidence_bundle)
    """

    def __init__(
        self,
        kg: KGLoader,
        z_active_thresh: float = 1.5,
        z_inactive_thresh: float = -0.5,
    ):
        self.kg = kg
        self.z_active_thresh = z_active_thresh
        self.z_inactive_thresh = z_inactive_thresh

    # ── public interface ──────────────────────────────────────────────────────

    def annotate_cluster(
        self,
        evidence: EvidenceBundle,
        tissue: str | None = None,
    ) -> dict[str, CellTypeVerdict]:
        """
        Correction 1 (FULL KG TEST MODE): test ALL cell types in the KG.
        Returns a dict of cell_type → CellTypeVerdict for every active cell type.

        Parameters
        ----------
        evidence : EvidenceBundle
        tissue : optional tissue filter (restricts KG rules to tissue-universal
                 and tissue-matching rules)
        """
        tissue = tissue or evidence.tissue or None
        cell_types = self.kg.active_cell_types()

        verdicts: dict[str, CellTypeVerdict] = {}
        for cell_type in cell_types:
            rules = self.kg.active_rules(cell_type=cell_type, tissue=tissue)
            firings = [self._fire_rule(rule, evidence) for rule in rules]
            verdicts[cell_type] = self._compute_verdict(cell_type, firings)

        return verdicts

    def get_best_verdict(
        self,
        evidence: EvidenceBundle,
        tissue: str | None = None,
    ) -> CellTypeVerdict:
        """
        Returns the single best CellTypeVerdict for this evidence bundle.

        Priority:
          1. PROVEN types (highest support_count first)
          2. UNCERTAIN types (highest support_count first)
          3. NOT_TESTABLE types (if nothing else)
          4. CONTRADICTED types (should never be best, returned only if all CONTRADICTED)
        """
        verdicts = self.annotate_cluster(evidence, tissue=tissue)
        return self._rank_verdicts(list(verdicts.values()))

    # ── rule firing ───────────────────────────────────────────────────────────

    def _fire_rule(self, rule: KGRule, evidence: EvidenceBundle) -> RuleFiring:
        if rule.rule_type == "circuit":
            return self._fire_circuit_rule(rule, evidence)
        elif rule.rule_type in ("positive", "negative"):
            return self._fire_unary_rule(rule, evidence)
        else:
            # spatial, trajectory — stub: NOT_TESTABLE until Tier 3 hooked in
            return RuleFiring(
                rule_id=rule.rule_id,
                cell_type=rule.cell_type,
                rule_type=rule.rule_type,
                verdict=RuleVerdict.NOT_TESTABLE,
                genes_tested=rule.gene_or_regulon,
                gene_values={g: None for g in rule.gene_or_regulon},
                mechanistic_basis=rule.mechanistic_basis,
                pmid=rule.pmid,
                notes=f"rule_type={rule.rule_type!r} not yet supported at Tier 2",
            )

    def _fire_unary_rule(self, rule: KGRule, evidence: EvidenceBundle) -> RuleFiring:
        """Fire a positive or negative rule (single gene or regulon)."""
        # For multi-gene positive/negative rules, all genes are evaluated individually;
        # any single gene's result is sufficient to trigger the rule.
        # (In practice, positive/negative rules have exactly 1 gene.)
        gene = rule.gene_or_regulon[0]
        value = evidence.lookup(gene, rule.evidence_source)

        if value is None:
            return RuleFiring(
                rule_id=rule.rule_id,
                cell_type=rule.cell_type,
                rule_type=rule.rule_type,
                verdict=RuleVerdict.NOT_TESTABLE,
                genes_tested=[gene],
                gene_values={gene: None},
                mechanistic_basis=rule.mechanistic_basis,
                pmid=rule.pmid,
                notes=f"{gene} not in evidence bundle ({rule.evidence_source})",
            )

        condition_met = self._satisfies_direction(value, rule.direction)

        if rule.rule_type == "positive":
            # PASS = gene is in the expected direction (supports this type)
            verdict = RuleVerdict.PASS if condition_met else RuleVerdict.FAIL
            notes = (
                f"{gene}={value:.2f} ≥ {self.z_active_thresh} ({rule.direction})"
                if condition_met
                else f"{gene}={value:.2f} does not satisfy {rule.direction} threshold"
            )
        else:
            # negative rule:
            # condition_met = the "bad" state IS present → FAIL → CONTRADICTED
            # condition NOT met → PASS (not contradicted)
            verdict = RuleVerdict.FAIL if condition_met else RuleVerdict.PASS
            notes = (
                f"CONTRADICTION: {gene}={value:.2f} satisfies {rule.direction} (should not for {rule.cell_type})"
                if condition_met
                else f"{gene}={value:.2f} does not trigger contradiction"
            )

        return RuleFiring(
            rule_id=rule.rule_id,
            cell_type=rule.cell_type,
            rule_type=rule.rule_type,
            verdict=verdict,
            genes_tested=[gene],
            gene_values={gene: value},
            mechanistic_basis=rule.mechanistic_basis,
            pmid=rule.pmid,
            notes=notes,
        )

    def _fire_circuit_rule(self, rule: KGRule, evidence: EvidenceBundle) -> RuleFiring:
        """
        Fire a circuit rule (co-expression requirement).

        Correction 2 (CIRCUIT COMPLETENESS): if ANY component gene is missing
        from the evidence bundle, the whole circuit returns NOT_TESTABLE.
        This prevents partial PASS verdicts from incomplete evidence.
        """
        gene_values: dict[str, float | None] = {}

        for gene in rule.gene_or_regulon:
            # Circuit rules may span both regulons and marker_genes;
            # lookup_any() checks regulons first, then marker_genes.
            val = evidence.lookup_any(gene)
            gene_values[gene] = val
            if val is None:
                # Correction 2: any missing gene → NOT_TESTABLE for whole circuit
                return RuleFiring(
                    rule_id=rule.rule_id,
                    cell_type=rule.cell_type,
                    rule_type=rule.rule_type,
                    verdict=RuleVerdict.NOT_TESTABLE,
                    genes_tested=list(gene_values.keys()),
                    gene_values=gene_values,
                    mechanistic_basis=rule.mechanistic_basis,
                    pmid=rule.pmid,
                    notes=(
                        f"Circuit NOT_TESTABLE: {gene!r} not in evidence "
                        f"(requires all of {rule.gene_or_regulon})"
                    ),
                )

        # All genes present — check all satisfy direction
        all_satisfied = all(
            self._satisfies_direction(v, rule.direction)
            for v in gene_values.values()
            if v is not None
        )

        verdict = RuleVerdict.PASS if all_satisfied else RuleVerdict.FAIL
        gene_summary = ", ".join(
            f"{g}={v:.2f}" for g, v in gene_values.items() if v is not None
        )
        notes = (
            f"Circuit PASS: all of [{gene_summary}] satisfy {rule.direction}"
            if all_satisfied
            else f"Circuit FAIL: not all of [{gene_summary}] satisfy {rule.direction}"
        )

        return RuleFiring(
            rule_id=rule.rule_id,
            cell_type=rule.cell_type,
            rule_type=rule.rule_type,
            verdict=verdict,
            genes_tested=list(gene_values.keys()),
            gene_values=gene_values,
            mechanistic_basis=rule.mechanistic_basis,
            pmid=rule.pmid,
            notes=notes,
        )

    # ── verdict aggregation ───────────────────────────────────────────────────

    def _compute_verdict(
        self, cell_type: str, firings: list[RuleFiring]
    ) -> CellTypeVerdict:
        """
        Aggregate individual rule firings into a single verdict for a cell type.

        Corrections applied:
          3. Only positive/circuit PASS counts as support (not negative PASS).
          4. PROVEN requires: circuit_passes AND (positive_passes OR no positive rules).
        """
        if not firings:
            return CellTypeVerdict(
                cell_type=cell_type,
                verdict=Verdict.NOT_TESTABLE,
                support_count=0,
                contradiction_count=0,
                rule_firings=[],
            )

        # Correction 4 bookkeeping
        has_positive_rules = any(f.rule_type == "positive" for f in firings)
        positive_passes = [
            f for f in firings
            if f.rule_type == "positive" and f.verdict == RuleVerdict.PASS
        ]
        circuit_passes = [
            f for f in firings
            if f.rule_type == "circuit" and f.verdict == RuleVerdict.PASS
        ]
        negative_fails = [
            f for f in firings
            if f.rule_type == "negative" and f.verdict == RuleVerdict.FAIL
        ]

        # Correction 3: support = positive PASS + circuit PASS only
        support_count = len(positive_passes) + len(circuit_passes)

        # One contradiction eliminates the candidate (proof-by-contradiction)
        if negative_fails:
            return CellTypeVerdict(
                cell_type=cell_type,
                verdict=Verdict.CONTRADICTED,
                support_count=support_count,
                contradiction_count=len(negative_fails),
                rule_firings=firings,
                confidence=0.0,
            )

        # Correction 4: PROVEN requires circuit PASS + (positive PASS OR no positive rules)
        if circuit_passes and (positive_passes or not has_positive_rules):
            confidence = min(1.0, support_count / max(1, len(firings)))
            return CellTypeVerdict(
                cell_type=cell_type,
                verdict=Verdict.PROVEN,
                support_count=support_count,
                contradiction_count=0,
                rule_firings=firings,
                confidence=confidence,
            )

        # Check if ALL rules were NOT_TESTABLE
        all_not_testable = all(f.verdict == RuleVerdict.NOT_TESTABLE for f in firings)
        if all_not_testable:
            return CellTypeVerdict(
                cell_type=cell_type,
                verdict=Verdict.NOT_TESTABLE,
                support_count=0,
                contradiction_count=0,
                rule_firings=firings,
            )

        confidence = min(1.0, support_count / max(1, len(firings))) if support_count else 0.0
        return CellTypeVerdict(
            cell_type=cell_type,
            verdict=Verdict.UNCERTAIN,
            support_count=support_count,
            contradiction_count=0,
            rule_firings=firings,
            confidence=confidence,
        )

    @staticmethod
    def _rank_verdicts(verdicts: list[CellTypeVerdict]) -> CellTypeVerdict:
        """Return best verdict: PROVEN > UNCERTAIN > NOT_TESTABLE > CONTRADICTED."""
        _priority = {
            Verdict.PROVEN: 0,
            Verdict.UNCERTAIN: 1,
            Verdict.NOT_TESTABLE: 2,
            Verdict.CONTRADICTED: 3,
        }
        return min(
            verdicts,
            key=lambda v: (_priority[v.verdict], -v.support_count),
        )

    # ── helpers ───────────────────────────────────────────────────────────────

    def _satisfies_direction(self, value: float, direction: str) -> bool:
        if direction in ("high", "present", "active"):
            return value >= self.z_active_thresh
        if direction in ("low", "absent", "inactive"):
            return value <= self.z_inactive_thresh
        raise ValueError(f"Unknown direction: {direction!r}")
