"""
Tier 5 novel attractor characterization.

Cells that fail all four annotation tiers are not discarded — they are characterized
as candidate novel cell states. This module extracts a structured CandidateState
from the available multi-omic evidence for a failing cluster.

The approach:
  1. Extract active TF circuits from SCENIC regulon z-scores (evidence.regulons)
  2. Identify nearest known velocity attractor from Tier 3 stream results
  3. Capture spatial niche context from spatial stream results
  4. Record cross-species conservation from OrthoFinder stream
  5. Query Cell Ontology by nearest marker distance (simple string heuristic)
  6. Generate falsifiable perturbation predictions from active TF circuits

Reference (attractor concept):
  Huang S. (2012) The molecular and mathematical basis of Waddington's epigenetic landscape.
  BioEssays 34:149. doi: 10.1002/bies.201100031

License: Apache 2.0
"""
from __future__ import annotations

from axiom_sc.tier2.evidence import EvidenceBundle, Z_ACTIVE_THRESH
from axiom_sc.tier5.candidate_state import CandidateState, PerturbationPrediction

# Minimum regulon z-score to call a TF circuit "active"
CIRCUIT_ACTIVE_Z: float = Z_ACTIVE_THRESH  # 1.5, same as Tier 2

# Heuristic Cell Ontology lookup — maps dominant active circuit pattern to ontology term.
# This is a first-pass approximation; Tier 4 LLM followup refines it.
_CIRCUIT_TO_ONTOLOGY: dict[frozenset, str] = {
    frozenset({"FOXP3"}):              "CL:0000792 — CD4-positive, CD25-positive, alpha-beta regulatory T cell",
    frozenset({"TBX21"}):              "CL:0000965 — Th1 cell",
    frozenset({"GATA3"}):              "CL:0000546 — Th2 cell",
    frozenset({"RORC"}):               "CL:0000899 — Th17 cell",
    frozenset({"BCL6"}):               "CL:0000844 — germinal center B cell",
    frozenset({"IRF7"}):               "CL:0000784 — plasmacytoid dendritic cell",
    frozenset({"AIRE"}):               "CL:0002363 — medullary thymic epithelial cell",
    frozenset({"FOXJ1"}):              "CL:0000067 — ciliated epithelial cell",
    frozenset({"TP63"}):               "CL:0000646 — basal cell",
    frozenset({"NEUROD1"}):            "CL:0000164 — enteroendocrine cell",
    frozenset({"SOX9"}):               "CL:1000497 — kidney collecting duct principal cell",
    frozenset({"TBX21", "EOMES"}):     "CL:0000625 — CD8-positive, alpha-beta T cell",
    frozenset({"RORC", "GATA3"}):      "CL:0001071 — group 2 innate lymphoid cell",
}


class AttractorDiscovery:
    """
    Characterizes novel cell attractors from clusters failing Tiers 1–4.

    Parameters
    ----------
    dataset_id : str
        Identifier for the dataset being analysed.
    """

    def __init__(self, dataset_id: str):
        self.dataset_id = dataset_id

    def characterize(
        self,
        evidence: EvidenceBundle,
        tier3_stream_results: list | None = None,
    ) -> CandidateState:
        """
        Build a CandidateState from the available evidence for one failing cluster.

        Parameters
        ----------
        evidence : EvidenceBundle
            Tier 2 evidence for the cluster.
        tier3_stream_results : list[StreamResult], optional
            Results from Tier 3 streams (velocity, chromatin, communication,
            spatial_niche, cross_species, spatial_type).
        """
        t3 = tier3_stream_results or []

        active_circuits = self._extract_active_circuits(evidence)
        velocity_sink = self._extract_velocity_sink(t3)
        spatial_context = self._extract_spatial_context(evidence, t3)
        cross_species = self._extract_cross_species(t3)
        ontology_term = self._nearest_ontology_term(active_circuits)
        perturbations = self._generate_perturbation_predictions(active_circuits)

        evidence_summary = {
            "n_markers": len(evidence.marker_genes),
            "n_regulons": len(evidence.regulons),
            "top_markers": sorted(
                evidence.marker_genes.items(), key=lambda x: -x[1]
            )[:10],
            "active_circuits_raw": active_circuits,
            "t3_streams_run": [r.stream for r in t3 if r.available],
        }

        return CandidateState(
            cluster_id=evidence.cluster_id,
            dataset_id=self.dataset_id,
            active_circuits=active_circuits,
            velocity_sink=velocity_sink,
            spatial_context=spatial_context,
            cross_species=cross_species,
            nearest_ontology_term=ontology_term,
            perturbation_predictions=perturbations,
            status="CANDIDATE",
            evidence_summary=evidence_summary,
        )

    # ── private helpers ───────────────────────────────────────────────────────

    def _extract_active_circuits(self, evidence: EvidenceBundle) -> list[str]:
        """Return TF names whose regulon z-score is ≥ CIRCUIT_ACTIVE_Z, sorted descending."""
        active = [
            tf for tf, z in evidence.regulons.items()
            if z >= CIRCUIT_ACTIVE_Z
        ]
        return sorted(active, key=lambda tf: -evidence.regulons[tf])

    def _extract_velocity_sink(self, t3: list) -> str | None:
        """Extract the nearest velocity attractor label from the velocity stream."""
        for r in t3:
            if r.stream == "velocity" and r.available and r.label:
                return r.label
        return None

    def _extract_spatial_context(
        self, evidence: EvidenceBundle, t3: list
    ) -> str:
        """Build spatial context string from evidence and spatial_niche stream."""
        parts = []
        if evidence.spatial_context:
            ctx = evidence.spatial_context
            if isinstance(ctx, dict):
                if "niche_label" in ctx:
                    parts.append(f"niche={ctx['niche_label']}")
                if "tissue" in ctx:
                    parts.append(f"tissue={ctx['tissue']}")
            else:
                parts.append(str(ctx))

        for r in t3:
            if r.stream == "spatial_niche" and r.available and r.label:
                parts.append(f"spatial_niche={r.label}(conf={r.confidence:.2f})")
                break

        return "; ".join(parts) if parts else "unknown"

    def _extract_cross_species(self, t3: list) -> dict[str, float]:
        """Extract species conservation scores from cross_species stream evidence."""
        for r in t3:
            if r.stream == "cross_species" and r.available:
                scores = r.evidence.get("conservation_scores", {})
                if scores:
                    return dict(scores)
        return {}

    def _nearest_ontology_term(self, active_circuits: list[str]) -> str:
        """Heuristic: match the active TF set to the closest Cell Ontology term."""
        circuit_set = frozenset(active_circuits)
        best_term = ""
        best_overlap = 0
        for ref_set, term in _CIRCUIT_TO_ONTOLOGY.items():
            overlap = len(circuit_set & ref_set)
            if overlap > best_overlap:
                best_overlap = overlap
                best_term = term
        return best_term or "CL:0000000 — cell (no match found)"

    def _generate_perturbation_predictions(
        self, active_circuits: list[str]
    ) -> list[PerturbationPrediction]:
        """
        Generate KO predictions for each active TF circuit.

        Heuristic: KO of an active master TF should dissolve the attractor state.
        All predictions are marked as low confidence until experimental validation.
        """
        predictions = []
        for tf in active_circuits[:5]:  # limit to top 5 circuits
            predictions.append(
                PerturbationPrediction(
                    gene=tf,
                    perturbation="KO",
                    predicted_effect=(
                        f"Knockout of {tf} is predicted to dissolve the novel attractor "
                        f"state if {tf} drives the transcriptional programme of this "
                        f"cluster. Expected outcome: cells re-route to the nearest "
                        f"known attractor."
                    ),
                    confidence=0.3,
                )
            )
        return predictions
