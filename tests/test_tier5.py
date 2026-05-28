"""
Tests for Tier 5: novel attractor discovery.

Covers:
  - CandidateState: dataclass, candidate_id hash, to_dict/from_dict, status validation
  - AttractorDiscovery: active circuit extraction, velocity/spatial/cross-species parsing,
    ontology matching, perturbation prediction generation
  - DiscoveryFeedback: generate_candidates rules, write_candidates JSON output

License: Apache 2.0
"""
from __future__ import annotations

import json
import tempfile
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from axiom_sc.tier5.candidate_state import CandidateState, PerturbationPrediction
from axiom_sc.tier5.attractor_discovery import AttractorDiscovery, CIRCUIT_ACTIVE_Z
from axiom_sc.tier5.feedback import DiscoveryFeedback, KGCandidate


# ── Helpers ──────────────────────────────────────────────────────────────────

def _make_evidence(
    cluster_id: str = "cluster_007",
    marker_genes: dict | None = None,
    regulons: dict | None = None,
    tissue: str = "thymus",
    spatial_context=None,
):
    """Build a minimal EvidenceBundle-like object for testing."""
    from axiom_sc.tier2.evidence import EvidenceBundle
    return EvidenceBundle(
        cluster_id=cluster_id,
        marker_genes=marker_genes or {"CD3E": 2.5, "FOXP3": 3.1, "CD4": 2.0},
        regulons=regulons or {"FOXP3": 2.8, "TBX21": 0.3, "RORC": 0.1},
        tissue=tissue,
        spatial_context=spatial_context,
    )


def _make_stream_result(stream: str, label: str = "", confidence: float = 0.8,
                        available: bool = True, evidence: dict | None = None):
    """Build a minimal StreamResult-like object."""
    from axiom_sc.tier3.stream_result import StreamResult
    return StreamResult(
        cluster_id="cluster_007",
        stream=stream,
        label=label,
        confidence=confidence,
        available=available,
        evidence=evidence or {},
    )


def _make_validated_state(
    cluster_id: str = "cluster_007",
    dataset_id: str = "thymus_2024",
    active_circuits: list | None = None,
    status: str = "VALIDATED",
    proposed_name: str = "Novel_Treg_Like",
    spatial_context: str = "niche=medulla",
) -> CandidateState:
    return CandidateState(
        cluster_id=cluster_id,
        dataset_id=dataset_id,
        active_circuits=active_circuits or ["FOXP3", "IKZF2"],
        velocity_sink="Treg",
        spatial_context=spatial_context,
        cross_species={"mouse": 0.85},
        nearest_ontology_term="CL:0000792 — CD4-positive, CD25-positive, alpha-beta regulatory T cell",
        perturbation_predictions=[],
        status=status,
        proposed_name=proposed_name,
    )


# ── CandidateState ────────────────────────────────────────────────────────────

class TestCandidateState:
    def test_candidate_id_is_deterministic(self):
        s1 = CandidateState(cluster_id="c1", dataset_id="ds1", status="CANDIDATE")
        s2 = CandidateState(cluster_id="c1", dataset_id="ds1", status="CANDIDATE")
        assert s1.candidate_id == s2.candidate_id

    def test_candidate_id_differs_for_different_inputs(self):
        s1 = CandidateState(cluster_id="c1", dataset_id="ds1", status="CANDIDATE")
        s2 = CandidateState(cluster_id="c2", dataset_id="ds1", status="CANDIDATE")
        assert s1.candidate_id != s2.candidate_id

    def test_candidate_id_is_16_chars(self):
        s = CandidateState(cluster_id="c1", dataset_id="ds1", status="CANDIDATE")
        assert len(s.candidate_id) == 16

    def test_invalid_status_raises(self):
        with pytest.raises(ValueError, match="Invalid status"):
            CandidateState(cluster_id="c1", dataset_id="ds1", status="UNKNOWN")  # type: ignore

    def test_valid_statuses(self):
        for status in ("CANDIDATE", "VALIDATED", "PUBLISHED"):
            s = CandidateState(cluster_id="c1", dataset_id="ds1", status=status)
            assert s.status == status

    def test_to_dict_has_required_keys(self):
        s = _make_validated_state()
        d = s.to_dict()
        for key in ("candidate_id", "cluster_id", "dataset_id", "active_circuits",
                    "velocity_sink", "spatial_context", "cross_species",
                    "nearest_ontology_term", "perturbation_predictions",
                    "status", "evidence_summary", "created_at", "proposed_name"):
            assert key in d, f"Missing key: {key}"

    def test_to_dict_active_circuits(self):
        s = _make_validated_state(active_circuits=["FOXP3", "IKZF2"])
        d = s.to_dict()
        assert d["active_circuits"] == ["FOXP3", "IKZF2"]

    def test_from_dict_roundtrip(self):
        original = _make_validated_state()
        original.perturbation_predictions = [
            PerturbationPrediction(
                gene="FOXP3", perturbation="KO",
                predicted_effect="Dissolves Treg attractor", confidence=0.3
            )
        ]
        d = original.to_dict()
        restored = CandidateState.from_dict(d)
        assert restored.cluster_id == original.cluster_id
        assert restored.dataset_id == original.dataset_id
        assert restored.active_circuits == original.active_circuits
        assert restored.status == original.status
        assert restored.proposed_name == original.proposed_name
        assert len(restored.perturbation_predictions) == 1
        assert restored.perturbation_predictions[0].gene == "FOXP3"
        assert restored.perturbation_predictions[0].perturbation == "KO"
        assert restored.perturbation_predictions[0].confidence == pytest.approx(0.3)

    def test_summary_line_format(self):
        s = _make_validated_state()
        line = s.summary_line()
        assert s.candidate_id in line
        assert "VALIDATED" in line
        assert "FOXP3" in line

    def test_created_at_is_utc_iso(self):
        s = CandidateState(cluster_id="c1", dataset_id="ds1", status="CANDIDATE")
        assert "T" in s.created_at
        assert "+" in s.created_at or "Z" in s.created_at or s.created_at.endswith("+00:00")

    def test_perturbation_prediction_defaults(self):
        p = PerturbationPrediction(gene="FOXP3", perturbation="KO", predicted_effect="test")
        assert p.confidence == 0.5


# ── AttractorDiscovery ────────────────────────────────────────────────────────

class TestAttractorDiscovery:
    def setup_method(self):
        self.discovery = AttractorDiscovery(dataset_id="thymus_2024")

    def test_characterize_returns_candidate_state(self):
        ev = _make_evidence()
        state = self.discovery.characterize(ev)
        assert isinstance(state, CandidateState)

    def test_characterize_status_is_candidate(self):
        ev = _make_evidence()
        state = self.discovery.characterize(ev)
        assert state.status == "CANDIDATE"

    def test_characterize_cluster_and_dataset_ids(self):
        ev = _make_evidence(cluster_id="clus_42")
        state = self.discovery.characterize(ev)
        assert state.cluster_id == "clus_42"
        assert state.dataset_id == "thymus_2024"

    def test_active_circuits_above_threshold(self):
        ev = _make_evidence(regulons={"FOXP3": 2.8, "TBX21": 0.3, "RORC": 1.6})
        state = self.discovery.characterize(ev)
        assert "FOXP3" in state.active_circuits
        assert "RORC" in state.active_circuits
        assert "TBX21" not in state.active_circuits

    def test_active_circuits_sorted_by_z_descending(self):
        ev = _make_evidence(regulons={"FOXP3": 3.5, "RORC": 2.0, "IRF7": 4.1})
        state = self.discovery.characterize(ev)
        assert state.active_circuits == ["IRF7", "FOXP3", "RORC"]

    def test_active_circuits_empty_when_none_active(self):
        ev = _make_evidence(regulons={"FOXP3": 0.5, "TBX21": -0.3})
        state = self.discovery.characterize(ev)
        assert state.active_circuits == []

    def test_velocity_sink_extracted_from_t3(self):
        t3 = [_make_stream_result("velocity", label="Treg")]
        ev = _make_evidence()
        state = self.discovery.characterize(ev, tier3_stream_results=t3)
        assert state.velocity_sink == "Treg"

    def test_velocity_sink_none_when_unavailable(self):
        t3 = [_make_stream_result("velocity", label="Treg", available=False)]
        ev = _make_evidence()
        state = self.discovery.characterize(ev, tier3_stream_results=t3)
        assert state.velocity_sink is None

    def test_velocity_sink_none_when_no_t3(self):
        ev = _make_evidence()
        state = self.discovery.characterize(ev, tier3_stream_results=None)
        assert state.velocity_sink is None

    def test_spatial_context_from_evidence_dict(self):
        ev = _make_evidence(spatial_context={"niche_label": "cortex", "tissue": "thymus"})
        state = self.discovery.characterize(ev)
        assert "niche=cortex" in state.spatial_context
        assert "tissue=thymus" in state.spatial_context

    def test_spatial_context_from_string(self):
        ev = _make_evidence(spatial_context="periportal_zone")
        state = self.discovery.characterize(ev)
        assert "periportal_zone" in state.spatial_context

    def test_spatial_context_from_t3_spatial_niche(self):
        t3 = [_make_stream_result("spatial_niche", label="germinal_center", confidence=0.9)]
        ev = _make_evidence(spatial_context=None)
        state = self.discovery.characterize(ev, tier3_stream_results=t3)
        assert "germinal_center" in state.spatial_context
        assert "0.90" in state.spatial_context

    def test_spatial_context_unknown_when_no_evidence(self):
        ev = _make_evidence(spatial_context=None)
        state = self.discovery.characterize(ev)
        assert state.spatial_context == "unknown"

    def test_cross_species_extracted_from_t3(self):
        t3 = [_make_stream_result(
            "cross_species",
            evidence={"conservation_scores": {"mouse": 0.92, "zebrafish": 0.55}}
        )]
        ev = _make_evidence()
        state = self.discovery.characterize(ev, tier3_stream_results=t3)
        assert state.cross_species == {"mouse": 0.92, "zebrafish": 0.55}

    def test_cross_species_empty_when_no_t3(self):
        ev = _make_evidence()
        state = self.discovery.characterize(ev)
        assert state.cross_species == {}

    def test_ontology_term_foxp3_match(self):
        ev = _make_evidence(regulons={"FOXP3": 3.0})
        state = self.discovery.characterize(ev)
        assert "CL:0000792" in state.nearest_ontology_term

    def test_ontology_term_tbx21_match(self):
        ev = _make_evidence(regulons={"TBX21": 2.5})
        state = self.discovery.characterize(ev)
        assert "CL:0000965" in state.nearest_ontology_term

    def test_ontology_term_no_match_fallback(self):
        ev = _make_evidence(regulons={"UNKNOWN_TF": 3.0})
        state = self.discovery.characterize(ev)
        assert "CL:0000000" in state.nearest_ontology_term

    def test_perturbation_predictions_generated(self):
        ev = _make_evidence(regulons={"FOXP3": 3.0, "IKZF2": 2.0})
        state = self.discovery.characterize(ev)
        assert len(state.perturbation_predictions) >= 1
        assert all(p.perturbation == "KO" for p in state.perturbation_predictions)

    def test_perturbation_predictions_capped_at_5(self):
        regulons = {f"TF{i}": float(3 - i * 0.01) for i in range(10)}
        ev = _make_evidence(regulons=regulons)
        state = self.discovery.characterize(ev)
        assert len(state.perturbation_predictions) <= 5

    def test_perturbation_prediction_confidence_is_low(self):
        ev = _make_evidence(regulons={"FOXP3": 3.0})
        state = self.discovery.characterize(ev)
        for pred in state.perturbation_predictions:
            assert pred.confidence == pytest.approx(0.3)

    def test_evidence_summary_populated(self):
        ev = _make_evidence()
        state = self.discovery.characterize(ev)
        assert "n_markers" in state.evidence_summary
        assert "n_regulons" in state.evidence_summary
        assert "active_circuits_raw" in state.evidence_summary


# ── DiscoveryFeedback ─────────────────────────────────────────────────────────

class TestDiscoveryFeedback:
    def setup_method(self):
        self.feedback = DiscoveryFeedback()

    def test_candidate_state_raises_value_error(self):
        state = _make_validated_state(status="CANDIDATE")
        with pytest.raises(ValueError, match="CANDIDATE"):
            self.feedback.generate_candidates(state)

    def test_validated_state_generates_candidates(self):
        state = _make_validated_state(status="VALIDATED")
        candidates = self.feedback.generate_candidates(state)
        assert len(candidates) > 0

    def test_published_state_generates_candidates(self):
        state = _make_validated_state(status="PUBLISHED")
        candidates = self.feedback.generate_candidates(state)
        assert len(candidates) > 0

    def test_positive_rules_generated_per_tf(self):
        state = _make_validated_state(
            active_circuits=["FOXP3", "IKZF2", "GATA3"],
            status="VALIDATED",
        )
        candidates = self.feedback.generate_candidates(state)
        positive = [c for c in candidates if c.rule_type == "positive"]
        assert len(positive) == 3

    def test_circuit_rule_generated_for_two_tfs(self):
        state = _make_validated_state(
            active_circuits=["FOXP3", "IKZF2"],
            status="VALIDATED",
        )
        candidates = self.feedback.generate_candidates(state)
        circuit = [c for c in candidates if c.rule_type == "circuit"]
        assert len(circuit) == 1

    def test_no_circuit_rule_for_single_tf(self):
        state = _make_validated_state(
            active_circuits=["FOXP3"],
            status="VALIDATED",
        )
        candidates = self.feedback.generate_candidates(state)
        circuit = [c for c in candidates if c.rule_type == "circuit"]
        assert len(circuit) == 0

    def test_circuit_rule_contains_top2_tfs(self):
        state = _make_validated_state(
            active_circuits=["FOXP3", "IKZF2", "RORC"],
            status="VALIDATED",
        )
        candidates = self.feedback.generate_candidates(state)
        circuit = [c for c in candidates if c.rule_type == "circuit"][0]
        assert "FOXP3" in circuit.gene_or_regulon
        assert "IKZF2" in circuit.gene_or_regulon
        assert "RORC" not in circuit.gene_or_regulon

    def test_positive_rules_capped_at_6(self):
        state = _make_validated_state(
            active_circuits=[f"TF{i}" for i in range(10)],
            status="VALIDATED",
        )
        candidates = self.feedback.generate_candidates(state)
        positive = [c for c in candidates if c.rule_type == "positive"]
        assert len(positive) == 6

    def test_all_candidates_have_pending_review_status(self):
        state = _make_validated_state(status="VALIDATED")
        candidates = self.feedback.generate_candidates(state)
        assert all(c.status == "PENDING_REVIEW" for c in candidates)

    def test_all_candidates_have_needs_review_pmid(self):
        state = _make_validated_state(status="VALIDATED")
        candidates = self.feedback.generate_candidates(state)
        assert all(c.pmid == "NEEDS_REVIEW" for c in candidates)

    def test_all_candidates_have_low_confidence(self):
        state = _make_validated_state(status="VALIDATED")
        candidates = self.feedback.generate_candidates(state)
        assert all(c.confidence == "low" for c in candidates)

    def test_candidate_links_back_to_candidate_state(self):
        state = _make_validated_state(status="VALIDATED")
        candidates = self.feedback.generate_candidates(state)
        for c in candidates:
            assert c.candidate_id == state.candidate_id

    def test_proposed_cell_type_overrides_state_name(self):
        state = _make_validated_state(status="VALIDATED", proposed_name="Default_Name")
        candidates = self.feedback.generate_candidates(
            state, proposed_cell_type="Override_Type"
        )
        assert all(c.cell_type == "Override_Type" for c in candidates)

    def test_state_proposed_name_used_when_no_override(self):
        state = _make_validated_state(status="VALIDATED", proposed_name="Novel_Treg_Like")
        candidates = self.feedback.generate_candidates(state)
        assert all(c.cell_type == "Novel_Treg_Like" for c in candidates)

    def test_fallback_cell_type_when_no_name(self):
        state = CandidateState(
            cluster_id="c1", dataset_id="ds1",
            active_circuits=["FOXP3", "IKZF2"],
            status="VALIDATED",
            proposed_name="",
        )
        candidates = self.feedback.generate_candidates(state)
        cid_prefix = state.candidate_id[:8]
        assert all(c.cell_type == f"Novel_{cid_prefix}" for c in candidates)

    def test_spatial_context_in_tissue_context(self):
        state = _make_validated_state(
            status="VALIDATED",
            spatial_context="cortex_niche",
        )
        candidates = self.feedback.generate_candidates(state)
        assert all("cortex_niche" in c.tissue_context for c in candidates)

    def test_empty_spatial_context_gives_empty_tissue_context(self):
        state = CandidateState(
            cluster_id="c1", dataset_id="ds1",
            active_circuits=["FOXP3", "IKZF2"],
            status="VALIDATED",
            spatial_context="",
        )
        candidates = self.feedback.generate_candidates(state)
        assert all(c.tissue_context == [] for c in candidates)

    def test_write_candidates_creates_json(self):
        state = _make_validated_state(status="VALIDATED")
        candidates = self.feedback.generate_candidates(state)
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "candidates.json"
            result = self.feedback.write_candidates(candidates, output)
            assert result == output
            assert output.exists()
            data = json.loads(output.read_text())
            assert "candidates" in data
            assert "metadata" in data
            assert data["metadata"]["n_candidates"] == len(candidates)

    def test_write_candidates_append_mode(self):
        state1 = _make_validated_state(
            cluster_id="c1", dataset_id="ds1",
            active_circuits=["FOXP3", "IKZF2"],
            status="VALIDATED",
            proposed_name="TypeA",
        )
        state2 = _make_validated_state(
            cluster_id="c2", dataset_id="ds1",
            active_circuits=["TBX21", "EOMES"],
            status="VALIDATED",
            proposed_name="TypeB",
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "candidates.json"
            c1 = self.feedback.generate_candidates(state1)
            c2 = self.feedback.generate_candidates(state2)
            self.feedback.write_candidates(c1, output, append=False)
            self.feedback.write_candidates(c2, output, append=True)
            data = json.loads(output.read_text())
            all_ids = {r["rule_id"] for r in data["candidates"]}
            # Both batches present
            assert len(all_ids) == len(c1) + len(c2)

    def test_write_candidates_deduplicates_by_rule_id(self):
        state = _make_validated_state(status="VALIDATED")
        candidates = self.feedback.generate_candidates(state)
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "candidates.json"
            self.feedback.write_candidates(candidates, output, append=False)
            # Write same candidates again — should deduplicate
            self.feedback.write_candidates(candidates, output, append=True)
            data = json.loads(output.read_text())
            assert data["metadata"]["n_candidates"] == len(candidates)

    def test_write_candidates_creates_parent_dirs(self):
        state = _make_validated_state(status="VALIDATED")
        candidates = self.feedback.generate_candidates(state)
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "deep" / "nested" / "candidates.json"
            self.feedback.write_candidates(candidates, output)
            assert output.exists()
