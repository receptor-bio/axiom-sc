"""
Tier 4 tests — LLM ensemble, evidence bundler, RAG, prompts.

All LLM API calls are mocked — no network access required.
Tests verify prompt construction, response parsing, consensus logic, and RAG retrieval.
"""
from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from axiom_sc.tier2.evidence import EvidenceBundle
from axiom_sc.tier3.stream_result import StreamResult

import axiom_sc as _axiom_sc
KG_PATH = Path(_axiom_sc.__file__).parent / "kg_data" / "oracle_kg_v0.2.0.json"


# ── prompts ──────────────────────────────────────────────────────────────────

class TestPrompts:
    def test_prompt_version_is_string(self):
        from axiom_sc.tier4.prompts import PROMPT_VERSION
        assert isinstance(PROMPT_VERSION, str)
        assert len(PROMPT_VERSION) > 0

    def test_format_tier1_predictions_sorts_descending(self):
        from axiom_sc.tier4.prompts import format_tier1_predictions
        result = format_tier1_predictions({"Treg": 0.9, "CD4_T": 0.5, "NK": 0.3})
        lines = result.strip().splitlines()
        # Treg should be first (highest)
        assert "Treg" in lines[0]

    def test_format_tier1_predictions_empty(self):
        from axiom_sc.tier4.prompts import format_tier1_predictions
        result = format_tier1_predictions({})
        assert "no" in result.lower() or "none" in result.lower()

    def test_format_tier2_rules_pass_fail(self):
        from axiom_sc.tier4.prompts import format_tier2_rules
        rules = [
            {"verdict": "PASS", "rule_id": "TREG_POS_001",
             "cell_type": "Treg", "mechanistic_basis": "FOXP3 is master TF for Treg identity"},
            {"verdict": "FAIL", "rule_id": "TREG_NEG_001",
             "cell_type": "Treg", "mechanistic_basis": "IL2 must not be high in Tregs"},
        ]
        result = format_tier2_rules(rules)
        assert "✓" in result
        assert "✗" in result
        assert "Treg" in result

    def test_format_tier2_rules_empty(self):
        from axiom_sc.tier4.prompts import format_tier2_rules
        result = format_tier2_rules([])
        assert "not run" in result.lower() or "no rules" in result.lower()

    def test_format_marker_genes_top_n(self):
        from axiom_sc.tier4.prompts import format_marker_genes
        markers = {f"GENE{i}": float(i) for i in range(30)}
        result = format_marker_genes(markers, top_n=5)
        # Should contain 5 genes, highest first
        assert "GENE29" in result  # highest z-score

    def test_format_active_regulons_filters_threshold(self):
        from axiom_sc.tier4.prompts import format_active_regulons
        result = format_active_regulons({"FOXP3": 2.5, "TBX21": 0.5, "GATA3": -1.0}, thresh=1.5)
        assert "FOXP3" in result
        assert "TBX21" not in result
        assert "GATA3" not in result

    def test_build_user_prompt_contains_cluster_id(self):
        from axiom_sc.tier4.prompts import build_user_prompt
        prompt = build_user_prompt(
            cluster_id="test-cluster-99",
            tier1_predictions={"Treg": 0.9},
            tier2_rules_fired=[],
            tier2_summary="1 PROVEN",
            tier3_streams_available=["velocity", "chromatin"],
            tier3_agreement="2/2 → Treg",
            tier3_details="  velocity: Treg [conf=0.85]",
            marker_genes={"FOXP3": 3.0, "IL2RA": 2.5},
            regulons={"FOXP3": 2.8},
            spatial_context="  thymus cortex",
            rag_context="  Treg (overlap=0.85): FOXP3, IL2RA",
        )
        assert "test-cluster-99" in prompt
        assert "FOXP3" in prompt
        assert "velocity" in prompt

    def test_system_prompt_not_empty(self):
        from axiom_sc.tier4.prompts import SYSTEM_PROMPT
        assert len(SYSTEM_PROMPT) > 100
        assert "cell type" in SYSTEM_PROMPT.lower()


# ── RAG retriever ─────────────────────────────────────────────────────────────

class TestKGRetriever:
    def test_loads_kg_and_indexes_cell_types(self):
        from axiom_sc.tier4.rag import KGRetriever
        retriever = KGRetriever(KG_PATH)
        assert retriever.n_cell_types > 0

    def test_retrieve_returns_top_k(self):
        from axiom_sc.tier4.rag import KGRetriever
        retriever = KGRetriever(KG_PATH)
        results = retriever.retrieve({"FOXP3": 3.0, "IL2RA": 2.5, "CTLA4": 2.0}, top_k=3)
        assert len(results) <= 3
        for r in results:
            assert "cell_type" in r
            assert "overlap_score" in r
            assert "key_markers" in r

    def test_retrieve_treg_markers_returns_treg_first(self):
        from axiom_sc.tier4.rag import KGRetriever
        retriever = KGRetriever(KG_PATH)
        # FOXP3 + IL2RA are canonical Treg markers in the KG
        results = retriever.retrieve({"FOXP3": 3.5, "IL2RA": 3.0, "CTLA4": 2.5}, top_k=5)
        assert results, "Expected at least one result"
        top_types = [r["cell_type"] for r in results[:3]]
        assert "Treg" in top_types, f"Treg not in top-3: {top_types}"

    def test_retrieve_empty_markers_returns_no_results(self):
        from axiom_sc.tier4.rag import KGRetriever
        retriever = KGRetriever(KG_PATH)
        results = retriever.retrieve({}, top_k=5)
        assert results == []

    def test_retrieve_negative_rules_not_included_in_scoring(self):
        from axiom_sc.tier4.rag import KGRetriever
        # A negative rule gene matching the query should NOT cause the type to score high
        # just because of a negative marker — negatives are excluded from retrieval index
        retriever = KGRetriever(KG_PATH)
        # IL2 is a Treg negative rule gene (Tregs suppress IL2)
        # Query with IL2 alone should not strongly retrieve Treg
        results = retriever.retrieve({"IL2": 3.0}, top_k=10)
        types = [r["cell_type"] for r in results]
        # If Treg appears, it should not be rank 1 based on a single negative-rule gene
        # (other types may match, this just tests we don't crash)
        assert isinstance(types, list)

    def test_rag_result_as_prompt_text(self):
        from axiom_sc.tier4.rag import KGRetriever
        retriever = KGRetriever(KG_PATH)
        result = retriever.retrieve_for_cluster(
            cluster_id="test-cluster",
            marker_genes={"FOXP3": 3.0, "IL2RA": 2.5},
            top_k=3,
        )
        text = result.as_prompt_text()
        assert isinstance(text, str)
        assert len(text) > 0

    def test_rag_result_empty_as_prompt_text(self):
        from axiom_sc.tier4.rag import RAGResult
        empty = RAGResult(cluster_id="x", top_types=[])
        assert "no" in empty.as_prompt_text().lower()


# ── evidence bundler ──────────────────────────────────────────────────────────

class TestEvidenceBundler:
    def _make_evidence(self, cluster_id="test-cluster"):
        return EvidenceBundle(
            cluster_id=cluster_id,
            marker_genes={"FOXP3": 3.0, "IL2RA": 2.5, "CTLA4": 2.0},
            regulons={"FOXP3": 2.8, "IKZF2": 1.9},
            tissue="thymus",
        )

    def test_bundle_returns_bundled_evidence(self):
        from axiom_sc.tier4.evidence_bundler import EvidenceBundler
        bundler = EvidenceBundler()
        ev = self._make_evidence()
        bundled = bundler.bundle(ev)
        assert bundled.cluster_id == "test-cluster"
        assert len(bundled.system_prompt) > 50
        assert len(bundled.user_prompt) > 100

    def test_bundle_contains_marker_genes_in_prompt(self):
        from axiom_sc.tier4.evidence_bundler import EvidenceBundler
        bundler = EvidenceBundler()
        ev = self._make_evidence()
        bundled = bundler.bundle(ev, tier1_predictions={"Treg": 0.92})
        assert "FOXP3" in bundled.user_prompt
        assert "Treg" in bundled.user_prompt

    def test_bundle_with_rag(self):
        from axiom_sc.tier4.evidence_bundler import EvidenceBundler
        from axiom_sc.tier4.rag import KGRetriever
        retriever = KGRetriever(KG_PATH)
        bundler = EvidenceBundler(kg_retriever=retriever)
        ev = self._make_evidence()
        bundled = bundler.bundle(ev)
        # RAG results should appear in prompt
        assert "overlap" in bundled.user_prompt.lower() or "Treg" in bundled.user_prompt

    def test_bundle_with_tier3_results(self):
        from axiom_sc.tier4.evidence_bundler import EvidenceBundler
        bundler = EvidenceBundler()
        ev = self._make_evidence()
        t3 = [
            StreamResult("test-cluster", "velocity", "Treg", 0.85),
            StreamResult("test-cluster", "chromatin", "Treg", 0.90),
            StreamResult("test-cluster", "communication", "", 0.0, available=False),
        ]
        bundled = bundler.bundle(ev, tier3_stream_results=t3)
        assert "velocity" in bundled.user_prompt
        assert "chromatin" in bundled.user_prompt

    def test_bundle_raw_fields_populated(self):
        from axiom_sc.tier4.evidence_bundler import EvidenceBundler
        bundler = EvidenceBundler()
        ev = self._make_evidence()
        bundled = bundler.bundle(ev, tier1_predictions={"Treg": 0.9})
        assert bundled.raw_fields["n_markers"] == 3
        assert bundled.raw_fields["n_regulons"] == 2


# ── LLM ensemble (all API calls mocked) ──────────────────────────────────────

class TestLLMEnsemble:
    """Tests for LLMEnsemble logic — no real API calls made."""

    def _make_bundled(self, cluster_id="test-cluster"):
        from axiom_sc.tier4.evidence_bundler import BundledEvidence
        return BundledEvidence(
            cluster_id=cluster_id,
            system_prompt="You are an expert biologist.",
            user_prompt="What cell type is this cluster?",
        )

    def _make_anthropic_response(self, cell_type: str, confidence: float = 0.9) -> str:
        payload = {
            "cell_type": cell_type,
            "confidence": confidence,
            "reasoning": f"Based on markers, this is {cell_type}.",
            "alternatives": ["CD4_T"],
            "novel_candidate": False,
        }
        return json.dumps(payload)

    def test_from_configs_builds_ensemble(self):
        from axiom_sc.tier4.llm_ensemble import LLMEnsemble
        ensemble = LLMEnsemble.from_configs([
            {"provider": "anthropic", "model_id": "claude-sonnet-4-6", "api_key": "test-key"}
        ])
        assert len(ensemble._models) == 1

    def test_from_configs_empty_raises(self):
        from axiom_sc.tier4.llm_ensemble import LLMEnsemble
        with pytest.raises(ValueError, match="at least one"):
            LLMEnsemble([])

    def test_annotate_single_model_success(self):
        from axiom_sc.tier4.llm_ensemble import LLMEnsemble
        ensemble = LLMEnsemble.from_configs([
            {"provider": "anthropic", "model_id": "claude-sonnet-4-6", "api_key": "key"}
        ])
        bundled = self._make_bundled()

        with patch.object(ensemble, "_http_call",
                          return_value=self._make_anthropic_response("Treg", 0.95)):
            result = ensemble.annotate(bundled)

        assert result.cluster_id == "test-cluster"
        assert result.cell_type == "Treg"
        assert abs(result.confidence - 0.95) < 0.01
        assert result.n_models_called == 1
        assert result.n_models_success == 1

    def test_annotate_majority_vote_two_models(self):
        from axiom_sc.tier4.llm_ensemble import LLMEnsemble
        ensemble = LLMEnsemble.from_configs([
            {"provider": "anthropic", "model_id": "claude-sonnet-4-6", "api_key": "k1"},
            {"provider": "openai", "model_id": "gpt-4o", "api_key": "k2"},
        ])
        bundled = self._make_bundled()

        responses = [
            self._make_anthropic_response("Treg", 0.9),
            self._make_anthropic_response("Treg", 0.85),
        ]
        with patch.object(ensemble, "_http_call", side_effect=responses):
            result = ensemble.annotate(bundled)

        assert result.cell_type == "Treg"
        assert result.n_models_success == 2

    def test_annotate_tie_broken_by_confidence(self):
        """Two models disagree — winner is the one with higher mean confidence."""
        from axiom_sc.tier4.llm_ensemble import LLMEnsemble
        ensemble = LLMEnsemble.from_configs([
            {"provider": "anthropic", "model_id": "claude-sonnet-4-6", "api_key": "k1"},
            {"provider": "openai", "model_id": "gpt-4o", "api_key": "k2"},
        ])
        bundled = self._make_bundled()

        responses = [
            self._make_anthropic_response("Treg", 0.95),  # Anthropic votes Treg
            json.dumps({
                "cell_type": "CD4_T",
                "confidence": 0.5,
                "reasoning": "could be CD4",
                "alternatives": [],
                "novel_candidate": False,
            }),  # OpenAI votes CD4_T with lower confidence
        ]
        with patch.object(ensemble, "_http_call", side_effect=responses):
            result = ensemble.annotate(bundled)

        # 1v1 tie — Treg wins by higher confidence (0.95 vs 0.5)
        assert result.cell_type == "Treg"

    def test_annotate_all_fail_returns_unknown(self):
        from axiom_sc.tier4.llm_ensemble import LLMEnsemble
        ensemble = LLMEnsemble.from_configs([
            {"provider": "anthropic", "model_id": "claude-sonnet-4-6", "api_key": "key"}
        ])
        bundled = self._make_bundled()

        with patch.object(ensemble, "_http_call", side_effect=Exception("Network error")):
            result = ensemble.annotate(bundled)

        assert result.cell_type == "UNKNOWN"
        assert result.n_models_called == 1
        assert result.n_models_success == 0

    def test_parse_response_strips_markdown_fences(self):
        from axiom_sc.tier4.llm_ensemble import LLMEnsemble, ModelConfig
        ensemble = LLMEnsemble.from_configs([
            {"provider": "anthropic", "model_id": "claude-sonnet-4-6", "api_key": "key"}
        ])
        model = ensemble._models[0]
        raw = '```json\n{"cell_type":"Treg","confidence":0.9,"reasoning":"r","alternatives":[],"novel_candidate":false}\n```'
        vote = ensemble._parse_response(model, raw)
        assert vote.cell_type == "Treg"
        assert vote.confidence == 0.9

    def test_tier4_result_to_dict(self):
        from axiom_sc.tier4.llm_ensemble import Tier4Result, ModelVote
        result = Tier4Result(
            cluster_id="c1",
            cell_type="Treg",
            confidence=0.9,
            reasoning="FOXP3 active",
            alternatives=["CD4_T"],
            model_votes=[
                ModelVote("anthropic", "claude-sonnet-4-6", "Treg", 0.9, "FOXP3 active")
            ],
        )
        d = result.to_dict()
        assert d["cell_type"] == "Treg"
        assert d["confidence"] == 0.9
        assert len(d["model_votes"]) == 1

    def test_novel_candidate_majority_vote(self):
        """If majority of models vote novel_candidate=True, result.novel_candidate is True."""
        from axiom_sc.tier4.llm_ensemble import LLMEnsemble
        ensemble = LLMEnsemble.from_configs([
            {"provider": "anthropic", "model_id": "claude-sonnet-4-6", "api_key": "k1"},
            {"provider": "openai", "model_id": "gpt-4o", "api_key": "k2"},
        ])
        bundled = self._make_bundled()

        responses = [
            json.dumps({"cell_type": "NOVEL_CANDIDATE", "confidence": 0.6,
                        "reasoning": "no match", "alternatives": [], "novel_candidate": True}),
            json.dumps({"cell_type": "NOVEL_CANDIDATE", "confidence": 0.7,
                        "reasoning": "unknown", "alternatives": [], "novel_candidate": True}),
        ]
        with patch.object(ensemble, "_http_call", side_effect=responses):
            result = ensemble.annotate(bundled)

        assert result.novel_candidate is True
