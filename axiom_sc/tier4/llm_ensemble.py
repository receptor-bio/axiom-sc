"""
Tier 4 LLM elite ensemble for cell type annotation.

Model selection based on DeepCellSeek benchmark (2025):
  Elite ensemble: Claude (Anthropic), GPT (OpenAI)
  Key finding: ensemble outperforms any single model on subtype annotation.

CRITICAL DIFFERENCE FROM CASSIA:
  CASSIA sends marker genes to LLM (single evidence stream).
  AXIOM-SC Tier 4 sends the FULL evidence bundle (Tier 1-3 results +
  velocity direction + chromatin state + spatial context + AXIOM rules fired).
  LLM reasons over all evidence simultaneously, not just marker genes.

API keys are user-provided and never logged or transmitted beyond the configured endpoint.
No single-cell expression data is sent to any LLM API — only the evidence summary text.

License: Apache 2.0
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from typing import Optional

from axiom_sc.tier4.evidence_bundler import BundledEvidence
from axiom_sc.tier4.prompts import PROMPT_VERSION


@dataclass
class ModelConfig:
    """Configuration for one LLM in the ensemble."""

    provider: str          # "anthropic" | "openai"
    model_id: str          # e.g. "claude-sonnet-4-6", "gpt-4o"
    api_key: str
    max_tokens: int = 512
    temperature: float = 0.0   # deterministic by default


@dataclass
class ModelVote:
    """Raw annotation result from one LLM model."""

    provider: str
    model_id: str
    cell_type: str
    confidence: float
    reasoning: str
    alternatives: list[str] = field(default_factory=list)
    novel_candidate: bool = False
    parse_error: str = ""


@dataclass
class Tier4Result:
    """
    Consensus annotation result from the Tier 4 LLM ensemble.

    cluster_id    : str     — which cluster was annotated
    cell_type     : str     — consensus cell type label (or 'NOVEL_CANDIDATE')
    confidence    : float   — mean confidence across voting models
    reasoning     : str     — reasoning from highest-confidence model
    alternatives  : list    — runner-up cell types
    novel_candidate: bool   — True if majority voted novel
    model_votes   : list    — per-model raw votes
    prompt_version: str     — PROMPT_VERSION used to build the prompt
    n_models_called: int    — how many models were queried
    n_models_success: int   — how many returned parseable responses
    """

    cluster_id: str
    cell_type: str
    confidence: float
    reasoning: str
    alternatives: list[str] = field(default_factory=list)
    novel_candidate: bool = False
    model_votes: list[ModelVote] = field(default_factory=list)
    prompt_version: str = PROMPT_VERSION
    n_models_called: int = 0
    n_models_success: int = 0

    def to_dict(self) -> dict:
        return {
            "cluster_id": self.cluster_id,
            "cell_type": self.cell_type,
            "confidence": self.confidence,
            "reasoning": self.reasoning,
            "alternatives": self.alternatives,
            "novel_candidate": self.novel_candidate,
            "prompt_version": self.prompt_version,
            "n_models_called": self.n_models_called,
            "n_models_success": self.n_models_success,
            "model_votes": [
                {
                    "provider": v.provider,
                    "model_id": v.model_id,
                    "cell_type": v.cell_type,
                    "confidence": v.confidence,
                }
                for v in self.model_votes
            ],
        }


class LLMEnsemble:
    """
    Multi-model LLM ensemble for hard-case cell type annotation.

    Calls each configured model, collects JSON responses, resolves consensus
    by majority vote on cell_type with mean confidence.

    Usage
    -----
    ensemble = LLMEnsemble.from_env()   # reads ANTHROPIC_API_KEY, OPENAI_API_KEY
    result = ensemble.annotate(bundled_evidence)
    """

    def __init__(self, models: list[ModelConfig]):
        if not models:
            raise ValueError("LLMEnsemble requires at least one ModelConfig")
        self._models = models

    # ── factory ───────────────────────────────────────────────────────────────

    @classmethod
    def from_env(
        cls,
        include_anthropic: bool = True,
        include_openai: bool = True,
    ) -> "LLMEnsemble":
        """Build ensemble from environment variables."""
        models: list[ModelConfig] = []
        if include_anthropic:
            key = os.environ.get("ANTHROPIC_API_KEY", "")
            if key:
                models.append(ModelConfig(
                    provider="anthropic",
                    model_id="claude-sonnet-4-6",
                    api_key=key,
                ))
        if include_openai:
            key = os.environ.get("OPENAI_API_KEY", "")
            if key:
                models.append(ModelConfig(
                    provider="openai",
                    model_id="gpt-4o",
                    api_key=key,
                ))
        if not models:
            raise RuntimeError(
                "No API keys found. Set ANTHROPIC_API_KEY and/or OPENAI_API_KEY."
            )
        return cls(models)

    @classmethod
    def from_configs(cls, configs: list[dict]) -> "LLMEnsemble":
        """Build from a list of config dicts (for programmatic construction)."""
        return cls([ModelConfig(**c) for c in configs])

    # ── public API ────────────────────────────────────────────────────────────

    def annotate(self, bundled: BundledEvidence) -> Tier4Result:
        """
        Call all models and return consensus annotation.

        All models receive the same prompt. Consensus:
          - cell_type: majority vote (ties broken by highest mean confidence)
          - confidence: mean of winning votes
          - reasoning: from the winning model with highest individual confidence
        """
        votes: list[ModelVote] = []
        for model in self._models:
            vote = self._call_model(model, bundled)
            votes.append(vote)

        successful = [v for v in votes if not v.parse_error]
        result = self._consensus(bundled.cluster_id, votes, successful)
        result.n_models_called = len(votes)
        result.n_models_success = len(successful)
        return result

    # ── private ───────────────────────────────────────────────────────────────

    def _call_model(self, model: ModelConfig, bundled: BundledEvidence) -> ModelVote:
        """Call one model and return its vote. Never raises — errors go into parse_error."""
        try:
            raw = self._http_call(model, bundled.system_prompt, bundled.user_prompt)
            return self._parse_response(model, raw)
        except Exception as exc:
            return ModelVote(
                provider=model.provider,
                model_id=model.model_id,
                cell_type="UNKNOWN",
                confidence=0.0,
                reasoning="",
                parse_error=str(exc),
            )

    def _http_call(self, model: ModelConfig, system: str, user: str) -> str:
        """Make the HTTP request and return the raw response text."""
        try:
            import httpx
        except ImportError as e:
            raise ImportError(
                "httpx is required for Tier 4. Install with: pip install httpx"
            ) from e

        if model.provider == "anthropic":
            return self._call_anthropic(model, system, user)
        if model.provider == "openai":
            return self._call_openai(model, system, user)
        raise ValueError(f"Unknown LLM provider: {model.provider!r}")

    def _call_anthropic(self, model: ModelConfig, system: str, user: str) -> str:
        import httpx

        payload = {
            "model": model.model_id,
            "max_tokens": model.max_tokens,
            "temperature": model.temperature,
            "system": system,
            "messages": [{"role": "user", "content": user}],
        }
        with httpx.Client(timeout=60.0) as client:
            resp = client.post(
                "https://api.anthropic.com/v1/messages",
                headers={
                    "x-api-key": model.api_key,
                    "anthropic-version": "2023-06-01",
                    "content-type": "application/json",
                },
                json=payload,
            )
        resp.raise_for_status()
        data = resp.json()
        return data["content"][0]["text"]

    def _call_openai(self, model: ModelConfig, system: str, user: str) -> str:
        import httpx

        payload = {
            "model": model.model_id,
            "temperature": model.temperature,
            "max_tokens": model.max_tokens,
            "messages": [
                {"role": "system", "content": system},
                {"role": "user", "content": user},
            ],
        }
        with httpx.Client(timeout=60.0) as client:
            resp = client.post(
                "https://api.openai.com/v1/chat/completions",
                headers={
                    "Authorization": f"Bearer {model.api_key}",
                    "content-type": "application/json",
                },
                json=payload,
            )
        resp.raise_for_status()
        data = resp.json()
        return data["choices"][0]["message"]["content"]

    def _parse_response(self, model: ModelConfig, raw: str) -> ModelVote:
        """Parse JSON from model response into a ModelVote."""
        text = raw.strip()
        # Strip markdown code fences if present
        if text.startswith("```"):
            lines = text.splitlines()
            text = "\n".join(
                l for l in lines if not l.startswith("```")
            ).strip()
        parsed = json.loads(text)
        return ModelVote(
            provider=model.provider,
            model_id=model.model_id,
            cell_type=str(parsed.get("cell_type", "UNKNOWN")),
            confidence=float(parsed.get("confidence", 0.0)),
            reasoning=str(parsed.get("reasoning", "")),
            alternatives=list(parsed.get("alternatives", [])),
            novel_candidate=bool(parsed.get("novel_candidate", False)),
        )

    def _consensus(
        self,
        cluster_id: str,
        all_votes: list[ModelVote],
        successful: list[ModelVote],
    ) -> Tier4Result:
        """Resolve majority vote across successful model responses."""
        if not successful:
            return Tier4Result(
                cluster_id=cluster_id,
                cell_type="UNKNOWN",
                confidence=0.0,
                reasoning="All LLM calls failed — check API keys and connectivity.",
                model_votes=all_votes,
            )

        # Majority vote on cell_type
        type_votes: dict[str, list[float]] = {}
        for v in successful:
            type_votes.setdefault(v.cell_type, []).append(v.confidence)

        winning_type = max(type_votes, key=lambda t: (len(type_votes[t]), sum(type_votes[t])))
        winning_confs = type_votes[winning_type]
        mean_conf = sum(winning_confs) / len(winning_confs)

        # Best reasoning from highest-confidence winner
        winners = [v for v in successful if v.cell_type == winning_type]
        best = max(winners, key=lambda v: v.confidence)

        # Alternatives: other cell types not chosen as winner
        alternatives = sorted(
            [t for t in type_votes if t != winning_type],
            key=lambda t: -sum(type_votes[t]) / len(type_votes[t]),
        )

        novel = sum(1 for v in successful if v.novel_candidate) > len(successful) / 2

        return Tier4Result(
            cluster_id=cluster_id,
            cell_type=winning_type,
            confidence=round(mean_conf, 4),
            reasoning=best.reasoning,
            alternatives=alternatives[:3],
            novel_candidate=novel,
            model_votes=all_votes,
        )
