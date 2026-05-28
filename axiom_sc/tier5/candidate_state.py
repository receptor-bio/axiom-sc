"""
CandidateState: dataclass representing a putative novel cell state from Tier 5.

Cells that fail all four annotation tiers are not discarded — they are characterized
as candidate novel states. This is the key Tier 5 output: a structured description
of an unknown attractor that can be validated experimentally.

License: Apache 2.0
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Literal

# Valid status values matching the KG status vocabulary
CandidateStatus = Literal["CANDIDATE", "VALIDATED", "PUBLISHED"]


@dataclass
class PerturbationPrediction:
    """One falsifiable experimental prediction for a candidate state."""

    gene: str
    perturbation: Literal["KO", "OE", "KD"]   # knockout / overexpression / knockdown
    predicted_effect: str                        # qualitative description
    confidence: float = 0.5                      # 0-1


@dataclass
class CandidateState:
    """
    Structured description of a putative novel cell attractor.

    Fields
    ------
    candidate_id       : unique hash of dataset_id + cluster_id
    cluster_id         : source cluster identifier
    dataset_id         : source dataset identifier
    active_circuits    : TF regulons with z-score ≥ 1.5 from SCENIC
    velocity_sink      : nearest known attractor in velocity space (None if unavailable)
    spatial_context    : tissue niche description from CELLama or manual annotation
    cross_species      : OrthoFinder conservation scores {species: score}
    nearest_ontology_term: closest Cell Ontology term by marker distance
    perturbation_predictions: falsifiable experimental predictions
    status             : "CANDIDATE" | "VALIDATED" | "PUBLISHED"
    evidence_summary   : raw evidence fields used to produce this characterization
    created_at         : ISO 8601 UTC timestamp of discovery
    proposed_name      : optional working name assigned by LLM or researcher
    """

    cluster_id: str
    dataset_id: str
    active_circuits: list[str] = field(default_factory=list)
    velocity_sink: str | None = None
    spatial_context: str = ""
    cross_species: dict[str, float] = field(default_factory=dict)
    nearest_ontology_term: str = ""
    perturbation_predictions: list[PerturbationPrediction] = field(default_factory=list)
    status: CandidateStatus = "CANDIDATE"
    evidence_summary: dict = field(default_factory=dict)
    created_at: str = field(
        default_factory=lambda: datetime.now(tz=timezone.utc).isoformat()
    )
    proposed_name: str = ""

    def __post_init__(self) -> None:
        # Validate status
        if self.status not in ("CANDIDATE", "VALIDATED", "PUBLISHED"):
            raise ValueError(
                f"Invalid status {self.status!r}. "
                "Must be 'CANDIDATE', 'VALIDATED', or 'PUBLISHED'."
            )

    @property
    def candidate_id(self) -> str:
        """Deterministic hash of dataset_id + cluster_id."""
        raw = f"{self.dataset_id}::{self.cluster_id}"
        return hashlib.sha256(raw.encode()).hexdigest()[:16]

    def to_dict(self) -> dict:
        """Serialise to a JSON-compatible dict."""
        return {
            "candidate_id": self.candidate_id,
            "cluster_id": self.cluster_id,
            "dataset_id": self.dataset_id,
            "active_circuits": self.active_circuits,
            "velocity_sink": self.velocity_sink,
            "spatial_context": self.spatial_context,
            "cross_species": self.cross_species,
            "nearest_ontology_term": self.nearest_ontology_term,
            "perturbation_predictions": [
                {
                    "gene": p.gene,
                    "perturbation": p.perturbation,
                    "predicted_effect": p.predicted_effect,
                    "confidence": p.confidence,
                }
                for p in self.perturbation_predictions
            ],
            "status": self.status,
            "evidence_summary": self.evidence_summary,
            "created_at": self.created_at,
            "proposed_name": self.proposed_name,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "CandidateState":
        """Deserialise from a dict (e.g. loaded from JSON)."""
        preds = [
            PerturbationPrediction(
                gene=p["gene"],
                perturbation=p["perturbation"],
                predicted_effect=p["predicted_effect"],
                confidence=p.get("confidence", 0.5),
            )
            for p in d.get("perturbation_predictions", [])
        ]
        return cls(
            cluster_id=d["cluster_id"],
            dataset_id=d["dataset_id"],
            active_circuits=d.get("active_circuits", []),
            velocity_sink=d.get("velocity_sink"),
            spatial_context=d.get("spatial_context", ""),
            cross_species=d.get("cross_species", {}),
            nearest_ontology_term=d.get("nearest_ontology_term", ""),
            perturbation_predictions=preds,
            status=d.get("status", "CANDIDATE"),
            evidence_summary=d.get("evidence_summary", {}),
            created_at=d.get("created_at", ""),
            proposed_name=d.get("proposed_name", ""),
        )

    def summary_line(self) -> str:
        circuits = ", ".join(self.active_circuits[:3]) or "none"
        return (
            f"[{self.status}] {self.candidate_id} | cluster={self.cluster_id} | "
            f"circuits={circuits} | sink={self.velocity_sink or 'n/a'} | "
            f"ontology={self.nearest_ontology_term or 'unknown'}"
        )
