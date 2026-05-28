"""
StreamResult: per-cluster result from a single Tier 3 evidence stream.

Shared type imported by velocity.py, chromatin.py, communication.py,
spatial_niche.py, cross_species.py, and convergence.py.

License: Apache-2.0
"""
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class StreamResult:
    """
    Annotation result produced by one Tier 3 stream for one cluster.

    Parameters
    ----------
    cluster_id : str
        Cluster identifier.
    stream : str
        Stream name: "velocity" | "chromatin" | "communication" |
        "spatial_niche" | "cross_species" | "spatial_type"
    label : str
        Best cell type label predicted by this stream.
        Empty string if stream could not assign a label.
    confidence : float
        Confidence in [0, 1].
    available : bool
        False if the stream could not run (missing data / missing deps).
    evidence : dict
        Stream-specific evidence details (gene values, scores, etc.).
    """

    cluster_id: str
    stream: str
    label: str
    confidence: float
    available: bool = True
    evidence: dict = field(default_factory=dict)

    def agrees_with(self, other: "StreamResult", tolerance: float = 0.0) -> bool:
        """True if both streams are available and returned the same label."""
        if not self.available or not other.available:
            return False
        return self.label == other.label

    def __repr__(self) -> str:
        status = "ok" if self.available else "n/a"
        return (
            f"StreamResult({self.stream}/{self.cluster_id}, "
            f"label={self.label!r}, conf={self.confidence:.2f}, {status})"
        )
