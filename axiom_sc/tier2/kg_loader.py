"""
KGLoader: loads and validates the AXIOM-SC knowledge graph from JSON.

Phase 1 corrections carried forward:
- Status filtering: only ACTIVE rules are returned for annotation.
- PMID enforcement: all rules must have non-empty pmid (validated at load).
- Tissue-context filtering: optional filter to restrict rules to a tissue.

License: Apache 2.0
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator

import jsonschema


_SCHEMA_PATH = Path(__file__).parent.parent / "kg" / "schema.json"


@dataclass
class KGRule:
    """
    Single mechanistic rule in the AXIOM-SC knowledge graph.
    Fields mirror kg/schema.json exactly.
    """

    rule_id: str
    cell_type: str
    rule_type: str           # positive | negative | circuit | spatial
    evidence_source: str     # marker_genes | regulon | trajectory | ...
    gene_or_regulon: list[str]
    direction: str           # high | low | absent | present | active | inactive
    mechanistic_basis: str
    pmid: str
    confidence: str          # high | medium | low
    tissue_context: list[str]
    source_db: str
    status: str              # ACTIVE | PENDING_REVIEW | DEPRECATED
    paired_with: list[str] = field(default_factory=list)
    incompatible_with: list[str] = field(default_factory=list)
    added_in_version: str = ""

    @classmethod
    def from_dict(cls, d: dict) -> "KGRule":
        return cls(
            rule_id=d["rule_id"],
            cell_type=d["cell_type"],
            rule_type=d["rule_type"],
            evidence_source=d["evidence_source"],
            gene_or_regulon=d["gene_or_regulon"],
            direction=d["direction"],
            mechanistic_basis=d["mechanistic_basis"],
            pmid=d["pmid"],
            confidence=d["confidence"],
            tissue_context=d.get("tissue_context", []),
            source_db=d["source_db"],
            status=d["status"],
            paired_with=d.get("paired_with", []),
            incompatible_with=d.get("incompatible_with", []),
            added_in_version=d.get("added_in_version", ""),
        )

    def to_dict(self) -> dict:
        return {
            "rule_id": self.rule_id,
            "cell_type": self.cell_type,
            "rule_type": self.rule_type,
            "evidence_source": self.evidence_source,
            "gene_or_regulon": self.gene_or_regulon,
            "direction": self.direction,
            "mechanistic_basis": self.mechanistic_basis,
            "pmid": self.pmid,
            "confidence": self.confidence,
            "tissue_context": self.tissue_context,
            "source_db": self.source_db,
            "status": self.status,
            "paired_with": self.paired_with,
            "incompatible_with": self.incompatible_with,
            "added_in_version": self.added_in_version,
        }


class KGLoader:
    """
    Loads, validates, and indexes an AXIOM-SC knowledge graph JSON file.

    Usage::

        kg = KGLoader.from_file("kg_data/oracle_kg_v0.2.0.json")
        rules = kg.active_rules(cell_type="Treg")
        cell_types = kg.cell_types()
    """

    def __init__(self, rules: list[KGRule]):
        self._rules = rules
        # index by cell_type for fast lookup
        self._by_cell_type: dict[str, list[KGRule]] = {}
        for rule in rules:
            self._by_cell_type.setdefault(rule.cell_type, []).append(rule)

    @classmethod
    def from_file(cls, path: str | Path) -> "KGLoader":
        """
        Load KG from a JSON file (list of rule dicts).
        Validates every rule against kg/schema.json.
        Raises jsonschema.ValidationError on first schema violation.
        """
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"KG file not found: {path}")

        with open(path) as f:
            raw: list[dict] = json.load(f)

        if not isinstance(raw, list):
            raise ValueError("KG JSON must be a top-level list of rule objects.")

        schema = cls._load_schema()
        rules: list[KGRule] = []
        for i, d in enumerate(raw):
            try:
                jsonschema.validate(d, schema)
            except jsonschema.ValidationError as e:
                raise jsonschema.ValidationError(
                    f"Rule #{i} (id={d.get('rule_id', '?')}) failed schema validation: {e.message}"
                ) from e
            rules.append(KGRule.from_dict(d))

        return cls(rules)

    @classmethod
    def from_list(cls, rules: list[dict]) -> "KGLoader":
        """Build a KGLoader from a list of rule dicts (useful in tests)."""
        schema = cls._load_schema()
        validated: list[KGRule] = []
        for d in rules:
            jsonschema.validate(d, schema)
            validated.append(KGRule.from_dict(d))
        return cls(validated)

    @staticmethod
    def _load_schema() -> dict:
        with open(_SCHEMA_PATH) as f:
            return json.load(f)

    # ── query interface ───────────────────────────────────────────────────────

    def active_rules(
        self,
        cell_type: str | None = None,
        tissue: str | None = None,
    ) -> list[KGRule]:
        """
        Returns ACTIVE rules, optionally filtered by cell_type and tissue.
        Tissue filter: returns rules whose tissue_context is empty (universal)
        OR includes the given tissue string.
        """
        rules = self._rules if cell_type is None else self._by_cell_type.get(cell_type, [])
        result = []
        for r in rules:
            if r.status != "ACTIVE":
                continue
            if tissue and r.tissue_context and tissue not in r.tissue_context:
                continue
            result.append(r)
        return result

    def cell_types(self) -> list[str]:
        """Returns the sorted list of unique cell types in the KG (all statuses)."""
        return sorted(self._by_cell_type.keys())

    def active_cell_types(self) -> list[str]:
        """Returns cell types that have at least one ACTIVE rule."""
        return sorted(
            ct for ct in self._by_cell_type
            if any(r.status == "ACTIVE" for r in self._by_cell_type[ct])
        )

    def rule_count(self, status: str = "ACTIVE") -> int:
        return sum(1 for r in self._rules if r.status == status)

    def __len__(self) -> int:
        return len(self._rules)

    def __iter__(self) -> Iterator[KGRule]:
        return iter(self._rules)

    def __repr__(self) -> str:
        active = self.rule_count("ACTIVE")
        return (
            f"KGLoader(total={len(self._rules)}, active={active}, "
            f"cell_types={len(self._by_cell_type)})"
        )
