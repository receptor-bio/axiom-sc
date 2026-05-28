"""
KG seeder from CellMarker 2.0.

Downloads the CellMarker 2.0 human cell marker Excel file and converts
marker entries into candidate KG rules with status="PENDING_REVIEW".
Candidates must be reviewed by an expert before promotion to ACTIVE.

License note:
  CellMarker 2.0 paper: Hu et al. (2023) Nucleic Acids Research 51:D870
  doi: 10.1093/nar/gkac947
  Data license: CC BY 4.0 (unrestricted including commercial).
  WARNING: CellMarker 1.0 was CC BY-NC 4.0 — this seeder uses v2.0 specifically.

DO NOT USE for seeding:
  - PanglaoDB: CC BY-NC 4.0 (non-commercial only)
  - scTypeDB: in GPL v3 repository

License: Apache 2.0
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import time
from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path
from typing import Optional

CELLMARKER2_URL = (
    "http://bio-bigdata.hrbmu.edu.cn/CellMarker/"
    "CellMarker_download_files/file/Cell_marker_Human.xlsx"
)

# Column names in the CellMarker 2.0 Excel file
_COL_SPECIES = "speciesType"
_COL_TISSUE_CLASS = "tissueType"
_COL_CELL_NAME = "cellName"
_COL_CELL_TYPE = "cellType"          # Normal / Cancer
_COL_GENE_SYMBOL = "geneSymbol"
_COL_PMID = "PMID"
_COL_MARKER_RESOURCE = "markerResource"

# Known transcription factors that should be suggested as regulon evidence
# (not exhaustive — just a heuristic for candidate generation)
_KNOWN_TFS = frozenset([
    "FOXP3", "TBX21", "GATA3", "RORC", "IRF7", "IRF8", "IRF4",
    "PAX5", "EBF1", "AIRE", "PRDM1", "XBP1", "BCL6", "TOX",
    "TCF7", "SPI1", "CEBPA", "IKZF2", "IKZF3", "BATF3", "ID2",
    "EOMES", "RUNX3", "NR4A1", "STAT1", "STAT3", "STAT4", "STAT6",
    "MYB", "E2F1", "NFKB1", "NFKB2", "RELA", "TP53", "MYC",
])

# Tissues supported in the default seeder run
DEFAULT_TISSUES = ["thymus", "blood", "lung", "liver", "gut"]


@dataclass
class CandidateRule:
    """A seeded candidate rule awaiting expert review."""
    cell_type: str
    rule_id: str
    rule_type: str
    evidence_source: str
    gene_or_regulon: list
    direction: str
    paired_with: list
    incompatible_with: list
    mechanistic_basis: str
    pmid: str
    confidence: str
    tissue_context: list
    source_db: str
    status: str
    added_in_version: str


def _normalize_cell_name(name: str) -> str:
    """Convert a CellMarker2 cell name to a KG-style identifier."""
    s = name.strip()
    s = re.sub(r"[^A-Za-z0-9 _]", "", s)
    s = re.sub(r"\s+", "_", s)
    return s


def _make_rule_id(cell_type: str, gene: str, rule_type: str, idx: int) -> str:
    """Generate a deterministic rule ID: CELLTYPE_RULETYPE_NNN."""
    prefix = re.sub(r"[^A-Z0-9]", "", cell_type.upper())[:8]
    rtype = {"positive": "POS", "negative": "NEG", "circuit": "CIR"}.get(rule_type, "POS")
    return f"{prefix}_{rtype}_{idx:03d}"


def _infer_evidence_source(gene: str) -> str:
    """Suggest regulon if gene is a known TF, else marker_genes."""
    return "regulon" if gene.upper() in _KNOWN_TFS else "marker_genes"


def _infer_direction(evidence_source: str) -> str:
    return "active" if evidence_source == "regulon" else "high"


def _pmid_clean(raw: str) -> str:
    """Extract first numeric PMID from a possibly multi-valued string."""
    if not raw or str(raw).strip().lower() in ("nan", "none", ""):
        return "NEEDS_REVIEW"
    nums = re.findall(r"\d{5,9}", str(raw))
    return nums[0] if nums else "NEEDS_REVIEW"


class CellMarker2Seeder:
    """
    Downloads CellMarker 2.0 and generates PENDING_REVIEW candidate KG rules.

    Usage::

        seeder = CellMarker2Seeder(cache_dir="kg_data/cellmarker2_cache")
        seeder.run(
            tissues=["thymus", "blood", "lung"],
            output_path="kg_data/candidates/cellmarker2_candidates.json",
        )
    """

    def __init__(
        self,
        cache_dir: str = "kg_data/cellmarker2_cache",
        url: str = CELLMARKER2_URL,
        kg_version: str = "0.2.0",
        timeout: int = 30,
    ):
        self.cache_dir = Path(cache_dir)
        self.url = url
        self.kg_version = kg_version
        self.timeout = timeout

    # ── public interface ──────────────────────────────────────────────────────

    def run(
        self,
        tissues: list[str] | None = None,
        output_path: str = "kg_data/candidates/cellmarker2_candidates.json",
        max_candidates: int | None = None,
    ) -> Path:
        """
        Full seeder pipeline: download → parse → generate candidates → save.

        Returns the output path.
        Raises RuntimeError if download fails and no cached file exists.
        """
        tissues = tissues or DEFAULT_TISSUES
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        xlsx_path = self._get_xlsx(download=True)
        df = self._parse_xlsx(xlsx_path, tissues=tissues)
        candidates = self._generate_candidates(df)

        if max_candidates:
            candidates = candidates[:max_candidates]

        output = Path(output_path)
        output.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "metadata": {
                "source": "CellMarker 2.0",
                "doi": "10.1093/nar/gkac947",
                "url": self.url,
                "download_date": datetime.utcnow().isoformat() + "Z",
                "tissues": tissues,
                "n_candidates": len(candidates),
                "axiom_sc_version": self.kg_version,
                "license": "CC BY 4.0",
            },
            "candidates": [asdict(c) for c in candidates],
        }
        with open(output, "w") as f:
            json.dump(payload, f, indent=2)

        print(f"[seeder] Generated {len(candidates)} candidates → {output}")
        return output

    def generate_from_dataframe(self, df: "pd.DataFrame") -> list[CandidateRule]:
        """Public method for testing — generates candidates from a pre-loaded DataFrame."""
        return self._generate_candidates(df)

    # ── internals ─────────────────────────────────────────────────────────────

    def _get_xlsx(self, download: bool = True) -> Path:
        """Return path to cached Excel file, downloading if necessary."""
        cached = self.cache_dir / "Cell_marker_Human.xlsx"
        if cached.exists():
            return cached
        if not download:
            raise FileNotFoundError(
                f"CellMarker 2.0 file not in cache ({cached}). "
                "Run seeder with download=True, or place the file manually."
            )
        return self._download_xlsx(cached)

    def _download_xlsx(self, dest: Path) -> Path:
        """Download the CellMarker 2.0 Excel file."""
        try:
            import requests
        except ImportError:
            raise ImportError("requests package required: pip install requests")

        print(f"[seeder] Downloading CellMarker 2.0 from {self.url} …")
        try:
            resp = requests.get(self.url, timeout=self.timeout, stream=True)
            resp.raise_for_status()
            with open(dest, "wb") as f:
                for chunk in resp.iter_content(chunk_size=65536):
                    f.write(chunk)
            print(f"[seeder] Downloaded → {dest} ({dest.stat().st_size / 1e6:.1f} MB)")
            return dest
        except Exception as e:
            raise RuntimeError(
                f"CellMarker 2.0 download failed: {e}\n"
                f"Manual download: {self.url}\n"
                f"Place at: {dest}"
            ) from e

    def _parse_xlsx(self, xlsx_path: Path, tissues: list[str]) -> "pd.DataFrame":
        """Parse the Excel file into a filtered DataFrame."""
        try:
            import pandas as pd
        except ImportError:
            raise ImportError("pandas required: pip install pandas openpyxl")

        df = pd.read_excel(xlsx_path, engine="openpyxl")

        # Normalise column names (CellMarker2 uses camelCase)
        df.columns = [c.strip() for c in df.columns]

        # Filter: human normal cells only
        if "speciesType" in df.columns:
            df = df[df["speciesType"].str.strip() == "Human"]
        if "cellType" in df.columns:
            df = df[df["cellType"].str.strip() == "Normal"]

        # Filter to requested tissues (case-insensitive substring match)
        if tissues and "tissueType" in df.columns:
            tissue_pattern = "|".join(re.escape(t) for t in tissues)
            mask = df["tissueType"].str.lower().str.contains(
                tissue_pattern.lower(), na=False
            )
            df = df[mask]

        # Drop rows with missing essential fields
        essential = [c for c in ["cellName", "geneSymbol"] if c in df.columns]
        df = df.dropna(subset=essential)
        df = df[df["geneSymbol"].str.strip() != ""]

        return df.reset_index(drop=True)

    def _generate_candidates(self, df: "pd.DataFrame") -> list[CandidateRule]:
        """Convert a filtered CellMarker2 DataFrame into candidate KG rules."""
        candidates: list[CandidateRule] = []
        seen: set[str] = set()   # deduplication key = (cell_type, gene)

        # Count rules per cell type for numbering
        type_counter: dict[str, int] = {}

        for _, row in df.iterrows():
            raw_cell = str(row.get("cellName", "")).strip()
            raw_gene = str(row.get("geneSymbol", "")).strip()

            if not raw_cell or not raw_gene or raw_gene.lower() == "nan":
                continue

            cell_type = _normalize_cell_name(raw_cell)
            gene = raw_gene.upper().strip()

            dedup_key = f"{cell_type}::{gene}"
            if dedup_key in seen:
                continue
            seen.add(dedup_key)

            pmid = _pmid_clean(str(row.get("PMID", "")))
            tissue_raw = str(row.get("tissueType", "")).strip().lower()

            evidence_source = _infer_evidence_source(gene)
            direction = _infer_direction(evidence_source)

            type_counter[cell_type] = type_counter.get(cell_type, 0) + 1
            rule_id = _make_rule_id(cell_type, gene, "positive", type_counter[cell_type])

            mechanistic_basis = (
                f"Marker gene {gene} is expressed in {raw_cell} cells "
                f"as catalogued in CellMarker 2.0 (Hu et al. 2023, doi:10.1093/nar/gkac947). "
                f"Expert review required to add mechanistic basis before activation."
            )

            candidates.append(CandidateRule(
                cell_type=cell_type,
                rule_id=rule_id,
                rule_type="positive",
                evidence_source=evidence_source,
                gene_or_regulon=[gene],
                direction=direction,
                paired_with=[],
                incompatible_with=[],
                mechanistic_basis=mechanistic_basis,
                pmid=pmid,
                confidence="low",   # bumped to medium/high after expert review
                tissue_context=[tissue_raw] if tissue_raw else [],
                source_db="CellMarker2",
                status="PENDING_REVIEW",
                added_in_version=self.kg_version,
            ))

        return candidates


def load_candidates(path: str) -> tuple[dict, list[dict]]:
    """Load a candidates JSON file. Returns (metadata, candidates_list)."""
    with open(path) as f:
        payload = json.load(f)
    return payload.get("metadata", {}), payload.get("candidates", [])


def save_candidates(metadata: dict, candidates: list[dict], path: str) -> None:
    """Save candidates back to JSON (used by review CLI after editing)."""
    with open(path, "w") as f:
        json.dump({"metadata": metadata, "candidates": candidates}, f, indent=2)
