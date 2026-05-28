"""
PMID lookup and REFERENCES.md generator for the AXIOM-SC knowledge graph.

Every KG rule must have a non-empty pmid field. This module:
  1. Fetches APA-format citations from PubMed via biopython Entrez.
  2. Generates a REFERENCES.md file from all unique PMIDs in the KG.

Run::
    axiom-kg references --kg kg_data/oracle_kg_v0.2.0.json --output REFERENCES.md

License: Apache 2.0
"""
from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Optional

_ENTREZ_EMAIL = "oss@receptor.bio"
_RATE_LIMIT_S = 0.34   # NCBI allows max 3 requests/second without API key


def fetch_citation(pmid: str, retries: int = 3) -> str:
    """
    Fetch a formatted APA-style citation string from PubMed for a given PMID.

    Parameters
    ----------
    pmid : str
        PubMed ID (numeric string).
    retries : int
        Number of retry attempts on network failure.

    Returns
    -------
    str — formatted citation, e.g.:
      "Hori S, et al. (2003) Control of regulatory T cell development by the
       transcription factor Foxp3. Science 299:1057-1061. PMID:12612578"

    Raises
    ------
    RuntimeError if all retries fail.
    """
    try:
        from Bio import Entrez
        from Bio import Medline
    except ImportError:
        raise ImportError(
            "biopython required: pip install biopython. "
            "Or install science extras: pip install 'axiom-sc[science]'"
        )

    Entrez.email = _ENTREZ_EMAIL
    pmid = str(pmid).strip()

    last_exc: Exception | None = None
    for attempt in range(retries):
        try:
            handle = Entrez.efetch(
                db="pubmed", id=pmid, rettype="medline", retmode="text"
            )
            records = list(Medline.parse(handle))
            handle.close()

            if not records:
                return f"PMID:{pmid} (no record found)"

            rec = records[0]
            authors = rec.get("AU", [])
            if len(authors) > 3:
                author_str = f"{authors[0]}, et al."
            elif authors:
                author_str = "; ".join(authors)
            else:
                author_str = "Unknown authors"

            year = rec.get("DP", "")[:4] or "n.d."
            title = rec.get("TI", "No title").rstrip(".")
            journal = rec.get("TA", rec.get("JT", "Unknown journal"))
            volume = rec.get("VI", "")
            pages = rec.get("PG", "")
            vol_pg = f"{volume}:{pages}" if volume and pages else ""

            citation = f"{author_str} ({year}) {title}. {journal}"
            if vol_pg:
                citation += f" {vol_pg}"
            citation += f". PMID:{pmid}"
            return citation

        except Exception as exc:
            last_exc = exc
            if attempt < retries - 1:
                time.sleep(_RATE_LIMIT_S * (attempt + 1))

    return f"PMID:{pmid} (fetch failed: {last_exc})"


def fetch_citations_batch(
    pmids: list[str],
    delay: float = _RATE_LIMIT_S,
) -> dict[str, str]:
    """
    Fetch citations for a list of PMIDs, respecting NCBI rate limits.

    Returns
    -------
    dict[pmid -> citation_string]
    """
    results: dict[str, str] = {}
    unique = list(dict.fromkeys(p for p in pmids if p and p != "NEEDS_REVIEW"))

    for i, pmid in enumerate(unique):
        results[pmid] = fetch_citation(pmid)
        if i < len(unique) - 1:
            time.sleep(delay)

    return results


def generate_references_md(
    kg_path: str,
    output_path: str = "REFERENCES.md",
    include_header: bool = True,
) -> Path:
    """
    Generate REFERENCES.md from all unique PMIDs in a KG JSON file.

    Fetches each PMID from PubMed and writes a sorted reference list.

    Parameters
    ----------
    kg_path : str
        Path to the KG JSON file.
    output_path : str
        Output path for REFERENCES.md.
    include_header : bool
        Whether to prepend a header explaining the file's origin.

    Returns
    -------
    Path to the generated file.
    """
    with open(kg_path) as f:
        rules = json.load(f)

    # Collect all unique PMIDs from ACTIVE rules
    pmids = sorted(set(
        r["pmid"] for r in rules
        if r.get("status") == "ACTIVE"
        and r.get("pmid")
        and r["pmid"] != "NEEDS_REVIEW"
    ))

    print(f"[references] Fetching {len(pmids)} unique PMIDs from PubMed …")
    citations = fetch_citations_batch(pmids)

    lines: list[str] = []
    if include_header:
        lines += [
            "# AXIOM-SC Knowledge Graph — References",
            "",
            "Auto-generated from KG rule PMIDs. Do not edit manually.",
            f"Source KG: {kg_path}",
            f"Generated: see git log",
            "",
            "---",
            "",
        ]

    for pmid in pmids:
        citation = citations.get(pmid, f"PMID:{pmid}")
        lines.append(f"- {citation}")

    out = Path(output_path)
    out.write_text("\n".join(lines) + "\n")
    print(f"[references] Written {len(pmids)} references → {out}")
    return out


def check_missing_pmids(kg_path: str) -> list[str]:
    """
    Returns list of rule_ids that are missing or have placeholder PMIDs.
    Used as a pre-release check.
    """
    with open(kg_path) as f:
        rules = json.load(f)

    missing = [
        r["rule_id"] for r in rules
        if r.get("status") == "ACTIVE"
        and (not r.get("pmid") or r.get("pmid") == "NEEDS_REVIEW")
    ]
    return missing
