"""
AXIOM-SC CLI entry points.

Commands:
  axiom-annotate        — annotate an h5ad file
  axiom-kg              — manage the knowledge graph
  axiom-benchmark       — compare AXIOM-SC vs baseline tools
  axiom-train-tier1     — train AXIOMTier1 on CELLxGENE Census
  axiom-sc download-resources — download all Tier 3 reference databases (Patch B)

License: Apache-2.0
"""
from __future__ import annotations

import sys


def annotate() -> None:
    """CLI entry point: axiom-annotate"""
    print("axiom-annotate: annotate a .h5ad file with AXIOM-SC.")
    print("Usage: axiom-annotate --input data.h5ad --tiers 1 2 --profile oss-apache")
    print("Full implementation in Phase 2 playground (Day 11+).")


def kg() -> None:
    """CLI entry point: axiom-kg"""
    # Route to review_cli interactive tool or references generator
    args = sys.argv[1:]
    if "review" in args:
        from axiom_sc.kg.review_cli import ReviewCLI
        ReviewCLI().run()
    elif "references" in args:
        from axiom_sc.kg.references import generate_references_md
        generate_references_md()
    elif "download-resources" in args or (len(args) == 0 and False):
        _download_resources_cmd()
    else:
        print("axiom-kg commands: review | references")
        print("  axiom-kg review --input <candidates.json>")
        print("  axiom-kg references --output REFERENCES.md")


def benchmark() -> None:
    """CLI entry point: axiom-benchmark"""
    print("axiom-benchmark: compare AXIOM-SC vs CASSIA / mLLMCelltype.")
    print("Full implementation in Phase 2 Day 17+.")


def download_resources() -> None:
    """CLI entry point: axiom-sc download-resources"""
    _download_resources_cmd()


def _download_resources_cmd(tier3_only: bool = False, scenic_only: bool = False) -> None:
    """Download all Tier 3 reference databases for full coverage."""
    import subprocess
    import importlib.util
    from pathlib import Path

    script = Path(__file__).parent.parent / "scripts" / "patch_b_download_resources.py"
    if not script.exists():
        # Fallback: run inline
        print("Downloading AXIOM-SC Tier 3 resources ...")
        _run_inline_download(tier3_only=tier3_only, scenic_only=scenic_only)
        return

    result = subprocess.run(
        [sys.executable, str(script)],
        check=False,
    )
    sys.exit(result.returncode)


def _run_inline_download(tier3_only: bool = False, scenic_only: bool = False) -> None:
    """Run the Patch B resource download inline (no subprocess)."""
    from pathlib import Path
    import json
    import urllib.request

    resources_dir = Path.home() / ".axiom_sc" / "resources"
    resources_dir.mkdir(parents=True, exist_ok=True)

    checks = {
        "SCENIC rankings DB": resources_dir / "scenic" / "hg38_rankings.feather",
        "SCENIC motif DB":    resources_dir / "scenic" / "motifs-v9.tbl",
        "SCENIC TF list":     resources_dir / "scenic" / "allTFs_hg38.txt",
        "CellPhoneDB L-R":    resources_dir / "communication" / "cellphonedb_interactions.csv",
        "L-R gene index":     resources_dir / "communication" / "lr_gene_index.json",
        "ENSEMBL orthologs":  resources_dir / "orthology" / "hs_mm_orthologs.tsv",
        "Conservation index": resources_dir / "orthology" / "conservation_index.json",
        "CELLama embeddings": resources_dir / "spatial" / "cellama_reference_embeddings.pkl",
    }

    print("\n=== AXIOM-SC Tier 3 Resource Status ===")
    all_ok = True
    for name, path in checks.items():
        if path.exists():
            print(f"  ✓  {name} ({path.stat().st_size / 1e6:.1f} MB)")
        else:
            print(f"  ✗  {name}  ← MISSING")
            all_ok = False

    if not all_ok:
        print("\nRun: python scripts/patch_b_download_resources.py")
    else:
        print("\n✓ All resources present.")
