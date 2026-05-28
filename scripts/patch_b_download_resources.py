"""
Patch B — download all Tier 3 reference databases.

After running this script, every RNA-only dataset gets 4/6 Tier 3 streams
(communication + cross_species + spatial_niche + SCENIC). Velocity and
chromatin activate automatically when the user's data has the right layers.

Usage:
    /Users/sivamoturi/miniconda/bin/python3.10 scripts/patch_b_download_resources.py

Resources (downloaded once, live in ~/.axiom_sc/resources/):
  1. pySCENIC cisTarget rankings DB (~1.3 GB) — enables full cisTarget pruning
  2. pySCENIC motif-to-TF annotations (~50 MB)
  3. pySCENIC human TF list
  4. CellPhoneDB v4 L-R interactions (~2,300+ pairs) → lr_gene_index.json
  5. ENSEMBL human-mouse ortholog table (~50 MB) → conservation_index.json
  6. CELLama reference embeddings (~2 MB) — enables embedding-mode in SpatialNicheStream

License notes:
  CellPhoneDB: CC BY 4.0 (Efremova et al. 2020 doi:10.1038/s41592-020-0906-9)
  ENSEMBL:     CC0 public domain (doi:10.1093/nar/gkac958)
  CELLama:     MIT (doi:10.1002/advs.202413514)
  pySCENIC:    GPL v3 (used via subprocess only in axiom_sc)
"""
from __future__ import annotations

import gzip
import json
import shutil
import sys
import urllib.request
from pathlib import Path

import pandas as pd

RESOURCES_DIR = Path.home() / ".axiom_sc" / "resources"
RESOURCES_DIR.mkdir(parents=True, exist_ok=True)


def download_file(url: str, dest: Path, desc: str = "") -> bool:
    """Download url to dest, skipping if dest already exists. Returns True if downloaded."""
    if dest.exists():
        print(f"  ✓ Already present: {dest.name} ({dest.stat().st_size / 1e6:.1f} MB)")
        return False
    print(f"  Downloading {desc or dest.name} ...")
    try:
        urllib.request.urlretrieve(url, dest)
        size_mb = dest.stat().st_size / 1e6
        print(f"  ✓ Saved: {dest.name} ({size_mb:.1f} MB)")
        return True
    except Exception as exc:
        print(f"  ✗ FAILED: {exc}")
        if dest.exists():
            dest.unlink()
        return False


# ── 1. pySCENIC resources ─────────────────────────────────────────────────────

SCENIC_DIR = RESOURCES_DIR / "scenic"
SCENIC_DIR.mkdir(exist_ok=True)

print("\n=== pySCENIC resources ===")
download_file(
    "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/"
    "mc9nr/gene_based/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather",
    SCENIC_DIR / "hg38_rankings.feather",
    "pySCENIC cisTarget rankings DB (1.3 GB)",
)
download_file(
    "https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl",
    SCENIC_DIR / "motifs-v9.tbl",
    "pySCENIC motif-to-TF annotations (50 MB)",
)
download_file(
    "https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt",
    SCENIC_DIR / "allTFs_hg38.txt",
    "pySCENIC human TF list",
)

# ── 2. CellPhoneDB L-R interactions ──────────────────────────────────────────

COMMOT_DIR = RESOURCES_DIR / "communication"
COMMOT_DIR.mkdir(exist_ok=True)

print("\n=== CellPhoneDB L-R interactions ===")

# Try CellPhoneDB v5 interactions CSV from GitHub release
CELLPHONEDB_URLS = [
    # CellPhoneDB v5 GitHub release artifact
    "https://github.com/ventolab/CellphoneDB/raw/master/cellphonedb/src/core/data/interaction_table.csv",
    # CellPhoneDB v4 from PyPI package data (alternative mirror)
    "https://raw.githubusercontent.com/ventolab/CellphoneDB/v4.0.0/cellphonedb/src/core/data/interaction_table.csv",
]

lr_csv_path = COMMOT_DIR / "cellphonedb_interactions.csv"
downloaded = False
for url in CELLPHONEDB_URLS:
    downloaded = download_file(url, lr_csv_path, "CellPhoneDB L-R interactions")
    if downloaded or lr_csv_path.exists():
        break

lr_index_path = COMMOT_DIR / "lr_gene_index.json"
if lr_csv_path.exists() and not lr_index_path.exists():
    print("  Building L-R gene index ...")
    try:
        df = pd.read_csv(lr_csv_path)
        lr_dict: dict = {}
        gene_cols = [c for c in ["partner_a", "partner_b", "gene_name_a", "gene_name_b"] if c in df.columns]
        for col in gene_cols:
            for gene in df[col].dropna().astype(str).str.strip():
                if gene and gene != "nan":
                    lr_dict[gene] = lr_dict.get(gene, [])
        with open(lr_index_path, "w") as f:
            json.dump(lr_dict, f)
        print(f"  ✓ L-R index: {len(lr_dict)} unique genes → {lr_index_path.name}")
    except Exception as exc:
        print(f"  ✗ Could not build L-R index: {exc}")

# ── 3. ENSEMBL human-mouse ortholog table ─────────────────────────────────────

ORTHO_DIR = RESOURCES_DIR / "orthology"
ORTHO_DIR.mkdir(exist_ok=True)

print("\n=== ENSEMBL ortholog table ===")

gz_path = ORTHO_DIR / "hs_orthologs.tsv.gz"
tsv_path = ORTHO_DIR / "hs_mm_orthologs.tsv"
conservation_path = ORTHO_DIR / "conservation_index.json"

if not tsv_path.exists():
    ENSEMBL_URL = (
        "https://ftp.ensembl.org/pub/release-110/tsv/homo_sapiens/"
        "Homo_sapiens.GRCh38.110.homologs.tsv.gz"
    )
    ok = download_file(ENSEMBL_URL, gz_path, "ENSEMBL human-mouse orthologs (50 MB)")
    if ok and gz_path.exists():
        print("  Decompressing ...")
        with gzip.open(gz_path, "rb") as gz_in, open(tsv_path, "wb") as tsv_out:
            shutil.copyfileobj(gz_in, tsv_out)
        gz_path.unlink()
        print(f"  ✓ Decompressed to {tsv_path.name}")
else:
    print(f"  ✓ Already present: {tsv_path.name}")

if tsv_path.exists() and not conservation_path.exists():
    print("  Building conservation index ...")
    try:
        df = pd.read_csv(tsv_path, sep="\t", usecols=lambda c: c in [
            "Gene name", "Mouse gene name", "Mouse homology type"
        ])
        df = df.dropna()
        one2one = set(df[df["Mouse homology type"] == "ortholog_one2one"]["Gene name"].unique())
        any_ortho = set(df["Gene name"].unique())
        index = {gene: ("one2one" if gene in one2one else "many") for gene in any_ortho}
        with open(conservation_path, "w") as f:
            json.dump(index, f)
        print(f"  ✓ Conservation index: {len(index)} genes, {len(one2one)} one-to-one orthologs")
    except Exception as exc:
        print(f"  ✗ Could not build conservation index: {exc}")

# ── 4. CELLama reference embeddings ──────────────────────────────────────────

SPATIAL_DIR = RESOURCES_DIR / "spatial"
SPATIAL_DIR.mkdir(exist_ok=True)

print("\n=== CELLama reference embeddings ===")
CELLAMA_URLS = [
    "https://github.com/portrai-io/CELLama/raw/main/data/reference_embeddings.pkl",
    "https://github.com/portrai-io/CELLama/raw/main/ref_embeddings/cellama_reference_embeddings.pkl",
]
embeddings_path = SPATIAL_DIR / "cellama_reference_embeddings.pkl"
if not embeddings_path.exists():
    for url in CELLAMA_URLS:
        ok = download_file(url, embeddings_path, "CELLama reference embeddings (~2 MB)")
        if embeddings_path.exists():
            break
else:
    print(f"  ✓ Already present: {embeddings_path.name}")

# ── 5. Resource status summary ────────────────────────────────────────────────

def resource_status() -> bool:
    checks = {
        "SCENIC rankings DB (1.3 GB)": SCENIC_DIR / "hg38_rankings.feather",
        "SCENIC motif DB":             SCENIC_DIR / "motifs-v9.tbl",
        "SCENIC TF list":              SCENIC_DIR / "allTFs_hg38.txt",
        "CellPhoneDB L-R pairs":       COMMOT_DIR / "cellphonedb_interactions.csv",
        "L-R gene index":              COMMOT_DIR / "lr_gene_index.json",
        "ENSEMBL orthologs":           ORTHO_DIR / "hs_mm_orthologs.tsv",
        "Conservation index":          ORTHO_DIR / "conservation_index.json",
        "CELLama embeddings":          SPATIAL_DIR / "cellama_reference_embeddings.pkl",
    }
    print("\n=== AXIOM-SC Tier 3 Resource Status ===")
    all_present = True
    for name, path in checks.items():
        if path.exists():
            size_mb = path.stat().st_size / 1e6
            print(f"  ✓  {name} ({size_mb:.1f} MB)")
        else:
            print(f"  ✗  {name}  ← MISSING")
            all_present = False

    if all_present:
        print("\n✓ All Tier 3 resources present — full 4-stream RNA-only coverage enabled")
    else:
        print("\n⚠ Some resources missing — check download errors above")
    return all_present


if __name__ == "__main__":
    resource_status()
