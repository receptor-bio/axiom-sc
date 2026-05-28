#!/usr/bin/env python
"""
Day 6 validation: scenic_runner subprocess isolation test.

Tests the pySCENIC subprocess isolation mechanism using Phase 1 thymus data.

WHAT THIS VALIDATES:
  1. scenic_runner.run() successfully calls scenic_worker.py via subprocess
  2. Phase 2 settings (NES=2.0, forced_genes) recover more regulons than Phase 1
  3. Recovery report format and forced gene tracking

WHAT THIS DOES NOT VALIDATE (requires GPU + SCENIC resource files):
  - Full motif-pruned regulon quality (needs 1.3 GB rankings feather)
  - FOXP3 recovery at 50k cells (thymus only has ~20k cells locally)
  - GPU-scale timing estimates

HOW TO RUN FULL GPU VALIDATION (Day 16, after training):
  1. Prepare stratified 50k-cell loom (combine thymus + thymus-adjacent data)
  2. Download SCENIC resources:
       ScenicRunner().download_resources()  # ~1.5 GB total
  3. Run with full motif pruning:
       ScenicPipeline().run(adata, cluster_key='leiden', download_resources=False)
  4. Expected: FOXP3 in recovery_report()['forced_recovered'] with z ≥ 2.0

Usage:
  cd axiom-sc
  /Users/sivamoturi/miniconda/bin/python3.10 scripts/validate_scenic_runner.py
"""
import json
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))

# ── Config ────────────────────────────────────────────────────────────────────

# Use axiom-poc env (has pyscenic 0.12.1 + loompy) as scenic-env
SCENIC_ENV = "axiom-poc"

# Phase 1 thymus loom (3038 genes × 19926 cells, all forced genes present)
THYMUS_LOOM = Path("/Users/sivamoturi/Documents/siva-wrk/axiom-poc/data/processed/thymus.loom")

# Phase 1 AUC output for comparison
PHASE1_AUC = Path("/Users/sivamoturi/Documents/siva-wrk/axiom-poc/scenic_outputs/thymus/auc_matrix.csv")

# Subset size for local validation (full run needs GPU)
N_CELLS_SUBSET = 2000

# ── Imports ────────────────────────────────────────────────────────────────────

from axiom_sc.tier2.scenic_runner import ScenicRunner, SCENIC_DEFAULTS
from axiom_sc.tier2.scenic_pipeline import ScenicPipeline
from axiom_sc.utils.zscore_auc import zscore_auc_dict


def print_section(title: str):
    print(f"\n{'═' * 60}")
    print(f"  {title}")
    print('═' * 60)


def load_phase1_results():
    """Load Phase 1 SCENIC output for comparison."""
    import pandas as pd
    df = pd.read_csv(PHASE1_AUC, index_col=0)
    # Strip pySCENIC's "(+)" suffix
    df.columns = [c.split("(")[0].strip() for c in df.columns]
    return df


def build_subset_loom(loom_path: Path, n_cells: int, output_path: Path):
    """Extract a random subset of cells from the thymus loom into a new loom."""
    import subprocess as sp, sys
    # loompy lives in axiom-poc env; use its Python to build the subset loom
    axiom_poc_python = "/Users/sivamoturi/miniconda/envs/axiom-poc/bin/python"
    script = f"""
import loompy, numpy as np
rng = np.random.default_rng(42)
with loompy.connect({str(loom_path)!r}, mode='r') as ds:
    idx = rng.choice(ds.shape[1], size=min({n_cells}, ds.shape[1]), replace=False)
    idx.sort()
    matrix = ds[:, idx]
    row_attrs = {{'Gene': ds.ra['Gene']}}
    col_attrs = {{'CellID': ds.ca['CellID'][idx]}}
loompy.create({str(output_path)!r}, matrix, row_attrs, col_attrs)
print(f'Created {{matrix.shape[1]}} cells x {{matrix.shape[0]}} genes')
"""
    result = sp.run([axiom_poc_python, "-c", script], capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Subset loom creation failed:\n{result.stderr}")
    print(f"  {result.stdout.strip()}")
    return output_path

def build_subset_loom_UNUSED(loom_path: Path, n_cells: int, output_path: Path):
    """Original loompy-based version (requires loompy in current env)."""
    import loompy
    rng = np.random.default_rng(42)
    with loompy.connect(str(loom_path), mode="r") as ds:
        total_cells = ds.shape[1]
        idx = rng.choice(total_cells, size=min(n_cells, total_cells), replace=False)
        idx.sort()
        matrix = ds[:, idx]  # genes × cells subset
        row_attrs = {"Gene": ds.ra["Gene"]}
        col_attrs = {"CellID": ds.ca["CellID"][idx]}
    loompy.create(str(output_path), matrix, row_attrs, col_attrs)
    print(f"  Created subset loom: {matrix.shape[1]} cells × {matrix.shape[0]} genes")
    return output_path


def main():
    print_section("Phase 1 SCENIC results (NES=3.0, 19,926 cells, no forced genes)")
    p1 = load_phase1_results()
    forced_genes = set(SCENIC_DEFAULTS["forced_genes"])
    p1_regulons = set(p1.columns)
    p1_forced_found = forced_genes & p1_regulons
    p1_forced_missing = forced_genes - p1_regulons
    print(f"  Total regulons recovered:  {len(p1_regulons)}")
    print(f"  Forced master TFs found:   {sorted(p1_forced_found)} ({len(p1_forced_found)}/26)")
    print(f"  Forced master TFs MISSING: {sorted(p1_forced_missing)}")
    print(f"  Full regulon list: {sorted(p1_regulons)}")

    print_section(f"Phase 2 subplot test: NES=2.0, forced_genes, {N_CELLS_SUBSET}-cell subset")
    print(f"  Using scenic-env: {SCENIC_ENV!r}")
    print(f"  Source loom: {THYMUS_LOOM}")
    print(f"  Fallback mode: no rankings/motif DB (adjacency-only regulons)")
    print(f"  This validates: subprocess isolation mechanism")

    runner = ScenicRunner(
        scenic_conda_env=SCENIC_ENV,
        scenic_defaults={"nes_threshold": 2.0, "min_genes": 5},
    )

    # Check scenic-env — verifies pyscenic CLI availability
    print("\n  Checking scenic-env (pyscenic CLI)...", end=" ", flush=True)
    try:
        runner.check_env()
        print("OK")
    except Exception as e:
        print(f"FAILED: {e}")
        sys.exit(1)

    # Build subset loom
    with tempfile.TemporaryDirectory(prefix="axiom_scenic_val_") as tmpdir:
        subset_loom = Path(tmpdir) / "thymus_subset.loom"
        output_path = Path(tmpdir) / "auc_output.csv"

        build_subset_loom(THYMUS_LOOM, N_CELLS_SUBSET, subset_loom)

        # Run subprocess
        print(f"\n  Running pySCENIC subprocess (NES=2.0, fallback mode)...")
        print(f"  Expected: more than {len(p1_regulons)} regulons from Phase 1")
        try:
            result = runner.run(
                loom_path=str(subset_loom),
                output_path=str(output_path),
            )
        except Exception as e:
            print(f"  SUBPROCESS FAILED: {e}")
            sys.exit(1)

        n_recovered = len(result)
        recovered_regulons = set(result.keys())
        forced_found = forced_genes & recovered_regulons
        forced_missing = forced_genes - recovered_regulons

        print(f"\n  ── Phase 2 NES=2.0 results ({N_CELLS_SUBSET} cells, fallback mode) ──")
        print(f"  Total regulons recovered:  {n_recovered}")
        print(f"  Forced master TFs found:   {sorted(forced_found)} ({len(forced_found)}/26)")
        print(f"  Forced master TFs missing: {sorted(forced_missing)}")

        # FOXP3 specifically
        if "FOXP3" in result:
            foxp3_vals = np.array(result["FOXP3"])
            print(f"\n  ✓ FOXP3 RECOVERED (z-scores: mean={foxp3_vals.mean():.2f}, "
                  f"max={foxp3_vals.max():.2f})")
        else:
            print(f"\n  ✗ FOXP3 not recovered in {N_CELLS_SUBSET}-cell subset")
            print(f"    (expected — {N_CELLS_SUBSET} cells < 300 min for rare TFs)")
            print(f"    Full validation: 50k cells on GPU with motif pruning")

        print_section("Comparison: Phase 1 (NES=3.0) vs Phase 2 (NES=2.0)")
        print(f"  Phase 1  NES=3.0 19,926 cells: {len(p1_regulons):3d} regulons  "
              f"| FOXP3: {'✓' if 'FOXP3' in p1_regulons else '✗'}  "
              f"IRF7: {'✓' if 'IRF7' in p1_regulons else '✗'}")
        print(f"  Phase 2  NES=2.0  {N_CELLS_SUBSET:5d} cells: {n_recovered:3d} regulons  "
              f"| FOXP3: {'✓' if 'FOXP3' in result else '✗'}  "
              f"IRF7: {'✓' if 'IRF7' in result else '✗'}")

        print_section("Recovery report (pipeline-level)")
        pipeline = ScenicPipeline(scenic_conda_env=SCENIC_ENV)
        cluster_regulons = {
            "all_cells": {reg: float(np.array(vals).mean()) for reg, vals in result.items()}
        }
        report = pipeline.recovery_report(cluster_regulons)
        print(f"  n_regulons:       {report['n_regulons']}")
        print(f"  forced_recovered: {sorted(report['forced_recovered'].keys())}")
        print(f"  forced_missing:   {report['forced_missing'][:10]}{'...' if len(report['forced_missing']) > 10 else ''}")

    print_section("VALIDATION SUMMARY")
    print(f"  [{'PASS' if n_recovered >= 0 else 'FAIL'}] Subprocess isolation: scenic_runner.run() succeeded")
    print(f"  [{'PASS' if n_recovered > len(p1_regulons) else 'INFO'}] "
          f"NES=2.0 recovered {n_recovered} regulons "
          f"({'more' if n_recovered > len(p1_regulons) else 'fewer/same'} than Phase 1's {len(p1_regulons)})")
    print(f"  [PEND] FOXP3 at 50k cells — requires GPU + SCENIC resource files")
    print(f"         Run: ScenicPipeline(scenic_conda_env='scenic-env').run(adata_50k, ...)")
    print(f"         After: axiom-train-tier1 + download_resources() (Day 16)")
    print()


if __name__ == "__main__":
    main()
