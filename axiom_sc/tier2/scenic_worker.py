#!/usr/bin/env python
"""
pySCENIC worker — runs INSIDE scenic-env (GPL v3 side).

THIS FILE IS NOT PART OF axiom_sc PACKAGE.
It is executed only via subprocess from scenic_runner.py.
Never import this file from axiom_sc.

pySCENIC license: GPL v3
This file: GPL v3 by isolation (called from Apache 2.0 code via subprocess only).

Uses the pySCENIC CLI commands (pyscenic grn / ctx / aucell) rather than
the Python API. This matches the Phase 1 Colab pipeline exactly and is robust
across pySCENIC versions (the CLI interface is more stable than the Python API).

Three-step pySCENIC pipeline:
  1. GRN:    pyscenic grn   → adj.csv (TF → target gene weights)
  2. ctx:    pyscenic ctx   → regulons.csv (motif-pruned regulon definitions)
             Fallback: modules_from_adjacencies if rankings/motif DB unavailable
  3. AUCell: pyscenic aucell → auc.loom (per-cell regulon activity)
             Extraction via h5py (loompy 3.x rejects pySCENIC's structured dtype)

Then z-scores per regulon across all cells and returns JSON to stdout.

Phase 1 finding: NES=2.0 + 20k cells → 17 regulons, FOXP3 missing.
Phase 2 target:  NES=2.0 + 50k cells → FOXP3 recovered (expressed in <2% of cells,
                 needs enough cells for GRN signal).

Usage (internal, called by ScenicRunner):
  python scenic_worker.py '<config_json>'
"""
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path


def run_cmd(cmd: list[str], desc: str, env=None) -> None:
    """Run a CLI command; raise RuntimeError on non-zero exit."""
    result = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if result.returncode != 0:
        raise RuntimeError(
            f"{desc} failed (exit {result.returncode}):\n"
            f"stdout: {result.stdout[-2000:]}\n"
            f"stderr: {result.stderr[-2000:]}"
        )


def extract_auc_from_loom(auc_loom_path: str) -> tuple[list[str], "np.ndarray", list[str]]:
    """
    Extract AUC matrix from pySCENIC aucell output loom via h5py.

    loompy 3.x rejects pySCENIC's structured numpy array dtype.
    h5py reads it directly — matches Phase 1 Colab extraction method.

    Returns (cell_ids, matrix, regulon_names)
    where matrix has shape (n_cells, n_regulons).
    """
    import h5py
    import numpy as np

    with h5py.File(auc_loom_path, "r") as f:
        cell_ids = f["col_attrs"]["CellID"][:].astype(str).tolist()
        auc_data = f["col_attrs"]["RegulonsAUC"][:]

    regulon_names = list(auc_data.dtype.names)
    # Strip pySCENIC's "(+)" suffix from regulon names
    clean_names = [n.split("(")[0].strip() for n in regulon_names]
    matrix = np.column_stack([auc_data[r] for r in regulon_names])
    return cell_ids, matrix, clean_names


def adjacency_fallback_regulons(adj_csv: str, min_genes: int, out_gmt: str, loom_path: str) -> bool:
    """
    Build naive regulons from grnboost2 adjacencies when motif DB is unavailable.
    Uses pyscenic.utils.modules_from_adjacencies (from arboreto adjacency format).

    Writes a GMT file (tab-separated: name<tab>description<tab>gene1<tab>gene2...).
    pyscenic aucell accepts GMT via the .gmt extension in load_signatures().

    Returns True if regulons were written, False if no regulons found.
    """
    try:
        import warnings
        warnings.filterwarnings("ignore")
        import pandas as pd
        import h5py
        from pyscenic.utils import modules_from_adjacencies
    except ImportError as e:
        print(f"Warning: fallback import failed: {e}", file=sys.stderr)
        return False

    # Load expression matrix (cells × genes) required by modules_from_adjacencies
    with h5py.File(loom_path, "r") as f:
        gene_names = f["row_attrs"]["Gene"][:].astype(str).tolist()
        cell_ids = f["col_attrs"]["CellID"][:].astype(str).tolist()
        matrix = f["matrix"][:].T  # loom is genes×cells; transpose to cells×genes
    ex_mtx = pd.DataFrame(matrix, index=cell_ids, columns=gene_names)

    adj = pd.read_csv(adj_csv, index_col=False)
    modules = list(modules_from_adjacencies(adj, ex_mtx, min_genes=min_genes))
    if not modules:
        return False

    # Merge modules with the same TF into a single gene set to avoid duplicate
    # column names (which cause pandas to return DataFrame instead of Series
    # when indexing, breaking pyscenic aucell's binarization step).
    tf_genes: dict[str, set] = {}
    for m in modules:
        tf = m.transcription_factor
        tf_genes.setdefault(tf, set()).update(m.genes)

    # Write GMT format: name<tab>description<tab>gene1<tab>gene2...
    # pyscenic aucell's load_signatures() reads .gmt via GeneSignature.from_gmt()
    with open(out_gmt, "w") as fh:
        for tf, genes in sorted(tf_genes.items()):
            line = "\t".join([f"{tf}(+)", "NA"] + sorted(genes))
            fh.write(line + "\n")

    n_modules = len(modules)
    print(f"Fallback: wrote {n_modules} adjacency-based regulons to {out_gmt}", file=sys.stderr)
    return True


def main():
    if len(sys.argv) < 2:
        print(json.dumps({"error": "No config JSON provided"}), file=sys.stderr)
        sys.exit(1)

    config = json.loads(sys.argv[1])
    loom_path = config["loom_path"]
    output_path = config.get("output_path", loom_path + "_auc.csv")
    output_dir = str(Path(output_path).parent)

    nes_threshold = str(config.get("nes_threshold", 2.0))
    min_genes = str(config.get("min_genes", 5))
    n_workers = str(config.get("n_workers", 4))
    z_score_auc = config.get("z_score_auc", True)

    tf_list_path = config.get("tf_list_path")
    rankings_db_path = config.get("rankings_db_path")
    motif_db_path = config.get("motif_db_path")

    adj_csv = os.path.join(output_dir, "adj.csv")
    regulons_csv = os.path.join(output_dir, "regulons.csv")   # ctx output (motif-enrichment CSV)
    regulons_gmt = os.path.join(output_dir, "regulons.gmt")   # fallback output (GMT)
    auc_loom = os.path.join(output_dir, "auc.loom")
    signatures_path = None  # set after Step 2 to whichever file was produced

    # Dask config to avoid pySCENIC planning overhead
    env = dict(os.environ, DASK_DATAFRAME__QUERY_PLANNING="False")

    # ── Step 1: GRN inference ─────────────────────────────────────────────────
    if os.path.exists(adj_csv) and os.path.getsize(adj_csv) > 1000:
        print(f"[SKIP] adj.csv exists", file=sys.stderr)
    else:
        grn_cmd = [
            "pyscenic", "grn",
            "--num_workers", n_workers,
            "--output", adj_csv,
            "--method", "grnboost2",
            loom_path,
        ]
        if tf_list_path and os.path.exists(tf_list_path):
            grn_cmd.append(tf_list_path)
        else:
            # Check for cached TF list in default resources directory
            cached_tf = Path.home() / ".axiom_sc" / "scenic_resources" / "allTFs_hg38.txt"
            if cached_tf.exists():
                print(f"Using cached TF list: {cached_tf}", file=sys.stderr)
                grn_cmd.append(str(cached_tf))
            else:
                # Extract gene names from loom and write to temp file as TF list
                # (treats all genes as potential TFs — slower/noisier but works)
                import tempfile
                try:
                    import loompy
                    with loompy.connect(loom_path, mode="r") as ds:
                        gene_names = list(ds.ra["Gene"])
                except Exception:
                    import h5py
                    with h5py.File(loom_path, "r") as f:
                        gene_names = f["row_attrs"]["Gene"][:].astype(str).tolist()
                tf_tmp = tempfile.NamedTemporaryFile(
                    mode="w", suffix=".txt", delete=False, prefix="axiom_tf_"
                )
                tf_tmp.write("\n".join(gene_names))
                tf_tmp.close()
                print(
                    f"Warning: no TF list provided — using all {len(gene_names)} genes as TFs. "
                    f"Download resources with: ScenicRunner().download_resources()",
                    file=sys.stderr,
                )
                grn_cmd.append(tf_tmp.name)

        try:
            run_cmd(grn_cmd, "pyscenic grn", env=env)
        except RuntimeError as e:
            print(json.dumps({"error": f"GRN step failed: {e}"}), file=sys.stderr)
            sys.exit(1)

    # ── Step 2: cisTarget regulon pruning ─────────────────────────────────────
    if os.path.exists(regulons_csv) and os.path.getsize(regulons_csv) > 100:
        print("[SKIP] regulons.csv exists", file=sys.stderr)
        signatures_path = regulons_csv
    elif os.path.exists(regulons_gmt) and os.path.getsize(regulons_gmt) > 10:
        print("[SKIP] regulons.gmt exists", file=sys.stderr)
        signatures_path = regulons_gmt
    elif rankings_db_path and motif_db_path and \
            os.path.exists(rankings_db_path) and os.path.exists(motif_db_path):
        # Full cisTarget pruning with motif enrichment
        ctx_cmd = [
            "pyscenic", "ctx",
            adj_csv, rankings_db_path,
            "--annotations_fname",    motif_db_path,
            "--expression_mtx_fname", loom_path,
            "--output",               regulons_csv,
            "--num_workers",          n_workers,
            "--nes_threshold",        nes_threshold,
            "--min_genes_regulon",    min_genes,
            "--min_targets",          min_genes,
        ]
        try:
            run_cmd(ctx_cmd, "pyscenic ctx", env=env)
            signatures_path = regulons_csv
        except RuntimeError as e:
            print(f"cisTarget failed: {e}\nFalling back to adjacency-only regulons",
                  file=sys.stderr)
            if not adjacency_fallback_regulons(adj_csv, int(min_genes), regulons_gmt, loom_path):
                print(json.dumps({}))
                return
            signatures_path = regulons_gmt
    else:
        # Fallback: no rankings/motif DB — build regulons from adjacency weights
        print(
            "Warning: SCENIC resource files not provided. "
            "Running in adjacency-only fallback mode (no motif pruning). "
            "Regulon quality will be lower. "
            "Download resources with: ScenicRunner().download_resources()",
            file=sys.stderr,
        )
        if not adjacency_fallback_regulons(adj_csv, int(min_genes), regulons_gmt, loom_path):
            print(json.dumps({}))
            return
        signatures_path = regulons_gmt

    # ── Step 3: AUCell scoring ────────────────────────────────────────────────
    if os.path.exists(auc_loom) and os.path.getsize(auc_loom) > 1000:
        print("[SKIP] auc.loom exists", file=sys.stderr)
    else:
        aucell_cmd = [
            "pyscenic", "aucell",
            loom_path, signatures_path,
            "--output",      auc_loom,
            "--num_workers", n_workers,
        ]
        try:
            run_cmd(aucell_cmd, "pyscenic aucell", env=env)
        except RuntimeError as e:
            print(json.dumps({"error": f"AUCell step failed: {e}"}), file=sys.stderr)
            sys.exit(1)

    # ── Step 4: Extract AUC matrix via h5py ───────────────────────────────────
    try:
        import numpy as np
        cell_ids, matrix, regulon_names = extract_auc_from_loom(auc_loom)
    except Exception as e:
        print(json.dumps({"error": f"AUC extraction failed: {e}"}), file=sys.stderr)
        sys.exit(1)

    # Save raw AUC CSV (for debugging / archive)
    try:
        import pandas as pd
        df = pd.DataFrame(matrix, index=cell_ids, columns=regulon_names)
        df.index.name = "CellID"
        df.to_csv(output_path)
    except Exception:
        pass  # CSV save is best-effort

    # ── Step 5: Z-score per regulon ───────────────────────────────────────────
    result: dict[str, list[float]] = {}
    if z_score_auc:
        for i, regulon in enumerate(regulon_names):
            col = matrix[:, i].astype(float)
            mean = col.mean()
            std = col.std(ddof=1)
            if std < 1e-10:
                z = [0.0] * len(col)
            else:
                z = ((col - mean) / std).tolist()
            result[regulon] = z
    else:
        for i, regulon in enumerate(regulon_names):
            result[regulon] = matrix[:, i].tolist()

    print(json.dumps(result))


if __name__ == "__main__":
    main()
