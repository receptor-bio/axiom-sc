"""
pySCENIC subprocess runner — Apache 2.0 side.

pySCENIC is GPL v3 and must NEVER be imported into axiom_sc.
This module calls scenic_worker.py inside a separate conda environment (scenic-env)
via subprocess, keeping the GPL code fully isolated.

Phase 1 validated defaults are in SCENIC_DEFAULTS below.
See CLAUDE.md Section 9 for the rationale behind each default.

References:
  Van de Sande et al. (2020) Nature Protocols 15:2247
  doi: 10.1038/s41596-020-0336-2

License: Apache 2.0 — this file only. scenic_worker.py is GPL-isolated.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import anndata

# Phase 1 validated SCENIC defaults (see CLAUDE.md Section 9 for corrections)
SCENIC_DEFAULTS: dict = {
    # Phase 1 lesson: NES 3.0 recovers ZERO immune master TFs at 20k cells
    # NES 2.0 recovers IRF7, GATA3, RORC
    "nes_threshold": 2.0,       # was 3.0 in pySCENIC default

    # Phase 1 lesson: default min_genes=30 too strict for rare TFs
    "min_genes": 5,

    # Phase 1 lesson: HVG selection excludes master TFs (low variance in mixed datasets)
    # These MUST be added to loom input regardless of HVG rank
    "forced_genes": [
        "FOXP3", "TBX21", "GATA3", "RORC", "IRF7", "PAX5", "EBF1",
        "AIRE", "DLL4", "PSMB11", "TOX", "TCF7", "NR4A1", "RUNX3",
        "IKZF2", "IKZF3", "BCL6", "PRDM1", "XBP1", "CEBPA", "SPI1",
        "IRF4", "IRF8", "BATF3", "ID2", "EOMES",
    ],

    # Phase 1 lesson: raw AUC values not comparable across runs/datasets
    "z_score_auc": True,
    "z_active_thresh": 1.5,
    "z_inactive_thresh": -0.5,

    # Phase 1 lesson: 20k cells insufficient for rare TF recovery
    "min_cells_per_type": 300,
    "recommended_total": 50000,

    # SCENIC resource URLs
    "tf_list_url": (
        "https://resources.aertslab.org/cistarget/tf_lists/"
        "allTFs_hg38.txt"
    ),
    "rankings_db_url": (
        "https://resources.aertslab.org/cistarget/databases/"
        "homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/"
        "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
    ),
    "motif_db_url": (
        "https://resources.aertslab.org/cistarget/motif2tf/"
        "motifs-v9-nr.hgnc-m0.001-o0.0.tbl"
    ),
}


class ScenicRunner:
    """
    Runs pySCENIC inside scenic-env via subprocess.
    Returns z-scored regulon AUC values as dict[regulon -> z_score_per_cell].

    Requires:
      - conda env scenic-env with pyscenic installed
      - scenic_worker.py accessible inside scenic-env
    """

    def __init__(
        self,
        scenic_conda_env: str | None = None,
        scenic_defaults: dict | None = None,
        resources_dir: str | None = None,
    ):
        self.conda_env = scenic_conda_env or os.environ.get("SCENIC_CONDA_ENV", "scenic-env")
        self.defaults = {**SCENIC_DEFAULTS, **(scenic_defaults or {})}
        self.resources_dir = Path(resources_dir or os.environ.get(
            "SCENIC_RESOURCES_DIR",
            Path.home() / ".axiom_sc" / "scenic_resources",
        ))

    # ── environment checks ────────────────────────────────────────────────────

    @staticmethod
    def _find_conda() -> str:
        """
        Return the conda executable path, searching PATH and common install locations.
        Raises ScenicEnvNotFoundError if conda cannot be found.
        """
        conda = shutil.which("conda")
        if conda:
            return conda
        # Common locations when conda is not on PATH (e.g. invoked via absolute Python)
        for candidate in [
            Path.home() / "miniconda" / "bin" / "conda",
            Path.home() / "miniconda3" / "bin" / "conda",
            Path.home() / "anaconda3" / "bin" / "conda",
            Path("/usr/local/miniconda") / "bin" / "conda",
            Path("/opt/conda") / "bin" / "conda",
        ]:
            if candidate.exists():
                return str(candidate)
        raise ScenicEnvNotFoundError(
            "conda not found on PATH or in common install locations "
            "(~/miniconda, ~/miniconda3, ~/anaconda3). "
            "Install conda/miniconda to use pySCENIC."
        )

    def check_env(self) -> bool:
        """
        Return True if the SCENIC conda env exists and the pyscenic CLI is available.
        Raises ScenicEnvNotFoundError if env or pyscenic CLI is missing.

        Checks for `pyscenic` CLI (used by scenic_worker.py via subprocess calls),
        not the Python import. This matches the Phase 1 pipeline approach.
        """
        conda = self._find_conda()
        # pyscenic has no --version; use grn --help (exit 0) to verify CLI works
        result = subprocess.run(
            [conda, "run", "-n", self.conda_env, "pyscenic", "grn", "--help"],
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            raise ScenicEnvNotFoundError(
                f"scenic-env conda environment not found or pyscenic CLI not installed.\n"
                f"Create it with:\n"
                f"  conda create -n {self.conda_env} python=3.10 -y\n"
                f"  conda activate {self.conda_env}\n"
                f"  pip install pyscenic loompy h5py\n"
                f"stderr: {result.stderr.strip()}"
            )
        version = (result.stdout + result.stderr).strip()
        return True

    def _get_conda_python(self) -> str:
        """Return the Python executable path inside the scenic conda env."""
        try:
            conda = self._find_conda()
        except ScenicEnvNotFoundError:
            return "python"
        result = subprocess.run(
            [conda, "run", "-n", self.conda_env, "python", "-c",
             "import sys; print(sys.executable)"],
            capture_output=True,
            text=True,
        )
        if result.returncode == 0:
            return result.stdout.strip()
        return "python"

    # ── resource management ───────────────────────────────────────────────────

    def download_resources(self, force: bool = False) -> dict[str, Path]:
        """
        Download SCENIC resource files (TF list, rankings DB, motif DB) to
        self.resources_dir if not already cached.

        Returns dict[name -> local_path].
        Large files (rankings DB ~1.3 GB): skips download if already present.
        """
        import requests

        self.resources_dir.mkdir(parents=True, exist_ok=True)
        resource_map = {
            "tf_list": ("tf_list_url", "allTFs_hg38.txt"),
            "rankings_db": ("rankings_db_url", "hg38_10kbp.mc9nr.genes.rankings.feather"),
            "motif_db": ("motif_db_url", "motifs-v9-nr.hgnc-m0.001-o0.0.tbl"),
        }
        paths: dict[str, Path] = {}
        for name, (url_key, filename) in resource_map.items():
            dest = self.resources_dir / filename
            paths[name] = dest
            if dest.exists() and not force:
                continue
            url = self.defaults[url_key]
            resp = requests.get(url, stream=True, timeout=300)
            resp.raise_for_status()
            with open(dest, "wb") as fh:
                for chunk in resp.iter_content(chunk_size=8192):
                    fh.write(chunk)
        return paths

    def resource_paths(self) -> dict[str, Path | None]:
        """
        Return paths to cached resource files (None if not yet downloaded).
        """
        resource_map = {
            "tf_list": "allTFs_hg38.txt",
            "rankings_db": "hg38_10kbp.mc9nr.genes.rankings.feather",
            "motif_db": "motifs-v9-nr.hgnc-m0.001-o0.0.tbl",
        }
        return {
            name: (self.resources_dir / fname if (self.resources_dir / fname).exists() else None)
            for name, fname in resource_map.items()
        }

    # ── core run methods ──────────────────────────────────────────────────────

    def run(
        self,
        loom_path: str,
        output_path: str,
        tf_list_path: str | None = None,
        rankings_db_path: str | None = None,
        motif_db_path: str | None = None,
    ) -> dict[str, list[float]]:
        """
        Run SCENIC on a loom file via subprocess.

        Returns z-scored AUC matrix: dict[regulon_name -> list[z_score_per_cell]].
        Raises RuntimeError if scenic-env is not found or subprocess fails.
        Raises ScenicEnvNotFoundError if the conda env is missing.
        """
        worker = Path(__file__).parent / "scenic_worker.py"
        config: dict = {
            **self.defaults,
            "loom_path": loom_path,
            "output_path": output_path,
        }
        if tf_list_path:
            config["tf_list_path"] = tf_list_path
        if rankings_db_path:
            config["rankings_db_path"] = rankings_db_path
        if motif_db_path:
            config["motif_db_path"] = motif_db_path

        # Fill in locally-cached resource paths if not explicitly provided
        cached = self.resource_paths()
        if not config.get("tf_list_path") and cached.get("tf_list"):
            config["tf_list_path"] = str(cached["tf_list"])
        if not config.get("rankings_db_path") and cached.get("rankings_db"):
            config["rankings_db_path"] = str(cached["rankings_db"])
        if not config.get("motif_db_path") and cached.get("motif_db"):
            config["motif_db_path"] = str(cached["motif_db"])

        config_json = json.dumps(config)
        conda = self._find_conda()
        cmd = [
            conda, "run", "-n", self.conda_env, "--no-capture-output",
            "python", str(worker), config_json,
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            if "No such environment" in result.stderr or "CondaEnvException" in result.stderr:
                raise ScenicEnvNotFoundError(
                    f"conda env '{self.conda_env}' not found.\n{result.stderr}"
                )
            raise RuntimeError(
                f"pySCENIC subprocess failed (conda env: {self.conda_env}):\n"
                f"{result.stderr}"
            )
        return json.loads(result.stdout)

    def run_from_adata(
        self,
        adata: "anndata.AnnData",
        output_dir: str | None = None,
        n_hvg: int = 2000,
        tf_list_path: str | None = None,
        rankings_db_path: str | None = None,
        motif_db_path: str | None = None,
    ) -> dict[str, list[float]]:
        """
        Convert AnnData to loom and run SCENIC in one call.

        Forces SCENIC_DEFAULTS["forced_genes"] into the gene set regardless
        of HVG rank (Phase 1 lesson: master TFs have low variance in mixed data).

        Parameters
        ----------
        adata : AnnData
            Normalized (10k counts, log1p) AnnData object.
        output_dir : str | None
            Directory for loom and SCENIC output files. Defaults to tempdir.
        n_hvg : int
            Number of highly variable genes to include (in addition to forced_genes).

        Returns
        -------
        dict[regulon_name -> list[z_score_per_cell]]
        """
        try:
            import scanpy as sc
            import loompy
            import numpy as np
        except ImportError as e:
            raise ImportError(
                f"scanpy and loompy required for run_from_adata. "
                f"Install with: pip install 'axiom-sc[science]' loompy\n{e}"
            )

        work_dir = Path(output_dir) if output_dir else Path(tempfile.mkdtemp(prefix="axiom_scenic_"))
        work_dir.mkdir(parents=True, exist_ok=True)
        loom_path = str(work_dir / "input.loom")
        output_path = str(work_dir / "auc_mtx.csv")

        # Select genes: top HVGs + forced genes
        sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, flavor="seurat")
        hvg_set = set(adata.var_names[adata.var["highly_variable"]])
        forced = set(self.defaults.get("forced_genes", []))
        # Only keep forced genes that are present in the dataset
        present_forced = forced & set(adata.var_names)
        gene_set = list(hvg_set | present_forced)

        adata_sub = adata[:, gene_set].copy()

        # Write to loom (pySCENIC expects cells × genes in loom format)
        X = adata_sub.X
        if hasattr(X, "toarray"):
            X = X.toarray()
        matrix = X.T.astype(np.float32)  # genes × cells for loompy convention

        row_attrs = {"Gene": np.array(adata_sub.var_names)}
        col_attrs = {"CellID": np.array(adata_sub.obs_names)}
        loompy.create(loom_path, matrix, row_attrs, col_attrs)

        return self.run(
            loom_path=loom_path,
            output_path=output_path,
            tf_list_path=tf_list_path,
            rankings_db_path=rankings_db_path,
            motif_db_path=motif_db_path,
        )


class ScenicEnvNotFoundError(RuntimeError):
    """Raised when the scenic-env conda environment is missing or incomplete."""
    pass
