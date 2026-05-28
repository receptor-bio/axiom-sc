"""
pySCENIC pipeline orchestrator for AXIOM-SC Tier 2.

Calls pySCENIC via subprocess isolation (GPL v3 — not imported directly).
Handles resource caching, stratified sampling, and per-cluster
EvidenceBundle population with z-scored regulon AUC values.

References:
  Van de Sande et al. (2020) Nature Protocols 15:2247
  doi: 10.1038/s41596-020-0336-2

License: Apache 2.0 — pySCENIC is GPL v3, called via subprocess only.
"""
from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from axiom_sc.tier2.scenic_runner import ScenicRunner, SCENIC_DEFAULTS
from axiom_sc.utils.zscore_auc import zscore_auc_dict, cluster_mean_zscore
from axiom_sc.utils.stratified_sample import stratified_sample, check_min_cells

if TYPE_CHECKING:
    import anndata
    from axiom_sc.tier2.evidence import EvidenceBundle

logger = logging.getLogger(__name__)


class ScenicPipeline:
    """
    Full pySCENIC pipeline for Tier 2 evidence generation.

    Typical usage:
        pipeline = ScenicPipeline()
        regulon_bundles = pipeline.run(adata, cluster_key="leiden")
        # regulon_bundles: dict[cluster_id -> dict[regulon -> mean_z_score]]
    """

    def __init__(
        self,
        scenic_conda_env: str | None = None,
        resources_dir: str | None = None,
        scenic_overrides: dict | None = None,
    ):
        self.runner = ScenicRunner(
            scenic_conda_env=scenic_conda_env,
            resources_dir=resources_dir,
            scenic_defaults=scenic_overrides,
        )
        self.defaults = {**SCENIC_DEFAULTS, **(scenic_overrides or {})}

    # ── public entry point ────────────────────────────────────────────────────

    def run(
        self,
        adata: "anndata.AnnData",
        cluster_key: str = "leiden",
        cell_type_key: str | None = None,
        output_dir: str | None = None,
        n_hvg: int = 2000,
        download_resources: bool = False,
        tf_list_path: str | None = None,
        rankings_db_path: str | None = None,
        motif_db_path: str | None = None,
    ) -> dict[str, dict[str, float]]:
        """
        Run the full SCENIC pipeline and return per-cluster mean z-scored regulon scores.

        Steps:
          1. Stratified subsample to SCENIC_DEFAULTS["recommended_total"] cells
          2. (Optionally) download resource files
          3. Call ScenicRunner.run_from_adata → z-scored regulon AUC dict
          4. Compute per-cluster mean z-scores for each regulon
          5. Return dict[cluster_id -> dict[regulon -> mean_z]]

        Parameters
        ----------
        adata : AnnData
            10k-normalized, log1p-transformed AnnData.
        cluster_key : str
            adata.obs key with cluster labels (used for per-cluster aggregation).
        cell_type_key : str | None
            adata.obs key with cell type labels for stratified sampling.
            If None, cluster_key is used for stratification.
        output_dir : str | None
            Directory for intermediate loom / AUC CSV files.
        n_hvg : int
            Number of highly variable genes (default 2000).
        download_resources : bool
            If True, download SCENIC resource files before running.
        tf_list_path, rankings_db_path, motif_db_path : str | None
            Paths to SCENIC resource files (overrides cached paths).

        Returns
        -------
        dict[cluster_id -> dict[regulon -> mean_z_score]]
            Per-cluster mean z-scores for each recovered regulon.
            These go directly into EvidenceBundle.regulons.
        """
        try:
            import numpy as np
        except ImportError:
            raise ImportError("numpy required. Install with: pip install 'axiom-sc[science]'")

        strat_key = cell_type_key or cluster_key
        recommended_total = self.defaults.get("recommended_total", 50_000)
        min_cells_per_type = self.defaults.get("min_cells_per_type", 300)

        # ── Step 1: stratified sampling ───────────────────────────────────────
        n_cells = adata.n_obs
        if n_cells > recommended_total:
            logger.info(
                f"Subsampling {n_cells} → {recommended_total} cells "
                f"(stratified by {strat_key!r})"
            )
            adata_sub = stratified_sample(
                adata,
                cell_type_key=strat_key,
                target_total=recommended_total,
                min_per_type=min_cells_per_type,
            )
        else:
            adata_sub = adata
            if n_cells < recommended_total:
                logger.warning(
                    f"Only {n_cells} cells available (recommended ≥ {recommended_total}). "
                    "Rare TF regulons (FOXP3, AIRE) may not be recovered. "
                    "Phase 1 found 20k cells insufficient for these TFs."
                )

        # Warn about cell types below minimum
        insufficient = check_min_cells(adata_sub, strat_key, min_cells_per_type)
        if insufficient:
            logger.warning(
                f"Cell types below {min_cells_per_type}-cell minimum after sampling: "
                f"{insufficient}. GRN signal may be weak for these types."
            )

        # ── Step 2: download resources if requested ───────────────────────────
        if download_resources:
            logger.info("Downloading SCENIC resource files...")
            self.runner.download_resources()

        # ── Step 2b: auto-detect cached Patch B resources ────────────────────
        from pathlib import Path as _Path
        _scenic_dir = _Path.home() / ".axiom_sc" / "resources" / "scenic"
        if rankings_db_path is None and (_scenic_dir / "hg38_rankings.feather").exists():
            rankings_db_path = str(_scenic_dir / "hg38_rankings.feather")
            logger.info("Using cached pySCENIC rankings DB: %s", rankings_db_path)
        if motif_db_path is None and (_scenic_dir / "motifs-v9.tbl").exists():
            motif_db_path = str(_scenic_dir / "motifs-v9.tbl")
            logger.info("Using cached pySCENIC motif DB: %s", motif_db_path)
        if tf_list_path is None and (_scenic_dir / "allTFs_hg38.txt").exists():
            tf_list_path = str(_scenic_dir / "allTFs_hg38.txt")
            logger.info("Using cached SCENIC TF list: %s", tf_list_path)

        # ── Step 3: run pySCENIC via subprocess ───────────────────────────────
        mode = "cisTarget" if rankings_db_path else "adjacency fallback"
        logger.info(
            f"Running pySCENIC on {adata_sub.n_obs} cells × {adata_sub.n_vars} genes "
            f"(conda env: {self.runner.conda_env!r}, mode: {mode})"
        )
        auc_zscores = self.runner.run_from_adata(
            adata_sub,
            output_dir=output_dir,
            n_hvg=n_hvg,
            tf_list_path=tf_list_path,
            rankings_db_path=rankings_db_path,
            motif_db_path=motif_db_path,
        )

        if not auc_zscores:
            logger.warning("pySCENIC returned zero regulons. Check SCENIC output and NES threshold.")
            return {}

        n_regulons = len(auc_zscores)
        logger.info(f"pySCENIC recovered {n_regulons} regulons.")

        # Check forced genes recovery
        forced = set(self.defaults.get("forced_genes", []))
        recovered = set(auc_zscores.keys())
        missing_forced = forced & set(adata_sub.var_names) - recovered
        if missing_forced:
            logger.warning(
                f"Forced master TFs not recovered as regulons: {sorted(missing_forced)}. "
                "Consider lowering nes_threshold (current: "
                f"{self.defaults.get('nes_threshold', 2.0)})."
            )

        # ── Step 4: per-cluster mean z-scores ─────────────────────────────────
        # auc_zscores already z-scored by scenic_worker.py (z_score_auc=True default)
        # Map cell indices in adata_sub back to cluster labels
        cluster_labels = adata_sub.obs[cluster_key]
        cluster_ids = cluster_labels.unique()

        result: dict[str, dict[str, float]] = {}
        for cluster_id in cluster_ids:
            cell_mask = cluster_labels == cluster_id
            cell_indices = list(np.where(cell_mask)[0])
            cluster_regulons: dict[str, float] = {}
            for regulon, zscores in auc_zscores.items():
                vals = [zscores[i] for i in cell_indices]
                cluster_regulons[regulon] = float(sum(vals) / len(vals)) if vals else 0.0
            result[str(cluster_id)] = cluster_regulons

        logger.info(
            f"Computed per-cluster regulon z-scores for "
            f"{len(result)} clusters, {n_regulons} regulons."
        )
        return result

    def _build_config(
        self,
        loom_path: str,
        output_dir: str | None = None,
        tf_list_path: str | None = None,
        rankings_db_path: str | None = None,
        motif_db_path: str | None = None,
    ) -> dict:
        """
        Build the resolved SCENIC config dict (auto-detecting cached Patch B resources).

        Used by run() internally and exposed for testing.
        """
        from pathlib import Path as _Path
        scenic_dir = _Path.home() / ".axiom_sc" / "resources" / "scenic"
        if rankings_db_path is None and (scenic_dir / "hg38_rankings.feather").exists():
            rankings_db_path = str(scenic_dir / "hg38_rankings.feather")
        if motif_db_path is None and (scenic_dir / "motifs-v9.tbl").exists():
            motif_db_path = str(scenic_dir / "motifs-v9.tbl")
        if tf_list_path is None and (scenic_dir / "allTFs_hg38.txt").exists():
            tf_list_path = str(scenic_dir / "allTFs_hg38.txt")
        return {
            **self.defaults,
            "loom_path": loom_path,
            "output_dir": output_dir,
            "tf_list_path": tf_list_path,
            "rankings_db_path": rankings_db_path,
            "motif_db_path": motif_db_path,
        }

    def populate_evidence_bundles(
        self,
        bundles: dict[str, "EvidenceBundle"],
        cluster_regulons: dict[str, dict[str, float]],
    ) -> None:
        """
        Merge SCENIC regulon z-scores into existing EvidenceBundle objects.

        Parameters
        ----------
        bundles : dict[cluster_id -> EvidenceBundle]
        cluster_regulons : dict[cluster_id -> dict[regulon -> mean_z]]
            Output of self.run().
        """
        for cluster_id, bundle in bundles.items():
            if cluster_id in cluster_regulons:
                bundle.regulons.update(cluster_regulons[cluster_id])
            else:
                logger.debug(f"No SCENIC regulons for cluster {cluster_id!r}")

    def recovery_report(
        self,
        cluster_regulons: dict[str, dict[str, float]],
        z_active_thresh: float | None = None,
    ) -> dict:
        """
        Generate a regulon recovery report for validation.

        Returns
        -------
        dict with keys:
          n_regulons: total regulons recovered
          forced_recovered: forced master TFs found as regulons (with max z across clusters)
          forced_missing: forced TFs not recovered
          active_per_cluster: per-cluster count of active regulons (z ≥ z_active_thresh)
        """
        thresh = z_active_thresh or self.defaults.get("z_active_thresh", 1.5)
        all_regulons = set(
            r for cluster in cluster_regulons.values() for r in cluster
        )
        forced = set(self.defaults.get("forced_genes", []))
        forced_recovered = {
            tf: max(cluster_regulons[c].get(tf, 0.0) for c in cluster_regulons)
            for tf in forced & all_regulons
        }
        forced_missing = forced - all_regulons

        active_per_cluster = {
            cluster_id: sum(1 for z in regulons.values() if z >= thresh)
            for cluster_id, regulons in cluster_regulons.items()
        }

        return {
            "n_regulons": len(all_regulons),
            "forced_recovered": forced_recovered,
            "forced_missing": sorted(forced_missing),
            "active_per_cluster": active_per_cluster,
        }
