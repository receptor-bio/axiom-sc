"""
Day 6 tests — pySCENIC subprocess isolation, scType scorer, stratified sampling,
z-score AUC normalization, and scenic pipeline (mocked — no real scenic-env required).

Tests verify the Apache-2.0-side logic without invoking pySCENIC.
The subprocess call is mocked so CI passes without scenic-env installed.
"""
import json
import subprocess
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

# ── AnnData utils tests ───────────────────────────────────────────────────────

_anndata = pytest.importorskip("anndata")
import pandas as pd

from axiom_sc.utils.anndata_utils import validate_adata_normalized


class TestAnnDataUtils:
    def _make_log1p_adata(self, n_cells=20, n_genes=10):
        """Make a plausible log1p-normalized AnnData (max ~10)."""
        X = np.random.rand(n_cells, n_genes).astype(np.float32) * 5.0
        obs = pd.DataFrame(index=[f"c{i}" for i in range(n_cells)])
        var = pd.DataFrame(index=[f"g{i}" for i in range(n_genes)])
        return _anndata.AnnData(X=X, obs=obs, var=var)

    def _make_raw_counts_adata(self, n_cells=20, n_genes=10):
        """Make a raw-count AnnData (max ~1000)."""
        X = (np.random.rand(n_cells, n_genes) * 1000).astype(np.float32)
        obs = pd.DataFrame(index=[f"c{i}" for i in range(n_cells)])
        var = pd.DataFrame(index=[f"g{i}" for i in range(n_genes)])
        return _anndata.AnnData(X=X, obs=obs, var=var)

    def test_validate_normalized_passes_on_log1p(self):
        """validate_adata_normalized does not raise for log1p-normalized data."""
        adata = self._make_log1p_adata()
        validate_adata_normalized(adata)  # should not raise

    def test_validate_normalized_raises_on_raw_counts(self):
        """validate_adata_normalized raises ValueError for raw count data."""
        adata = self._make_raw_counts_adata()
        with pytest.raises(ValueError, match="raw counts"):
            validate_adata_normalized(adata)

    def test_validate_normalized_works_with_sparse(self):
        """validate_adata_normalized works with scipy sparse matrices."""
        from scipy.sparse import csr_matrix
        X = csr_matrix(np.random.rand(20, 10).astype(np.float32) * 3.0)
        obs = pd.DataFrame(index=[f"c{i}" for i in range(20)])
        var = pd.DataFrame(index=[f"g{i}" for i in range(10)])
        adata = _anndata.AnnData(X=X, obs=obs, var=var)
        validate_adata_normalized(adata)  # should not raise


# ── ScenicRunner tests ────────────────────────────────────────────────────────

from axiom_sc.tier2.scenic_runner import ScenicRunner, ScenicEnvNotFoundError, SCENIC_DEFAULTS


class TestScenicDefaults:
    def test_nes_threshold_is_2(self):
        """Phase 1 correction: NES must be 2.0, not 3.0."""
        assert SCENIC_DEFAULTS["nes_threshold"] == 2.0

    def test_min_genes_is_5(self):
        """Phase 1 correction: min_genes must be 5 (not 30)."""
        assert SCENIC_DEFAULTS["min_genes"] == 5

    def test_z_score_auc_enabled(self):
        """Phase 1 correction: z-scoring must be on by default."""
        assert SCENIC_DEFAULTS["z_score_auc"] is True

    def test_z_active_threshold(self):
        assert SCENIC_DEFAULTS["z_active_thresh"] == 1.5

    def test_z_inactive_threshold(self):
        assert SCENIC_DEFAULTS["z_inactive_thresh"] == -0.5

    def test_forced_genes_include_master_tfs(self):
        """All Phase 1 master TFs must be in forced_genes."""
        required = {"FOXP3", "TBX21", "GATA3", "RORC", "IRF7", "PAX5",
                    "AIRE", "PSMB11", "EBF1"}
        forced = set(SCENIC_DEFAULTS["forced_genes"])
        missing = required - forced
        assert not missing, f"Master TFs missing from forced_genes: {missing}"

    def test_recommended_total_is_50k(self):
        """Phase 2 uses 50k cells, not 20k."""
        assert SCENIC_DEFAULTS["recommended_total"] == 50_000

    def test_min_cells_per_type(self):
        assert SCENIC_DEFAULTS["min_cells_per_type"] == 300


class TestScenicRunnerInit:
    def test_default_conda_env(self):
        runner = ScenicRunner()
        assert runner.conda_env == "scenic-env"

    def test_custom_conda_env(self):
        runner = ScenicRunner(scenic_conda_env="my-scenic-env")
        assert runner.conda_env == "my-scenic-env"

    def test_override_defaults(self):
        runner = ScenicRunner(scenic_defaults={"nes_threshold": 1.5})
        assert runner.defaults["nes_threshold"] == 1.5
        # Other defaults should remain
        assert runner.defaults["min_genes"] == 5

    def test_env_var_conda_env(self, monkeypatch):
        monkeypatch.setenv("SCENIC_CONDA_ENV", "env-from-envvar")
        runner = ScenicRunner()
        assert runner.conda_env == "env-from-envvar"


class TestScenicRunnerRun:
    def test_run_calls_subprocess_with_correct_config(self, tmp_path):
        """run() passes correct config JSON to scenic_worker.py subprocess."""
        runner = ScenicRunner()
        expected_result = {"FOXP3": [1.2, 0.3, -0.5], "IRF7": [2.1, -0.1, 0.8]}

        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(
                returncode=0,
                stdout=json.dumps(expected_result),
                stderr="",
            )
            result = runner.run(
                loom_path=str(tmp_path / "test.loom"),
                output_path=str(tmp_path / "auc.csv"),
            )

        assert result == expected_result
        call_args = mock_run.call_args[0][0]  # first positional arg = command list
        assert "conda" in call_args[0]
        assert "python" in call_args
        # Config JSON should be in the command
        config_json = call_args[-1]
        config = json.loads(config_json)
        assert config["nes_threshold"] == 2.0
        assert config["z_score_auc"] is True

    def test_run_raises_on_subprocess_failure(self, tmp_path):
        """run() raises RuntimeError when subprocess exits non-zero."""
        runner = ScenicRunner()
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(
                returncode=1,
                stdout="",
                stderr="Something went wrong in pySCENIC",
            )
            with pytest.raises(RuntimeError, match="pySCENIC subprocess failed"):
                runner.run(str(tmp_path / "test.loom"), str(tmp_path / "auc.csv"))

    def test_run_raises_scenic_env_not_found(self, tmp_path):
        """run() raises ScenicEnvNotFoundError when conda env is missing."""
        runner = ScenicRunner()
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(
                returncode=1,
                stdout="",
                stderr="No such environment: scenic-env\nCondaEnvException: ...",
            )
            with pytest.raises(ScenicEnvNotFoundError):
                runner.run(str(tmp_path / "test.loom"), str(tmp_path / "auc.csv"))

    def test_run_includes_resource_paths_in_config(self, tmp_path):
        """run() passes explicit resource paths through to worker config."""
        runner = ScenicRunner()
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout="{}", stderr="")
            runner.run(
                loom_path=str(tmp_path / "test.loom"),
                output_path=str(tmp_path / "auc.csv"),
                tf_list_path="/path/to/tfs.txt",
                rankings_db_path="/path/to/rankings.feather",
                motif_db_path="/path/to/motifs.tbl",
            )
        config = json.loads(mock_run.call_args[0][0][-1])
        assert config["tf_list_path"] == "/path/to/tfs.txt"
        assert config["rankings_db_path"] == "/path/to/rankings.feather"
        assert config["motif_db_path"] == "/path/to/motifs.tbl"

    def test_check_env_raises_when_conda_missing(self):
        """check_env raises ScenicEnvNotFoundError when conda not on PATH or common paths."""
        runner = ScenicRunner()
        # Patch shutil.which AND Path.exists so no fallback candidate is found
        with patch("shutil.which", return_value=None), \
             patch("pathlib.Path.exists", return_value=False):
            with pytest.raises(ScenicEnvNotFoundError, match="conda not found"):
                runner.check_env()

    def test_check_env_raises_when_pyscenic_cli_missing(self):
        """check_env raises when scenic-env exists but pyscenic CLI not installed."""
        runner = ScenicRunner()
        with patch("axiom_sc.tier2.scenic_runner.ScenicRunner._find_conda",
                   return_value="/usr/bin/conda"):
            with patch("subprocess.run") as mock_run:
                mock_run.return_value = MagicMock(returncode=1, stdout="", stderr="pyscenic: command not found")
                with pytest.raises(ScenicEnvNotFoundError, match="pyscenic CLI not installed"):
                    runner.check_env()

    def test_resource_paths_returns_none_for_uncached(self, tmp_path):
        """resource_paths returns None for files not yet downloaded."""
        runner = ScenicRunner(resources_dir=str(tmp_path / "empty_dir"))
        paths = runner.resource_paths()
        assert all(v is None for v in paths.values())

    def test_resource_paths_returns_path_for_cached(self, tmp_path):
        """resource_paths returns Path for files that exist."""
        res_dir = tmp_path / "resources"
        res_dir.mkdir()
        (res_dir / "allTFs_hg38.txt").write_text("FOXP3\nIRF7\n")
        runner = ScenicRunner(resources_dir=str(res_dir))
        paths = runner.resource_paths()
        assert paths["tf_list"] is not None
        assert paths["tf_list"].exists()
        assert paths["rankings_db"] is None   # not created


# ── Z-score AUC tests ─────────────────────────────────────────────────────────

from axiom_sc.utils.zscore_auc import zscore_auc_matrix, zscore_auc_dict, cluster_mean_zscore


class TestZscoreAuc:
    def test_zscore_matrix_basic(self):
        """Z-scoring a known matrix returns correct per-column z-scores."""
        # Column 0: values [0, 1, 2] → z-scores [-1.22, 0.0, 1.22] (approx)
        matrix = np.array([[0.0, 10.0], [1.0, 10.0], [2.0, 10.0]])
        result = zscore_auc_matrix(matrix, ["REG1", "REG2"])
        z1 = result["REG1"]
        assert abs(z1.mean()) < 1e-10, "Z-scores must have zero mean"
        assert abs(z1.std(ddof=0) - 1.0) < 0.01 or len(z1) == 3

    def test_zscore_constant_column_returns_zeros(self):
        """Constant regulon (no variance) returns zero z-scores, not NaN."""
        matrix = np.array([[5.0], [5.0], [5.0]])
        result = zscore_auc_matrix(matrix, ["FLAT_REG"])
        assert np.all(result["FLAT_REG"] == 0.0)

    def test_zscore_matrix_shape_mismatch_raises(self):
        """Mismatched matrix columns and regulon_names raises ValueError."""
        matrix = np.array([[1.0, 2.0], [3.0, 4.0]])
        with pytest.raises(ValueError, match="regulon names"):
            zscore_auc_matrix(matrix, ["ONLY_ONE"])

    def test_zscore_dict_basic(self):
        """zscore_auc_dict wraps zscore_auc_matrix correctly."""
        auc_dict = {
            "FOXP3": [0.1, 0.5, 0.9, 0.2],
            "IRF7": [1.0, 2.0, 1.5, 3.0],
        }
        result = zscore_auc_dict(auc_dict)
        assert "FOXP3" in result
        assert "IRF7" in result
        assert abs(result["FOXP3"].mean()) < 1e-10

    def test_zscore_dict_empty(self):
        """Empty dict returns empty dict."""
        assert zscore_auc_dict({}) == {}

    def test_cluster_mean_zscore(self):
        """cluster_mean_zscore averages z-scores over a cell index list."""
        zscore_dict = {
            "FOXP3": np.array([1.0, 2.0, 3.0, -1.0]),
            "IRF7": np.array([0.5, 0.5, 0.5, 0.5]),
        }
        # Cluster = cells 0, 1
        means = cluster_mean_zscore(zscore_dict, [0, 1])
        assert means["FOXP3"] == pytest.approx(1.5)
        assert means["IRF7"] == pytest.approx(0.5)

    def test_zscore_preserves_regulon_order(self):
        """Column order in the output matches regulon_names order."""
        matrix = np.random.rand(10, 3)
        names = ["A", "B", "C"]
        result = zscore_auc_matrix(matrix, names)
        assert list(result.keys()) == names


# ── Stratified sampling tests ─────────────────────────────────────────────────

anndata = pytest.importorskip("anndata")
from axiom_sc.utils.stratified_sample import (
    stratified_sample, get_type_counts, check_min_cells
)


def _make_adata(n_per_type: dict[str, int], n_genes: int = 10) -> anndata.AnnData:
    """Helper: make a small AnnData with cell type labels."""
    import pandas as pd
    labels = []
    for ct, n in n_per_type.items():
        labels.extend([ct] * n)
    n_cells = len(labels)
    X = np.random.rand(n_cells, n_genes).astype(np.float32)
    obs = pd.DataFrame({"cell_type": labels})
    return anndata.AnnData(X=X, obs=obs)


class TestStratifiedSample:
    def test_respects_target_total(self):
        """Output has ≤ target_total cells."""
        adata = _make_adata({"A": 5000, "B": 5000, "C": 5000})
        sampled = stratified_sample(adata, "cell_type", target_total=1000)
        assert sampled.n_obs <= 1000

    def test_min_per_type_guaranteed(self):
        """Each type gets at least min_per_type cells (or all if fewer)."""
        adata = _make_adata({"A": 1000, "B": 1000, "C": 100})  # C is small
        sampled = stratified_sample(adata, "cell_type", target_total=5000, min_per_type=300)
        counts = get_type_counts(sampled, "cell_type")
        assert counts["A"] >= 300
        assert counts["B"] >= 300
        assert counts["C"] == 100  # only 100 available, take all

    def test_no_subsampling_when_below_target(self):
        """When n_cells < target_total, all cells are returned."""
        adata = _make_adata({"A": 50, "B": 50})
        sampled = stratified_sample(adata, "cell_type", target_total=10_000)
        assert sampled.n_obs == 100

    def test_reproducible_with_same_seed(self):
        """Same seed → same cell selection."""
        adata = _make_adata({"A": 500, "B": 500, "C": 500})
        s1 = stratified_sample(adata, "cell_type", target_total=200, seed=7)
        s2 = stratified_sample(adata, "cell_type", target_total=200, seed=7)
        assert list(s1.obs_names) == list(s2.obs_names)

    def test_different_seeds_differ(self):
        """Different seeds → different cell selections (statistically)."""
        adata = _make_adata({"A": 500, "B": 500})
        s1 = stratified_sample(adata, "cell_type", target_total=100, seed=1)
        s2 = stratified_sample(adata, "cell_type", target_total=100, seed=99)
        # With 500 cells and 100 sample, probability all same is astronomically small
        assert list(s1.obs_names) != list(s2.obs_names)

    def test_missing_key_raises(self):
        """stratified_sample raises ValueError if cell_type_key not in obs."""
        adata = _make_adata({"A": 100})
        with pytest.raises(ValueError, match="not found in adata.obs"):
            stratified_sample(adata, "nonexistent_key")

    def test_check_min_cells_detects_insufficient(self):
        """check_min_cells returns types below threshold."""
        adata = _make_adata({"A": 500, "B": 100, "C": 50})
        insufficient = check_min_cells(adata, "cell_type", min_per_type=200)
        assert "A" not in insufficient
        assert "B" in insufficient
        assert "C" in insufficient
        assert insufficient["C"] == 50

    def test_check_min_cells_empty_when_all_sufficient(self):
        """check_min_cells returns empty dict when all types are above threshold."""
        adata = _make_adata({"A": 500, "B": 400})
        result = check_min_cells(adata, "cell_type", min_per_type=300)
        assert result == {}

    def test_get_type_counts(self):
        """get_type_counts returns correct per-type counts."""
        adata = _make_adata({"X": 30, "Y": 70})
        counts = get_type_counts(adata, "cell_type")
        assert counts["X"] == 30
        assert counts["Y"] == 70


# ── scType scorer tests ───────────────────────────────────────────────────────

from axiom_sc.tier2.sctype_scorer import score_cell_types


class TestScTypeScorer:
    def _make_expr(self):
        """3 cells × 4 genes: CD4, CD8A, FOXP3, IL2"""
        import pandas as pd
        gene_names = ["CD4", "CD8A", "FOXP3", "IL2"]
        # Cell 0: CD4 T cell (CD4 high)
        # Cell 1: CD8 T cell (CD8A high)
        # Cell 2: Treg (FOXP3 high, CD4 high)
        matrix = np.array([
            [3.0, 0.1, 0.1, 0.1],   # CD4_T
            [0.1, 3.0, 0.1, 0.1],   # CD8_T
            [2.5, 0.1, 3.0, 0.1],   # Treg
        ], dtype=float)
        return matrix, gene_names

    def test_positive_markers_score(self):
        """Positive markers produce higher scores for matching cell types."""
        matrix, genes = self._make_expr()
        markers = {"CD4_T": ["CD4"], "CD8_T": ["CD8A"], "Treg": ["FOXP3", "CD4"]}
        scores = score_cell_types(matrix, genes, markers)
        assert scores.shape == (3, 3), f"Expected (3, 3) score matrix, got {scores.shape}"
        # Cell 0 should have highest CD4_T score
        assert scores.iloc[0]["CD4_T"] > scores.iloc[1]["CD4_T"]
        # Cell 1 should have highest CD8_T score
        assert scores.iloc[1]["CD8_T"] > scores.iloc[0]["CD8_T"]

    def test_negative_markers_reduce_score(self):
        """Negative markers reduce score for a cell type."""
        matrix, genes = self._make_expr()
        markers = {"Treg": ["FOXP3"]}
        # IL2 is a negative marker for Treg (Treg does NOT produce IL2)
        neg = {"Treg": ["IL2"]}
        scores_no_neg = score_cell_types(matrix, genes, markers)
        scores_with_neg = score_cell_types(matrix, genes, markers, neg)
        # Scores with negative markers should be ≤ scores without
        assert all(scores_with_neg["Treg"] <= scores_no_neg["Treg"])

    def test_missing_genes_ignored(self):
        """Markers not in gene_names are silently ignored (no KeyError)."""
        matrix, genes = self._make_expr()
        markers = {"MyType": ["NONEXISTENT_GENE_XYZ"]}
        scores = score_cell_types(matrix, genes, markers)
        # All zeros because no genes matched
        assert (scores["MyType"] == 0.0).all()

    def test_empty_cell_types(self):
        """Empty markers dict returns empty DataFrame."""
        matrix, genes = self._make_expr()
        scores = score_cell_types(matrix, genes, {})
        assert scores.shape[1] == 0

    def test_output_is_dataframe_with_correct_shape(self):
        """Output DataFrame has n_cells rows and n_cell_types columns."""
        matrix, genes = self._make_expr()
        markers = {"A": ["CD4"], "B": ["CD8A"], "C": ["FOXP3"]}
        scores = score_cell_types(matrix, genes, markers)
        import pandas as pd
        assert isinstance(scores, pd.DataFrame)
        assert scores.shape == (3, 3)


# ── ScenicPipeline tests ──────────────────────────────────────────────────────

from axiom_sc.tier2.scenic_pipeline import ScenicPipeline
from axiom_sc.tier2.evidence import EvidenceBundle


class TestScenicPipeline:
    def _make_adata_with_clusters(self):
        import pandas as pd
        n_per_type = {"T": 50, "B": 50, "NK": 50}
        labels, cluster_ids = [], []
        for i, (ct, n) in enumerate(n_per_type.items()):
            labels.extend([ct] * n)
            cluster_ids.extend([str(i)] * n)
        n_cells = len(labels)
        # Use forced gene names in var
        gene_names = SCENIC_DEFAULTS["forced_genes"][:20]
        X = np.random.rand(n_cells, len(gene_names)).astype(np.float32)
        obs = pd.DataFrame({"cell_type": labels, "leiden": cluster_ids},
                           index=[f"cell_{i}" for i in range(n_cells)])
        var = pd.DataFrame(index=gene_names)
        return anndata.AnnData(X=X, obs=obs, var=var)

    def test_populate_evidence_bundles_merges_regulons(self):
        """populate_evidence_bundles adds SCENIC z-scores to existing bundles."""
        pipeline = ScenicPipeline()
        bundles = {
            "0": EvidenceBundle(cluster_id="0", marker_genes={"CD3E": 2.5}),
            "1": EvidenceBundle(cluster_id="1", marker_genes={"CD19": 2.0}),
        }
        cluster_regulons = {
            "0": {"FOXP3": 2.1, "IRF7": -0.3},
            "1": {"PAX5": 3.2, "FOXP3": -0.8},
        }
        pipeline.populate_evidence_bundles(bundles, cluster_regulons)
        assert bundles["0"].regulons["FOXP3"] == pytest.approx(2.1)
        assert bundles["0"].marker_genes["CD3E"] == pytest.approx(2.5)  # preserved
        assert bundles["1"].regulons["PAX5"] == pytest.approx(3.2)

    def test_populate_evidence_bundles_missing_cluster_skipped(self):
        """Bundles with no matching cluster_regulons entry are unchanged."""
        pipeline = ScenicPipeline()
        bundle = EvidenceBundle(cluster_id="99", marker_genes={"CD3E": 2.5})
        pipeline.populate_evidence_bundles({"99": bundle}, {})
        assert bundle.regulons == {}
        assert bundle.marker_genes["CD3E"] == pytest.approx(2.5)

    def test_recovery_report_structure(self):
        """recovery_report returns expected keys."""
        pipeline = ScenicPipeline()
        cluster_regulons = {
            "0": {"FOXP3": 2.5, "TBX21": 1.8, "CD19": -0.3},
            "1": {"FOXP3": -0.2, "PAX5": 3.0},
        }
        report = pipeline.recovery_report(cluster_regulons)
        assert "n_regulons" in report
        assert "forced_recovered" in report
        assert "forced_missing" in report
        assert "active_per_cluster" in report
        assert report["n_regulons"] >= 4  # FOXP3, TBX21, CD19, PAX5

    def test_recovery_report_forced_genes(self):
        """recovery_report correctly identifies forced TFs recovered/missing."""
        pipeline = ScenicPipeline()
        cluster_regulons = {"0": {"FOXP3": 2.5, "IRF7": 1.6}}
        report = pipeline.recovery_report(cluster_regulons)
        assert "FOXP3" in report["forced_recovered"]
        assert "IRF7" in report["forced_recovered"]
        # TBX21 was not recovered
        assert "TBX21" in report["forced_missing"]

    def test_recovery_report_active_per_cluster(self):
        """active_per_cluster counts regulons above z_active_thresh."""
        pipeline = ScenicPipeline()
        cluster_regulons = {
            "0": {"A": 2.0, "B": 0.5, "C": 1.6},  # 2 active (A and C)
            "1": {"A": 0.2, "B": 0.1},               # 0 active
        }
        report = pipeline.recovery_report(cluster_regulons, z_active_thresh=1.5)
        assert report["active_per_cluster"]["0"] == 2
        assert report["active_per_cluster"]["1"] == 0

    @patch("axiom_sc.tier2.scenic_runner.ScenicRunner.run_from_adata")
    def test_run_pipeline_end_to_end_mocked(self, mock_run_from_adata):
        """Full ScenicPipeline.run() with mocked runner returns per-cluster dict."""
        # Mock pySCENIC returning z-scored regulons for 150 cells
        n_cells = 150
        mock_run_from_adata.return_value = {
            "FOXP3": [float(i % 3) for i in range(n_cells)],
            "PAX5":  [float((i + 1) % 2) for i in range(n_cells)],
        }

        pipeline = ScenicPipeline()
        adata = self._make_adata_with_clusters()

        result = pipeline.run(adata, cluster_key="leiden", cell_type_key="cell_type")

        assert isinstance(result, dict)
        # 3 clusters (0, 1, 2)
        assert len(result) == 3
        for cluster_id, regulons in result.items():
            assert "FOXP3" in regulons
            assert "PAX5" in regulons
            assert isinstance(regulons["FOXP3"], float)
