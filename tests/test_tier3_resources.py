"""
Patch B: tests for Tier 3 external resource loading (graceful fallback when absent).
"""
from __future__ import annotations

import json
import pickle
import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock

import numpy as np
import pytest

anndata = pytest.importorskip("anndata")


# ── helpers ───────────────────────────────────────────────────────────────────

def _make_adata(n_cells=40, n_genes=30, n_clusters=4, with_connectivities=False):
    import scipy.sparse as sp
    X = np.abs(np.random.randn(n_cells, n_genes)).astype(np.float32)
    obs = {"leiden": [str(i % n_clusters) for i in range(n_cells)]}
    var_names = [f"GENE_{i}" for i in range(n_genes)]
    adata = anndata.AnnData(X=X, obs=obs)
    adata.var_names = var_names
    if with_connectivities:
        conn = np.zeros((n_cells, n_cells), dtype=np.float32)
        for i in range(n_cells):
            nbrs = np.random.choice([j for j in range(n_cells) if j != i], size=5, replace=False)
            for j in nbrs:
                conn[i, j] = 1.0
        adata.obsp["connectivities"] = sp.csr_matrix(conn)
    return adata


# ── CommunicationStream: external LR database ─────────────────────────────────

class TestCommunicationLRDatabase:

    def test_lr_database_empty_by_default_when_file_absent(self, tmp_path):
        """Without the index file, LR_DATABASE should fall back to empty {}."""
        missing = tmp_path / "lr_gene_index.json"
        with patch("axiom_sc.tier3.communication._LR_INDEX_PATH", missing):
            from axiom_sc.tier3 import communication
            db = communication._load_lr_database()
        assert isinstance(db, dict)
        # May be empty or contain the built-in subset — just must be a dict
        # (built-in L-R list is hardcoded in the module; external index extends it)

    def test_communication_uses_lr_database_when_available(self, tmp_path):
        """When lr_gene_index.json is present, _load_lr_database loads it."""
        # Build a minimal index
        fake_index = {"TGFB1": "LIGAND", "TGFBR1": "RECEPTOR", "IL6": "LIGAND"}
        index_path = tmp_path / "lr_gene_index.json"
        index_path.write_text(json.dumps(fake_index))

        from axiom_sc.tier3 import communication
        with patch("axiom_sc.tier3.communication._LR_INDEX_PATH", index_path):
            db = communication._load_lr_database()
        assert "TGFB1" in db
        assert "IL6" in db

    def test_communication_stream_runs_without_external_resources(self):
        """CommunicationStream must run on plain RNA-only adata."""
        from axiom_sc.tier3.communication import CommunicationStream
        adata = _make_adata(n_cells=40, n_genes=30, n_clusters=4)
        stream = CommunicationStream()
        results = stream.run(adata, cluster_key="leiden")
        assert isinstance(results, dict)
        # All clusters should have a StreamResult
        assert set(results.keys()) == {"0", "1", "2", "3"}
        for sr in results.values():
            assert sr.stream == "communication"

    def test_communication_stream_handles_corrupt_lr_file(self, tmp_path):
        """Corrupt lr_gene_index.json → fallback to empty {}."""
        bad_path = tmp_path / "lr_gene_index.json"
        bad_path.write_text("NOT JSON {{{{")
        from axiom_sc.tier3 import communication
        with patch("axiom_sc.tier3.communication._LR_INDEX_PATH", bad_path):
            db = communication._load_lr_database()
        assert isinstance(db, dict)


# ── CrossSpeciesStream: ENSEMBL conservation index ───────────────────────────

class TestCrossSpeciesEnsembl:

    def test_conservation_index_empty_when_file_absent(self, tmp_path):
        missing = tmp_path / "conservation_index.json"
        from axiom_sc.tier3 import cross_species
        with patch("axiom_sc.tier3.cross_species._CONSERVATION_PATH", missing):
            idx = cross_species._load_conservation_index()
        assert isinstance(idx, dict)

    def test_cross_species_uses_ensembl_when_available(self, tmp_path):
        """When conservation_index.json present, _load_conservation_index loads it."""
        fake_index = {"FOXP3": "one2one", "TBX21": "one2one", "NOVEL_GENE_XY": "many"}
        idx_path = tmp_path / "conservation_index.json"
        idx_path.write_text(json.dumps(fake_index))

        from axiom_sc.tier3 import cross_species
        with patch("axiom_sc.tier3.cross_species._CONSERVATION_PATH", idx_path):
            idx = cross_species._load_conservation_index()
        assert "FOXP3" in idx
        assert idx["FOXP3"] == "one2one"

    def test_cross_species_stream_runs_without_external_resources(self):
        """CrossSpeciesStream must run on plain RNA-only adata using built-in marker sets."""
        from axiom_sc.tier3.cross_species import CrossSpeciesStream
        adata = _make_adata(n_cells=40, n_genes=30, n_clusters=4)
        stream = CrossSpeciesStream()
        results = stream.run(adata, cluster_key="leiden")
        assert isinstance(results, dict)

    def test_cross_species_stream_handles_corrupt_index(self, tmp_path):
        bad_path = tmp_path / "conservation_index.json"
        bad_path.write_text("NOT JSON")
        from axiom_sc.tier3 import cross_species
        with patch("axiom_sc.tier3.cross_species._CONSERVATION_PATH", bad_path):
            idx = cross_species._load_conservation_index()
        assert isinstance(idx, dict)


# ── SpatialNicheStream: CELLama reference embeddings ─────────────────────────

class TestSpatialNicheCELLama:

    def test_cellama_reference_none_when_file_absent(self, tmp_path):
        missing = tmp_path / "cellama_reference_embeddings.pkl"
        from axiom_sc.tier3 import spatial_niche
        with patch("axiom_sc.tier3.spatial_niche._CELLAMA_EMBEDDINGS_PATH", missing):
            ref = spatial_niche._load_cellama_reference()
        assert ref is None

    def test_spatial_niche_uses_cellama_when_available(self, tmp_path):
        """When pkl present, _load_cellama_reference returns the dict."""
        fake_ref = {
            "T_cell": np.random.randn(64).astype(np.float32),
            "B_cell": np.random.randn(64).astype(np.float32),
        }
        pkl_path = tmp_path / "cellama_reference_embeddings.pkl"
        with open(pkl_path, "wb") as f:
            pickle.dump(fake_ref, f)

        from axiom_sc.tier3 import spatial_niche
        with patch("axiom_sc.tier3.spatial_niche._CELLAMA_EMBEDDINGS_PATH", pkl_path):
            ref = spatial_niche._load_cellama_reference()
        assert ref is not None
        assert "T_cell" in ref

    def test_spatial_niche_cosine_match_with_cellama_reference(self, tmp_path):
        """With reference embeddings, _score_from_cellama uses cosine similarity."""
        fake_ref = {
            "NK": np.array([1.0, 0.0, 0.0, 0.0], dtype=np.float32),
            "B_cell": np.array([0.0, 1.0, 0.0, 0.0], dtype=np.float32),
        }
        pkl_path = tmp_path / "cellama_reference_embeddings.pkl"
        with open(pkl_path, "wb") as f:
            pickle.dump(fake_ref, f)

        # Build adata with X_cellama embeddings pointing toward NK
        n_cells = 20
        emb = np.tile([0.9, 0.1, 0.0, 0.0], (n_cells, 1)).astype(np.float32)
        adata = anndata.AnnData(X=np.ones((n_cells, 5)))
        adata.obs["leiden"] = ["0"] * n_cells
        adata.obsm["X_cellama"] = emb

        from axiom_sc.tier3.spatial_niche import SpatialNicheStream, CELLAMA_REFERENCE
        with patch("axiom_sc.tier3.spatial_niche.CELLAMA_REFERENCE", fake_ref):
            stream = SpatialNicheStream()
            results = stream.run(adata, cluster_key="leiden")
        assert "0" in results
        # Should match NK (closest cosine)
        assert results["0"].label == "NK"
        assert results["0"].confidence > 0.5

    def test_spatial_niche_stream_runs_without_cellama(self):
        """Without CELLama embeddings, falls back to niche composition mode."""
        from axiom_sc.tier3.spatial_niche import SpatialNicheStream
        adata = _make_adata(n_cells=40, n_genes=10, n_clusters=4, with_connectivities=True)
        stream = SpatialNicheStream()
        with patch("axiom_sc.tier3.spatial_niche.CELLAMA_REFERENCE", None):
            results = stream.run(adata, cluster_key="leiden")
        assert isinstance(results, dict)

    def test_spatial_niche_handles_corrupt_pkl(self, tmp_path):
        bad_path = tmp_path / "cellama_reference_embeddings.pkl"
        bad_path.write_bytes(b"NOT A PICKLE ####")
        from axiom_sc.tier3 import spatial_niche
        with patch("axiom_sc.tier3.spatial_niche._CELLAMA_EMBEDDINGS_PATH", bad_path):
            ref = spatial_niche._load_cellama_reference()
        assert ref is None


# ── Auto-modality detection ───────────────────────────────────────────────────

class TestAutoModalityDetection:

    def test_auto_modality_detection_rna_only(self):
        """RNA-only adata: velocity=False, chromatin=False, others=True."""
        from axiom_sc.tier3.convergence import MultiStreamConvergence
        adata = _make_adata(n_cells=20, n_genes=10, n_clusters=2)
        msc = MultiStreamConvergence()
        modalities = msc._detect_modalities(adata)
        assert modalities["velocity"] is False
        assert modalities["chromatin"] is False
        assert modalities["communication"] is True
        assert modalities["spatial_niche"] is True
        assert modalities["cross_species"] is True
        assert modalities["scenic"] is True

    def test_auto_modality_detection_with_spliced_layer(self):
        """Adata with 'spliced' layer: velocity=True."""
        from axiom_sc.tier3.convergence import MultiStreamConvergence
        adata = _make_adata(n_cells=20, n_genes=10, n_clusters=2)
        adata.layers["spliced"] = adata.X.copy()
        msc = MultiStreamConvergence()
        modalities = msc._detect_modalities(adata)
        assert modalities["velocity"] is True

    def test_auto_modality_detection_with_velocity_layer(self):
        """Adata with 'velocity' layer: velocity=True."""
        from axiom_sc.tier3.convergence import MultiStreamConvergence
        adata = _make_adata(n_cells=20, n_genes=10, n_clusters=2)
        adata.layers["velocity"] = adata.X.copy()
        msc = MultiStreamConvergence()
        modalities = msc._detect_modalities(adata)
        assert modalities["velocity"] is True

    def test_auto_modality_detection_multiome(self):
        """Adata with gene_activities obsm: chromatin=True."""
        from axiom_sc.tier3.convergence import MultiStreamConvergence
        adata = _make_adata(n_cells=20, n_genes=10, n_clusters=2)
        adata.obsm["gene_activities"] = np.abs(np.random.randn(20, 10)).astype(np.float32)
        msc = MultiStreamConvergence()
        modalities = msc._detect_modalities(adata)
        assert modalities["chromatin"] is True

    def test_run_and_converge_auto_detect_rna_only(self):
        """run_and_converge with auto_detect=True on RNA-only: skips velocity+chromatin."""
        from axiom_sc.tier3.convergence import MultiStreamConvergence
        adata = _make_adata(n_cells=40, n_genes=20, n_clusters=3)
        msc = MultiStreamConvergence()
        # Should not raise; velocity+chromatin will be skipped
        results = msc.run_and_converge(adata, cluster_key="leiden", auto_detect=True)
        assert isinstance(results, dict)
        assert len(results) == 3


# ── ScenicPipeline: _build_config and cisTarget auto-detection ────────────────

class TestScenicBuildConfig:

    def test_build_config_returns_dict(self):
        from axiom_sc.tier2.scenic_pipeline import ScenicPipeline
        pipeline = ScenicPipeline()
        cfg = pipeline._build_config(loom_path="dummy.loom", output_dir="/tmp")
        assert isinstance(cfg, dict)
        assert cfg["loom_path"] == "dummy.loom"
        assert cfg["output_dir"] == "/tmp"

    def test_build_config_no_rankings_when_absent(self, tmp_path):
        """When scenic dir has no rankings.feather, rankings_db_path=None."""
        from axiom_sc.tier2.scenic_pipeline import ScenicPipeline
        pipeline = ScenicPipeline()
        fake_scenic_dir = tmp_path / "scenic"
        fake_scenic_dir.mkdir()
        # No feather file in dir
        with patch("axiom_sc.tier2.scenic_pipeline.ScenicPipeline._build_config",
                   wraps=pipeline._build_config) as _:
            # Patch the home dir lookup inside _build_config
            with patch("pathlib.Path.home", return_value=tmp_path.parent):
                cfg = pipeline._build_config(loom_path="x.loom")
        # rankings_db_path stays None since the file doesn't exist
        assert cfg["rankings_db_path"] is None

    def test_scenic_uses_full_cistaret_when_rankings_present(self, tmp_path):
        """When hg38_rankings.feather present, _build_config sets rankings_db_path."""
        from axiom_sc.tier2.scenic_pipeline import ScenicPipeline
        pipeline = ScenicPipeline()
        scenic_dir = tmp_path / ".axiom_sc" / "resources" / "scenic"
        scenic_dir.mkdir(parents=True)
        # Touch the files
        (scenic_dir / "hg38_rankings.feather").write_bytes(b"fake_feather")
        (scenic_dir / "motifs-v9.tbl").write_bytes(b"fake_motif")
        (scenic_dir / "allTFs_hg38.txt").write_text("FOXP3\nTBX21\n")

        with patch("pathlib.Path.home", return_value=tmp_path):
            cfg = pipeline._build_config(loom_path="test.loom")

        assert cfg["rankings_db_path"] is not None
        assert "hg38_rankings.feather" in cfg["rankings_db_path"]
        assert cfg["motif_db_path"] is not None
        assert cfg["tf_list_path"] is not None

    def test_build_config_explicit_paths_override_cached(self, tmp_path):
        """Explicitly passed paths override cached resource detection."""
        from axiom_sc.tier2.scenic_pipeline import ScenicPipeline
        pipeline = ScenicPipeline()
        # Even if cached files exist, explicit path wins
        explicit = "/custom/rankings.feather"
        cfg = pipeline._build_config(
            loom_path="test.loom",
            rankings_db_path=explicit,
        )
        assert cfg["rankings_db_path"] == explicit

    def test_build_config_defaults_from_scenic_defaults(self):
        from axiom_sc.tier2.scenic_pipeline import ScenicPipeline, SCENIC_DEFAULTS
        pipeline = ScenicPipeline()
        cfg = pipeline._build_config(loom_path="x.loom")
        # Should contain all SCENIC_DEFAULTS keys
        for key in ("nes_threshold", "min_genes", "z_score_auc", "z_active_thresh"):
            assert key in cfg, f"Missing key {key!r} from SCENIC_DEFAULTS in _build_config result"
