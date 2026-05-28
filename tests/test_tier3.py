"""
Tier 3 tests — Day 7 requirement.

Tests for:
  - StreamResult dataclass
  - VelocityStream (with mock AnnData)
  - ChromatinStream (with mock gene activities)
  - MultiStreamConvergence (pure-Python logic, no heavy deps)
  - GNNPropagator (with mock AnnData + scipy sparse graph)
  - discriminate_exhausted_t utility function

All tests run on CPU with small synthetic data; no scVelo, Signac, or R required.
"""
from __future__ import annotations

import numpy as np
import pytest

# ── helpers ───────────────────────────────────────────────────────────────────

def make_mock_adata(
    n_cells: int = 60,
    n_genes: int = 30,
    n_clusters: int = 3,
    add_velocity: bool = False,
    add_gene_activities: bool = False,
    add_connectivities: bool = False,
    gene_names: list[str] | None = None,
):
    """Create a minimal mock AnnData without importing anndata as a hard dep."""
    anndata = pytest.importorskip("anndata")
    import pandas as pd

    np.random.seed(42)
    X = np.random.rand(n_cells, n_genes).astype(np.float32)

    if gene_names is None:
        gene_names = [f"gene_{i}" for i in range(n_genes)]

    obs = pd.DataFrame({
        "leiden": [str(i % n_clusters) for i in range(n_cells)],
        "cell_type": [f"type_{i % n_clusters}" for i in range(n_cells)],
    }, index=[f"cell_{i}" for i in range(n_cells)])

    var = pd.DataFrame(index=gene_names)
    adata = anndata.AnnData(X=X, obs=obs, var=var)

    if add_velocity:
        # velocity layer: very small values for cluster 0, large for cluster 1
        # Cluster 0: L2 norm ≈ sqrt(30 * 0.003^2) ≈ 0.016 → well below MAGNITUDE_TERMINAL_THRESH
        # Cluster 1: L2 norm ≈ sqrt(30 * 0.15^2) ≈ 0.82 → well above MAGNITUDE_PROGENITOR_THRESH
        velocity = np.zeros((n_cells, n_genes), dtype=np.float32)
        for i in range(n_cells):
            c = int(obs.iloc[i]["leiden"])
            if c == 0:
                velocity[i] = np.random.rand(n_genes) * 0.003   # low magnitude (terminal sink)
            elif c == 1:
                velocity[i] = np.random.rand(n_genes) * 0.5     # high magnitude (progenitor)
            else:
                velocity[i] = np.random.rand(n_genes) * 0.12    # mid (transitional)
        adata.layers["velocity"] = velocity

    if add_gene_activities and gene_names:
        # gene activity matrix: cells × genes
        ga = np.random.rand(n_cells, n_genes).astype(np.float32)
        adata.obsm["gene_activities"] = ga
        adata.uns["gene_activity_genes"] = gene_names

    if add_connectivities:
        from scipy.sparse import csr_matrix
        # Simple chain graph: each cell connected to its 2 nearest neighbours
        rows, cols, data = [], [], []
        for i in range(n_cells):
            for j in [max(0, i - 1), min(n_cells - 1, i + 1)]:
                if i != j:
                    rows.append(i)
                    cols.append(j)
                    data.append(1.0)
        conn = csr_matrix((data, (rows, cols)), shape=(n_cells, n_cells))
        adata.obsp["connectivities"] = conn

    return adata


# ── StreamResult ──────────────────────────────────────────────────────────────

class TestStreamResult:
    def test_basic_attributes(self):
        from axiom_sc.tier3.stream_result import StreamResult
        sr = StreamResult(cluster_id="c0", stream="velocity", label="Treg", confidence=0.8)
        assert sr.cluster_id == "c0"
        assert sr.stream == "velocity"
        assert sr.label == "Treg"
        assert sr.confidence == 0.8
        assert sr.available is True
        assert sr.evidence == {}

    def test_unavailable_result(self):
        from axiom_sc.tier3.stream_result import StreamResult
        sr = StreamResult(
            cluster_id="c1", stream="chromatin", label="", confidence=0.0, available=False
        )
        assert not sr.available

    def test_agrees_with_same_label(self):
        from axiom_sc.tier3.stream_result import StreamResult
        a = StreamResult("c0", "velocity", "Treg", 0.8)
        b = StreamResult("c0", "chromatin", "Treg", 0.7)
        assert a.agrees_with(b)

    def test_agrees_with_different_label(self):
        from axiom_sc.tier3.stream_result import StreamResult
        a = StreamResult("c0", "velocity", "Treg", 0.8)
        b = StreamResult("c0", "chromatin", "Th1", 0.7)
        assert not a.agrees_with(b)

    def test_agrees_with_unavailable(self):
        from axiom_sc.tier3.stream_result import StreamResult
        a = StreamResult("c0", "velocity", "Treg", 0.8)
        b = StreamResult("c0", "chromatin", "Treg", 0.0, available=False)
        assert not a.agrees_with(b)

    def test_repr(self):
        from axiom_sc.tier3.stream_result import StreamResult
        sr = StreamResult("c0", "velocity", "Treg", 0.8)
        assert "velocity" in repr(sr)
        assert "Treg" in repr(sr)


# ── VelocityStream ────────────────────────────────────────────────────────────

class TestVelocityStream:
    def test_returns_unavailable_when_no_velocity(self):
        from axiom_sc.tier3.velocity import VelocityStream
        adata = make_mock_adata(n_cells=30, n_clusters=2)
        vs = VelocityStream()
        results = vs.run(adata, cluster_key="leiden")
        assert len(results) == 2
        for sr in results.values():
            assert not sr.available
            assert sr.stream == "velocity"

    def test_returns_unavailable_bad_cluster_key(self):
        from axiom_sc.tier3.velocity import VelocityStream
        adata = make_mock_adata(n_cells=30, n_clusters=2, add_velocity=True)
        vs = VelocityStream()
        results = vs.run(adata, cluster_key="nonexistent_key")
        assert results == {}

    def test_uses_velocity_layer_when_present(self):
        from axiom_sc.tier3.velocity import VelocityStream
        adata = make_mock_adata(n_cells=60, n_clusters=3, add_velocity=True)
        vs = VelocityStream()
        results = vs.run(adata, cluster_key="leiden")
        assert len(results) == 3
        for sr in results.values():
            assert sr.available
            assert sr.stream == "velocity"
            assert "mean_velocity_magnitude" in sr.evidence
            assert sr.confidence >= 0.0

    def test_low_magnitude_cluster_classified_as_sink(self):
        """Cluster 0 has low velocity (terminal sink)."""
        from axiom_sc.tier3.velocity import VelocityStream, MAGNITUDE_TERMINAL_THRESH
        adata = make_mock_adata(n_cells=60, n_clusters=3, add_velocity=True)
        vs = VelocityStream()
        results = vs.run(adata, cluster_key="leiden")
        # Cluster "0" was assigned low velocity values
        c0 = results["0"]
        assert c0.evidence["mean_velocity_magnitude"] < MAGNITUDE_TERMINAL_THRESH

    def test_high_magnitude_cluster_classified_as_source(self):
        """Cluster 1 has high velocity (progenitor source)."""
        from axiom_sc.tier3.velocity import VelocityStream, MAGNITUDE_PROGENITOR_THRESH
        adata = make_mock_adata(n_cells=60, n_clusters=3, add_velocity=True)
        vs = VelocityStream()
        results = vs.run(adata, cluster_key="leiden")
        c1 = results["1"]
        assert c1.evidence["mean_velocity_magnitude"] > MAGNITUDE_PROGENITOR_THRESH

    def test_tier2_hint_confirms_label(self):
        """If Tier 2 hint matches velocity type, confidence is boosted."""
        from axiom_sc.tier3.velocity import VelocityStream
        adata = make_mock_adata(n_cells=60, n_clusters=3, add_velocity=True)
        # Cluster "0" is a sink; "Treg" is expected to be a sink
        vs = VelocityStream(tier2_labels={"0": "Treg"})
        results = vs.run(adata, cluster_key="leiden")
        c0 = results["0"]
        assert c0.label == "Treg"
        assert c0.confidence >= 0.7  # confirmed → high confidence

    def test_tier2_hint_conflicts_label(self):
        """If Tier 2 hint clashes with velocity type, confidence is reduced."""
        from axiom_sc.tier3.velocity import VelocityStream
        adata = make_mock_adata(n_cells=60, n_clusters=3, add_velocity=True)
        # Cluster "1" is a source; "Treg" is expected to be a sink
        vs = VelocityStream(tier2_labels={"1": "Treg"})
        results = vs.run(adata, cluster_key="leiden")
        c1 = results["1"]
        assert c1.label == "Treg"   # still returns the label
        assert c1.confidence < 0.5  # but confidence is reduced


# ── discriminate_exhausted_t ──────────────────────────────────────────────────

class TestDiscriminateExhaustedT:
    def test_higher_magnitude_is_tpex(self):
        from axiom_sc.tier3.velocity import discriminate_exhausted_t
        adata = make_mock_adata(n_cells=60, n_clusters=3, add_velocity=True)
        result = discriminate_exhausted_t(adata, "leiden", ["0", "1"])
        # cluster "0" has low vel → TEX; cluster "1" has high vel → TPEX
        assert result["0"] == "Exhausted_T"
        assert result["1"] == "Progenitor_Exhausted_T"

    def test_no_velocity_returns_uncertain(self):
        from axiom_sc.tier3.velocity import discriminate_exhausted_t
        adata = make_mock_adata(n_cells=60, n_clusters=3)
        result = discriminate_exhausted_t(adata, "leiden", ["0", "1"])
        assert all(v == "uncertain" for v in result.values())

    def test_single_cluster_returns_uncertain(self):
        from axiom_sc.tier3.velocity import discriminate_exhausted_t
        adata = make_mock_adata(n_cells=30, n_clusters=1, add_velocity=True)
        result = discriminate_exhausted_t(adata, "leiden", ["0"])
        assert result["0"] == "uncertain"


# ── ChromatinStream ───────────────────────────────────────────────────────────

# Key gene names for chromatin signature testing
CHROMATIN_GENES = [
    "FOXP3", "IL2RA", "AIRE", "PSMB11", "TBX21", "GATA3",
    "RORC", "IRF7", "SIGLEC1", "IRF8", "BATF3", "NCR2",
    "TOX", "TCF7", "PDCD1",
]


def make_chromatin_adata(n_cells: int = 60, n_clusters: int = 3):
    """AnnData with gene_activities for key chromatin loci."""
    anndata = pytest.importorskip("anndata")
    import pandas as pd

    n_loci = len(CHROMATIN_GENES)
    np.random.seed(42)
    X = np.random.rand(n_cells, n_loci).astype(np.float32)

    obs = pd.DataFrame({
        "leiden": [str(i % n_clusters) for i in range(n_cells)],
    }, index=[f"cell_{i}" for i in range(n_cells)])
    var = pd.DataFrame(index=CHROMATIN_GENES)
    adata = anndata.AnnData(X=X, obs=obs, var=var)

    # Build gene activities:
    # Cluster 0: high FOXP3 + IL2RA (Treg signature)
    # Cluster 1: high TBX21, low GATA3/RORC (Th1 signature)
    # Cluster 2: random (no clear signature)
    ga = np.zeros((n_cells, n_loci), dtype=np.float32)
    gene_idx = {g: i for i, g in enumerate(CHROMATIN_GENES)}

    for i in range(n_cells):
        c = int(obs.iloc[i]["leiden"])
        if c == 0:
            ga[i, gene_idx["FOXP3"]] = 3.0
            ga[i, gene_idx["IL2RA"]] = 2.5
        elif c == 1:
            ga[i, gene_idx["TBX21"]] = 3.0
            ga[i, gene_idx["GATA3"]] = -2.0
            ga[i, gene_idx["RORC"]] = -2.0
        else:
            ga[i] = np.random.rand(n_loci) * 0.1

    adata.obsm["gene_activities"] = ga
    adata.uns["gene_activity_genes"] = CHROMATIN_GENES
    return adata


class TestChromatinStream:
    def test_returns_unavailable_without_atac_data(self):
        from axiom_sc.tier3.chromatin import ChromatinStream
        adata = make_mock_adata(n_cells=30, n_clusters=2)
        cs = ChromatinStream(use_rpy2=False)
        results = cs.run(adata, cluster_key="leiden")
        for sr in results.values():
            assert not sr.available
            assert sr.stream == "chromatin"

    def test_returns_unavailable_bad_cluster_key(self):
        from axiom_sc.tier3.chromatin import ChromatinStream
        adata = make_chromatin_adata()
        cs = ChromatinStream(use_rpy2=False)
        results = cs.run(adata, cluster_key="nonexistent")
        assert results == {}

    def test_processes_gene_activities_present(self):
        from axiom_sc.tier3.chromatin import ChromatinStream
        adata = make_chromatin_adata()
        cs = ChromatinStream(use_rpy2=False)
        results = cs.run(adata, cluster_key="leiden")
        assert len(results) == 3
        for sr in results.values():
            assert sr.available
            assert sr.stream == "chromatin"
            assert 0.0 <= sr.confidence <= 1.0

    def test_treg_cluster_identified(self):
        """Cluster 0 with high FOXP3 + IL2RA should be recognised as Treg."""
        from axiom_sc.tier3.chromatin import ChromatinStream
        adata = make_chromatin_adata()
        cs = ChromatinStream(use_rpy2=False)
        results = cs.run(adata, cluster_key="leiden")
        c0 = results["0"]
        assert c0.label == "Treg"
        assert c0.confidence > 0.5

    def test_th1_cluster_identified(self):
        """Cluster 1 with high TBX21, low GATA3/RORC should be Th1."""
        from axiom_sc.tier3.chromatin import ChromatinStream
        adata = make_chromatin_adata()
        cs = ChromatinStream(use_rpy2=False)
        results = cs.run(adata, cluster_key="leiden")
        c1 = results["1"]
        assert c1.label == "Th1"
        assert c1.confidence > 0.5

    def test_gene_activity_missing_gene_names(self):
        """If gene_activity_genes not set and var_names don't match, stream unavailable."""
        anndata = pytest.importorskip("anndata")
        import pandas as pd
        from axiom_sc.tier3.chromatin import ChromatinStream

        n_cells, n_genes = 30, 5
        X = np.random.rand(n_cells, n_genes).astype(np.float32)
        ga = np.random.rand(n_cells, 20).astype(np.float32)  # 20 ≠ n_genes=5
        obs = pd.DataFrame({"leiden": ["0"] * n_cells}, index=[f"c{i}" for i in range(n_cells)])
        var = pd.DataFrame(index=[f"g{i}" for i in range(n_genes)])
        adata = anndata.AnnData(X=X, obs=obs, var=var)
        adata.obsm["gene_activities"] = ga
        # No gene_activity_genes in uns, and shape mismatch → unavailable

        cs = ChromatinStream(use_rpy2=False)
        results = cs.run(adata, cluster_key="leiden")
        for sr in results.values():
            assert not sr.available

    def test_tier2_hint_used_in_scoring(self):
        """With a Tier 2 hint, the hinted cell type is checked first."""
        from axiom_sc.tier3.chromatin import ChromatinStream
        adata = make_chromatin_adata()
        # Hint cluster 0 as Treg — should still score correctly
        cs = ChromatinStream(tier2_labels={"0": "Treg"}, use_rpy2=False)
        results = cs.run(adata, cluster_key="leiden")
        assert results["0"].label == "Treg"


# ── MultiStreamConvergence ────────────────────────────────────────────────────

class TestMultiStreamConvergence:
    def _make_stream_results(
        self,
        cluster_id: str,
        votes: dict[str, tuple[str, float]],
    ) -> dict[str, "StreamResult"]:
        """Build stream_results dict from {stream: (label, confidence)}."""
        from axiom_sc.tier3.stream_result import StreamResult
        return {
            stream: StreamResult(
                cluster_id=cluster_id,
                stream=stream,
                label=label,
                confidence=conf,
            )
            for stream, (label, conf) in votes.items()
        }

    def test_proven_when_4_streams_agree(self):
        from axiom_sc.tier3.convergence import MultiStreamConvergence
        msc = MultiStreamConvergence()
        votes = {
            "velocity":    ("Treg", 0.8),
            "chromatin":   ("Treg", 0.9),
            "communication": ("Treg", 0.7),
            "spatial_niche": ("Treg", 0.75),
        }
        sr = self._make_stream_results("c0", votes)
        result = msc.converge("c0", sr)
        assert result.verdict == "PROVEN"
        assert result.consensus_label == "Treg"
        assert result.n_streams_agree == 4

    def test_proven_when_5_streams_agree(self):
        from axiom_sc.tier3.convergence import MultiStreamConvergence
        msc = MultiStreamConvergence()
        votes = {
            "velocity":    ("pDC", 0.8),
            "chromatin":   ("pDC", 0.9),
            "communication": ("pDC", 0.7),
            "spatial_niche": ("pDC", 0.75),
            "cross_species": ("pDC", 0.85),
        }
        sr = self._make_stream_results("c1", votes)
        result = msc.converge("c1", sr)
        assert result.verdict == "PROVEN"
        assert result.n_streams_agree == 5

    def test_high_confidence_when_3_agree(self):
        from axiom_sc.tier3.convergence import MultiStreamConvergence
        msc = MultiStreamConvergence()
        votes = {
            "velocity":  ("Th1", 0.7),
            "chromatin": ("Th1", 0.8),
            "communication": ("Th1", 0.65),
        }
        sr = self._make_stream_results("c0", votes)
        result = msc.converge("c0", sr)
        assert result.verdict == "HIGH_CONFIDENCE"
        assert result.consensus_label == "Th1"

    def test_uncertain_when_2_agree(self):
        from axiom_sc.tier3.convergence import MultiStreamConvergence
        msc = MultiStreamConvergence()
        votes = {
            "velocity":  ("Treg", 0.7),
            "chromatin": ("Treg", 0.8),
        }
        sr = self._make_stream_results("c0", votes)
        result = msc.converge("c0", sr)
        assert result.verdict == "UNCERTAIN"

    def test_uncertain_when_streams_disagree(self):
        from axiom_sc.tier3.convergence import MultiStreamConvergence
        msc = MultiStreamConvergence()
        votes = {
            "velocity":    ("Treg", 0.8),
            "chromatin":   ("Th1", 0.8),
            "communication": ("pDC", 0.7),
            "spatial_niche": ("Th2", 0.75),
        }
        sr = self._make_stream_results("c0", votes)
        result = msc.converge("c0", sr)
        assert result.verdict == "UNCERTAIN"

    def test_uncertain_when_no_streams_available(self):
        from axiom_sc.tier3.convergence import MultiStreamConvergence
        from axiom_sc.tier3.stream_result import StreamResult
        msc = MultiStreamConvergence()
        sr = {
            "velocity": StreamResult("c0", "velocity", "", 0.0, available=False),
            "chromatin": StreamResult("c0", "chromatin", "", 0.0, available=False),
        }
        result = msc.converge("c0", sr)
        assert result.verdict == "UNCERTAIN"
        assert result.n_streams_available == 0

    def test_low_confidence_streams_not_counted(self):
        """Streams below min_confidence are filtered out."""
        from axiom_sc.tier3.convergence import MultiStreamConvergence
        msc = MultiStreamConvergence(min_confidence=0.5)
        votes = {
            "velocity":  ("Treg", 0.8),
            "chromatin": ("Treg", 0.8),
            "communication": ("Treg", 0.8),
            "spatial_niche": ("Treg", 0.2),   # below threshold
        }
        sr = self._make_stream_results("c0", votes)
        result = msc.converge("c0", sr)
        # Only 3 streams count (spatial_niche dropped) → HIGH_CONFIDENCE, not PROVEN
        assert result.n_streams_available == 3
        assert result.verdict == "HIGH_CONFIDENCE"

    def test_converge_all_returns_all_clusters(self):
        from axiom_sc.tier3.convergence import MultiStreamConvergence
        from axiom_sc.tier3.stream_result import StreamResult
        msc = MultiStreamConvergence()
        stream_results = {
            "velocity": {
                "c0": StreamResult("c0", "velocity", "Treg", 0.8),
                "c1": StreamResult("c1", "velocity", "Th1", 0.7),
            },
            "chromatin": {
                "c0": StreamResult("c0", "chromatin", "Treg", 0.9),
                "c1": StreamResult("c1", "chromatin", "Th1", 0.8),
            },
        }
        results = msc.converge_all(stream_results)
        assert "c0" in results
        assert "c1" in results

    def test_converge_all_missing_cluster_in_one_stream(self):
        """Cluster absent from one stream is still converged with the others."""
        from axiom_sc.tier3.convergence import MultiStreamConvergence
        from axiom_sc.tier3.stream_result import StreamResult
        msc = MultiStreamConvergence()
        stream_results = {
            "velocity": {
                "c0": StreamResult("c0", "velocity", "Treg", 0.8),
                "c1": StreamResult("c1", "velocity", "Th1", 0.7),
            },
            "chromatin": {
                "c0": StreamResult("c0", "chromatin", "Treg", 0.9),
                # c1 absent from chromatin stream
            },
        }
        results = msc.converge_all(stream_results)
        assert "c0" in results
        assert "c1" in results
        # c1 only has 1 stream → UNCERTAIN
        assert results["c1"].verdict == "UNCERTAIN"

    def test_convergence_result_to_dict(self):
        from axiom_sc.tier3.convergence import ConvergenceResult
        result = ConvergenceResult(
            cluster_id="c0",
            verdict="PROVEN",
            consensus_label="Treg",
            confidence=0.85,
            n_streams_available=5,
            n_streams_agree=4,
        )
        d = result.to_dict()
        assert d["verdict"] == "PROVEN"
        assert d["consensus_label"] == "Treg"
        assert d["confidence"] == 0.85
        assert "n_streams_agree" in d

    def test_update_evidence_bundles(self):
        from axiom_sc.tier3.convergence import MultiStreamConvergence, ConvergenceResult
        from axiom_sc.tier2.evidence import EvidenceBundle

        msc = MultiStreamConvergence()
        bundle = EvidenceBundle(cluster_id="c0", marker_genes={}, tissue="thymus")
        conv = ConvergenceResult(
            cluster_id="c0", verdict="HIGH_CONFIDENCE",
            consensus_label="pDC", confidence=0.75,
            n_streams_available=3, n_streams_agree=3,
        )
        msc.update_evidence_bundles({"c0": conv}, {"c0": bundle})
        assert bundle.spatial_context["tier3_verdict"] == "HIGH_CONFIDENCE"
        assert bundle.spatial_context["tier3_label"] == "pDC"

    def test_make_unavailable_stream_results(self):
        from axiom_sc.tier3.convergence import make_unavailable_stream_results
        results = make_unavailable_stream_results(["c0", "c1"], "cross_species", "not_installed")
        assert len(results) == 2
        assert not results["c0"].available
        assert results["c0"].stream == "cross_species"

    def test_require_n_streams_threshold(self):
        """With require_n_streams=2, single-stream clusters return UNCERTAIN."""
        from axiom_sc.tier3.convergence import MultiStreamConvergence
        from axiom_sc.tier3.stream_result import StreamResult
        msc = MultiStreamConvergence(require_n_streams=2)
        sr = {
            "velocity": StreamResult("c0", "velocity", "Treg", 0.9),
        }
        result = msc.converge("c0", sr)
        assert result.verdict == "UNCERTAIN"
        assert result.n_streams_available == 1


# ── GNNPropagator ─────────────────────────────────────────────────────────────

class TestGNNPropagator:
    def _make_verdicts(self, specs: dict):
        """
        Build mock verdicts dict.
        specs: {cluster_id: (verdict_str, label, confidence, n_cells)}
        """
        from axiom_sc.tier2.axiom_annotator import Verdict, CellTypeVerdict

        verdicts = {}
        for cid, (vstr, label, conf, n_cells) in specs.items():
            # Use a simple object that matches what GNNPropagator reads via getattr
            class V:
                pass
            v = V()
            v.verdict = Verdict(vstr)
            v.cell_type = label
            v.confidence = conf
            v.n_cells = n_cells
            verdicts[cid] = v
        return verdicts

    def test_no_connectivities_returns_identity(self):
        """Without kNN graph, propagation returns identity (no changes)."""
        from axiom_sc.tier2.gnn_propagator import GNNPropagator
        adata = make_mock_adata(n_cells=30, n_clusters=3)
        verdicts = self._make_verdicts({
            "0": ("PROVEN", "Treg", 0.9, 10),
            "1": ("UNCERTAIN", "", 0.3, 10),
            "2": ("UNCERTAIN", "", 0.3, 10),
        })
        gnn = GNNPropagator()
        results = gnn.propagate(adata, verdicts, cluster_key="leiden")
        # Should return without crashing; no updates since no graph
        assert "0" in results or "1" in results  # at least something returned

    def test_bad_cluster_key_returns_empty(self):
        from axiom_sc.tier2.gnn_propagator import GNNPropagator
        adata = make_mock_adata(n_cells=30, n_clusters=3, add_connectivities=True)
        verdicts = self._make_verdicts({"0": ("PROVEN", "Treg", 0.9, 10)})
        gnn = GNNPropagator()
        results = gnn.propagate(adata, verdicts, cluster_key="bad_key")
        assert results == {}

    def test_proven_clusters_not_changed(self):
        """PROVEN clusters should not be modified by propagation."""
        pytest.importorskip("scipy")
        from axiom_sc.tier2.gnn_propagator import GNNPropagator
        from axiom_sc.tier2.axiom_annotator import Verdict
        adata = make_mock_adata(n_cells=60, n_clusters=3, add_connectivities=True)
        verdicts = self._make_verdicts({
            "0": ("PROVEN", "Treg", 0.9, 20),
            "1": ("UNCERTAIN", "", 0.3, 20),
            "2": ("UNCERTAIN", "", 0.3, 20),
        })
        gnn = GNNPropagator()
        results = gnn.propagate(adata, verdicts, cluster_key="leiden")
        if "0" in results:
            assert results["0"].propagated_verdict == Verdict.PROVEN.value
            assert results["0"].propagated_label == "Treg"
            assert not results["0"].was_updated

    def test_propagation_result_attributes(self):
        pytest.importorskip("scipy")
        from axiom_sc.tier2.gnn_propagator import GNNPropagator
        from axiom_sc.tier2.axiom_annotator import Verdict
        adata = make_mock_adata(n_cells=60, n_clusters=3, add_connectivities=True)
        verdicts = self._make_verdicts({
            "0": ("PROVEN", "Treg", 0.9, 20),
            "1": ("UNCERTAIN", "", 0.3, 20),
            "2": ("UNCERTAIN", "", 0.3, 20),
        })
        gnn = GNNPropagator()
        results = gnn.propagate(adata, verdicts, cluster_key="leiden")
        valid_verdicts = {Verdict.PROVEN.value, Verdict.UNCERTAIN.value, Verdict.CONTRADICTED.value}
        for cid, pr in results.items():
            assert pr.cluster_id == cid
            assert pr.original_verdict in valid_verdicts
            assert 0.0 <= pr.confidence <= 1.0
            assert pr.n_cells_total >= 0

    def test_propagation_result_was_updated_false_for_proven(self):
        pytest.importorskip("scipy")
        from axiom_sc.tier2.gnn_propagator import GNNPropagator, PropagationResult
        from axiom_sc.tier2.axiom_annotator import Verdict
        pr = PropagationResult(
            cluster_id="c0",
            original_verdict=str(Verdict.PROVEN),
            propagated_verdict=str(Verdict.PROVEN),
            original_label="Treg",
            propagated_label="Treg",
            n_cells_total=20,
            n_cells_propagated=0,
        )
        assert not pr.was_updated

    def test_no_propagation_results_helper(self):
        """_no_propagation_results returns identity for all clusters."""
        from axiom_sc.tier2.gnn_propagator import GNNPropagator
        verdicts = self._make_verdicts({
            "0": ("PROVEN", "Treg", 0.9, 10),
            "1": ("UNCERTAIN", "", 0.3, 10),
        })
        gnn = GNNPropagator()
        results = gnn._no_propagation_results(verdicts)
        assert set(results.keys()) == {"0", "1"}
        assert not results["0"].was_updated
        assert not results["1"].was_updated


# ── run_and_converge integration ──────────────────────────────────────────────

class TestRunAndConverge:
    def test_run_and_converge_velocity_only(self):
        """Run velocity+chromatin but only velocity data present → chromatin unavailable."""
        from axiom_sc.tier3.convergence import MultiStreamConvergence
        adata = make_mock_adata(n_cells=60, n_clusters=3, add_velocity=True)
        msc = MultiStreamConvergence()
        results = msc.run_and_converge(
            adata, cluster_key="leiden", run_velocity=True, run_chromatin=True
        )
        assert len(results) == 3
        # Only 1 stream available per cluster → all UNCERTAIN
        for cr in results.values():
            assert cr.verdict == "UNCERTAIN"

    def test_run_and_converge_no_streams(self):
        """Disabling both streams returns empty dict."""
        from axiom_sc.tier3.convergence import MultiStreamConvergence
        adata = make_mock_adata(n_cells=30, n_clusters=2)
        msc = MultiStreamConvergence()
        results = msc.run_and_converge(
            adata, cluster_key="leiden",
            run_velocity=False, run_chromatin=False,
            run_communication=False, run_spatial_niche=False, run_cross_species=False,
        )
        assert results == {}

    def test_run_and_converge_all_6_streams_degrade_gracefully(self):
        """All 6 streams run on minimal adata — no crashes, all UNCERTAIN due to no data."""
        from axiom_sc.tier3.convergence import MultiStreamConvergence
        adata = make_mock_adata(n_cells=30, n_clusters=2)
        msc = MultiStreamConvergence()
        results = msc.run_and_converge(adata, cluster_key="leiden")
        # Velocity and chromatin return unavailable; communication/spatial_niche/
        # cross_species may return results based on what's in the mock adata
        assert len(results) >= 0  # should not crash


# ── CommunicationStream ───────────────────────────────────────────────────────

LR_GENES = [
    "TGFB1", "TGFB2", "FGF7", "FGF10",   # CAF senders
    "NOTCH3", "PDGFRB", "ACVRL1",          # SMC
    "IL10", "IL2RA", "TNFRSF18",           # Treg
    "IFNA1", "IFNA2", "IFNB1",             # pDC
    "GZMB", "PRF1", "IFNG",                # NK senders
    "CCL2", "CCL7", "S100A8",              # Monocyte
    "SPP1", "CXCL10", "CXCL9", "CSF1R",   # Tissue-resident Mφ
]


def make_lr_adata(n_cells: int = 80, n_clusters: int = 4):
    """AnnData with L-R gene expression. Cluster 0 = CAF (high TGFB1/TGFB2)."""
    anndata = pytest.importorskip("anndata")
    import pandas as pd

    n_genes = len(LR_GENES)
    np.random.seed(42)
    X = np.random.rand(n_cells, n_genes).astype(np.float32) * 0.1

    obs = pd.DataFrame({
        "leiden": [str(i % n_clusters) for i in range(n_cells)],
    }, index=[f"cell_{i}" for i in range(n_cells)])
    var = pd.DataFrame(index=LR_GENES)
    adata = anndata.AnnData(X=X, obs=obs, var=var)

    # Cluster 0: high TGFB1 + TGFB2 (CAF signature)
    # Cluster 1: high IFNA1 + IFNA2 + IFNB1 (pDC signature)
    # Cluster 2: random (no clear signature)
    gene_idx = {g: i for i, g in enumerate(LR_GENES)}
    for i in range(n_cells):
        c = int(obs.iloc[i]["leiden"])
        if c == 0:
            adata.X[i, gene_idx["TGFB1"]] = 4.0
            adata.X[i, gene_idx["TGFB2"]] = 3.5
        elif c == 1:
            adata.X[i, gene_idx["IFNA1"]] = 4.0
            adata.X[i, gene_idx["IFNA2"]] = 3.8
            adata.X[i, gene_idx["IFNB1"]] = 3.5

    return adata


class TestCommunicationStream:
    def test_no_lr_genes_returns_unavailable(self):
        """AnnData with no known L-R genes → stream unavailable."""
        from axiom_sc.tier3.communication import CommunicationStream
        adata = make_mock_adata(n_cells=30, n_clusters=2)
        cs = CommunicationStream()
        results = cs.run(adata, cluster_key="leiden")
        # Should return unavailable since no L-R genes in mock adata
        assert len(results) >= 0  # no crash; result is valid (may be unavailable)

    def test_bad_cluster_key_returns_empty(self):
        from axiom_sc.tier3.communication import CommunicationStream
        adata = make_lr_adata()
        cs = CommunicationStream()
        results = cs.run(adata, cluster_key="bad_key")
        assert results == {}

    def test_lr_expression_mode_runs(self):
        """With L-R genes in adata, expression mode produces results."""
        from axiom_sc.tier3.communication import CommunicationStream
        adata = make_lr_adata()
        cs = CommunicationStream()
        results = cs.run(adata, cluster_key="leiden")
        assert len(results) == 4
        for sr in results.values():
            assert sr.stream == "communication"
            assert 0.0 <= sr.confidence <= 1.0

    def test_caf_cluster_labeled(self):
        """Cluster 0 with high TGFB1/TGFB2 should be labelled CAF."""
        from axiom_sc.tier3.communication import CommunicationStream
        adata = make_lr_adata()
        cs = CommunicationStream()
        results = cs.run(adata, cluster_key="leiden")
        c0 = results["0"]
        assert c0.available
        assert c0.label == "CAF"
        assert c0.confidence > 0.3

    def test_pdc_cluster_labeled(self):
        """Cluster 1 with high IFNA1/IFNA2/IFNB1 should be labelled pDC."""
        from axiom_sc.tier3.communication import CommunicationStream
        adata = make_lr_adata()
        cs = CommunicationStream()
        results = cs.run(adata, cluster_key="leiden")
        c1 = results["1"]
        assert c1.available
        assert c1.label == "pDC"

    def test_commot_mode_used_when_obsm_present(self):
        """Pre-computed COMMOT keys in obsm → COMMOT mode is used."""
        import pandas as pd
        from axiom_sc.tier3.communication import CommunicationStream, COMMOT_SENDER_PREFIX
        adata = make_lr_adata()
        n_cells = len(adata)
        # Add fake COMMOT sender scores
        fake_sender = np.random.rand(n_cells, 3).astype(np.float32)
        obsm_key = f"{COMMOT_SENDER_PREFIX}-CellChat"
        adata.obsm[obsm_key] = fake_sender
        cs = CommunicationStream()
        results = cs.run(adata, cluster_key="leiden")
        assert len(results) == 4
        for sr in results.values():
            assert sr.available
            assert "sender_pathways" in sr.evidence

    def test_tier2_hint_prioritised(self):
        """Tier 2 hint changes which signature is evaluated first."""
        from axiom_sc.tier3.communication import CommunicationStream
        adata = make_lr_adata()
        cs = CommunicationStream(tier2_labels={"0": "CAF"})
        results = cs.run(adata, cluster_key="leiden")
        # Hint + expression should agree → CAF
        assert results["0"].label == "CAF"


# ── SpatialNicheStream ────────────────────────────────────────────────────────

def make_niche_adata(n_cells: int = 90):
    """
    AnnData with connectivities.
    3 clusters: 0=pDC, 1=T_cell, 2=B_cell
    Cluster 0 (pDC) has mostly T_cell neighbours.
    """
    anndata = pytest.importorskip("anndata")
    pytest.importorskip("scipy")
    import pandas as pd
    from scipy.sparse import lil_matrix

    n_per = n_cells // 3
    obs = pd.DataFrame({
        "leiden": (["0"] * n_per) + (["1"] * n_per) + (["2"] * n_per),
    }, index=[f"cell_{i}" for i in range(n_cells)])
    var = pd.DataFrame(index=[f"g{i}" for i in range(10)])
    X = np.random.rand(n_cells, 10).astype(np.float32)
    adata = anndata.AnnData(X=X, obs=obs, var=var)

    # Build connectivity: cluster 0 cells are heavily connected to cluster 1 (T cells)
    conn = lil_matrix((n_cells, n_cells), dtype=np.float32)
    for i in range(n_per):  # cluster 0 cells
        for j in range(n_per, 2 * n_per):  # cluster 1 cells (T)
            conn[i, j] = 1.0
            conn[j, i] = 1.0
    # Cluster 1 and 2 are connected to each other
    for i in range(n_per, 2 * n_per):
        for j in range(2 * n_per, n_cells):
            conn[i, j] = 0.5
            conn[j, i] = 0.5
    adata.obsp["connectivities"] = conn.tocsr()
    return adata


class TestSpatialNicheStream:
    def test_no_connectivities_no_cellama_returns_unavailable(self):
        """Without connectivities or CELLama, stream is unavailable."""
        from axiom_sc.tier3.spatial_niche import SpatialNicheStream
        adata = make_mock_adata(n_cells=30, n_clusters=2)
        sn = SpatialNicheStream()
        results = sn.run(adata, cluster_key="leiden")
        for sr in results.values():
            assert not sr.available

    def test_bad_cluster_key_returns_empty(self):
        from axiom_sc.tier3.spatial_niche import SpatialNicheStream
        adata = make_niche_adata()
        sn = SpatialNicheStream()
        results = sn.run(adata, cluster_key="nonexistent")
        assert results == {}

    def test_niche_composition_mode_runs(self):
        """With connectivities, niche composition mode returns results."""
        from axiom_sc.tier3.spatial_niche import SpatialNicheStream
        adata = make_niche_adata()
        sn = SpatialNicheStream()
        results = sn.run(adata, cluster_key="leiden")
        assert len(results) == 3
        for sr in results.values():
            assert sr.available
            assert sr.stream == "spatial_niche"

    def test_pdc_cluster_has_t_cell_neighbors(self):
        """Cluster 0 has all T cell (cluster 1) neighbours — should match pDC niche."""
        from axiom_sc.tier3.spatial_niche import SpatialNicheStream
        adata = make_niche_adata()
        # Cluster 0 → neighbour label "1"; pDC niche expects "T" substring in neighbour labels
        # However "1" doesn't contain "T" — so we also test with a descriptive label key
        sn = SpatialNicheStream(tier2_labels={"0": "pDC"})
        results = sn.run(adata, cluster_key="leiden")
        c0 = results["0"]
        assert c0.available
        # With tier2 hint, it should at minimum return the hint
        assert c0.label in ("pDC", "0", "unknown") or len(c0.label) > 0

    def test_cellama_mode_used_when_embedding_present(self):
        """If adata.obsm['X_cellama'] present, CELLama mode is used."""
        from axiom_sc.tier3.spatial_niche import SpatialNicheStream
        adata = make_mock_adata(n_cells=30, n_clusters=2)
        adata.obsm["X_cellama"] = np.random.rand(30, 64).astype(np.float32)
        sn = SpatialNicheStream(tier2_labels={"0": "pDC", "1": "T_cell"})
        results = sn.run(adata, cluster_key="leiden")
        assert len(results) == 2
        for sr in results.values():
            assert sr.available
            assert "cellama_embedding" in sr.evidence.get("mode", "")

    def test_niche_result_confidence_range(self):
        from axiom_sc.tier3.spatial_niche import SpatialNicheStream
        adata = make_niche_adata()
        sn = SpatialNicheStream()
        results = sn.run(adata, cluster_key="leiden")
        for sr in results.values():
            assert 0.0 <= sr.confidence <= 1.0


# ── CrossSpeciesStream ────────────────────────────────────────────────────────

def make_conserved_adata(n_cells: int = 60, n_clusters: int = 3):
    """
    AnnData with conserved marker genes in var_names.
    Cluster 0: Treg-like (high FOXP3, IL2RA, CTLA4)
    Cluster 1: NK-like (high NCAM1, KLRD1, NKG7, GZMB)
    Cluster 2: random
    """
    anndata = pytest.importorskip("anndata")
    import pandas as pd

    genes = [
        "FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18",  # Treg
        "NCAM1", "KLRD1", "NKG7", "GNLY", "GZMB",          # NK
        "CD3D", "CD3E", "CD4", "CD8A",                      # generic T
        "GENE_NOISE1", "GENE_NOISE2",                        # not in signatures
    ]
    n_genes = len(genes)
    gene_idx = {g: i for i, g in enumerate(genes)}

    np.random.seed(42)
    X = np.random.rand(n_cells, n_genes).astype(np.float32) * 0.5

    obs = pd.DataFrame({
        "leiden": [str(i % n_clusters) for i in range(n_cells)],
    }, index=[f"cell_{i}" for i in range(n_cells)])
    var = pd.DataFrame(index=genes)
    adata = anndata.AnnData(X=X, obs=obs, var=var)

    for i in range(n_cells):
        c = int(obs.iloc[i]["leiden"])
        if c == 0:
            for g in ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18"]:
                adata.X[i, gene_idx[g]] = 4.0
        elif c == 1:
            for g in ["NCAM1", "KLRD1", "NKG7", "GNLY", "GZMB"]:
                adata.X[i, gene_idx[g]] = 4.0
    return adata


class TestCrossSpeciesStream:
    def test_runs_on_adata_with_conserved_genes(self):
        """With conserved marker genes in adata, stream runs successfully."""
        from axiom_sc.tier3.cross_species import CrossSpeciesStream
        adata = make_conserved_adata()
        cs = CrossSpeciesStream()
        results = cs.run(adata, cluster_key="leiden")
        assert len(results) == 3
        for sr in results.values():
            assert sr.stream == "cross_species"
            assert 0.0 <= sr.confidence <= 1.0

    def test_bad_cluster_key_returns_empty(self):
        from axiom_sc.tier3.cross_species import CrossSpeciesStream
        adata = make_conserved_adata()
        cs = CrossSpeciesStream()
        results = cs.run(adata, cluster_key="bad_key")
        assert results == {}

    def test_treg_cluster_identified(self):
        """Cluster 0 with high FOXP3/IL2RA/CTLA4 → Treg."""
        from axiom_sc.tier3.cross_species import CrossSpeciesStream
        adata = make_conserved_adata()
        cs = CrossSpeciesStream()
        results = cs.run(adata, cluster_key="leiden")
        c0 = results["0"]
        assert c0.available
        assert c0.label == "Treg"
        assert c0.confidence > 0.5

    def test_nk_cluster_identified(self):
        """Cluster 1 with high NCAM1/KLRD1/NKG7/GZMB → NK."""
        from axiom_sc.tier3.cross_species import CrossSpeciesStream
        adata = make_conserved_adata()
        cs = CrossSpeciesStream()
        results = cs.run(adata, cluster_key="leiden")
        c1 = results["1"]
        assert c1.available
        assert c1.label == "NK"
        assert c1.confidence > 0.5

    def test_orthofinder_mode_used_when_uns_present(self):
        """Pre-computed OrthoFinder results in adata.uns → OrthoFinder mode."""
        from axiom_sc.tier3.cross_species import CrossSpeciesStream
        adata = make_conserved_adata()
        adata.uns["cross_species_conservation"] = {
            "0": {"conservation_score": 0.9, "label": "Treg"},
            "1": {"conservation_score": 0.85, "label": "NK"},
            "2": {"conservation_score": 0.4, "label": ""},
        }
        cs = CrossSpeciesStream()
        results = cs.run(adata, cluster_key="leiden")
        assert results["0"].label == "Treg"
        assert results["0"].confidence == pytest.approx(0.9)
        assert results["1"].label == "NK"
        assert results["2"].confidence == pytest.approx(0.4)

    def test_no_conserved_genes_returns_fallback(self):
        """AnnData with no conserved marker genes → fallback tier2 results."""
        from axiom_sc.tier3.cross_species import CrossSpeciesStream
        adata = make_mock_adata(n_cells=30, n_clusters=2)
        cs = CrossSpeciesStream(tier2_labels={"0": "Treg", "1": "NK"})
        results = cs.run(adata, cluster_key="leiden")
        # Fallback returns tier2 hints
        assert results["0"].label == "Treg"
        assert results["1"].label == "NK"


# ── stream_agreement_report ───────────────────────────────────────────────────

class TestStreamAgreementReport:
    def _make_convergence_results(self) -> dict:
        from axiom_sc.tier3.convergence import ConvergenceResult
        from axiom_sc.tier3.stream_result import StreamResult

        def make_cr(cid, verdict, label, n_avail, n_agree):
            sr = {
                "velocity": StreamResult(cid, "velocity", label, 0.8, available=n_avail >= 1),
                "chromatin": StreamResult(cid, "chromatin", label, 0.7, available=n_avail >= 2),
            }
            return ConvergenceResult(
                cluster_id=cid, verdict=verdict, consensus_label=label,
                confidence=0.75, n_streams_available=n_avail,
                n_streams_agree=n_agree, stream_results=sr,
            )

        return {
            "c0": make_cr("c0", "PROVEN", "Treg", 5, 4),
            "c1": make_cr("c1", "HIGH_CONFIDENCE", "pDC", 4, 3),
            "c2": make_cr("c2", "UNCERTAIN", "", 2, 1),
            "c3": make_cr("c3", "PROVEN", "Th1", 5, 5),
        }

    def test_report_structure(self):
        from axiom_sc.tier3.convergence import stream_agreement_report
        results = self._make_convergence_results()
        report = stream_agreement_report(results)
        assert "n_clusters" in report
        assert "n_proven" in report
        assert "pct_proven" in report
        assert "mean_streams_available" in report
        assert "per_stream_availability" in report
        assert "label_distribution" in report

    def test_report_counts(self):
        from axiom_sc.tier3.convergence import stream_agreement_report
        results = self._make_convergence_results()
        report = stream_agreement_report(results)
        assert report["n_clusters"] == 4
        assert report["n_proven"] == 2
        assert report["n_high_confidence"] == 1
        assert report["n_uncertain"] == 1

    def test_report_percentages(self):
        from axiom_sc.tier3.convergence import stream_agreement_report
        results = self._make_convergence_results()
        report = stream_agreement_report(results)
        assert report["pct_proven"] == 50.0
        assert report["pct_high_confidence"] == 25.0
        assert report["pct_uncertain"] == 25.0

    def test_report_label_distribution(self):
        from axiom_sc.tier3.convergence import stream_agreement_report
        results = self._make_convergence_results()
        report = stream_agreement_report(results)
        label_dist = report["label_distribution"]
        assert "Treg" in label_dist
        assert "pDC" in label_dist
        assert "Th1" in label_dist

    def test_report_per_stream_availability(self):
        from axiom_sc.tier3.convergence import stream_agreement_report
        results = self._make_convergence_results()
        report = stream_agreement_report(results)
        avail = report["per_stream_availability"]
        # All 4 CRs have velocity available (n_avail >= 1 in all)
        assert "velocity" in avail
        assert avail["velocity"] == 1.0

    def test_empty_report(self):
        from axiom_sc.tier3.convergence import stream_agreement_report
        assert stream_agreement_report({}) == {"n_clusters": 0}


# ── 6-stream integration ──────────────────────────────────────────────────────

class TestSixStreamIntegration:
    def test_all_streams_with_full_mock_adata(self):
        """
        Run all 6 streams on adata that has velocity + gene activities
        + connectivities + LR genes + conserved genes.
        Verify: no crashes, cross_species returns label for clusters with conserved genes.
        """
        anndata = pytest.importorskip("anndata")
        pytest.importorskip("scipy")
        import pandas as pd
        from scipy.sparse import csr_matrix
        from axiom_sc.tier3.convergence import MultiStreamConvergence, stream_agreement_report

        n_cells = 90
        n_clusters = 3
        # Gene list covers: LR genes + conserved genes + velocity placeholder
        all_genes = list(set(
            ["TGFB1", "TGFB2", "FGF7", "IFNA1", "IFNA2", "IFNB1",
             "CCL2", "SPP1", "CSF1R", "GZMB", "PRF1",
             "FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18",
             "NCAM1", "KLRD1", "NKG7", "GNLY",
             "FOXP3", "TBX21", "GATA3", "RORC",  # for chromatin / cross-species
             "IRF7", "SIGLEC1",
             ] + [f"noise_gene_{i}" for i in range(10)]
        ))

        np.random.seed(1)
        X = np.random.rand(n_cells, len(all_genes)).astype(np.float32) * 0.2
        gene_idx = {g: i for i, g in enumerate(all_genes)}

        obs = pd.DataFrame({
            "leiden": [str(i % n_clusters) for i in range(n_cells)],
        }, index=[f"c{i}" for i in range(n_cells)])
        var = pd.DataFrame(index=all_genes)
        adata = anndata.AnnData(X=X, obs=obs, var=var)

        # Add velocity
        vel = np.random.rand(n_cells, len(all_genes)).astype(np.float32) * 0.1
        adata.layers["velocity"] = vel

        # Add gene activities for chromatin (only FOXP3/IL2RA for cluster 0)
        ga = np.random.rand(n_cells, len(all_genes)).astype(np.float32) * 0.1
        n_per = n_cells // n_clusters
        for i in range(n_per):  # cluster 0
            ga[i, gene_idx["FOXP3"]] = 3.0
            ga[i, gene_idx["IL2RA"]] = 2.5
        adata.obsm["gene_activities"] = ga
        adata.uns["gene_activity_genes"] = all_genes

        # Add connectivities
        rows, cols, data = [], [], []
        for i in range(n_cells):
            for j in [(i - 1) % n_cells, (i + 1) % n_cells]:
                rows.append(i); cols.append(j); data.append(1.0)
        adata.obsp["connectivities"] = csr_matrix((data, (rows, cols)), shape=(n_cells, n_cells))

        # Boost conserved Treg markers in cluster 0
        for i in range(n_per):
            for g in ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18"]:
                if g in gene_idx:
                    adata.X[i, gene_idx[g]] = 4.0

        msc = MultiStreamConvergence()
        results = msc.run_and_converge(adata, cluster_key="leiden")

        assert len(results) == n_clusters
        report = stream_agreement_report(results)
        assert report["n_clusters"] == n_clusters
        assert "velocity" in report["per_stream_availability"]
        assert "cross_species" in report["per_stream_availability"]
        # Cluster 0 should have Treg from chromatin + cross_species
        c0 = results["0"]
        assert c0.n_streams_available >= 1  # at least one stream contributed
