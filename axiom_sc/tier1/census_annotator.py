"""
CELLxGENE Census kNN fallback annotator for Tier 1.

Used as the primary Tier 1 method while AXIOMTier1 weights are training,
and as a fallback / sanity check alongside the MLP ensemble.

References:
  CELLxGENE Census: Chanzuckerberg Initiative (2023)
    doi: 10.1126/science.abl4896
  Tabula Sapiens: Tabula Sapiens Consortium (2022) Science 376:eabl4896
    doi: 10.1126/science.abl4896

License note:
  cellxgene-census Python API: MIT license
  Underlying data: CC BY 4.0 per dataset

License: Apache 2.0
"""
from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    import pandas as pd

logger = logging.getLogger(__name__)

# Confidence thresholds — consistent with AXIOMTier1Ensemble
_ACCEPT_THRESH = 0.85
_VERIFY_THRESH = 0.50

# Census organism and measurement name
_CENSUS_ORGANISM = "Homo sapiens"
_CENSUS_MEASUREMENT = "RNA"

# Number of nearest neighbours for kNN annotation
_N_NEIGHBORS = 50

# Minimum vote fraction to assign a label (otherwise → 'tier2_full')
_MIN_VOTE_FRACTION = 0.30


class CensusAnnotator:
    """
    Nearest-neighbour cell type annotator using CELLxGENE Census embeddings.

    Retrieves a reference set of cells from Census, builds a kNN index over
    their embeddings, and annotates query cells by majority vote over the k
    nearest Census neighbours.

    Parameters
    ----------
    tissue_filter : list[str] | None
        Limit Census reference to specific tissues (e.g. ["blood", "lung"]).
        None = all tissues (slower but broader coverage).
    n_neighbors : int
        k for kNN search (default 50).
    embedding_name : str
        Which Census embedding to use.
        Default "scvi" — the Census-wide scVI latent embedding.
    use_existing_obsm : str | None
        If the query adata already has an embedding (e.g. adata.obsm["X_scvi"]),
        use it directly instead of recomputing. Key name to look up in .obsm.
    """

    def __init__(
        self,
        tissue_filter: list[str] | None = None,
        n_neighbors: int = _N_NEIGHBORS,
        embedding_name: str = "scvi",
        use_existing_obsm: str | None = None,
    ):
        self.tissue_filter = tissue_filter
        self.n_neighbors = n_neighbors
        self.embedding_name = embedding_name
        self.use_existing_obsm = use_existing_obsm

        self._ref_embeddings: np.ndarray | None = None
        self._ref_labels: list[str] | None = None
        self._index = None   # annoy / sklearn NearestNeighbors index

    # ── build reference ───────────────────────────────────────────────────────

    def build_reference(
        self,
        max_cells_per_type: int = 2000,
        cell_type_column: str = "cell_type",
        census_version: str = "stable",
    ) -> None:
        """
        Download reference embeddings from CELLxGENE Census and build kNN index.

        Parameters
        ----------
        max_cells_per_type : int
            Maximum cells per cell type to include (stratified cap).
        cell_type_column : str
            obs column name for cell type labels.
        census_version : str
            Census version tag ("stable" or a date string like "2023-12-15").
        """
        try:
            import cellxgene_census
            import tiledbsoma
        except ImportError:
            raise ImportError(
                "cellxgene-census and tiledbsoma required: "
                "pip install 'axiom-sc[science]'"
            )

        logger.info("CensusAnnotator: opening Census version=%s …", census_version)

        with cellxgene_census.open_soma(census_version=census_version) as census:
            obs_value_filter = f'is_primary_data == True'
            if self.tissue_filter:
                tissue_str = " or ".join(
                    f'tissue_general == "{t}"' for t in self.tissue_filter
                )
                obs_value_filter += f" and ({tissue_str})"

            # Fetch embeddings
            logger.info(
                "CensusAnnotator: fetching %s embeddings …", self.embedding_name
            )
            emb_uri = (
                census["census_data"][_CENSUS_ORGANISM]
                .ms[_CENSUS_MEASUREMENT]
                .obsm[self.embedding_name]
            )

            # Get obs metadata for cell type labels
            obs_df = (
                census["census_data"][_CENSUS_ORGANISM]
                .obs.read(
                    value_filter=obs_value_filter,
                    column_names=["soma_joinid", cell_type_column],
                )
                .concat()
                .to_pandas()
            )

            # Stratified subsample
            obs_df = self._stratified_sample(obs_df, cell_type_column, max_cells_per_type)

            # Fetch the embeddings for the selected joinids
            join_ids = obs_df["soma_joinid"].tolist()
            embeddings = emb_uri.read(coords=(join_ids,)).concat().to_pandas()

            # Align
            emb_matrix = embeddings.drop(columns=["soma_joinid"]).values.astype(np.float32)
            labels = obs_df[cell_type_column].tolist()

        self._ref_embeddings = emb_matrix
        self._ref_labels = labels
        self._build_knn_index(emb_matrix)

        logger.info(
            "CensusAnnotator: reference built — %d cells, %d cell types",
            len(labels), len(set(labels)),
        )

    @staticmethod
    def _stratified_sample(df, col: str, max_per_type: int):
        """Return a stratified subsample with at most max_per_type rows per group."""
        import pandas as pd
        return (
            df.groupby(col, group_keys=False)
            .apply(lambda g: g.sample(min(len(g), max_per_type), random_state=42))
            .reset_index(drop=True)
        )

    def _build_knn_index(self, embeddings: np.ndarray) -> None:
        """Build sklearn NearestNeighbors index over the reference embeddings."""
        try:
            from sklearn.neighbors import NearestNeighbors
        except ImportError:
            raise ImportError("scikit-learn required: pip install 'axiom-sc[science]'")

        logger.info("CensusAnnotator: building kNN index (k=%d) …", self.n_neighbors)
        self._index = NearestNeighbors(
            n_neighbors=self.n_neighbors,
            metric="euclidean",
            algorithm="ball_tree",
            n_jobs=-1,
        )
        self._index.fit(embeddings)

    # ── load / save reference ─────────────────────────────────────────────────

    def save_reference(self, path: str) -> None:
        """Save built reference embeddings and labels for reuse (no re-download needed)."""
        import pickle
        with open(path, "wb") as f:
            pickle.dump({
                "embeddings": self._ref_embeddings,
                "labels": self._ref_labels,
            }, f)
        logger.info("CensusAnnotator: saved reference → %s", path)

    def load_reference(self, path: str) -> None:
        """Load pre-built reference and rebuild kNN index (fast, no Census needed)."""
        import pickle
        with open(path, "rb") as f:
            data = pickle.load(f)
        self._ref_embeddings = data["embeddings"]
        self._ref_labels = data["labels"]
        self._build_knn_index(self._ref_embeddings)
        logger.info(
            "CensusAnnotator: loaded reference from %s (%d cells)", path,
            len(self._ref_labels),
        )

    # ── inference ─────────────────────────────────────────────────────────────

    def predict(self, adata) -> dict:
        """
        Annotate cells in adata via kNN majority vote over Census reference.

        Parameters
        ----------
        adata : anndata.AnnData
            Must have .obsm[self.use_existing_obsm] if use_existing_obsm is set,
            OR expression matrix compatible with the Census embedding space.

        Returns
        -------
        dict with keys:
            'label'      : pd.Series[str]
            'confidence' : pd.Series[float]  — vote fraction for winning label
            'route'      : pd.Series[str]    — 'accept' | 'tier2_verify' | 'tier2_full'
        """
        import pandas as pd

        if self._index is None:
            raise RuntimeError(
                "CensusAnnotator: reference not built. "
                "Call build_reference() or load_reference() first."
            )

        query_emb = self._get_query_embedding(adata)
        _, indices = self._index.kneighbors(query_emb)   # (n_cells, k)

        labels = []
        confidences = []

        for cell_nbrs in indices:
            nbr_labels = [self._ref_labels[i] for i in cell_nbrs]
            # Majority vote
            from collections import Counter
            counts = Counter(nbr_labels)
            top_label, top_count = counts.most_common(1)[0]
            vote_frac = top_count / len(cell_nbrs)

            if vote_frac < _MIN_VOTE_FRACTION:
                labels.append("Unresolved")
                confidences.append(0.0)
            else:
                labels.append(top_label)
                confidences.append(vote_frac)

        confidence_arr = np.array(confidences)
        obs_names = adata.obs_names

        return {
            "label":      pd.Series(labels, index=obs_names, name="tier1_label"),
            "confidence": pd.Series(confidence_arr.tolist(), index=obs_names, name="tier1_confidence"),
            "route":      pd.Series(_route_from_confidence(confidence_arr), index=obs_names, name="tier1_route"),
        }

    def _get_query_embedding(self, adata) -> np.ndarray:
        """Extract query embedding from adata or raise informative error."""
        if self.use_existing_obsm:
            if self.use_existing_obsm not in adata.obsm:
                raise ValueError(
                    f"CensusAnnotator: use_existing_obsm={self.use_existing_obsm!r} "
                    f"not in adata.obsm. Available: {list(adata.obsm.keys())}"
                )
            return adata.obsm[self.use_existing_obsm].astype(np.float32)
        raise RuntimeError(
            "CensusAnnotator: no embedding available. "
            "Set use_existing_obsm to an obsm key, "
            "or call build_reference() with matching Census version."
        )

    def __repr__(self) -> str:
        n_ref = len(self._ref_labels) if self._ref_labels else 0
        return (
            f"CensusAnnotator(k={self.n_neighbors}, "
            f"embedding={self.embedding_name!r}, "
            f"ref_cells={n_ref})"
        )


def _route_from_confidence(confidence: np.ndarray) -> list[str]:
    routes = []
    for c in confidence:
        if c >= _ACCEPT_THRESH:
            routes.append("accept")
        elif c >= _VERIFY_THRESH:
            routes.append("tier2_verify")
        else:
            routes.append("tier2_full")
    return routes
