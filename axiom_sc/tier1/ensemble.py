"""
Tier 1 ensemble: fuses AXIOMTier1 MLP predictions with Census kNN fallback.

Fusion strategy
---------------
- Primary: AXIOMTier1Ensemble (MLP, Apache 2.0 weights)
- Fallback: CensusAnnotator (kNN, MIT API, CC BY 4.0 data)

When MLP weights are available: MLP is primary, Census used for sanity check.
When MLP weights unavailable (e.g. during initial training): Census is primary.

Conflict handling:
  MLP and Census agree   → use MLP confidence
  MLP and Census disagree on high-confidence calls → route to tier2_verify
  MLP unavailable        → use Census directly

License: Apache 2.0
"""
from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    import pandas as pd

logger = logging.getLogger(__name__)

# Agreement threshold: if MLP top-1 label matches Census top-1, treat as confirmed
_AGREEMENT_CONFIDENCE_BOOST = 0.05   # small boost when both agree


class Tier1Ensemble:
    """
    Fuses AXIOMTier1 MLP and Census kNN fallback for robust Tier 1 annotation.

    Parameters
    ----------
    mlp : AXIOMTier1Ensemble | None
        MLP ensemble. If None, Census-only mode is used.
    census : CensusAnnotator | None
        Census kNN annotator. If None, MLP-only mode is used.
    """

    def __init__(self, mlp=None, census=None):
        if mlp is None and census is None:
            raise ValueError("At least one of mlp or census must be provided.")
        self.mlp = mlp
        self.census = census

    def predict(self, adata) -> dict:
        """
        Predict cell types using available annotators.

        Returns the same dict schema as AXIOMTier1Ensemble.predict().
        """
        import pandas as pd

        if self.mlp is not None and self.census is not None:
            return self._fused_predict(adata)
        elif self.mlp is not None:
            return self.mlp.predict(adata)
        else:
            return self.census.predict(adata)

    def _fused_predict(self, adata) -> dict:
        """Run both annotators and fuse results."""
        import pandas as pd

        mlp_result = self.mlp.predict(adata)
        census_result = self.census.predict(adata)

        mlp_labels = mlp_result["label"]
        census_labels = census_result["label"]
        mlp_conf = mlp_result["confidence"].values
        census_conf = census_result["confidence"].values

        fused_labels = []
        fused_conf = []
        fused_routes = []

        for i, obs_name in enumerate(adata.obs_names):
            ml = mlp_labels.iloc[i]
            cl = census_labels.iloc[i]
            mc = mlp_conf[i]
            cc = census_conf[i]

            if ml == cl:
                # Both agree — boost confidence slightly
                conf = min(1.0, mc + _AGREEMENT_CONFIDENCE_BOOST)
                fused_labels.append(ml)
                fused_conf.append(conf)
            elif mc >= 0.85 and cc < 0.50:
                # MLP is confident, Census uncertain → trust MLP
                fused_labels.append(ml)
                fused_conf.append(mc)
            elif cc >= 0.85 and mc < 0.50:
                # Census is confident, MLP uncertain → use Census but downgrade
                fused_labels.append(cl)
                fused_conf.append(min(cc, 0.80))   # cap at verify threshold
            else:
                # Conflict between two reasonable predictions → demote to verify
                fused_labels.append(ml)            # MLP label as candidate
                fused_conf.append(min(mc, 0.70))   # force into tier2_verify

            # Routing
            c = fused_conf[-1]
            if c >= 0.85:
                fused_routes.append("accept")
            elif c >= 0.50:
                fused_routes.append("tier2_verify")
            else:
                fused_routes.append("tier2_full")

        obs_names = adata.obs_names
        return {
            "label":      pd.Series(fused_labels, index=obs_names, name="tier1_label"),
            "confidence": pd.Series(fused_conf, index=obs_names, name="tier1_confidence"),
            "route":      pd.Series(fused_routes, index=obs_names, name="tier1_route"),
            "probs":      mlp_result.get("probs"),   # MLP softmax if available
        }

    def acceptance_rate(self, predict_result: dict) -> float:
        """Fraction of cells routed as 'accept' (Tier 1 terminal)."""
        routes = predict_result["route"]
        return (routes == "accept").mean()

    def __repr__(self) -> str:
        parts = []
        if self.mlp is not None:
            parts.append(f"mlp={self.mlp!r}")
        if self.census is not None:
            parts.append(f"census={self.census!r}")
        return f"Tier1Ensemble({', '.join(parts)})"
