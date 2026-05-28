# License: Apache 2.0

from axiom_sc.profiles.registry import ProfileType


class Annotator:
    """
    Core AXIOM-SC annotator. Instantiated via AXIOMAnnotator() factory.
    Routes cells through enabled tiers in order.
    """

    def __init__(
        self,
        profile: ProfileType | str,
        components: set[str],
        tiers: list[int],
        **kwargs,
    ):
        self.profile = profile
        self.components = components
        self.tiers = sorted(tiers)

    def annotate(self, adata_or_path, progress_callback=None):
        raise NotImplementedError("Annotation pipeline implemented in Phase 2 tiers.")

    def __repr__(self) -> str:
        return (
            f"Annotator(profile={self.profile!r}, "
            f"tiers={self.tiers}, "
            f"n_components={len(self.components)})"
        )
