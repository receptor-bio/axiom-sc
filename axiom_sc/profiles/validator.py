# License: Apache 2.0

from axiom_sc.profiles.registry import COMPONENT_REGISTRY, ComponentMeta


class ProfileValidationError(Exception):
    pass


def validate_components(component_ids: set[str]) -> list[str]:
    """
    Validates a set of component IDs.
    Returns a list of warning strings (empty = clean).
    Raises ProfileValidationError on hard violations.
    """
    warnings: list[str] = []

    unknown = component_ids - set(COMPONENT_REGISTRY.keys())
    if unknown:
        raise ProfileValidationError(f"Unknown component IDs: {unknown}")

    non_commercial = [
        cid for cid in component_ids
        if not COMPONENT_REGISTRY[cid].commercial_ok
    ]
    if non_commercial:
        names = [COMPONENT_REGISTRY[c].display_name for c in non_commercial]
        warnings.append(
            f"Non-commercial components enabled ({len(non_commercial)}): "
            + ", ".join(names)
        )

    gpl = [
        cid for cid in component_ids
        if COMPONENT_REGISTRY[cid].license == "GPL-3.0"
    ]
    if gpl:
        names = [COMPONENT_REGISTRY[c].display_name for c in gpl]
        warnings.append(
            f"GPL-3.0 components enabled (subprocess isolation required): "
            + ", ".join(names)
        )

    nc = [
        cid for cid in component_ids
        if "NC" in COMPONENT_REGISTRY[cid].license
    ]
    if nc:
        names = [COMPONENT_REGISTRY[c].display_name for c in nc]
        warnings.append(
            f"Non-commercial-licensed data sources enabled: " + ", ".join(names)
        )

    return warnings


def check_tier_coverage(component_ids: set[str]) -> dict[int, list[str]]:
    """Returns a mapping of tier → list of enabled component IDs for that tier."""
    coverage: dict[int, list[str]] = {1: [], 2: [], 3: [], 4: [], 5: []}
    for cid in component_ids:
        meta = COMPONENT_REGISTRY[cid]
        if meta.tier in coverage:
            coverage[meta.tier].append(cid)
    return coverage
