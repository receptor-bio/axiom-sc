# License: Apache 2.0

from axiom_sc.version import __version__
from axiom_sc.profiles.registry import ProfileType, PROFILE_COMPONENT_MAP
from axiom_sc.profiles.loader import load_profile

__all__ = [
    "AXIOMAnnotator", "__version__", "kg_version", "kg_rule_count",
    "list_kg_rules", "add_kg_rule", "deprecate_kg_rule",
    "list_profiles", "get_profile_components", "create_custom_profile",
]

# KG metadata — updated when kg_data/oracle_kg_v*.json is loaded
kg_version: str = "0.2.0"


def kg_rule_count() -> int:
    """Returns the number of ACTIVE rules in the bundled KG snapshot."""
    import json
    from pathlib import Path

    kg_path = Path(__file__).parent.parent / "kg_data" / "oracle_kg_v0.2.0.json"
    if not kg_path.exists():
        return 0
    with open(kg_path) as f:
        rules = json.load(f)
    if isinstance(rules, list):
        return sum(1 for r in rules if r.get("status") == "ACTIVE")
    return 0


def AXIOMAnnotator(
    profile: str | ProfileType = "oss-apache",
    profile_path: str | None = None,
    tiers: list[int] | None = None,
    **kwargs,
):
    """
    Main entry point for AXIOM-SC annotation.

    profile: "oss-mit" | "oss-apache" (default) | "commercial"
             OR path to a custom axiom_profile.json
    tiers:   subset of [1,2,3,4,5] to run. None = all tiers in profile.

    Environment override: AXIOM_PROFILE=oss-apache
    """
    from axiom_sc.annotator import Annotator

    resolved_profile, components = load_profile(profile, profile_path)
    return Annotator(
        profile=resolved_profile,
        components=components,
        tiers=tiers or [1, 2, 3, 4, 5],
        **kwargs,
    )


# ─── KG API ──────────────────────────────────────────────────────────────────

def _load_kg_rules() -> list[dict]:
    import json
    from pathlib import Path
    kg_path = Path(__file__).parent.parent / "kg_data" / "oracle_kg_v0.2.0.json"
    if not kg_path.exists():
        return []
    with open(kg_path) as f:
        return json.load(f)


def _save_kg_rules(rules: list[dict]) -> None:
    import json
    from pathlib import Path
    kg_path = Path(__file__).parent.parent / "kg_data" / "oracle_kg_v0.2.0.json"
    with open(kg_path, "w") as f:
        json.dump(rules, f, indent=2)


def list_kg_rules(
    cell_type: str | None = None,
    tissue: str | None = None,
    rule_type: str | None = None,
    status: str = "ACTIVE",
) -> list[dict]:
    """Returns KG rules matching the given filters."""
    rules = _load_kg_rules()
    if status:
        rules = [r for r in rules if r.get("status") == status]
    if cell_type:
        rules = [r for r in rules if r.get("cell_type") == cell_type]
    if tissue:
        rules = [r for r in rules
                 if tissue in (r.get("tissue_context") or [])]
    if rule_type:
        rules = [r for r in rules if r.get("rule_type") == rule_type]
    return rules


def add_kg_rule(rule: dict) -> dict:
    """
    Validates and appends a new rule to the KG JSON.
    Status is set to PENDING_REVIEW regardless of what's passed.
    Raises ValueError if required fields are missing.
    """
    import json
    from pathlib import Path
    import jsonschema

    schema_path = Path(__file__).parent / "kg" / "schema.json"
    with open(schema_path) as f:
        schema = json.load(f)

    rule["status"] = "PENDING_REVIEW"
    jsonschema.validate(rule, schema)

    rules = _load_kg_rules()
    if any(r["rule_id"] == rule["rule_id"] for r in rules):
        raise ValueError(f"Duplicate rule_id: {rule['rule_id']}")
    rules.append(rule)
    _save_kg_rules(rules)
    return rule


def deprecate_kg_rule(rule_id: str) -> None:
    """Soft-deletes a KG rule by setting status to DEPRECATED."""
    rules = _load_kg_rules()
    found = False
    for r in rules:
        if r["rule_id"] == rule_id:
            r["status"] = "DEPRECATED"
            found = True
            break
    if not found:
        raise ValueError(f"Rule not found: {rule_id}")
    _save_kg_rules(rules)


# ─── Profile API ─────────────────────────────────────────────────────────────

def list_profiles() -> list[dict]:
    """Returns metadata for all built-in profiles."""
    from axiom_sc.profiles.registry import COMPONENT_REGISTRY, PROFILE_COMPONENT_MAP, ProfileType

    result = []
    for pt in ProfileType:
        enabled = PROFILE_COMPONENT_MAP[pt]
        result.append({
            "name": pt.value,
            "display_name": pt.value.replace("-", " ").title(),
            "component_count": len(enabled),
            "tiers_available": sorted({COMPONENT_REGISTRY[c].tier
                                       for c in enabled if c in COMPONENT_REGISTRY}),
            "commercial_ok": all(COMPONENT_REGISTRY[c].commercial_ok
                                 for c in enabled if c in COMPONENT_REGISTRY),
        })
    return result


def get_profile_components(name: str) -> list[dict]:
    """
    Returns detailed component list for a named profile.
    name: "oss-mit" | "oss-apache" | "commercial"
    """
    from axiom_sc.profiles.loader import _load_builtin_profile
    from axiom_sc.profiles.registry import COMPONENT_REGISTRY, PROFILE_COMPONENT_MAP, ProfileType

    try:
        pt = ProfileType(name)
    except ValueError:
        raise ValueError(f"Unknown profile: {name}. Valid: {[p.value for p in ProfileType]}")

    enabled = PROFILE_COMPONENT_MAP[pt]
    result = []
    for comp_id, meta in COMPONENT_REGISTRY.items():
        result.append({
            "component_id": comp_id,
            "display_name": meta.display_name,
            "license": meta.license,
            "commercial_ok": meta.commercial_ok,
            "tier": meta.tier,
            "enabled": comp_id in enabled,
            "bundled": meta.bundled,
            "notes": meta.notes,
            "paper_doi": meta.paper_doi,
        })
    return sorted(result, key=lambda x: (x["tier"], x["component_id"]))


def create_custom_profile(profile: dict) -> dict:
    """
    Validates a custom profile dict and returns it with validation_warnings.
    Does NOT persist to disk — caller saves the JSON themselves.
    """
    from axiom_sc.profiles.registry import COMPONENT_REGISTRY, PROFILE_COMPONENT_MAP, ProfileType
    from axiom_sc.profiles.loader import _load_builtin_profile

    warnings: list[str] = []

    base_name = profile.get("base_profile", "oss-apache")
    try:
        pt = ProfileType(base_name)
    except ValueError:
        raise ValueError(f"Unknown base_profile: {base_name}")

    additional = set(profile.get("additional_components", []))
    disabled = set(profile.get("disabled_components", []))

    unknown = (additional | disabled) - set(COMPONENT_REGISTRY.keys())
    if unknown:
        raise ValueError(f"Unknown component IDs: {unknown}")

    base_components = PROFILE_COMPONENT_MAP[pt]
    enabled = (base_components | additional) - disabled

    non_commercial = [c for c in enabled
                      if c in COMPONENT_REGISTRY and not COMPONENT_REGISTRY[c].commercial_ok]
    if non_commercial:
        warnings.append(
            f"Non-commercial components enabled: {non_commercial}. "
            "Commercial use would require separate licensing."
        )

    missing_ack = [c for c in enabled
                   if c in COMPONENT_REGISTRY
                   and COMPONENT_REGISTRY[c].license in ("GPL-3.0",)
                   and c not in profile.get("license_acknowledgements", {})]
    if missing_ack:
        warnings.append(
            f"GPL components missing license_acknowledgements: {missing_ack}"
        )

    profile_id = profile.get("name", "custom-" + str(hash(str(sorted(enabled))))[:8])
    return {
        "profile_id": profile_id,
        "enabled_components": sorted(enabled),
        "component_count": len(enabled),
        "validation_warnings": warnings,
    }
