"""
Profile smoke tests — Day 1 requirement.
Tests: built-in profile loading, component inclusion/exclusion, custom profile.
"""
import json
import tempfile
from pathlib import Path

import pytest

from axiom_sc.profiles.registry import (
    COMPONENT_REGISTRY,
    OSS_MIT_COMPONENTS,
    OSS_APACHE_COMPONENTS,
    COMMERCIAL_COMPONENTS,
    ProfileType,
    PROFILE_COMPONENT_MAP,
)
from axiom_sc.profiles.loader import load_profile
from axiom_sc.profiles.validator import validate_components, ProfileValidationError


# ── Test 1: oss-apache profile includes tier2_axiom_kg ────────────────────────

def test_oss_apache_includes_axiom_kg():
    profile, components = load_profile("oss-apache")
    assert profile == ProfileType.OSS_APACHE
    assert "tier2_axiom_kg" in components, (
        "oss-apache must include tier2_axiom_kg (Apache-2.0, commercial_ok=True)"
    )


def test_oss_apache_excludes_gpl():
    _, components = load_profile("oss-apache")
    gpl_components = {
        k for k, v in COMPONENT_REGISTRY.items() if v.license == "GPL-3.0"
    }
    overlap = components & gpl_components
    assert not overlap, (
        f"oss-apache must not include GPL-3.0 components: {overlap}"
    )


def test_oss_apache_excludes_panglao():
    _, components = load_profile("oss-apache")
    assert "kg_panglao" not in components, (
        "oss-apache must not include PanglaoDB (CC BY-NC, non-commercial)"
    )


# ── Test 2: oss-mit profile only includes MIT + Apache-2.0 components ─────────

def test_oss_mit_only_permissive():
    profile, components = load_profile("oss-mit")
    assert profile == ProfileType.OSS_MIT
    for cid in components:
        meta = COMPONENT_REGISTRY[cid]
        assert meta.license in ("MIT", "Apache-2.0"), (
            f"oss-mit must only include MIT/Apache-2.0 components, "
            f"but found {cid} with license {meta.license}"
        )
        assert meta.commercial_ok, (
            f"oss-mit must only include commercial_ok components, "
            f"but {cid} has commercial_ok=False"
        )


def test_oss_mit_includes_census():
    _, components = load_profile("oss-mit")
    assert "tier1_census" in components, (
        "oss-mit must include tier1_census (MIT license)"
    )


def test_oss_mit_excludes_bsd_scvelo():
    # scVelo is BSD-3-Clause, not MIT or Apache — excluded from oss-mit
    _, components = load_profile("oss-mit")
    assert "tier3_scvelo" not in components, (
        "oss-mit must exclude tier3_scvelo (BSD-3-Clause, not MIT/Apache-2.0)"
    )


# ── Test 3: custom profile validation ────────────────────────────────────────

def test_custom_profile_with_panglao():
    """Custom profile can add kg_panglao for academic use."""
    custom = {
        "name": "academic-test",
        "base_profile": "oss-apache",
        "additional_components": ["kg_panglao"],
        "disabled_components": [],
        "license_acknowledgements": {
            "kg_panglao": "Academic non-commercial research only per CC BY-NC 4.0"
        },
    }
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".json", delete=False
    ) as f:
        json.dump(custom, f)
        tmp_path = f.name

    profile_name, components = load_profile(profile_path=tmp_path)
    assert profile_name == "academic-test"
    assert "kg_panglao" in components, "Custom profile must include kg_panglao when added"
    assert "tier2_axiom_kg" in components, "Custom profile must inherit base oss-apache components"

    Path(tmp_path).unlink()


def test_custom_profile_disable_openai():
    """Custom profile can disable Tier 4 OpenAI."""
    custom = {
        "name": "no-openai",
        "base_profile": "oss-apache",
        "additional_components": [],
        "disabled_components": ["tier4_openai"],
    }
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".json", delete=False
    ) as f:
        json.dump(custom, f)
        tmp_path = f.name

    _, components = load_profile(profile_path=tmp_path)
    assert "tier4_openai" not in components
    Path(tmp_path).unlink()


def test_custom_profile_unknown_component_raises():
    """Custom profile with an unknown component ID must raise ValueError."""
    custom = {
        "name": "bad-profile",
        "base_profile": "oss-apache",
        "additional_components": ["nonexistent_component_xyz"],
        "disabled_components": [],
    }
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".json", delete=False
    ) as f:
        json.dump(custom, f)
        tmp_path = f.name

    with pytest.raises(ValueError, match="Unknown component IDs"):
        load_profile(profile_path=tmp_path)

    Path(tmp_path).unlink()


# ── Test 4: validator warnings ────────────────────────────────────────────────

def test_validator_warns_non_commercial():
    components_with_nc = OSS_APACHE_COMPONENTS | {"kg_panglao"}
    warnings = validate_components(components_with_nc)
    assert any("Non-commercial" in w for w in warnings)


def test_validator_clean_for_oss_apache():
    _, components = load_profile("oss-apache")
    warnings = validate_components(components)
    # oss-apache should produce no non-commercial warnings
    nc_warnings = [w for w in warnings if "Non-commercial" in w]
    assert not nc_warnings, f"oss-apache should have no NC warnings, got: {nc_warnings}"


# ── Test 5: commercial profile covers all non-panglao components ──────────────

def test_commercial_excludes_panglao():
    _, components = load_profile("commercial")
    assert "kg_panglao" not in components


def test_commercial_includes_scenic_subprocess():
    _, components = load_profile("commercial")
    # pySCENIC subprocess is allowed (GPL isolated, not imported) even commercially
    assert "tier2_scenic" in components


# ── Test 6: unknown profile string raises ────────────────────────────────────

def test_unknown_profile_raises():
    with pytest.raises(ValueError, match="Unknown profile"):
        load_profile("bogus-profile")


# ── Test 7: AXIOMAnnotator factory + kg_rule_count (covers __init__.py) ──────

def test_axiom_annotator_factory_instantiation():
    """AXIOMAnnotator() factory returns an Annotator with correct profile/tiers."""
    import axiom_sc
    annotator = axiom_sc.AXIOMAnnotator(profile="oss-apache", tiers=[1, 2])
    assert repr(annotator).startswith("Annotator(")
    assert "oss-apache" in repr(annotator)


def test_axiom_annotator_tiers_subset():
    """AXIOMAnnotator with tiers=[2] stores only tier 2."""
    import axiom_sc
    annotator = axiom_sc.AXIOMAnnotator(profile="oss-apache", tiers=[2])
    assert annotator.tiers == [2]


def test_kg_rule_count_returns_int():
    """kg_rule_count() returns the number of ACTIVE rules in the bundled KG."""
    import axiom_sc
    count = axiom_sc.kg_rule_count()
    assert isinstance(count, int)
    assert count >= 220, f"Expected ≥ 220 ACTIVE rules after Day 4, got {count}"


# ── Test 8: check_tier_coverage covers validator.py lines 59-64 ──────────────

def test_check_tier_coverage_returns_dict():
    """check_tier_coverage maps tier numbers to component IDs."""
    from axiom_sc.profiles.validator import check_tier_coverage
    coverage = check_tier_coverage(OSS_APACHE_COMPONENTS)
    assert isinstance(coverage, dict)
    assert set(coverage.keys()) == {1, 2, 3, 4, 5}
    # Tier 2 must include the axiom KG engine
    assert "tier2_axiom_kg" in coverage[2]


def test_validate_components_raises_on_unknown():
    """validate_components raises ProfileValidationError for unknown component IDs."""
    with pytest.raises(ProfileValidationError):
        validate_components({"completely_unknown_component_xyz"})
