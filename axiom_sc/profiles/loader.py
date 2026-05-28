# License: Apache 2.0

import json
import os
from pathlib import Path
from typing import Optional

from axiom_sc.profiles.registry import (
    COMPONENT_REGISTRY,
    PROFILE_COMPONENT_MAP,
    ProfileType,
)


def load_profile(
    profile: str | ProfileType = "oss-apache",
    profile_path: Optional[str] = None,
) -> tuple[ProfileType | str, set[str]]:
    """
    Returns (profile_name, enabled_component_ids).

    If profile_path is given, loads a custom axiom_profile.json.
    Otherwise resolves profile string → built-in ProfileType.
    Environment variable AXIOM_PROFILE overrides profile argument.
    """
    env_override = os.environ.get("AXIOM_PROFILE")
    if env_override:
        profile = env_override

    if profile_path is not None:
        return _load_custom_profile(profile_path)

    # Try to parse as a file path first
    if isinstance(profile, str) and profile.endswith(".json"):
        return _load_custom_profile(profile)

    return _load_builtin_profile(profile)


def _load_builtin_profile(profile: str | ProfileType) -> tuple[ProfileType, set[str]]:
    if isinstance(profile, str):
        try:
            profile = ProfileType(profile)
        except ValueError:
            raise ValueError(
                f"Unknown profile '{profile}'. "
                f"Valid built-in profiles: {[p.value for p in ProfileType]}"
            )
    components = PROFILE_COMPONENT_MAP[profile]
    return profile, components


def _load_custom_profile(path: str) -> tuple[str, set[str]]:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Custom profile not found: {path}")

    with open(p) as f:
        data = json.load(f)

    base_name = data.get("base_profile", "oss-apache")
    _, base_components = _load_builtin_profile(base_name)

    additional = set(data.get("additional_components", []))
    disabled = set(data.get("disabled_components", []))

    unknown = (additional | disabled) - set(COMPONENT_REGISTRY.keys())
    if unknown:
        raise ValueError(f"Unknown component IDs in custom profile: {unknown}")

    components = (base_components | additional) - disabled
    profile_name = data.get("name", p.stem)
    return profile_name, components


def load_profile_from_env() -> tuple[ProfileType | str, set[str]]:
    """Loads profile from AXIOM_PROFILE env var, defaulting to oss-apache."""
    profile_name = os.environ.get("AXIOM_PROFILE", "oss-apache")
    profile_path = os.environ.get("AXIOM_PROFILE_PATH")
    return load_profile(profile_name, profile_path)
