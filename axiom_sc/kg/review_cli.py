"""
Interactive KG rule review tool.

Walks a researcher through PENDING_REVIEW candidate rules, allowing:
  - Accept (promote to ACTIVE) after editing mechanistic_basis + pmid
  - Skip (leave as PENDING_REVIEW)
  - Reject (set DEPRECATED)
  - Batch approve a list of rule IDs

Usage::
    # Interactive mode (researcher reviews one at a time)
    axiom-kg review --candidates kg_data/candidates/cellmarker2_candidates.json \
                    --kg-out kg_data/oracle_kg_v0.2.0.json

    # Batch mode (approve pre-validated list)
    from axiom_sc.kg.review_cli import batch_approve
    approved = batch_approve(candidates_path, approved_ids, kg_path)

License: Apache 2.0
"""
from __future__ import annotations

import json
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional

from axiom_sc.kg.seeder import load_candidates, save_candidates


class RuleReviewCLI:
    """
    Interactive terminal-based rule review tool.
    Loads PENDING_REVIEW candidates and walks the researcher through each one.
    """

    def __init__(
        self,
        candidates_path: str,
        kg_out_path: str,
        kg_version: str = "0.2.0",
    ):
        self.candidates_path = Path(candidates_path)
        self.kg_out_path = Path(kg_out_path)
        self.kg_version = kg_version

    def run_interactive(
        self,
        cell_type_filter: str | None = None,
        max_rules: int | None = None,
    ) -> dict[str, int]:
        """
        Run interactive review session.
        Returns stats: {accepted, skipped, rejected}.
        """
        metadata, candidates = load_candidates(str(self.candidates_path))
        pending = [
            c for c in candidates
            if c.get("status") == "PENDING_REVIEW"
            and (cell_type_filter is None or c.get("cell_type") == cell_type_filter)
        ]

        if max_rules:
            pending = pending[:max_rules]

        print(f"\n{'='*60}")
        print(f"AXIOM-SC KG Review — {len(pending)} rules to review")
        print(f"Candidates: {self.candidates_path}")
        print(f"Output KG:  {self.kg_out_path}")
        print("="*60)
        print("Commands: [a]ccept  [s]kip  [r]eject  [q]uit")
        print()

        stats = {"accepted": 0, "skipped": 0, "rejected": 0}
        accepted_rules: list[dict] = []

        for i, rule in enumerate(pending):
            self._print_rule(rule, i + 1, len(pending))
            action = self._prompt_action(rule)

            if action == "a":
                edited = self._edit_rule(rule)
                edited["status"] = "ACTIVE"
                edited["added_in_version"] = self.kg_version
                accepted_rules.append(edited)
                # Update in candidates list
                rule["status"] = "ACTIVE"
                stats["accepted"] += 1
                print(f"  ✓ ACCEPTED: {rule['rule_id']}")

            elif action == "r":
                rule["status"] = "DEPRECATED"
                stats["rejected"] += 1
                print(f"  ✗ REJECTED: {rule['rule_id']}")

            elif action == "q":
                print("\nSession ended by user.")
                break

            else:
                stats["skipped"] += 1
                print(f"  → SKIPPED: {rule['rule_id']}")

        # Append accepted rules to KG
        if accepted_rules:
            self._append_to_kg(accepted_rules)

        # Save updated candidates (with status changes)
        save_candidates(metadata, candidates, str(self.candidates_path))

        print(f"\n{'='*60}")
        print(f"Session complete: {stats['accepted']} accepted, "
              f"{stats['skipped']} skipped, {stats['rejected']} rejected")
        print(f"KG updated: {self.kg_out_path}")
        return stats

    # ── helpers ───────────────────────────────────────────────────────────────

    def _print_rule(self, rule: dict, n: int, total: int) -> None:
        print(f"\n[{n}/{total}] {rule.get('rule_id')} — {rule.get('cell_type')}")
        print(f"  Gene:      {rule.get('gene_or_regulon')} ({rule.get('direction')})")
        print(f"  Source:    {rule.get('evidence_source')}")
        print(f"  PMID:      {rule.get('pmid')}")
        print(f"  Tissue:    {rule.get('tissue_context')}")
        basis = rule.get("mechanistic_basis", "")
        if len(basis) > 120:
            basis = basis[:117] + "..."
        print(f"  Basis:     {basis}")

    def _prompt_action(self, rule: dict) -> str:
        try:
            return input("  Action [a/s/r/q]: ").strip().lower() or "s"
        except (EOFError, KeyboardInterrupt):
            return "q"

    def _edit_rule(self, rule: dict) -> dict:
        """Allow researcher to edit mechanistic_basis and pmid before accepting."""
        edited = dict(rule)
        print(f"\n  Current mechanistic_basis:\n  {edited.get('mechanistic_basis', '')}")
        try:
            new_basis = input("  New mechanistic_basis (Enter to keep): ").strip()
            if new_basis:
                edited["mechanistic_basis"] = new_basis

            new_pmid = input(f"  PMID [{edited.get('pmid')}] (Enter to keep): ").strip()
            if new_pmid:
                edited["pmid"] = new_pmid

            new_conf = input(f"  Confidence [low/medium/high] "
                             f"[{edited.get('confidence', 'low')}] (Enter to keep): ").strip()
            if new_conf in ("low", "medium", "high"):
                edited["confidence"] = new_conf
        except (EOFError, KeyboardInterrupt):
            pass
        return edited

    def _append_to_kg(self, rules: list[dict]) -> None:
        """Append accepted rules to the main KG JSON file."""
        kg_path = self.kg_out_path
        if kg_path.exists():
            with open(kg_path) as f:
                existing = json.load(f)
        else:
            existing = []

        # Dedup: don't add if rule_id already in KG
        existing_ids = {r.get("rule_id") for r in existing}
        new_rules = [r for r in rules if r.get("rule_id") not in existing_ids]

        if not new_rules:
            print("  (no new rules to append — all already in KG)")
            return

        existing.extend(new_rules)
        with open(kg_path, "w") as f:
            json.dump(existing, f, indent=2)
        print(f"  Appended {len(new_rules)} rules to {kg_path}")


def batch_approve(
    candidates_path: str,
    rule_ids: list[str],
    kg_out_path: str,
    mechanistic_overrides: dict[str, str] | None = None,
    pmid_overrides: dict[str, str] | None = None,
    confidence_overrides: dict[str, str] | None = None,
    kg_version: str = "0.2.0",
) -> list[dict]:
    """
    Non-interactive batch approval: promotes a list of rule_ids to ACTIVE.

    Parameters
    ----------
    candidates_path : str
    rule_ids : list[str]
        Rule IDs to approve (must exist in candidates with PENDING_REVIEW).
    kg_out_path : str
        KG JSON file to append approved rules to.
    mechanistic_overrides : dict[rule_id -> new_basis]
        Override mechanistic_basis for specific rules.
    pmid_overrides : dict[rule_id -> pmid]
        Override pmid for specific rules.
    confidence_overrides : dict[rule_id -> confidence]
    kg_version : str

    Returns
    -------
    list[dict] — the approved rule dicts.
    """
    metadata, candidates = load_candidates(candidates_path)
    rule_id_set = set(rule_ids)
    mechanistic_overrides = mechanistic_overrides or {}
    pmid_overrides = pmid_overrides or {}
    confidence_overrides = confidence_overrides or {}

    approved: list[dict] = []
    for rule in candidates:
        if rule.get("rule_id") not in rule_id_set:
            continue
        if rule.get("status") != "PENDING_REVIEW":
            continue

        rule = dict(rule)
        rule["status"] = "ACTIVE"
        rule["added_in_version"] = kg_version

        rid = rule["rule_id"]
        if rid in mechanistic_overrides:
            rule["mechanistic_basis"] = mechanistic_overrides[rid]
        if rid in pmid_overrides:
            rule["pmid"] = pmid_overrides[rid]
        if rid in confidence_overrides:
            rule["confidence"] = confidence_overrides[rid]

        approved.append(rule)

    # Append to KG
    kg_path = Path(kg_out_path)
    if kg_path.exists():
        with open(kg_path) as f:
            existing = json.load(f)
    else:
        existing = []

    existing_ids = {r.get("rule_id") for r in existing}
    new_rules = [r for r in approved if r.get("rule_id") not in existing_ids]
    existing.extend(new_rules)

    with open(kg_path, "w") as f:
        json.dump(existing, f, indent=2)

    # Update candidates status
    for rule in candidates:
        if rule.get("rule_id") in rule_id_set:
            rule["status"] = "ACTIVE"
    save_candidates(metadata, candidates, candidates_path)

    print(f"[batch_approve] Approved {len(approved)} rules → {kg_out_path}")
    return approved
