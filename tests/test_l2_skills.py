from __future__ import annotations

import pytest

from tools.eval.l2_validator import L2Validator, SKILLS_DIR

_EXPECTED_SKILLS = {
    "virtual_screen.md",
    "hit_gen_sbdd.md",
    "lead_opt.md",
    "fep_campaign.md",
    "antibody_design.md",
}


def test_all_expected_l2_skills_present():
    found = {f.name for f in SKILLS_DIR.rglob("*.md")}
    missing = _EXPECTED_SKILLS - found
    assert not missing, f"Missing L2 skill files: {missing}"


def test_l2_skills_no_unknown_tool_params():
    results = L2Validator().validate_all(SKILLS_DIR)
    errors = [e for r in results for e in r.errors]
    if errors:
        lines = [
            f"  {e.skill_file.name} [{e.step_hint}]: {e.message}"
            for e in errors
        ]
        pytest.fail("L2 skill validation errors:\n" + "\n".join(lines))


def test_l2_skills_required_sections(capsys):
    results = L2Validator().validate_all(SKILLS_DIR)
    warnings = [w for r in results for w in r.warnings if "Missing" in w.message]
    for w in warnings:
        print(f"WARN: {w.skill_file.name}: {w.message}")
