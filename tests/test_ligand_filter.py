from __future__ import annotations

import pytest

from tools.filtering.ligand_filter import run_ligand_filter


_IBUPROFEN = "CC(C)Cc1ccc(C(C)C(=O)O)cc1"   # clean — passes PAINS + BRENK
_CATECHOL  = "Oc1ccccc1O"                   # canonical PAINS_A hit
_NITRO     = "O=[N+]([O-])c1ccccc1"


# ── clean / positive control ──────────────────────────────────────────────────

def test_clean_molecule_passes():
    result = run_ligand_filter(_IBUPROFEN)
    assert result["passed"] is True
    assert result["n_alerts"] == 0
    assert result["failures"] == []


def test_pains_positive_control():
    result = run_ligand_filter(_CATECHOL, filter_sets=["PAINS"])
    assert result["passed"] is False
    assert result["n_alerts"] >= 1
    assert any(f["catalog"].startswith("PAINS") for f in result["failures"])


# ── invalid SMILES ────────────────────────────────────────────────────────────

def test_invalid_smiles_raises_single():
    with pytest.raises(ValueError, match="Invalid SMILES"):
        run_ligand_filter("not!!valid")


def test_invalid_smiles_in_batch_captured():
    result = run_ligand_filter(["CCO", "not!!valid"])
    assert result["n_failed"] >= 1
    bad = next(r for r in result["results"] if "error" in r)
    assert bad["error"] == "Invalid SMILES"


# ── batch behaviour ───────────────────────────────────────────────────────────

def test_batch_imatinib_set_partition(all_smiles):
    result = run_ligand_filter(all_smiles, filter_sets="default")
    assert result["n_passed"] + result["n_failed"] == len(all_smiles)
    assert set(result["passed_smiles"]).isdisjoint(result["failed_smiles"])


def test_batch_list_returns_batch_schema():
    result = run_ligand_filter([_IBUPROFEN])
    expected = {"passed_smiles", "failed_smiles", "n_passed", "n_failed",
                "results", "filter_sets_used"}
    assert set(result.keys()) == expected


# ── output schemas ────────────────────────────────────────────────────────────

def test_single_str_returns_single_schema():
    result = run_ligand_filter(_IBUPROFEN)
    expected = {"passed", "n_alerts", "failures", "smiles_canonical", "filter_sets_used"}
    assert set(result.keys()) == expected


def test_failure_entry_schema():
    result = run_ligand_filter(_CATECHOL, filter_sets=["PAINS"])
    assert result["failures"], "expected at least one failure for catechol"
    for f in result["failures"]:
        assert set(f.keys()) == {"catalog", "alert_name", "description", "smarts"}
        assert isinstance(f["catalog"], str)
        assert isinstance(f["alert_name"], str)
        assert isinstance(f["description"], str)
        assert isinstance(f["smarts"], str)


# ── presets ───────────────────────────────────────────────────────────────────

def test_preset_default_resolves_to_pains_brenk():
    result = run_ligand_filter(_IBUPROFEN, filter_sets="default")
    assert result["filter_sets_used"] == ["PAINS", "BRENK"]


def test_preset_strict_includes_nih_zinc():
    result = run_ligand_filter(_IBUPROFEN, filter_sets="strict")
    assert "NIH" in result["filter_sets_used"]
    assert "ZINC" in result["filter_sets_used"]


def test_explicit_catalog_list_only_applies_listed():
    # When only BRENK is requested, no PAINS-attributed failures should appear
    # in the output regardless of whether the molecule independently hits BRENK
    result = run_ligand_filter(_CATECHOL, filter_sets=["BRENK"])
    assert all(not f["catalog"].startswith("PAINS") for f in result["failures"])


# ── custom SMARTS ─────────────────────────────────────────────────────────────

def test_custom_smarts_flags_match():
    result = run_ligand_filter(
        _NITRO,
        filter_sets=["BRENK"],
        custom_smarts=["[N+](=O)[O-]"],
    )
    assert any(f["catalog"] == "custom" for f in result["failures"])


# ── lilly preset ──────────────────────────────────────────────────────────────

def test_lilly_preset_not_implemented():
    with pytest.raises(NotImplementedError):
        run_ligand_filter(_IBUPROFEN, filter_sets="lilly")
