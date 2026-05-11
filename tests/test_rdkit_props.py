from __future__ import annotations

import pytest

from tools.admet.rdkit_props import run_rdkit_props


# ── single-mode basics ────────────────────────────────────────────────────────

def test_imatinib_passes_lipinski(imatinib):
    result = run_rdkit_props(imatinib)
    # Imatinib MW ~493, logP ~3, HBD=2, HBA=6 — passes Ro5
    assert result["lipinski_pass"] is True
    assert 480 < result["mw"] < 510


def test_single_returns_required_keys(imatinib):
    result = run_rdkit_props(imatinib)
    expected = {
        "smiles", "mw", "logp", "hbd", "hba", "tpsa", "rotatable_bonds",
        "aromatic_rings", "qed", "sa_score",
        "lipinski_pass", "veber_pass", "egan_pass",
    }
    assert expected.issubset(result.keys())


def test_qed_in_unit_interval(imatinib):
    result = run_rdkit_props(imatinib)
    assert 0.0 <= result["qed"] <= 1.0


def test_invalid_smiles_raises_single():
    with pytest.raises(ValueError, match="Invalid SMILES"):
        run_rdkit_props("not!!valid")


# ── batch mode ────────────────────────────────────────────────────────────────

def test_batch_returns_results_list(all_smiles):
    result = run_rdkit_props(all_smiles)
    assert result["n_compounds"] == len(all_smiles)
    assert len(result["results"]) == len(all_smiles)


def test_batch_invalid_smiles_captured():
    result = run_rdkit_props(["CCO", "not!!valid", "CCN"])
    assert result["n_compounds"] == 3
    bad = [r for r in result["results"] if "error" in r]
    assert len(bad) == 1
    assert bad[0]["error"] == "Invalid SMILES"


def test_all_inhibitors_have_qed(all_smiles):
    result = run_rdkit_props(all_smiles)
    for r in result["results"]:
        assert "qed" in r
        assert isinstance(r["qed"], float)


# ── rule flags are bools ──────────────────────────────────────────────────────

def test_rule_flags_are_bools(imatinib):
    result = run_rdkit_props(imatinib)
    for flag in ("lipinski_pass", "veber_pass", "egan_pass"):
        assert isinstance(result[flag], bool)


# ── fingerprint optional ──────────────────────────────────────────────────────

def test_no_fingerprint_by_default(imatinib):
    result = run_rdkit_props(imatinib)
    assert "fingerprint" not in result


def test_morgan2_fingerprint_returns_bits(imatinib):
    result = run_rdkit_props(imatinib, fingerprint="morgan2", fp_bits=1024)
    assert "fingerprint" in result
    assert len(result["fingerprint"]) == 1024
    assert result["fingerprint_type"] == "morgan2"
    assert result["fingerprint_bits"] == 1024
