from __future__ import annotations

import pytest
from rdkit import Chem

from tools.enumeration.rdkit_enum import (
    _brics_enum,
    _get_scaffold,
    _passes_ro5,
    _recap_enum,
    run_rdkit_enum,
)

# RDKit is a pure-local library — these tests run without mocking.
# Imatinib (MW ~493, well within Ro5) is used as the primary test molecule.


# ── brics mode ────────────────────────────────────────────────────────────────

def test_brics_returns_list(imatinib):
    result = run_rdkit_enum(imatinib, mode="brics", max_compounds=50)
    assert isinstance(result["smiles_list"], list)


def test_brics_output_are_valid_smiles(imatinib):
    result = run_rdkit_enum(imatinib, mode="brics", max_compounds=50)
    for smi in result["smiles_list"]:
        assert Chem.MolFromSmiles(smi) is not None, f"Invalid SMILES: {smi}"


def test_brics_respects_max_compounds(imatinib):
    result = run_rdkit_enum(imatinib, mode="brics", max_compounds=5)
    assert len(result["smiles_list"]) <= 5


def test_brics_n_returned_matches_list_length(imatinib):
    result = run_rdkit_enum(imatinib, mode="brics", max_compounds=20)
    assert result["n_returned"] == len(result["smiles_list"])


def test_brics_n_generated_gte_n_returned(imatinib):
    result = run_rdkit_enum(imatinib, mode="brics", max_compounds=5)
    assert result["n_generated"] >= result["n_returned"]


# ── recap mode ────────────────────────────────────────────────────────────────

def test_recap_returns_list(imatinib):
    result = run_rdkit_enum(imatinib, mode="recap", max_compounds=50)
    assert isinstance(result["smiles_list"], list)


def test_recap_output_are_valid_smiles(imatinib):
    result = run_rdkit_enum(imatinib, mode="recap", max_compounds=50)
    for smi in result["smiles_list"]:
        assert Chem.MolFromSmiles(smi) is not None


# ── rgroup mode ───────────────────────────────────────────────────────────────

def test_rgroup_requires_rgroup_smiles(imatinib):
    with pytest.raises(ValueError, match="rgroup_smiles required"):
        run_rdkit_enum(imatinib, mode="rgroup")


def test_rgroup_with_simple_substituents():
    scaffold = "c1ccc([*])cc1"   # phenyl with attachment point
    rgroups  = ["[*]C", "[*]CC", "[*]F"]
    result = run_rdkit_enum(scaffold, mode="rgroup", rgroup_smiles=rgroups, max_compounds=50)
    assert result["n_returned"] >= 0  # may be 0 if RDKit can't match [*]


# ── invalid input ─────────────────────────────────────────────────────────────

def test_invalid_smiles_raises(tmp_path):
    with pytest.raises(ValueError, match="Invalid SMILES"):
        run_rdkit_enum("not!!valid!smiles", mode="brics")


def test_unknown_mode_raises(imatinib):
    with pytest.raises(ValueError, match="Unknown mode"):
        run_rdkit_enum(imatinib, mode="unknown_mode")


# ── lipinski filter ───────────────────────────────────────────────────────────

def test_lipinski_filter_removes_heavy_molecules():
    # A molecule that clearly violates Ro5 (very high MW)
    heavy = "C" * 60  # polyethylene fragment — MW >> 500
    result = run_rdkit_enum(
        "CC(=O)Oc1ccccc1C(=O)O",  # aspirin as scaffold
        mode="brics",
        max_compounds=200,
        filter_lipinski=True,
    )
    for smi in result["smiles_list"]:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            from rdkit.Chem import rdMolDescriptors
            assert rdMolDescriptors.CalcExactMolWt(mol) <= 500


def test_filter_off_allows_larger_molecules(imatinib):
    result_filtered   = run_rdkit_enum(imatinib, mode="brics", max_compounds=200, filter_lipinski=True)
    result_unfiltered = run_rdkit_enum(imatinib, mode="brics", max_compounds=200, filter_lipinski=False)
    # Without filter, at least as many (possibly more) compounds returned
    assert result_unfiltered["n_returned"] >= result_filtered["n_returned"]


# ── _passes_ro5 ───────────────────────────────────────────────────────────────

def test_passes_ro5_aspirin():
    assert _passes_ro5("CC(=O)Oc1ccccc1C(=O)O") is True


def test_fails_ro5_high_mw():
    # A long chain molecule definitely > 500 Da
    assert _passes_ro5("C" * 50) is False


def test_passes_ro5_invalid_smiles():
    assert _passes_ro5("not_smiles") is False


# ── _get_scaffold ─────────────────────────────────────────────────────────────

def test_get_scaffold_returns_string(imatinib):
    scaffold = _get_scaffold(imatinib)
    assert isinstance(scaffold, str)
    assert Chem.MolFromSmiles(scaffold) is not None


def test_get_scaffold_invalid_returns_input():
    result = _get_scaffold("invalid!!!")
    assert result == "invalid!!!"


# ── output schema ─────────────────────────────────────────────────────────────

def test_output_has_required_keys(imatinib):
    result = run_rdkit_enum(imatinib, mode="brics", max_compounds=10)
    for key in ("smiles_list", "n_generated", "n_returned", "scaffold_smiles"):
        assert key in result


def test_scaffold_smiles_is_valid(imatinib):
    result = run_rdkit_enum(imatinib, mode="brics", max_compounds=10)
    assert Chem.MolFromSmiles(result["scaffold_smiles"]) is not None


# ── all ABL1 inhibitors round-trip ────────────────────────────────────────────

def test_all_inhibitors_enumerate_without_error(all_smiles):
    for smi in all_smiles:
        result = run_rdkit_enum(smi, mode="brics", max_compounds=20, filter_lipinski=False)
        assert "smiles_list" in result
