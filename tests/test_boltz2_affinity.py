from __future__ import annotations

from unittest.mock import patch

import pytest

from tools.binding_affinity.boltz2_affinity import _predict_single, run_boltz2_affinity

_SINGLE_RESULT = {
    "pdb_path":            "/tmp/fake.pdb",
    "all_pdb_paths":       ["/tmp/fake.pdb"],
    "affinity_kcal_mol":   -9.2,
    "affinity_confidence": 0.88,
    "ptm_score":           0.91,
    "iptm_score":          0.82,
    "runtime_s":           21.0,
}


def _patch_boltz2(result=None):
    return patch(
        "tools.binding_affinity.boltz2_affinity.run_boltz2",
        return_value=result or _SINGLE_RESULT,
    )


# ── single-ligand mode ────────────────────────────────────────────────────────

def test_single_string_calls_boltz2_once(imatinib, abl1_fasta):
    with _patch_boltz2() as mock_b2:
        run_boltz2_affinity(abl1_fasta, imatinib, nim_api_key="k")
    mock_b2.assert_called_once()
    assert mock_b2.call_args.kwargs["ligand_smiles"] == imatinib


def test_single_returns_affinity_fields(imatinib, abl1_fasta):
    with _patch_boltz2():
        result = run_boltz2_affinity(abl1_fasta, imatinib, nim_api_key="k")

    assert result["affinity_kcal_mol"] == pytest.approx(-9.2)
    assert result["affinity_confidence"] == pytest.approx(0.88)
    assert result["iptm_score"] == pytest.approx(0.82)


def test_single_strips_pdb_path(imatinib, abl1_fasta):
    with _patch_boltz2():
        result = run_boltz2_affinity(abl1_fasta, imatinib, nim_api_key="k")
    assert "pdb_path" not in result
    assert "all_pdb_paths" not in result


# ── batch mode ────────────────────────────────────────────────────────────────

def test_batch_calls_boltz2_once_per_ligand(all_smiles, abl1_fasta):
    with _patch_boltz2() as mock_b2:
        run_boltz2_affinity(abl1_fasta, all_smiles, nim_api_key="k")
    assert mock_b2.call_count == len(all_smiles)


def test_batch_returns_all_predictions(all_smiles, abl1_fasta):
    with _patch_boltz2():
        result = run_boltz2_affinity(abl1_fasta, all_smiles, nim_api_key="k")
    assert result["n_ligands"] == len(all_smiles)
    assert len(result["predictions"]) == len(all_smiles)


def test_batch_sorted_ascending_by_affinity(abl1_fasta):
    affinities = [-7.0, -11.0, -9.5, -8.1]
    smiles_list = ["C", "CC", "CCC", "CCCC"]

    results = iter([
        {**_SINGLE_RESULT, "affinity_kcal_mol": a} for a in affinities
    ])
    with patch("tools.binding_affinity.boltz2_affinity.run_boltz2",
               side_effect=results):
        result = run_boltz2_affinity(abl1_fasta, smiles_list, nim_api_key="k")

    returned_affs = [p["affinity_kcal_mol"] for p in result["predictions"]]
    assert returned_affs == sorted(returned_affs)


def test_batch_best_smiles_is_tightest_binder(abl1_fasta):
    affinities = [-7.0, -11.0, -9.5]
    smiles_list = ["C", "CC", "CCC"]
    results = iter([{**_SINGLE_RESULT, "affinity_kcal_mol": a} for a in affinities])
    with patch("tools.binding_affinity.boltz2_affinity.run_boltz2", side_effect=results):
        result = run_boltz2_affinity(abl1_fasta, smiles_list, nim_api_key="k")

    assert result["best_smiles"] == "CC"  # -11.0 is lowest/tightest
    assert result["best_affinity_kcal_mol"] == pytest.approx(-11.0)


def test_batch_failed_ligand_counted(abl1_fasta):
    def _raise_on_second(*args, **kwargs):
        calls = getattr(_raise_on_second, "_calls", 0)
        _raise_on_second._calls = calls + 1
        if calls == 1:
            raise RuntimeError("NIM error")
        return _SINGLE_RESULT

    with patch("tools.binding_affinity.boltz2_affinity.run_boltz2",
               side_effect=_raise_on_second):
        result = run_boltz2_affinity(abl1_fasta, ["C", "CC", "CCC"], nim_api_key="k")

    assert result["n_failed"] == 1
    assert result["n_ligands"] == 3


def test_batch_output_has_smiles_field(all_smiles, abl1_fasta):
    with _patch_boltz2():
        result = run_boltz2_affinity(abl1_fasta, all_smiles, nim_api_key="k")
    for pred in result["predictions"]:
        assert "smiles" in pred
        assert pred["smiles"] in all_smiles
