from __future__ import annotations

import subprocess
from unittest.mock import MagicMock, patch

import pytest

from tools.docking.gnina import _parse_gnina_output, _smiles_to_sdf, run_gnina

CENTER = dict(center_x=10.0, center_y=20.0, center_z=30.0)


# ── dispatch ──────────────────────────────────────────────────────────────────

def test_dispatches_to_gnina_binary(imatinib, mock_receptor_pdb, tmp_path):
    fake_sdf = tmp_path / "out_0.sdf"
    fake_sdf.write_text("")

    with patch("tools.docking.gnina._smiles_to_sdf") as mock_sdf, \
         patch("subprocess.run") as mock_run:
        mock_sdf.return_value = tmp_path / "lig_0.sdf"
        mock_run.return_value = MagicMock(returncode=0)
        run_gnina(mock_receptor_pdb, imatinib, **CENTER, gnina_bin="gnina")

    mock_run.assert_called_once()
    cmd = mock_run.call_args[0][0]
    assert "gnina" in cmd
    assert "--receptor" in cmd
    assert mock_receptor_pdb in cmd


def test_passes_box_params(imatinib, mock_receptor_pdb, tmp_path):
    with patch("tools.docking.gnina._smiles_to_sdf") as mock_sdf, \
         patch("subprocess.run") as mock_run:
        mock_sdf.return_value = tmp_path / "lig.sdf"
        mock_run.return_value = MagicMock(returncode=0)
        run_gnina(mock_receptor_pdb, imatinib, **CENTER, box_size=25.0)

    cmd = mock_run.call_args[0][0]
    assert "--size_x" in cmd
    assert "25.0" in cmd


def test_passes_cnn_scoring(imatinib, mock_receptor_pdb, tmp_path):
    with patch("tools.docking.gnina._smiles_to_sdf") as mock_sdf, \
         patch("subprocess.run") as mock_run:
        mock_sdf.return_value = tmp_path / "lig.sdf"
        mock_run.return_value = MagicMock(returncode=0)
        run_gnina(mock_receptor_pdb, imatinib, **CENTER, cnn_scoring="refinement")

    cmd = mock_run.call_args[0][0]
    assert "--cnn_scoring" in cmd
    assert "refinement" in cmd


def test_covalent_flag_appended_when_provided(imatinib, mock_receptor_pdb, tmp_path):
    with patch("tools.docking.gnina._smiles_to_sdf") as mock_sdf, \
         patch("subprocess.run") as mock_run:
        mock_sdf.return_value = tmp_path / "lig.sdf"
        mock_run.return_value = MagicMock(returncode=0)
        run_gnina(mock_receptor_pdb, imatinib, **CENTER, covalent_res="CYS145")

    cmd = mock_run.call_args[0][0]
    assert "--covalent_rec_atom" in cmd
    assert "CYS145" in cmd


def test_batch_input_runs_once_per_ligand(all_smiles, mock_receptor_pdb, tmp_path):
    with patch("tools.docking.gnina._smiles_to_sdf") as mock_sdf, \
         patch("subprocess.run") as mock_run:
        mock_sdf.side_effect = [tmp_path / f"lig_{i}.sdf" for i in range(len(all_smiles))]
        mock_run.return_value = MagicMock(returncode=0)
        result = run_gnina(mock_receptor_pdb, all_smiles, **CENTER)

    assert mock_run.call_count == len(all_smiles)
    assert len(result["poses"]) == len(all_smiles)


# ── output schema ─────────────────────────────────────────────────────────────

def test_output_has_required_keys(imatinib, mock_receptor_pdb, tmp_path):
    with patch("tools.docking.gnina._smiles_to_sdf") as mock_sdf, \
         patch("subprocess.run") as mock_run:
        mock_sdf.return_value = tmp_path / "lig.sdf"
        mock_run.return_value = MagicMock(returncode=0)
        result = run_gnina(mock_receptor_pdb, imatinib, **CENTER)

    for key in ("poses", "best_pose_sdf", "cnn_affinity", "vina_score", "runtime_s"):
        assert key in result


def test_string_smiles_wrapped_in_list(imatinib, mock_receptor_pdb, tmp_path):
    with patch("tools.docking.gnina._smiles_to_sdf") as mock_sdf, \
         patch("subprocess.run") as mock_run:
        mock_sdf.return_value = tmp_path / "lig.sdf"
        mock_run.return_value = MagicMock(returncode=0)
        result = run_gnina(mock_receptor_pdb, imatinib, **CENTER)

    # String SMILES should be handled same as single-item list
    assert isinstance(result["poses"], list)
    assert len(result["poses"]) == 1


# ── _parse_gnina_output ───────────────────────────────────────────────────────

def test_parse_returns_empty_for_missing_sdf(tmp_path):
    poses = _parse_gnina_output(str(tmp_path / "nonexistent.sdf"), "C")
    assert poses == []


def test_parse_gnina_output_valid_sdf(tmp_path, imatinib):
    from rdkit import Chem
    from rdkit.Chem import AllChem

    sdf_path = tmp_path / "out.sdf"
    mol = Chem.MolFromSmiles(imatinib)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    mol.SetDoubleProp("CNNaffinity", 7.5)
    mol.SetDoubleProp("minimizedAffinity", -9.1)

    w = Chem.SDWriter(str(sdf_path))
    w.write(mol)
    w.close()

    poses = _parse_gnina_output(str(sdf_path), imatinib)
    assert len(poses) == 1
    assert poses[0]["cnn_affinity"] == pytest.approx(7.5)
    assert poses[0]["vina_score"] == pytest.approx(-9.1)
    assert poses[0]["smiles"] == imatinib


# ── _smiles_to_sdf ────────────────────────────────────────────────────────────

def test_smiles_to_sdf_creates_file(tmp_path, imatinib):
    from pathlib import Path
    out = tmp_path / "lig.sdf"
    _smiles_to_sdf(imatinib, out)
    assert out.exists()
    assert out.stat().st_size > 0


def test_smiles_to_sdf_invalid_smiles_raises(tmp_path):
    from rdkit.Chem.rdchem import AtomValenceException
    with pytest.raises(Exception):
        _smiles_to_sdf("not_a_smiles!!!", tmp_path / "bad.sdf")
