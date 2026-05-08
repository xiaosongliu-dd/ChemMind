from __future__ import annotations

import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from tools.docking.vina import _parse_vina_pdbqt, run_vina


def _write_pdbqt(path: Path, scores: list[float]) -> None:
    lines = []
    for i, score in enumerate(scores):
        lines.append(f"MODEL {i + 1}\n")
        lines.append(f"REMARK VINA RESULT:   {score}    0.000    0.000\n")
        lines.append("ATOM      1  C   LIG     1       0.000   0.000   0.000  1.00  0.00    +0.000 C\n")
        lines.append("ENDMDL\n")
    path.write_text("".join(lines))


# ── dispatch ──────────────────────────────────────────────────────────────────

def test_uses_vina_gpu_by_default(imatinib, tmp_path):
    with patch("tools.docking.vina._smiles_to_pdbqt", return_value=tmp_path / "lig.pdbqt"), \
         patch("tools.docking.vina._parse_vina_pdbqt", return_value=[]), \
         patch("subprocess.run") as mock_run:
        run_vina("receptor.pdbqt", imatinib, 10.0, 20.0, 30.0)
    cmd = mock_run.call_args[0][0]
    assert cmd[0] == "vina_gpu"


def test_uses_cpu_vina_when_gpu_false(imatinib, tmp_path):
    with patch("tools.docking.vina._smiles_to_pdbqt", return_value=tmp_path / "lig.pdbqt"), \
         patch("tools.docking.vina._parse_vina_pdbqt", return_value=[]), \
         patch("subprocess.run") as mock_run:
        run_vina("receptor.pdbqt", imatinib, 10.0, 20.0, 30.0, use_gpu=False)
    cmd = mock_run.call_args[0][0]
    assert cmd[0] == "vina"


def test_custom_vina_bin(imatinib, tmp_path):
    with patch("tools.docking.vina._smiles_to_pdbqt", return_value=tmp_path / "lig.pdbqt"), \
         patch("tools.docking.vina._parse_vina_pdbqt", return_value=[]), \
         patch("subprocess.run") as mock_run:
        run_vina("receptor.pdbqt", imatinib, 10.0, 20.0, 30.0, vina_bin="/opt/vina/vina_gpu2")
    cmd = mock_run.call_args[0][0]
    assert cmd[0] == "/opt/vina/vina_gpu2"


# ── box parameters ────────────────────────────────────────────────────────────

def test_passes_center_coords(imatinib, tmp_path):
    with patch("tools.docking.vina._smiles_to_pdbqt", return_value=tmp_path / "lig.pdbqt"), \
         patch("tools.docking.vina._parse_vina_pdbqt", return_value=[]), \
         patch("subprocess.run") as mock_run:
        run_vina("receptor.pdbqt", imatinib, 12.5, -3.7, 8.1)
    cmd = mock_run.call_args[0][0]
    assert "--center_x" in cmd
    idx = cmd.index("--center_x")
    assert cmd[idx + 1] == "12.5"
    assert cmd[cmd.index("--center_y") + 1] == "-3.7"
    assert cmd[cmd.index("--center_z") + 1] == "8.1"


def test_passes_box_size(imatinib, tmp_path):
    with patch("tools.docking.vina._smiles_to_pdbqt", return_value=tmp_path / "lig.pdbqt"), \
         patch("tools.docking.vina._parse_vina_pdbqt", return_value=[]), \
         patch("subprocess.run") as mock_run:
        run_vina("receptor.pdbqt", imatinib, 0.0, 0.0, 0.0, box_size=30.0)
    cmd = mock_run.call_args[0][0]
    assert cmd[cmd.index("--size_x") + 1] == "30.0"


def test_passes_exhaustiveness(imatinib, tmp_path):
    with patch("tools.docking.vina._smiles_to_pdbqt", return_value=tmp_path / "lig.pdbqt"), \
         patch("tools.docking.vina._parse_vina_pdbqt", return_value=[]), \
         patch("subprocess.run") as mock_run:
        run_vina("receptor.pdbqt", imatinib, 0.0, 0.0, 0.0, exhaustiveness=32)
    cmd = mock_run.call_args[0][0]
    assert cmd[cmd.index("--exhaustiveness") + 1] == "32"


def test_passes_flex_residues_when_set(imatinib, tmp_path):
    with patch("tools.docking.vina._smiles_to_pdbqt", return_value=tmp_path / "lig.pdbqt"), \
         patch("tools.docking.vina._parse_vina_pdbqt", return_value=[]), \
         patch("subprocess.run") as mock_run:
        run_vina("receptor.pdbqt", imatinib, 0.0, 0.0, 0.0,
                 flex_residues_pdbqt="flex.pdbqt")
    cmd = mock_run.call_args[0][0]
    assert "--flex" in cmd
    assert cmd[cmd.index("--flex") + 1] == "flex.pdbqt"


def test_no_flex_flag_by_default(imatinib, tmp_path):
    with patch("tools.docking.vina._smiles_to_pdbqt", return_value=tmp_path / "lig.pdbqt"), \
         patch("tools.docking.vina._parse_vina_pdbqt", return_value=[]), \
         patch("subprocess.run") as mock_run:
        run_vina("receptor.pdbqt", imatinib, 0.0, 0.0, 0.0)
    cmd = mock_run.call_args[0][0]
    assert "--flex" not in cmd


# ── batch mode ────────────────────────────────────────────────────────────────

def test_batch_calls_subprocess_once_per_ligand(all_smiles, tmp_path):
    with patch("tools.docking.vina._smiles_to_pdbqt", return_value=tmp_path / "lig.pdbqt"), \
         patch("tools.docking.vina._parse_vina_pdbqt", return_value=[]), \
         patch("subprocess.run") as mock_run:
        run_vina("receptor.pdbqt", all_smiles, 0.0, 0.0, 0.0)
    assert mock_run.call_count == len(all_smiles)


# ── PDBQT parser ──────────────────────────────────────────────────────────────

def test_parse_vina_pdbqt_reads_scores(tmp_path):
    pdbqt = tmp_path / "out.pdbqt"
    _write_pdbqt(pdbqt, [-7.5, -6.2, -5.9])
    poses = _parse_vina_pdbqt(str(pdbqt), "C")
    assert len(poses) == 3
    assert poses[0]["vina_score"] == pytest.approx(-7.5)


def test_parse_vina_pdbqt_rank_increments(tmp_path):
    pdbqt = tmp_path / "out.pdbqt"
    _write_pdbqt(pdbqt, [-8.0, -7.0])
    poses = _parse_vina_pdbqt(str(pdbqt), "CC")
    assert poses[0]["rank"] == 1
    assert poses[1]["rank"] == 2


def test_parse_vina_pdbqt_smiles_propagated(tmp_path):
    pdbqt = tmp_path / "out.pdbqt"
    _write_pdbqt(pdbqt, [-6.0])
    poses = _parse_vina_pdbqt(str(pdbqt), "CCO")
    assert poses[0]["smiles"] == "CCO"


def test_parse_vina_pdbqt_missing_file_returns_empty(tmp_path):
    poses = _parse_vina_pdbqt(str(tmp_path / "nonexistent.pdbqt"), "C")
    assert poses == []


# ── output schema ─────────────────────────────────────────────────────────────

def test_output_has_required_keys(imatinib, tmp_path):
    with patch("tools.docking.vina._smiles_to_pdbqt", return_value=tmp_path / "lig.pdbqt"), \
         patch("tools.docking.vina._parse_vina_pdbqt", return_value=[{"pose_pdbqt": "x.pdbqt", "vina_score": -7.0, "rank": 1, "smiles": imatinib}]), \
         patch("subprocess.run"):
        result = run_vina("receptor.pdbqt", imatinib, 0.0, 0.0, 0.0)
    for key in ("poses", "best_pose_pdbqt", "vina_score", "runtime_s"):
        assert key in result
