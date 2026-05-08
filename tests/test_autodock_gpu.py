from __future__ import annotations

import math
import tempfile
from pathlib import Path
from unittest.mock import patch

import pytest

from tools.docking.autodock_gpu import _parse_dlg, run_autodock_gpu


def _write_dlg(path: Path, energies: list[float]) -> None:
    lines = []
    for e in energies:
        lines.append(
            f"DOCKED: USER    Estimated Free Energy of Binding    = {e} kcal/mol\n"
        )
    path.write_text("".join(lines))


def _make_ligand_dir(tmp_path: Path, n: int = 2) -> Path:
    lig_dir = tmp_path / "ligands"
    lig_dir.mkdir()
    for i in range(n):
        (lig_dir / f"lig_{i}.pdbqt").write_text(f"ATOM  {i}\n")
    return lig_dir


# ── argument building ─────────────────────────────────────────────────────────

def test_passes_fld_path(tmp_path):
    lig_dir = _make_ligand_dir(tmp_path)
    with patch("subprocess.run") as mock_run, \
         patch("tools.docking.autodock_gpu.Path.glob", return_value=[]):
        mock_run.return_value = None
        try:
            run_autodock_gpu("receptor.pdbqt", str(lig_dir), "/data/receptor.gpf")
        except Exception:
            pass
    cmd = mock_run.call_args[0][0]
    assert "--ffile" in cmd
    assert cmd[cmd.index("--ffile") + 1] == "/data/receptor.maps.fld"


def test_passes_nrun(tmp_path):
    lig_dir = _make_ligand_dir(tmp_path)
    with patch("subprocess.run") as mock_run, \
         patch("tools.docking.autodock_gpu.Path.glob", return_value=[]):
        mock_run.return_value = None
        try:
            run_autodock_gpu("receptor.pdbqt", str(lig_dir), "r.gpf", n_runs=50)
        except Exception:
            pass
    cmd = mock_run.call_args[0][0]
    assert cmd[cmd.index("--nrun") + 1] == "50"


def test_passes_heuristics_flag(tmp_path):
    lig_dir = _make_ligand_dir(tmp_path)
    with patch("subprocess.run") as mock_run, \
         patch("tools.docking.autodock_gpu.Path.glob", return_value=[]):
        mock_run.return_value = None
        try:
            run_autodock_gpu("receptor.pdbqt", str(lig_dir), "r.gpf", heuristics=True)
        except Exception:
            pass
    cmd = mock_run.call_args[0][0]
    assert "--heuristics" in cmd
    assert cmd[cmd.index("--heuristics") + 1] == "1"


def test_no_heuristics_flag_when_false(tmp_path):
    lig_dir = _make_ligand_dir(tmp_path)
    with patch("subprocess.run") as mock_run, \
         patch("tools.docking.autodock_gpu.Path.glob", return_value=[]):
        mock_run.return_value = None
        try:
            run_autodock_gpu("receptor.pdbqt", str(lig_dir), "r.gpf", heuristics=False)
        except Exception:
            pass
    cmd = mock_run.call_args[0][0]
    assert "--heuristics" not in cmd


def test_custom_binary(tmp_path):
    lig_dir = _make_ligand_dir(tmp_path)
    with patch("subprocess.run") as mock_run, \
         patch("tools.docking.autodock_gpu.Path.glob", return_value=[]):
        mock_run.return_value = None
        try:
            run_autodock_gpu("r.pdbqt", str(lig_dir), "r.gpf",
                             autodock_gpu_bin="/usr/local/bin/autodock_gpu_64wi")
        except Exception:
            pass
    cmd = mock_run.call_args[0][0]
    assert cmd[0] == "/usr/local/bin/autodock_gpu_64wi"


# ── no PDBQT files ────────────────────────────────────────────────────────────

def test_raises_if_no_pdbqt_files(tmp_path):
    empty_dir = tmp_path / "empty"
    empty_dir.mkdir()
    with pytest.raises(FileNotFoundError):
        run_autodock_gpu("r.pdbqt", str(empty_dir), "r.gpf")


# ── DLG parser ────────────────────────────────────────────────────────────────

def test_parse_dlg_extracts_best_energy(tmp_path):
    dlg = tmp_path / "results_lig0.dlg"
    _write_dlg(dlg, [-8.5, -7.2, -6.9])
    result = _parse_dlg(dlg)
    assert result["best_energy"] == pytest.approx(-8.5)


def test_parse_dlg_returns_none_for_empty_file(tmp_path):
    dlg = tmp_path / "empty.dlg"
    dlg.write_text("NO USEFUL CONTENT\n")
    assert _parse_dlg(dlg) is None


def test_parse_dlg_ki_estimate_formula(tmp_path):
    dlg = tmp_path / "results_lig1.dlg"
    energy = -9.0
    _write_dlg(dlg, [energy])
    result = _parse_dlg(dlg)
    expected_ki = math.exp(energy * 1000 / (1.987 * 298.15))
    assert result["ki_estimate"] == pytest.approx(expected_ki)


def test_parse_dlg_sets_ligand_id_from_stem(tmp_path):
    dlg = tmp_path / "results_mylig.dlg"
    _write_dlg(dlg, [-7.0])
    result = _parse_dlg(dlg)
    assert result["ligand_id"] == "results_mylig"


# ── output schema ─────────────────────────────────────────────────────────────

def test_output_has_required_keys(tmp_path):
    lig_dir = _make_ligand_dir(tmp_path, 1)
    workdir = tempfile.mkdtemp()
    dlg_path = Path(workdir) / "results_lig_0.dlg"
    _write_dlg(dlg_path, [-9.5, -8.0])

    with patch("subprocess.run"), \
         patch("tools.docking.autodock_gpu.Path.glob", return_value=[dlg_path]):
        result = run_autodock_gpu("r.pdbqt", str(lig_dir), "r.gpf")

    for key in ("results", "best_energy", "n_docked", "failed", "runtime_s"):
        assert key in result


def test_results_sorted_by_energy_ascending(tmp_path):
    lig_dir = _make_ligand_dir(tmp_path, 2)
    workdir = tempfile.mkdtemp()
    dlg1 = Path(workdir) / "results_a.dlg"
    dlg2 = Path(workdir) / "results_b.dlg"
    _write_dlg(dlg1, [-6.0])
    _write_dlg(dlg2, [-9.0])

    with patch("subprocess.run"), \
         patch("tools.docking.autodock_gpu.Path.glob", return_value=[dlg1, dlg2]):
        result = run_autodock_gpu("r.pdbqt", str(lig_dir), "r.gpf")

    energies = [r["best_energy"] for r in result["results"]]
    assert energies == sorted(energies)


def test_best_energy_is_minimum(tmp_path):
    lig_dir = _make_ligand_dir(tmp_path, 2)
    workdir = tempfile.mkdtemp()
    dlg1 = Path(workdir) / "results_a.dlg"
    dlg2 = Path(workdir) / "results_b.dlg"
    _write_dlg(dlg1, [-6.0])
    _write_dlg(dlg2, [-9.0])

    with patch("subprocess.run"), \
         patch("tools.docking.autodock_gpu.Path.glob", return_value=[dlg1, dlg2]):
        result = run_autodock_gpu("r.pdbqt", str(lig_dir), "r.gpf")

    assert result["best_energy"] == pytest.approx(-9.0)
