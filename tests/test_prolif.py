from __future__ import annotations

from unittest.mock import MagicMock, call, patch

import pytest

from tools.analysis.prolif import run_prolif


# ── helpers ───────────────────────────────────────────────────────────────────

def _make_prolif_mocks():
    """Build a minimal mock for prolif + MDAnalysis."""
    import pandas as pd

    # Fingerprint DataFrame with 3 interaction columns
    cols = pd.MultiIndex.from_tuples(
        [("LIG", "ASP189", "HBAcceptor"),
         ("LIG", "PHE382", "Hydrophobic"),
         ("LIG", "LYS67",  "Ionic")],
        names=["ligand", "protein", "interaction"],
    )
    df = pd.DataFrame([[1, 1, 0]], columns=cols)

    fp_mock = MagicMock()
    fp_mock.to_dataframe.return_value = df

    plf_mock = MagicMock()
    plf_mock.Fingerprint.return_value = fp_mock
    plf_mock.Molecule.from_mda.return_value = MagicMock()

    # Trajectory: two frames (two poses)
    traj = [MagicMock(), MagicMock()]
    u_mock = MagicMock()
    u_mock.trajectory.__iter__ = MagicMock(return_value=iter(traj))
    u_mock.atoms = MagicMock()

    mda_mock = MagicMock()
    mda_mock.Universe.return_value = u_mock

    return plf_mock, mda_mock, df


def _patch_imports(plf_mock, mda_mock):
    return patch.dict("sys.modules", {
        "prolif":        plf_mock,
        "MDAnalysis":    mda_mock,
    })


# ── import fallback ────────────────────────────────────────────────────────────

def test_raises_when_prolif_unavailable(tmp_path, mock_receptor_pdb):
    sdf = tmp_path / "lig.sdf"
    sdf.write_text("dummy")

    with patch.dict("sys.modules", {"prolif": None, "MDAnalysis": None}):
        with patch("tools.analysis.prolif._prolif_subprocess") as mock_sub:
            mock_sub.return_value = {
                "interactions": {}, "summary": {"n_contacts": 0},
                "fingerprint": [], "n_poses": 0, "csv_path": None,
            }
            run_prolif(mock_receptor_pdb, str(sdf))
            mock_sub.assert_called_once()


# ── python path (mocked prolif + MDAnalysis) ─────────────────────────────────

def test_returns_required_keys(tmp_path, mock_receptor_pdb):
    import pandas as pd
    plf_mock, mda_mock, _ = _make_prolif_mocks()
    sdf = tmp_path / "lig.sdf"
    sdf.write_text("dummy")

    with _patch_imports(plf_mock, mda_mock):
        result = run_prolif(mock_receptor_pdb, str(sdf), output_dir=str(tmp_path))

    for key in ("interactions", "summary", "fingerprint", "n_poses", "csv_path"):
        assert key in result


def test_summary_has_required_keys(tmp_path, mock_receptor_pdb):
    plf_mock, mda_mock, _ = _make_prolif_mocks()
    sdf = tmp_path / "lig.sdf"
    sdf.write_text("dummy")

    with _patch_imports(plf_mock, mda_mock):
        result = run_prolif(mock_receptor_pdb, str(sdf), output_dir=str(tmp_path))

    for key in ("n_hbond", "n_hydrophobic", "n_ionic", "n_pistack", "n_contacts"):
        assert key in result["summary"]


def test_interactions_parses_residues(tmp_path, mock_receptor_pdb):
    plf_mock, mda_mock, _ = _make_prolif_mocks()
    sdf = tmp_path / "lig.sdf"
    sdf.write_text("dummy")

    with _patch_imports(plf_mock, mda_mock):
        result = run_prolif(mock_receptor_pdb, str(sdf), output_dir=str(tmp_path))

    # Two columns have value 1 (HBAcceptor + Hydrophobic), one has 0 (Ionic)
    assert result["summary"]["n_contacts"] == 2


def test_hbond_counted_in_summary(tmp_path, mock_receptor_pdb):
    plf_mock, mda_mock, _ = _make_prolif_mocks()
    sdf = tmp_path / "lig.sdf"
    sdf.write_text("dummy")

    with _patch_imports(plf_mock, mda_mock):
        result = run_prolif(mock_receptor_pdb, str(sdf), output_dir=str(tmp_path))

    assert result["summary"]["n_hbond"] == 1
    assert result["summary"]["n_hydrophobic"] == 1


def test_fingerprint_is_list_of_ints(tmp_path, mock_receptor_pdb):
    plf_mock, mda_mock, _ = _make_prolif_mocks()
    sdf = tmp_path / "lig.sdf"
    sdf.write_text("dummy")

    with _patch_imports(plf_mock, mda_mock):
        result = run_prolif(mock_receptor_pdb, str(sdf), output_dir=str(tmp_path))

    assert isinstance(result["fingerprint"], list)
    assert all(isinstance(b, int) for b in result["fingerprint"])


def test_csv_written_to_output_dir(tmp_path, mock_receptor_pdb):
    plf_mock, mda_mock, _ = _make_prolif_mocks()
    sdf = tmp_path / "lig.sdf"
    sdf.write_text("dummy")

    with _patch_imports(plf_mock, mda_mock):
        result = run_prolif(mock_receptor_pdb, str(sdf), output_dir=str(tmp_path))

    assert result["csv_path"] is not None
    assert "prolif_fingerprint" in result["csv_path"]


def test_protein_universe_loaded_from_pdb(tmp_path, mock_receptor_pdb):
    plf_mock, mda_mock, _ = _make_prolif_mocks()
    sdf = tmp_path / "lig.sdf"
    sdf.write_text("dummy")

    with _patch_imports(plf_mock, mda_mock):
        run_prolif(mock_receptor_pdb, str(sdf), output_dir=str(tmp_path))

    first_call_arg = mda_mock.Universe.call_args_list[0][0][0]
    assert first_call_arg == mock_receptor_pdb
