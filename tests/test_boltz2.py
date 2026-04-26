from __future__ import annotations

import base64
import json
import subprocess
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from tools.structure.boltz2 import _run_local, _run_nim, run_boltz2

FASTA = "MTEYKLVVVGAGGVGKSALTIQLIQNHFV"
SMILES = "CC1=CC=CC=C1"

# ── helpers ──────────────────────────────────────────────────────────────────

def _fake_structure(affinity=-8.5, confidence=0.9, ptm=0.85, iptm=0.80):
    pdb_bytes = b"ATOM      1  CA  ALA A   1       1.000   2.000   3.000\n"
    return {
        "pdb": base64.b64encode(pdb_bytes).decode(),
        "affinity_kcal_mol": affinity,
        "affinity_confidence": confidence,
        "ptm": ptm,
        "iptm": iptm,
    }


def _nim_response(n=2):
    return {"structures": [_fake_structure() for _ in range(n)], "runtime_s": 42.0}


# ── run_boltz2 dispatch ───────────────────────────────────────────────────────

def test_run_boltz2_dispatches_to_nim():
    with patch("tools.structure.boltz2._run_nim") as mock_nim:
        mock_nim.return_value = {"pdb_path": "/tmp/x.pdb"}
        result = run_boltz2(FASTA, use_nim=True, nim_api_key="test-key")
    mock_nim.assert_called_once()
    assert result["pdb_path"] == "/tmp/x.pdb"


def test_run_boltz2_dispatches_to_local():
    with patch("tools.structure.boltz2._run_local") as mock_local:
        mock_local.return_value = {"pdb_path": "/tmp/x.pdb"}
        result = run_boltz2(FASTA, use_nim=False)
    mock_local.assert_called_once()
    assert result["pdb_path"] == "/tmp/x.pdb"


def test_run_boltz2_reads_api_key_from_env(monkeypatch):
    monkeypatch.setenv("NVIDIA_API_KEY", "env-key")
    with patch("tools.structure.boltz2._run_nim") as mock_nim:
        mock_nim.return_value = {}
        run_boltz2(FASTA, use_nim=True)
    _, kwargs = mock_nim.call_args
    # api_key is the 5th positional arg
    assert mock_nim.call_args[0][4] == "env-key"


# ── _run_nim ─────────────────────────────────────────────────────────────────

def _mock_post(payload, n_structures=2):
    mock_resp = MagicMock()
    mock_resp.json.return_value = _nim_response(n_structures)
    mock_resp.raise_for_status = MagicMock()
    return mock_resp


def test_nim_protein_only():
    with patch("requests.post", return_value=_mock_post(None)) as mock_post:
        result = _run_nim(FASTA, None, None, 2, "key")

    payload = mock_post.call_args.kwargs["json"]
    sequences = payload["sequences"]
    assert len(sequences) == 1
    assert sequences[0]["protein"]["sequence"] == FASTA
    assert result["pdb_path"] is not None
    assert Path(result["pdb_path"]).exists()


def test_nim_protein_plus_ligand():
    with patch("requests.post", return_value=_mock_post(None)) as mock_post:
        _run_nim(FASTA, SMILES, None, 2, "key")

    sequences = mock_post.call_args.kwargs["json"]["sequences"]
    assert any("ligand" in s for s in sequences)
    assert any(s.get("ligand", {}).get("smiles") == SMILES for s in sequences)


def test_nim_protein_plus_antibody():
    antibody = "EVQLVESGGGLVQPGGSLRLSCAASGFTFS"
    with patch("requests.post", return_value=_mock_post(None)) as mock_post:
        _run_nim(FASTA, None, antibody, 2, "key")

    sequences = mock_post.call_args.kwargs["json"]["sequences"]
    ids = [s.get("protein", {}).get("id") for s in sequences if "protein" in s]
    assert "B" in ids


def test_nim_returns_all_pdb_paths():
    with patch("requests.post", return_value=_mock_post(None, n_structures=3)):
        result = _run_nim(FASTA, None, None, 3, "key")

    assert len(result["all_pdb_paths"]) == 3
    for p in result["all_pdb_paths"]:
        assert Path(p).exists()


def test_nim_returns_scores():
    with patch("requests.post", return_value=_mock_post(None)):
        result = _run_nim(FASTA, None, None, 2, "key")

    assert result["affinity_kcal_mol"] == pytest.approx(-8.5)
    assert result["affinity_confidence"] == pytest.approx(0.9)
    assert result["ptm_score"] == pytest.approx(0.85)
    assert result["iptm_score"] == pytest.approx(0.80)
    assert result["runtime_s"] == pytest.approx(42.0)


def test_nim_empty_structures():
    mock_resp = MagicMock()
    mock_resp.json.return_value = {"structures": [], "runtime_s": 1.0}
    mock_resp.raise_for_status = MagicMock()
    with patch("requests.post", return_value=mock_resp):
        result = _run_nim(FASTA, None, None, 1, "key")

    assert result["pdb_path"] is None
    assert result["all_pdb_paths"] == []


def test_nim_http_error_propagates():
    import requests

    mock_resp = MagicMock()
    mock_resp.raise_for_status.side_effect = requests.HTTPError("401 Unauthorized")
    with patch("requests.post", return_value=mock_resp):
        with pytest.raises(requests.HTTPError):
            _run_nim(FASTA, None, None, 1, "bad-key")


def test_nim_authorization_header():
    with patch("requests.post", return_value=_mock_post(None)) as mock_post:
        _run_nim(FASTA, None, None, 1, "my-secret")

    headers = mock_post.call_args.kwargs["headers"]
    assert headers["Authorization"] == "Bearer my-secret"


# ── _run_local ────────────────────────────────────────────────────────────────

def _make_local_output(workdir: Path, n_pdbs=2, with_summary=True):
    out = workdir / "out"
    out.mkdir(parents=True, exist_ok=True)
    for i in range(n_pdbs):
        (out / f"sample_{i}.pdb").write_text(f"ATOM {i}\n")
    if with_summary:
        summary = {
            "affinity_kcal_mol": -7.2,
            "affinity_confidence": 0.88,
            "ptm": 0.91,
            "iptm": 0.78,
            "runtime_s": 120.0,
        }
        (out / "summary.json").write_text(json.dumps(summary))


def _fake_subprocess(workdir):
    def _run(cmd, **kwargs):
        _make_local_output(workdir)
        return MagicMock(returncode=0)
    return _run


def test_local_writes_fasta_and_calls_boltz():
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        with patch("tempfile.mkdtemp", return_value=str(tmp_path)):
            with patch("subprocess.run", side_effect=_fake_subprocess(tmp_path)) as mock_run:
                result = _run_local(FASTA, None, 2)

        cmd = mock_run.call_args[0][0]
        assert "boltz" in cmd
        assert "predict" in cmd
        assert str(tmp_path / "input.fasta") in cmd
        assert "--num_diffusion_samples" in cmd
        assert "2" in cmd


def test_local_passes_ligand_flag():
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        with patch("tempfile.mkdtemp", return_value=str(tmp_path)):
            with patch("subprocess.run", side_effect=_fake_subprocess(tmp_path)) as mock_run:
                _run_local(FASTA, SMILES, 1)

        cmd = mock_run.call_args[0][0]
        assert "--ligand" in cmd


def test_local_no_ligand_flag_when_none():
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        with patch("tempfile.mkdtemp", return_value=str(tmp_path)):
            with patch("subprocess.run", side_effect=_fake_subprocess(tmp_path)) as mock_run:
                _run_local(FASTA, None, 1)

        cmd = mock_run.call_args[0][0]
        assert "--ligand" not in cmd


def test_local_returns_scores_from_summary():
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        with patch("tempfile.mkdtemp", return_value=str(tmp_path)):
            with patch("subprocess.run", side_effect=_fake_subprocess(tmp_path)):
                result = _run_local(FASTA, None, 2)

    assert result["affinity_kcal_mol"] == pytest.approx(-7.2)
    assert result["ptm_score"] == pytest.approx(0.91)
    assert result["runtime_s"] == pytest.approx(120.0)


def test_local_missing_summary_returns_none_scores():
    def _run_no_summary(workdir):
        def _run(cmd, **kwargs):
            _make_local_output(workdir, with_summary=False)
            return MagicMock(returncode=0)
        return _run

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        with patch("tempfile.mkdtemp", return_value=str(tmp_path)):
            with patch("subprocess.run", side_effect=_run_no_summary(tmp_path)):
                result = _run_local(FASTA, None, 1)

    assert result["affinity_kcal_mol"] is None
    assert result["ptm_score"] is None


def test_local_subprocess_failure_raises():
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        with patch("tempfile.mkdtemp", return_value=str(tmp_path)):
            with patch("subprocess.run", side_effect=subprocess.CalledProcessError(1, "boltz")):
                with pytest.raises(subprocess.CalledProcessError):
                    _run_local(FASTA, None, 1)
