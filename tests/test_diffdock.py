from __future__ import annotations

import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
import requests

from tools.docking.diffdock import _run_local, _run_nim, run_diffdock


def _nim_payload(n=3):
    return {
        "ligand_positions": [f"ATOM {i}  CA  ALA\n" for i in range(n)],
        "position_confidence": [0.9 - i * 0.1 for i in range(n)],
        "runtime_s": 28.5,
    }


# ── dispatch ──────────────────────────────────────────────────────────────────

def test_dispatches_to_nim_by_default(imatinib, mock_receptor_pdb, mock_nim_response):
    with patch("requests.post", return_value=mock_nim_response(_nim_payload())) as mock_post:
        run_diffdock(mock_receptor_pdb, imatinib, use_nim=True, nim_api_key="key")
    mock_post.assert_called_once()


def test_dispatches_to_local_when_nim_false(imatinib, mock_receptor_pdb):
    with patch("tools.docking.diffdock._run_local") as mock_local:
        mock_local.return_value = {"poses": [], "best_pose_pdb": None,
                                   "confidence_score": None, "n_poses": 0, "runtime_s": None}
        run_diffdock(mock_receptor_pdb, imatinib, use_nim=False)
    mock_local.assert_called_once()


def test_reads_api_key_from_env(imatinib, mock_receptor_pdb, mock_nim_response, monkeypatch):
    monkeypatch.setenv("NVIDIA_API_KEY", "env-key")
    with patch("requests.post", return_value=mock_nim_response(_nim_payload())) as mock_post:
        run_diffdock(mock_receptor_pdb, imatinib, use_nim=True)
    headers = mock_post.call_args.kwargs["headers"]
    assert headers["Authorization"] == "Bearer env-key"


# ── NIM payload ───────────────────────────────────────────────────────────────

def test_nim_sends_pdb_content(imatinib, mock_receptor_pdb, mock_nim_response):
    with patch("requests.post", return_value=mock_nim_response(_nim_payload())) as mock_post:
        _run_nim(mock_receptor_pdb, imatinib, 5, 20, "key")

    payload = mock_post.call_args.kwargs["json"]
    assert "protein" in payload
    assert "ATOM" in payload["protein"]   # actual PDB content was sent


def test_nim_sends_smiles(imatinib, mock_receptor_pdb, mock_nim_response):
    with patch("requests.post", return_value=mock_nim_response(_nim_payload())) as mock_post:
        _run_nim(mock_receptor_pdb, imatinib, 5, 20, "key")

    payload = mock_post.call_args.kwargs["json"]
    assert payload["ligand"] == imatinib
    assert payload["ligand_file_type"] == "smi"


def test_nim_passes_n_samples(imatinib, mock_receptor_pdb, mock_nim_response):
    with patch("requests.post", return_value=mock_nim_response(_nim_payload())) as mock_post:
        _run_nim(mock_receptor_pdb, imatinib, 7, 20, "key")

    payload = mock_post.call_args.kwargs["json"]
    assert payload["num_poses"] == 7


# ── NIM output ────────────────────────────────────────────────────────────────

def test_nim_writes_pose_pdbs(imatinib, mock_receptor_pdb, mock_nim_response):
    with patch("requests.post", return_value=mock_nim_response(_nim_payload(3))):
        result = _run_nim(mock_receptor_pdb, imatinib, 3, 20, "key")

    assert result["n_poses"] == 3
    for d in result["poses"]:
        assert Path(d["pose_pdb"]).exists()


def test_nim_poses_sorted_by_confidence(imatinib, mock_receptor_pdb, mock_nim_response):
    with patch("requests.post", return_value=mock_nim_response(_nim_payload(3))):
        result = _run_nim(mock_receptor_pdb, imatinib, 3, 20, "key")

    confidences = [d["confidence_score"] for d in result["poses"]]
    assert confidences == sorted(confidences, reverse=True)


def test_nim_best_pose_is_highest_confidence(imatinib, mock_receptor_pdb, mock_nim_response):
    with patch("requests.post", return_value=mock_nim_response(_nim_payload(3))):
        result = _run_nim(mock_receptor_pdb, imatinib, 3, 20, "key")

    best_conf = result["poses"][0]["confidence_score"]
    assert result["confidence_score"] == pytest.approx(best_conf)


def test_nim_returns_runtime(imatinib, mock_receptor_pdb, mock_nim_response):
    with patch("requests.post", return_value=mock_nim_response(_nim_payload())):
        result = _run_nim(mock_receptor_pdb, imatinib, 3, 20, "key")
    assert result["runtime_s"] == pytest.approx(28.5)


def test_nim_empty_positions(imatinib, mock_receptor_pdb, mock_nim_response):
    with patch("requests.post",
               return_value=mock_nim_response({"ligand_positions": [], "position_confidence": []})):
        result = _run_nim(mock_receptor_pdb, imatinib, 3, 20, "key")

    assert result["n_poses"] == 0
    assert result["best_pose_pdb"] is None


def test_nim_http_error_propagates(imatinib, mock_receptor_pdb):
    resp = MagicMock()
    resp.raise_for_status.side_effect = requests.HTTPError("403")
    with patch("requests.post", return_value=resp):
        with pytest.raises(requests.HTTPError):
            _run_nim(mock_receptor_pdb, imatinib, 3, 20, "bad-key")


# ── output schema ─────────────────────────────────────────────────────────────

def test_output_keys_present(imatinib, mock_receptor_pdb, mock_nim_response):
    with patch("requests.post", return_value=mock_nim_response(_nim_payload())):
        result = run_diffdock(mock_receptor_pdb, imatinib, use_nim=True, nim_api_key="k")

    for key in ("poses", "best_pose_pdb", "confidence_score", "n_poses", "runtime_s"):
        assert key in result
