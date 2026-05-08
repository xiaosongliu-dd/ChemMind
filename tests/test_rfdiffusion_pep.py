from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
import requests

from tools.biologics.rfdiffusion_pep import run_rfdiffusion_pep


_HOTSPOTS = ["A25", "A30", "A35"]


def _nim_payload(n=3, with_plddt=True):
    pdbs = [f"ATOM  {i}  CA  ALA A   1       0.000   0.000   0.000\nEND\n" for i in range(n)]
    payload = {"pdbs": pdbs}
    if with_plddt:
        payload["plddt_array"] = [[0.9 - i * 0.05] * 10 for i in range(n)]
    return payload


def _mock_nim(data):
    resp = MagicMock()
    resp.json.return_value = data
    resp.raise_for_status = MagicMock()
    return resp


# ── dispatch ──────────────────────────────────────────────────────────────────

def test_dispatches_to_nim_by_default(mock_receptor_pdb):
    with patch("requests.post", return_value=_mock_nim(_nim_payload())):
        run_rfdiffusion_pep(mock_receptor_pdb, _HOTSPOTS, nim_api_key="key")
    # no assertion needed — would raise if local was called without RFdiffusion installed


def test_dispatches_to_local_when_nim_false(mock_receptor_pdb):
    with patch("tools.biologics.rfdiffusion_pep._run_local") as mock_local:
        mock_local.return_value = {"designs": [], "best_pdb_path": None,
                                   "n_designs": 0, "hotspot_res": _HOTSPOTS,
                                   "binder_length": (10, 20)}
        run_rfdiffusion_pep(mock_receptor_pdb, _HOTSPOTS, use_nim=False)
    mock_local.assert_called_once()


def test_reads_api_key_from_env(mock_receptor_pdb, monkeypatch):
    monkeypatch.setenv("NVIDIA_API_KEY", "env-nim-key")
    with patch("requests.post", return_value=_mock_nim(_nim_payload())) as mock_post:
        run_rfdiffusion_pep(mock_receptor_pdb, _HOTSPOTS)
    headers = mock_post.call_args.kwargs["headers"]
    assert headers["Authorization"] == "Bearer env-nim-key"


# ── NIM payload ───────────────────────────────────────────────────────────────

def test_nim_sends_pdb_content(mock_receptor_pdb):
    with patch("requests.post", return_value=_mock_nim(_nim_payload())) as mock_post:
        run_rfdiffusion_pep(mock_receptor_pdb, _HOTSPOTS, nim_api_key="k")
    payload = mock_post.call_args.kwargs["json"]
    assert "ATOM" in payload["pdb_str"]


def test_nim_sends_hotspot_residues(mock_receptor_pdb):
    with patch("requests.post", return_value=_mock_nim(_nim_payload())) as mock_post:
        run_rfdiffusion_pep(mock_receptor_pdb, _HOTSPOTS, nim_api_key="k")
    payload = mock_post.call_args.kwargs["json"]
    assert payload["hotspot_res"] == _HOTSPOTS


def test_nim_contig_string_from_binder_length(mock_receptor_pdb):
    with patch("requests.post", return_value=_mock_nim(_nim_payload())) as mock_post:
        run_rfdiffusion_pep(mock_receptor_pdb, _HOTSPOTS,
                            binder_length=(15, 25), nim_api_key="k")
    payload = mock_post.call_args.kwargs["json"]
    assert payload["contigs"] == ["15-25"]


def test_nim_passes_num_designs(mock_receptor_pdb):
    with patch("requests.post", return_value=_mock_nim(_nim_payload(5))) as mock_post:
        run_rfdiffusion_pep(mock_receptor_pdb, _HOTSPOTS, n_designs=5, nim_api_key="k")
    payload = mock_post.call_args.kwargs["json"]
    assert payload["num_designs"] == 5


def test_nim_passes_noise_scale(mock_receptor_pdb):
    with patch("requests.post", return_value=_mock_nim(_nim_payload())) as mock_post:
        run_rfdiffusion_pep(mock_receptor_pdb, _HOTSPOTS, noise_scale=0.5, nim_api_key="k")
    payload = mock_post.call_args.kwargs["json"]
    assert payload["noise_scale_ca"] == pytest.approx(0.5)
    assert payload["noise_scale_frame"] == pytest.approx(0.5)


# ── NIM output ────────────────────────────────────────────────────────────────

def test_returns_correct_n_designs(mock_receptor_pdb):
    with patch("requests.post", return_value=_mock_nim(_nim_payload(4))):
        result = run_rfdiffusion_pep(mock_receptor_pdb, _HOTSPOTS, nim_api_key="k")
    assert result["n_designs"] == 4


def test_designs_sorted_by_plddt_descending(mock_receptor_pdb):
    payload = {
        "pdbs": ["ATOM  0  CA  ALA\nEND\n"] * 3,
        "plddt_array": [[0.6] * 5, [0.9] * 5, [0.75] * 5],
    }
    with patch("requests.post", return_value=_mock_nim(payload)):
        result = run_rfdiffusion_pep(mock_receptor_pdb, _HOTSPOTS, nim_api_key="k")

    plddts = [d["mean_plddt"] for d in result["designs"]]
    assert plddts == sorted(plddts, reverse=True)


def test_best_pdb_path_is_highest_plddt(mock_receptor_pdb):
    payload = {
        "pdbs": ["ATOM  0  CA  ALA\nEND\n"] * 3,
        "plddt_array": [[0.6] * 5, [0.9] * 5, [0.75] * 5],
    }
    with patch("requests.post", return_value=_mock_nim(payload)):
        result = run_rfdiffusion_pep(mock_receptor_pdb, _HOTSPOTS, nim_api_key="k")

    assert result["best_pdb_path"] == result["designs"][0]["pdb_path"]
    assert result["designs"][0]["mean_plddt"] == pytest.approx(0.9)


def test_pdb_files_are_written(mock_receptor_pdb):
    with patch("requests.post", return_value=_mock_nim(_nim_payload(3))):
        result = run_rfdiffusion_pep(mock_receptor_pdb, _HOTSPOTS, nim_api_key="k")

    for d in result["designs"]:
        assert Path(d["pdb_path"]).exists()


def test_hotspot_res_preserved_in_output(mock_receptor_pdb):
    with patch("requests.post", return_value=_mock_nim(_nim_payload())):
        result = run_rfdiffusion_pep(mock_receptor_pdb, _HOTSPOTS, nim_api_key="k")
    assert result["hotspot_res"] == _HOTSPOTS


def test_binder_length_preserved_in_output(mock_receptor_pdb):
    with patch("requests.post", return_value=_mock_nim(_nim_payload())):
        result = run_rfdiffusion_pep(mock_receptor_pdb, _HOTSPOTS,
                                     binder_length=(12, 18), nim_api_key="k")
    assert result["binder_length"] == (12, 18)


def test_no_plddt_array_in_response(mock_receptor_pdb):
    payload = {"pdbs": ["ATOM  0  CA  ALA\nEND\n"] * 2}
    with patch("requests.post", return_value=_mock_nim(payload)):
        result = run_rfdiffusion_pep(mock_receptor_pdb, _HOTSPOTS, nim_api_key="k")
    assert result["n_designs"] == 2
    assert all(d["mean_plddt"] is None for d in result["designs"])


def test_empty_designs_response(mock_receptor_pdb):
    with patch("requests.post", return_value=_mock_nim({"pdbs": []})):
        result = run_rfdiffusion_pep(mock_receptor_pdb, _HOTSPOTS, nim_api_key="k")
    assert result["n_designs"] == 0
    assert result["best_pdb_path"] is None


# ── HTTP error ────────────────────────────────────────────────────────────────

def test_http_error_propagates(mock_receptor_pdb):
    resp = MagicMock()
    resp.raise_for_status.side_effect = requests.HTTPError("403")
    with patch("requests.post", return_value=resp):
        with pytest.raises(requests.HTTPError):
            run_rfdiffusion_pep(mock_receptor_pdb, _HOTSPOTS, nim_api_key="bad")
