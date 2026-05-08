from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest
import requests

from tools.biologics.proteinmpnn import run_proteinmpnn


_NIM_RESPONSE = {
    "sequences": [
        {"sequence": "ACDEFGHIKLM", "global_score": 0.92},
        {"sequence": "ACDEFGHIKLN", "global_score": 0.88},
        {"sequence": "ACDEFGHIKLY", "global_score": 0.85},
    ]
}


def _mock_nim(data=None):
    resp = MagicMock()
    resp.json.return_value = data or _NIM_RESPONSE
    resp.raise_for_status = MagicMock()
    return resp


# ── dispatch ──────────────────────────────────────────────────────────────────

def test_dispatches_to_nim_by_default(mock_receptor_pdb):
    with patch("requests.post", return_value=_mock_nim()) as mock_post:
        run_proteinmpnn(mock_receptor_pdb, ["A"], nim_api_key="key")
    mock_post.assert_called_once()


def test_dispatches_to_local_when_nim_false(mock_receptor_pdb):
    with patch("tools.biologics.proteinmpnn._run_local") as mock_local:
        mock_local.return_value = {"sequences": [], "scores": [], "n_returned": 0,
                                   "chains_designed": ["A"], "sampling_temp": 0.1}
        run_proteinmpnn(mock_receptor_pdb, ["A"], use_nim=False)
    mock_local.assert_called_once()


def test_reads_api_key_from_env(mock_receptor_pdb, monkeypatch):
    monkeypatch.setenv("NVIDIA_API_KEY", "env-key")
    with patch("requests.post", return_value=_mock_nim()) as mock_post:
        run_proteinmpnn(mock_receptor_pdb, ["A"])
    headers = mock_post.call_args.kwargs["headers"]
    assert headers["Authorization"] == "Bearer env-key"


# ── NIM payload ───────────────────────────────────────────────────────────────

def test_nim_sends_pdb_content(mock_receptor_pdb):
    with patch("requests.post", return_value=_mock_nim()) as mock_post:
        run_proteinmpnn(mock_receptor_pdb, ["A"], nim_api_key="k")
    payload = mock_post.call_args.kwargs["json"]
    assert "pdb_str" in payload
    assert "ATOM" in payload["pdb_str"]


def test_nim_sends_chains_comma_joined(mock_receptor_pdb):
    with patch("requests.post", return_value=_mock_nim()) as mock_post:
        run_proteinmpnn(mock_receptor_pdb, ["H", "L"], nim_api_key="k")
    payload = mock_post.call_args.kwargs["json"]
    assert payload["chains_to_design"] == "H,L"


def test_nim_passes_n_sequences(mock_receptor_pdb):
    with patch("requests.post", return_value=_mock_nim()) as mock_post:
        run_proteinmpnn(mock_receptor_pdb, ["A"], n_sequences=16, nim_api_key="k")
    payload = mock_post.call_args.kwargs["json"]
    assert payload["num_seq_per_target"] == 16


def test_nim_passes_sampling_temp(mock_receptor_pdb):
    with patch("requests.post", return_value=_mock_nim()) as mock_post:
        run_proteinmpnn(mock_receptor_pdb, ["A"], sampling_temp=0.5, nim_api_key="k")
    payload = mock_post.call_args.kwargs["json"]
    assert payload["sampling_temp"] == "0.5"


def test_nim_sends_fixed_positions_when_set(mock_receptor_pdb):
    fixed = {"A": [10, 20, 30]}
    with patch("requests.post", return_value=_mock_nim()) as mock_post:
        run_proteinmpnn(mock_receptor_pdb, ["A"], fixed_positions=fixed, nim_api_key="k")
    payload = mock_post.call_args.kwargs["json"]
    assert "fixed_positions_jsonl" in payload


def test_nim_no_fixed_positions_key_by_default(mock_receptor_pdb):
    with patch("requests.post", return_value=_mock_nim()) as mock_post:
        run_proteinmpnn(mock_receptor_pdb, ["A"], nim_api_key="k")
    payload = mock_post.call_args.kwargs["json"]
    assert "fixed_positions_jsonl" not in payload


# ── NIM output ────────────────────────────────────────────────────────────────

def test_returns_correct_n_returned(mock_receptor_pdb):
    with patch("requests.post", return_value=_mock_nim()):
        result = run_proteinmpnn(mock_receptor_pdb, ["A"], nim_api_key="k")
    assert result["n_returned"] == 3


def test_sequences_list_matches_n_returned(mock_receptor_pdb):
    with patch("requests.post", return_value=_mock_nim()):
        result = run_proteinmpnn(mock_receptor_pdb, ["A"], nim_api_key="k")
    assert len(result["sequences"]) == result["n_returned"]


def test_scores_list_matches_n_returned(mock_receptor_pdb):
    with patch("requests.post", return_value=_mock_nim()):
        result = run_proteinmpnn(mock_receptor_pdb, ["A"], nim_api_key="k")
    assert len(result["scores"]) == result["n_returned"]


def test_chains_designed_preserved(mock_receptor_pdb):
    with patch("requests.post", return_value=_mock_nim()):
        result = run_proteinmpnn(mock_receptor_pdb, ["H", "L"], nim_api_key="k")
    assert result["chains_designed"] == ["H", "L"]


def test_sampling_temp_in_output(mock_receptor_pdb):
    with patch("requests.post", return_value=_mock_nim()):
        result = run_proteinmpnn(mock_receptor_pdb, ["A"], sampling_temp=0.3, nim_api_key="k")
    assert result["sampling_temp"] == pytest.approx(0.3)


# ── output schema ─────────────────────────────────────────────────────────────

def test_output_has_required_keys(mock_receptor_pdb):
    with patch("requests.post", return_value=_mock_nim()):
        result = run_proteinmpnn(mock_receptor_pdb, ["A"], nim_api_key="k")
    for key in ("sequences", "scores", "n_returned", "chains_designed", "sampling_temp"):
        assert key in result


def test_empty_sequences_response(mock_receptor_pdb):
    with patch("requests.post", return_value=_mock_nim({"sequences": []})):
        result = run_proteinmpnn(mock_receptor_pdb, ["A"], nim_api_key="k")
    assert result["n_returned"] == 0
    assert result["sequences"] == []


# ── HTTP errors ───────────────────────────────────────────────────────────────

def test_http_error_propagates(mock_receptor_pdb):
    resp = MagicMock()
    resp.raise_for_status.side_effect = requests.HTTPError("401")
    with patch("requests.post", return_value=resp):
        with pytest.raises(requests.HTTPError):
            run_proteinmpnn(mock_receptor_pdb, ["A"], nim_api_key="bad-key")
