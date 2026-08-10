from __future__ import annotations

import os
from unittest.mock import MagicMock, patch

import pytest

from tools.data.askcos import run_askcos


def _mock_post(payload: dict, status: int = 200):
    resp = MagicMock()
    resp.json.return_value = payload
    resp.raise_for_status = MagicMock()
    resp.status_code = status
    return MagicMock(return_value=resp)


_ROUTES_RESPONSE = {
    "result": [
        {
            "steps": [
                {"smiles": "CCO", "reaction": "esterification"},
                {"smiles": "CC(=O)Cl", "reaction": "acyl_chloride_formation"},
            ],
            "buyable_leaves": ["CCO", "CC(=O)Cl"],
            "overall_score": 0.85,
        },
        {
            "steps": [
                {"smiles": "CCN", "reaction": "amination"},
            ],
            "buyable_leaves": ["CCN"],
            "overall_score": 0.72,
        },
    ]
}


# ── basic dispatch ────────────────────────────────────────────────────────────

def test_returns_routes(imatinib):
    with patch.dict(os.environ, {"ASKCOS_API_URL": "http://askcos.local"}):
        with patch("requests.post", _mock_post(_ROUTES_RESPONSE)):
            result = run_askcos(imatinib)
    assert "routes" in result
    assert len(result["routes"]) == 2


def test_output_schema(imatinib):
    with patch.dict(os.environ, {"ASKCOS_API_URL": "http://askcos.local"}):
        with patch("requests.post", _mock_post(_ROUTES_RESPONSE)):
            result = run_askcos(imatinib)
    for key in ("routes", "n_routes", "target_smiles", "runtime_s"):
        assert key in result


def test_n_routes_matches(imatinib):
    with patch.dict(os.environ, {"ASKCOS_API_URL": "http://askcos.local"}):
        with patch("requests.post", _mock_post(_ROUTES_RESPONSE)):
            result = run_askcos(imatinib)
    assert result["n_routes"] == len(result["routes"])


def test_target_smiles_echoed(imatinib):
    with patch.dict(os.environ, {"ASKCOS_API_URL": "http://askcos.local"}):
        with patch("requests.post", _mock_post(_ROUTES_RESPONSE)):
            result = run_askcos(imatinib)
    assert result["target_smiles"] == imatinib


def test_routes_have_required_fields(imatinib):
    with patch.dict(os.environ, {"ASKCOS_API_URL": "http://askcos.local"}):
        with patch("requests.post", _mock_post(_ROUTES_RESPONSE)):
            result = run_askcos(imatinib)
    for route in result["routes"]:
        assert "steps" in route
        assert "overall_score" in route


# ── n_steps parameter ─────────────────────────────────────────────────────────

def test_n_steps_passed_in_payload(imatinib):
    with patch.dict(os.environ, {"ASKCOS_API_URL": "http://askcos.local"}):
        with patch("requests.post", _mock_post(_ROUTES_RESPONSE)) as mock_post:
            run_askcos(imatinib, n_steps=5)
    payload = mock_post.call_args.kwargs.get("json") or mock_post.call_args[1].get("json", {})
    assert payload.get("max_depth") == 5


# ── api_key ───────────────────────────────────────────────────────────────────

def test_api_key_added_to_header(imatinib):
    with patch.dict(os.environ, {"ASKCOS_API_URL": "http://askcos.local"}):
        with patch("requests.post", _mock_post(_ROUTES_RESPONSE)) as mock_post:
            run_askcos(imatinib, api_key="mysecret")
    headers = mock_post.call_args.kwargs.get("headers", {})
    assert "mysecret" in headers.get("Authorization", "")


# ── error handling ────────────────────────────────────────────────────────────

def test_no_api_url_raises(imatinib):
    env = {k: v for k, v in os.environ.items() if k != "ASKCOS_API_URL"}
    with patch.dict(os.environ, env, clear=True):
        with pytest.raises(ValueError, match="No ASKCOS API URL"):
            run_askcos(imatinib)


def test_raises_if_requests_not_installed(imatinib):
    with patch.dict("sys.modules", {"requests": None}):
        with pytest.raises(RuntimeError, match="requests not installed"):
            run_askcos(imatinib, api_url="http://askcos.local")


def test_network_error_raises_runtime(imatinib):
    import requests as req
    mock = MagicMock()
    mock.return_value.raise_for_status.side_effect = req.RequestException("timeout")
    with patch("requests.post", mock):
        with pytest.raises(RuntimeError, match="ASKCOS API unreachable"):
            run_askcos(imatinib, api_url="http://askcos.local")
