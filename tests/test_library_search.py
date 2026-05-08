from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest
import requests

from tools.enumeration.library_search import (
    _enamine_search,
    _zinc_search,
    run_library_search,
)


def _enamine_hits(n=3):
    return {"data": [{"smiles": f"C{'C'*i}", "similarity": 0.9 - i*0.05, "id": f"EN{i}"} for i in range(n)]}


def _zinc_hits(n=3):
    return {"substances": [{"smiles": f"C{'C'*i}", "tanimoto": 0.85 - i*0.05, "zinc_id": f"ZINC{i}"} for i in range(n)]}


# ── dispatch ──────────────────────────────────────────────────────────────────

def test_dispatches_to_enamine_by_default(imatinib, mock_nim_response):
    with patch("tools.enumeration.library_search._enamine_search") as mock_en:
        mock_en.return_value = {"smiles_list": [], "similarities": [],
                                "catalog_ids": [], "n_returned": 0, "library_searched": "enamine_real"}
        run_library_search(imatinib)
    mock_en.assert_called_once()


def test_dispatches_to_zinc(imatinib):
    with patch("tools.enumeration.library_search._zinc_search") as mock_zinc:
        mock_zinc.return_value = {"smiles_list": [], "similarities": [],
                                  "catalog_ids": [], "n_returned": 0, "library_searched": "zinc22"}
        run_library_search(imatinib, library="zinc22")
    mock_zinc.assert_called_once()


def test_unknown_library_raises(imatinib):
    with pytest.raises(ValueError, match="Unknown library"):
        run_library_search(imatinib, library="unknown_db")


def test_invalid_smiles_raises():
    with pytest.raises(ValueError, match="Invalid SMILES"):
        run_library_search("not_a_smiles!!!")


# ── Enamine API ───────────────────────────────────────────────────────────────

def test_enamine_similarity_uses_correct_endpoint(imatinib):
    resp = MagicMock()
    resp.json.return_value = _enamine_hits()
    resp.raise_for_status = MagicMock()
    with patch("requests.get", return_value=resp) as mock_get:
        _enamine_search(imatinib, "similarity", 0.7, 10, "ecfp4")
    url = mock_get.call_args[0][0]
    assert "search/sim" in url


def test_enamine_substructure_uses_correct_endpoint(imatinib):
    resp = MagicMock()
    resp.json.return_value = _enamine_hits()
    resp.raise_for_status = MagicMock()
    with patch("requests.get", return_value=resp) as mock_get:
        _enamine_search(imatinib, "substructure", 0.7, 10, "ecfp4")
    url = mock_get.call_args[0][0]
    assert "search/sub" in url


def test_enamine_passes_threshold(imatinib):
    resp = MagicMock()
    resp.json.return_value = _enamine_hits()
    resp.raise_for_status = MagicMock()
    with patch("requests.get", return_value=resp) as mock_get:
        _enamine_search(imatinib, "similarity", 0.85, 10, "ecfp4")
    params = mock_get.call_args.kwargs["params"]
    assert params["threshold"] == 0.85


def test_enamine_maps_ecfp4_fingerprint(imatinib):
    resp = MagicMock()
    resp.json.return_value = _enamine_hits()
    resp.raise_for_status = MagicMock()
    with patch("requests.get", return_value=resp) as mock_get:
        _enamine_search(imatinib, "similarity", 0.7, 10, "ecfp4")
    params = mock_get.call_args.kwargs["params"]
    assert params["fp"] == "ECFP4"


def test_enamine_returns_correct_structure(imatinib):
    resp = MagicMock()
    resp.json.return_value = _enamine_hits(4)
    resp.raise_for_status = MagicMock()
    with patch("requests.get", return_value=resp):
        result = _enamine_search(imatinib, "similarity", 0.7, 100, "ecfp4")

    assert result["n_returned"] == 4
    assert result["library_searched"] == "enamine_real"
    assert len(result["smiles_list"]) == 4
    assert len(result["catalog_ids"]) == 4


def test_enamine_http_error_propagates(imatinib):
    resp = MagicMock()
    resp.raise_for_status.side_effect = requests.HTTPError("500")
    with patch("requests.get", return_value=resp):
        with pytest.raises(requests.HTTPError):
            _enamine_search(imatinib, "similarity", 0.7, 10, "ecfp4")


# ── ZINC API ──────────────────────────────────────────────────────────────────

def test_zinc_passes_smiles(imatinib):
    resp = MagicMock()
    resp.json.return_value = _zinc_hits()
    resp.raise_for_status = MagicMock()
    with patch("requests.get", return_value=resp) as mock_get:
        _zinc_search(imatinib, "similarity", 0.7, 10)
    params = mock_get.call_args.kwargs["params"]
    assert params["smiles"] == imatinib


def test_zinc_substructure_sets_search_type(imatinib):
    resp = MagicMock()
    resp.json.return_value = _zinc_hits()
    resp.raise_for_status = MagicMock()
    with patch("requests.get", return_value=resp) as mock_get:
        _zinc_search(imatinib, "substructure", 0.7, 10)
    params = mock_get.call_args.kwargs["params"]
    assert params.get("search_type") == "substructure"


def test_zinc_similarity_no_search_type(imatinib):
    resp = MagicMock()
    resp.json.return_value = _zinc_hits()
    resp.raise_for_status = MagicMock()
    with patch("requests.get", return_value=resp) as mock_get:
        _zinc_search(imatinib, "similarity", 0.7, 10)
    params = mock_get.call_args.kwargs["params"]
    assert "search_type" not in params


def test_zinc_returns_correct_structure(imatinib):
    resp = MagicMock()
    resp.json.return_value = _zinc_hits(5)
    resp.raise_for_status = MagicMock()
    with patch("requests.get", return_value=resp):
        result = _zinc_search(imatinib, "similarity", 0.7, 100)

    assert result["n_returned"] == 5
    assert result["library_searched"] == "zinc22"
    assert all(z.startswith("ZINC") for z in result["catalog_ids"])
