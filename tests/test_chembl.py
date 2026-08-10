from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from tools.data.chembl import _safe_float, run_chembl


def _mock_get(payload: dict):
    resp = MagicMock()
    resp.json.return_value = payload
    resp.raise_for_status = MagicMock()
    return MagicMock(return_value=resp)


_BIOACTIVITY_RESPONSE = {
    "activities": [
        {
            "molecule_chembl_id": "CHEMBL1234",
            "canonical_smiles":   "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1",
            "standard_value":     "50.0",
            "standard_units":     "nM",
            "standard_type":      "IC50",
            "pchembl_value":      "7.3",
            "assay_chembl_id":    "CHEMBL99",
            "assay_description":  "Inhibition of CDK2",
        },
        {
            "molecule_chembl_id": "CHEMBL5678",
            "canonical_smiles":   "CCOc1ccc2nccc(Oc3ccc(NC4=CC(=O)N(C)C4)cc3)c2c1",
            "standard_value":     "120.0",
            "standard_units":     "nM",
            "standard_type":      "IC50",
            "pchembl_value":      "6.9",
            "assay_chembl_id":    "CHEMBL99",
            "assay_description":  "Inhibition of CDK2",
        },
    ]
}

_TARGET_RESPONSE = {
    "target_chembl_id": "CHEMBL301",
    "pref_name":        "Cyclin-dependent kinase 2",
    "organism":         "Homo sapiens",
    "target_type":      "SINGLE PROTEIN",
    "target_components": [],
}


# ── _safe_float ───────────────────────────────────────────────────────────────

def test_safe_float_string():
    assert _safe_float("7.3") == pytest.approx(7.3)

def test_safe_float_none():
    assert _safe_float(None) is None

def test_safe_float_bad_string():
    assert _safe_float("N/A") is None


# ── bioactivity query ─────────────────────────────────────────────────────────

def test_bioactivity_returns_records():
    with patch("requests.get", _mock_get(_BIOACTIVITY_RESPONSE)):
        result = run_chembl(target_id="CHEMBL301", standard_type="IC50")
    assert result["n_records"] == 2


def test_bioactivity_schema():
    with patch("requests.get", _mock_get(_BIOACTIVITY_RESPONSE)):
        result = run_chembl(target_id="CHEMBL301")
    for key in ("records", "n_records", "query_type", "target_id", "runtime_s"):
        assert key in result


def test_bioactivity_record_has_smiles():
    with patch("requests.get", _mock_get(_BIOACTIVITY_RESPONSE)):
        result = run_chembl(target_id="CHEMBL301")
    assert result["records"][0]["smiles"] is not None


def test_bioactivity_pchembl_value_is_float():
    with patch("requests.get", _mock_get(_BIOACTIVITY_RESPONSE)):
        result = run_chembl(target_id="CHEMBL301")
    assert isinstance(result["records"][0]["pchembl_value"], float)


def test_bioactivity_passes_target_id():
    with patch("requests.get", _mock_get(_BIOACTIVITY_RESPONSE)) as mock_get:
        run_chembl(target_id="CHEMBL301", standard_type="IC50")
    params = mock_get.return_value.json.call_args  # irrelevant — check call args
    called_params = mock_get.call_args.kwargs.get("params", {})
    assert called_params.get("target_chembl_id") == "CHEMBL301"


# ── target query ──────────────────────────────────────────────────────────────

def test_target_query_returns_pref_name():
    with patch("requests.get", _mock_get(_TARGET_RESPONSE)):
        result = run_chembl(target_id="CHEMBL301", query_type="target")
    assert result["records"][0]["pref_name"] == "Cyclin-dependent kinase 2"


def test_target_query_returns_organism():
    with patch("requests.get", _mock_get(_TARGET_RESPONSE)):
        result = run_chembl(target_id="CHEMBL301", query_type="target")
    assert result["records"][0]["organism"] == "Homo sapiens"


# ── compound query ────────────────────────────────────────────────────────────

def test_compound_query_schema(imatinib):
    compound_resp = {"molecules": [
        {
            "molecule_chembl_id": "CHEMBL941",
            "molecule_structures": {"canonical_smiles": imatinib},
            "molecule_properties": {"mw_freebase": "493.6", "alogp": "3.7",
                                    "hbd": "3", "hba": "9"},
            "pref_name": "IMATINIB",
        }
    ]}
    with patch("requests.get", _mock_get(compound_resp)):
        result = run_chembl(smiles=imatinib, query_type="compound")
    assert result["records"][0]["smiles"] == imatinib


# ── error handling ────────────────────────────────────────────────────────────

def test_unknown_query_type_raises():
    with pytest.raises(ValueError, match="Unknown query_type"):
        run_chembl(target_id="CHEMBL301", query_type="magic")


def test_bioactivity_requires_target_id():
    with pytest.raises(ValueError, match="target_id is required"):
        run_chembl(query_type="bioactivity")


def test_raises_if_requests_not_installed():
    with patch.dict("sys.modules", {"requests": None}):
        with pytest.raises(RuntimeError, match="requests not installed"):
            run_chembl(target_id="CHEMBL301")


def test_network_error_raises_runtime_error():
    import requests as req
    mock = MagicMock()
    mock.return_value.raise_for_status.side_effect = req.RequestException("timeout")
    with patch("requests.get", mock):
        with pytest.raises(RuntimeError, match="ChEMBL API unreachable"):
            run_chembl(target_id="CHEMBL301")
