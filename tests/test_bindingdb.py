from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from tools.data.bindingdb import _parse, _safe_float, run_bindingdb


# ── helpers ───────────────────────────────────────────────────────────────────

def _bdb_response(n: int = 3) -> dict:
    return {
        "affinities": [
            {
                "ligandSmiles": f"C{'C'*i}N",
                "target":       "CDK2",
                "affinity":     str(100 * (i + 1)),
                "affinity_type":"IC50",
                "uniprotID":    "P24941",
                "pmid":         f"2023{i}",
                "ligandName":   f"compound_{i}",
            }
            for i in range(n)
        ]
    }


def _mock_get(payload: dict):
    resp = MagicMock()
    resp.json.return_value = payload
    resp.raise_for_status = MagicMock()
    return MagicMock(return_value=resp)


# ── _safe_float ───────────────────────────────────────────────────────────────

def test_safe_float_numeric_string():
    assert _safe_float("123.4") == pytest.approx(123.4)

def test_safe_float_none_returns_none():
    assert _safe_float(None) is None

def test_safe_float_non_numeric_returns_none():
    assert _safe_float("not_a_number") is None


# ── _parse ────────────────────────────────────────────────────────────────────

def test_parse_returns_correct_count():
    records = _parse(_bdb_response(4), assay_type="IC50", max_results=100)
    assert len(records) == 4

def test_parse_filters_by_assay_type():
    data = {"affinities": [
        {"affinity_type": "IC50", "affinity": "50", "ligandSmiles": "CCO",
         "target": "T", "uniprotID": None, "pmid": None, "ligandName": None},
        {"affinity_type": "Ki",   "affinity": "30", "ligandSmiles": "CCC",
         "target": "T", "uniprotID": None, "pmid": None, "ligandName": None},
    ]}
    records = _parse(data, assay_type="IC50", max_results=100)
    assert len(records) == 1
    assert records[0]["affinity_type"] == "IC50"

def test_parse_respects_max_results():
    records = _parse(_bdb_response(10), assay_type=None, max_results=3)
    assert len(records) == 3

def test_parse_empty_affinities():
    records = _parse({}, assay_type="IC50", max_results=100)
    assert records == []


# ── run_bindingdb ─────────────────────────────────────────────────────────────

def test_by_target_returns_records():
    with patch("requests.get", _mock_get(_bdb_response(3))):
        result = run_bindingdb(target_name="CDK2")
    assert result["n_records"] == 3


def test_by_uniprot_returns_records():
    with patch("requests.get", _mock_get(_bdb_response(2))):
        result = run_bindingdb(uniprot_id="P24941")
    assert result["n_records"] == 2


def test_by_smiles_returns_records():
    with patch("requests.get", _mock_get(_bdb_response(1))):
        result = run_bindingdb(smiles="CCO")
    assert result["n_records"] == 1


def test_output_schema():
    with patch("requests.get", _mock_get(_bdb_response(3))):
        result = run_bindingdb(target_name="CDK2")
    for key in ("records", "n_records", "query", "runtime_s"):
        assert key in result


def test_record_schema():
    with patch("requests.get", _mock_get(_bdb_response(1))):
        result = run_bindingdb(target_name="CDK2")
    rec = result["records"][0]
    for key in ("ligand_smiles", "target_name", "affinity_nm", "affinity_type",
                "uniprot_id", "pubmed_id", "ligand_name"):
        assert key in rec


def test_sorted_most_potent_first():
    with patch("requests.get", _mock_get(_bdb_response(3))):
        result = run_bindingdb(target_name="CDK2")
    affs = [r["affinity_nm"] for r in result["records"]]
    assert affs == sorted(affs)


def test_min_affinity_filter():
    # only keep records with affinity_nm <= 150 nM
    with patch("requests.get", _mock_get(_bdb_response(3))):
        result = run_bindingdb(target_name="CDK2", min_affinity_nm=150)
    # bdb_response generates: 100, 200, 300 nM
    for r in result["records"]:
        assert r["affinity_nm"] <= 150


def test_no_input_raises():
    with pytest.raises(ValueError, match="Provide at least one"):
        run_bindingdb()


def test_raises_if_requests_not_installed():
    with patch.dict("sys.modules", {"requests": None}):
        with pytest.raises(RuntimeError, match="requests not installed"):
            run_bindingdb(target_name="CDK2")


def test_query_dict_echoed():
    with patch("requests.get", _mock_get(_bdb_response(1))):
        result = run_bindingdb(target_name="CDK2", assay_type="Ki")
    assert result["query"]["target_name"] == "CDK2"
    assert result["query"]["assay_type"] == "Ki"
