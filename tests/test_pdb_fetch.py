from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from tools.data.pdb_fetch import run_pdb_fetch


def _mock_get(json_payload=None, text_payload="", status=200):
    resp = MagicMock()
    resp.json.return_value = json_payload or {}
    resp.text = text_payload
    resp.content = text_payload.encode() if isinstance(text_payload, str) else text_payload
    resp.raise_for_status = MagicMock()
    return MagicMock(return_value=resp)


_GRAPHQL_METADATA = {
    "data": {
        "entry": {
            "struct": {"title": "CDK2 complex"},
            "refine": [{"ls_d_res_high": 1.8}],
            "exptl": [{"method": "X-RAY DIFFRACTION"}],
            "polymer_entities": [
                {"entity_poly": {"pdbx_strand_id": "A,B"}}
            ],
            "nonpolymer_entities": [
                {"nonpolymer_comp": {"chem_comp": {"id": "ATP", "name": "Adenosine triphosphate",
                                                   "formula": "C10H16N5O13P3"}}}
            ],
        }
    }
}


# ── structure download ────────────────────────────────────────────────────────

def test_structure_download_returns_pdb_path(tmp_path):
    pdb_content = "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n"
    with patch("requests.get", _mock_get(text_payload=pdb_content)):
        result = run_pdb_fetch(pdb_id="1ATP", query_type="structure", output_dir=str(tmp_path))
    assert "pdb_path" in result
    assert Path(result["pdb_path"]).exists()


def test_structure_download_file_contains_pdb_content(tmp_path):
    pdb_content = "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n"
    with patch("requests.get", _mock_get(text_payload=pdb_content)):
        result = run_pdb_fetch(pdb_id="1ATP", query_type="structure", output_dir=str(tmp_path))
    content = Path(result["pdb_path"]).read_text()
    assert "ATOM" in content


def test_structure_pdb_id_echoed(tmp_path):
    with patch("requests.get", _mock_get(text_payload="END\n")):
        result = run_pdb_fetch(pdb_id="2HHB", query_type="structure", output_dir=str(tmp_path))
    assert result.get("pdb_id", "").upper() == "2HHB"


# ── metadata query ────────────────────────────────────────────────────────────

def test_metadata_returns_title():
    with patch("requests.post", _mock_get(_GRAPHQL_METADATA)):
        result = run_pdb_fetch(pdb_id="1ATP", query_type="metadata")
    assert "title" in result
    assert result["title"] == "CDK2 complex"


def test_metadata_returns_resolution():
    with patch("requests.post", _mock_get(_GRAPHQL_METADATA)):
        result = run_pdb_fetch(pdb_id="1ATP", query_type="metadata")
    assert result.get("resolution") == pytest.approx(1.8)


def test_metadata_returns_ligands():
    with patch("requests.post", _mock_get(_GRAPHQL_METADATA)):
        result = run_pdb_fetch(pdb_id="1ATP", query_type="metadata")
    assert "ligands" in result
    assert any(lig["id"] == "ATP" for lig in result["ligands"])


# ── error handling ────────────────────────────────────────────────────────────

def test_unknown_query_type_raises():
    with pytest.raises(ValueError, match="Unknown query_type"):
        run_pdb_fetch(pdb_id="1ATP", query_type="invalid_mode")


def test_raises_if_requests_not_installed(tmp_path):
    with patch.dict("sys.modules", {"requests": None}):
        with pytest.raises(RuntimeError, match="requests not installed"):
            run_pdb_fetch(pdb_id="1ATP", output_dir=str(tmp_path))


def test_network_error_raises_runtime_error(tmp_path):
    import requests as req
    mock_get = MagicMock()
    mock_get.return_value.raise_for_status.side_effect = req.HTTPError("404")
    with patch("requests.get", mock_get):
        with pytest.raises(RuntimeError, match="RCSB API unreachable"):
            run_pdb_fetch(pdb_id="XXXX", query_type="structure", output_dir=str(tmp_path))
