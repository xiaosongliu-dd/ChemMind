from __future__ import annotations

import time


def run_chembl(
    target_id: str | None = None,
    smiles: str | None = None,
    query_type: str = "bioactivity",
    standard_type: str = "IC50",
    limit: int = 100,
    api_url: str = "https://www.ebi.ac.uk/chembl/api/data",
) -> dict:
    """
    Query the ChEMBL REST API for bioactivity data, target information, or
    compound properties.

    `query_type`:
        "bioactivity" — fetch activity records for a target (requires target_id)
        "target"      — fetch target metadata (requires target_id)
        "compound"    — fetch compound data by SMILES similarity (requires smiles)

    Returns:
        {records: list[dict], n_records, query_type, target_id, runtime_s}

    Raises RuntimeError("ChEMBL API unreachable: ...") on network failure.
    """
    try:
        import requests
    except ImportError as e:
        raise RuntimeError(f"requests not installed (required for ChEMBL): {e}")

    t0 = time.perf_counter()

    if query_type == "bioactivity":
        records = _fetch_bioactivity(requests, api_url, target_id, standard_type, limit)
    elif query_type == "target":
        records = _fetch_target(requests, api_url, target_id)
    elif query_type == "compound":
        records = _fetch_compound(requests, api_url, smiles, limit)
    else:
        raise ValueError(f"Unknown query_type: {query_type!r}. "
                         "Choose 'bioactivity', 'target', or 'compound'.")

    runtime = round(time.perf_counter() - t0, 1)
    return {
        "records":    records,
        "n_records":  len(records),
        "query_type": query_type,
        "target_id":  target_id,
        "runtime_s":  runtime,
    }


def _get_json(requests, url: str, params: dict) -> dict:
    try:
        resp = requests.get(url, params=params, timeout=30)
        resp.raise_for_status()
        return resp.json()
    except requests.RequestException as e:
        raise RuntimeError(f"ChEMBL API unreachable: {e}")


def _fetch_bioactivity(requests, api_url, target_id, standard_type, limit) -> list[dict]:
    if not target_id:
        raise ValueError("target_id is required for query_type='bioactivity'")

    params = {
        "target_chembl_id": target_id,
        "standard_type":    standard_type,
        "limit":            min(limit, 1000),
        "format":           "json",
    }
    data = _get_json(requests, f"{api_url}/activity.json", params)
    activities = data.get("activities", [])
    return [
        {
            "molecule_chembl_id": r.get("molecule_chembl_id"),
            "smiles":             r.get("canonical_smiles"),
            "standard_value":     _safe_float(r.get("standard_value")),
            "standard_units":     r.get("standard_units"),
            "standard_type":      r.get("standard_type"),
            "pchembl_value":      _safe_float(r.get("pchembl_value")),
            "assay_chembl_id":    r.get("assay_chembl_id"),
            "assay_description":  r.get("assay_description"),
        }
        for r in activities
    ]


def _fetch_target(requests, api_url, target_id) -> list[dict]:
    if not target_id:
        raise ValueError("target_id is required for query_type='target'")

    data = _get_json(requests, f"{api_url}/target/{target_id}.json", {})
    t = data
    return [{
        "target_chembl_id":  t.get("target_chembl_id"),
        "pref_name":         t.get("pref_name"),
        "organism":          t.get("organism"),
        "target_type":       t.get("target_type"),
        "components":        t.get("target_components", []),
    }]


def _fetch_compound(requests, api_url, smiles, limit) -> list[dict]:
    if not smiles:
        raise ValueError("smiles is required for query_type='compound'")

    try:
        from urllib.parse import quote
        encoded = quote(smiles, safe="")
    except Exception:
        encoded = smiles

    params = {"smiles": smiles, "limit": min(limit, 1000), "format": "json"}
    data = _get_json(requests, f"{api_url}/molecule.json", params)
    molecules = data.get("molecules", [])
    return [
        {
            "molecule_chembl_id": m.get("molecule_chembl_id"),
            "smiles":             (m.get("molecule_structures") or {}).get("canonical_smiles"),
            "mw":                 _safe_float((m.get("molecule_properties") or {}).get("mw_freebase")),
            "alogp":              _safe_float((m.get("molecule_properties") or {}).get("alogp")),
            "hbd":                _safe_float((m.get("molecule_properties") or {}).get("hbd")),
            "hba":                _safe_float((m.get("molecule_properties") or {}).get("hba")),
            "pref_name":          m.get("pref_name"),
        }
        for m in molecules
    ]


def _safe_float(val) -> float | None:
    try:
        return float(val)
    except (TypeError, ValueError):
        return None
