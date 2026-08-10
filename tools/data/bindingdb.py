from __future__ import annotations

import time
from typing import Any


_BDB_BASE = "https://bindingdb.org/axis2/services/BDBService"


def run_bindingdb(
    target_name: str | None = None,
    smiles: str | None = None,
    uniprot_id: str | None = None,
    *,
    assay_type: str = "IC50",
    max_results: int = 100,
    min_affinity_nm: float | None = None,
    max_affinity_nm: float | None = None,
) -> dict[str, Any]:
    """
    Query the BindingDB REST API for binding affinity data.

    Exactly one of target_name, uniprot_id, or smiles must be provided.

    Parameters
    ----------
    target_name : str | None
        Target protein name (e.g. "CDK2", "EGFR").
    smiles : str | None
        Query compound SMILES for similarity-based lookup.
    uniprot_id : str | None
        UniProt accession (e.g. "P24941" for human CDK2).
    assay_type : str
        Filter by assay type: "IC50", "Ki", "Kd", "EC50".
    max_results : int
        Maximum entries to return (capped at 500).
    min_affinity_nm : float | None
        Only keep records with affinity ≤ this value in nM (e.g. 1000 = ≤ 1 µM).
    max_affinity_nm : float | None
        Only keep records with affinity ≥ this value in nM (upper bound filter).

    Returns
    -------
    {
        records   : list[{ligand_smiles, target_name, affinity_nm, affinity_type,
                          uniprot_id, pubmed_id, ligand_name}]
        n_records : int
        query     : {target_name, smiles, uniprot_id, assay_type}
        runtime_s : float
    }
    """
    try:
        import requests
    except ImportError:
        raise RuntimeError("requests not installed: pip install requests")

    if not any([target_name, smiles, uniprot_id]):
        raise ValueError("Provide at least one of: target_name, smiles, uniprot_id.")

    t0 = time.perf_counter()
    max_results = min(max_results, 500)

    if uniprot_id:
        records = _by_uniprot(requests, uniprot_id, assay_type, max_results)
    elif target_name:
        records = _by_target(requests, target_name, assay_type, max_results)
    else:
        records = _by_smiles(requests, smiles, max_results)  # type: ignore[arg-type]

    # Affinity filters
    if min_affinity_nm is not None:
        records = [r for r in records
                   if r["affinity_nm"] is not None and r["affinity_nm"] <= min_affinity_nm]
    if max_affinity_nm is not None:
        records = [r for r in records
                   if r["affinity_nm"] is not None and r["affinity_nm"] >= max_affinity_nm]

    # Sort most-potent first
    records.sort(key=lambda r: (r["affinity_nm"] is None, r["affinity_nm"] or 0))

    return {
        "records":   records,
        "n_records": len(records),
        "query": {
            "target_name": target_name,
            "smiles":      smiles,
            "uniprot_id":  uniprot_id,
            "assay_type":  assay_type,
        },
        "runtime_s": round(time.perf_counter() - t0, 2),
    }


def _by_uniprot(requests, uniprot_id: str, assay_type: str, max_results: int) -> list[dict]:
    params = {"uniprot": uniprot_id, "cutoff": "10000000", "response": "json"}
    data = _get(requests, f"{_BDB_BASE}/getLigandsByUniprots", params)
    return _parse(data, assay_type, max_results)


def _by_target(requests, target_name: str, assay_type: str, max_results: int) -> list[dict]:
    params = {"targetname": target_name, "cutoff": "10000000", "response": "json"}
    data = _get(requests, f"{_BDB_BASE}/getLigandsByTargetName", params)
    return _parse(data, assay_type, max_results)


def _by_smiles(requests, smiles: str, max_results: int) -> list[dict]:
    params = {"smiles": smiles, "cutoff": "10000000", "response": "json"}
    data = _get(requests, f"{_BDB_BASE}/getLigandsBySmiles", params)
    return _parse(data, assay_type=None, max_results=max_results)


def _get(requests, url: str, params: dict) -> dict:
    try:
        resp = requests.get(url, params=params, timeout=30)
        resp.raise_for_status()
        return resp.json()
    except Exception as e:
        raise RuntimeError(f"BindingDB API unreachable: {e}")


def _parse(data: dict, assay_type: str | None, max_results: int) -> list[dict]:
    affinities = data.get("affinities") or []
    records = []
    for entry in affinities:
        atype = entry.get("affinity_type", "")
        if assay_type and assay_type.upper() not in atype.upper():
            continue
        records.append({
            "ligand_smiles": entry.get("ligandSmiles"),
            "target_name":   entry.get("target"),
            "affinity_nm":   _safe_float(entry.get("affinity")),
            "affinity_type": atype or assay_type,
            "uniprot_id":    entry.get("uniprotID"),
            "pubmed_id":     entry.get("pmid"),
            "ligand_name":   entry.get("ligandName"),
        })
        if len(records) >= max_results:
            break
    return records


def _safe_float(val: Any) -> float | None:
    try:
        return float(val)
    except (TypeError, ValueError):
        return None
