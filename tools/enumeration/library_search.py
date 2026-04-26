from __future__ import annotations

import requests

ENAMINE_API = "https://new.enaminestore.com/api/v1/"
ZINC_API    = "https://zinc22.docking.org/substances/"


def run_library_search(
    smiles: str,
    library: str = "enamine_real",
    mode: str = "similarity",
    tanimoto_threshold: float = 0.7,
    max_results: int = 500,
    fingerprint: str = "ecfp4",
) -> dict:
    from rdkit.Chem import MolFromSmiles

    if MolFromSmiles(smiles) is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    if library == "enamine_real":
        return _enamine_search(smiles, mode, tanimoto_threshold, max_results, fingerprint)
    elif library == "zinc22":
        return _zinc_search(smiles, mode, tanimoto_threshold, max_results)
    else:
        raise ValueError(f"Unknown library: {library}")


def _enamine_search(
    smiles: str,
    mode: str,
    threshold: float,
    max_results: int,
    fp_type: str,
) -> dict:
    fp_map = {"ecfp4": "ECFP4", "fcfp4": "FCFP4"}
    endpoint = "search/sim" if mode == "similarity" else "search/sub"
    params = {
        "query":     smiles,
        "threshold": threshold,
        "limit":     max_results,
        "fp":        fp_map.get(fp_type, "ECFP4"),
        "db":        "REAL",
    }
    resp = requests.get(ENAMINE_API + endpoint, params=params, timeout=60)
    resp.raise_for_status()
    hits = resp.json().get("data", [])
    return {
        "smiles_list":      [h["smiles"] for h in hits],
        "similarities":     [h.get("similarity") for h in hits],
        "catalog_ids":      [h.get("id", "") for h in hits],
        "n_returned":       len(hits),
        "library_searched": "enamine_real",
    }


def _zinc_search(
    smiles: str,
    mode: str,
    threshold: float,
    max_results: int,
) -> dict:
    params: dict = {"smiles": smiles, "count": max_results, "Tanimoto": threshold}
    if mode == "substructure":
        params["search_type"] = "substructure"
    resp = requests.get(ZINC_API + "search.json", params=params, timeout=60)
    resp.raise_for_status()
    hits = resp.json().get("substances", [])
    return {
        "smiles_list":      [h["smiles"] for h in hits],
        "similarities":     [h.get("tanimoto") for h in hits],
        "catalog_ids":      [h.get("zinc_id", "") for h in hits],
        "n_returned":       len(hits),
        "library_searched": "zinc22",
    }
