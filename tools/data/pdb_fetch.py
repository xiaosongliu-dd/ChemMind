from __future__ import annotations

import json
import tempfile
import time
from pathlib import Path


_RCSB_DOWNLOAD = "https://files.rcsb.org/download"
_RCSB_SEARCH   = "https://search.rcsb.org/rcsbsearch/v2/query"
_RCSB_GRAPHQL  = "https://data.rcsb.org/graphql"


def run_pdb_fetch(
    pdb_id: str | None = None,
    uniprot_id: str | None = None,
    query_type: str = "structure",
    download_format: str = "pdb",
    output_dir: str | None = None,
) -> dict:
    """
    Fetch PDB structures and metadata from the RCSB Protein Data Bank.

    `query_type`:
        "structure" — download a PDB/mmCIF file by PDB ID (requires pdb_id)
        "search"    — resolve UniProt ID → list of PDB IDs (requires uniprot_id)
        "metadata"  — fetch resolution, experiment, ligands for a PDB ID

    `download_format`: "pdb" (default) or "cif"

    Returns:
        {pdb_path?, pdb_id?, title?, resolution?, experiment_type?,
         chain_ids?, ligands?, uniprot_ids?, pdb_ids?}

    Raises RuntimeError("RCSB API unreachable: ...") on network failure.
    """
    try:
        import requests
    except ImportError as e:
        raise RuntimeError(f"requests not installed (required for pdb_fetch): {e}")

    t0 = time.perf_counter()
    if query_type == "structure":
        result = _fetch_structure(requests, pdb_id, download_format, output_dir)
    elif query_type == "search":
        result = _search_by_uniprot(requests, uniprot_id or pdb_id)
    elif query_type == "metadata":
        result = _fetch_metadata(requests, pdb_id)
    else:
        raise ValueError(f"Unknown query_type: {query_type!r}. "
                         "Choose 'structure', 'search', or 'metadata'.")

    result["runtime_s"] = round(time.perf_counter() - t0, 1)
    return result


def _fetch_structure(requests, pdb_id: str | None, fmt: str, output_dir) -> dict:
    if not pdb_id:
        raise ValueError("pdb_id is required for query_type='structure'")

    pdb_id = pdb_id.upper()
    ext = "cif" if fmt == "cif" else "pdb"
    url = f"{_RCSB_DOWNLOAD}/{pdb_id}.{ext}"

    try:
        resp = requests.get(url, timeout=60)
        resp.raise_for_status()
    except requests.RequestException as e:
        raise RuntimeError(f"RCSB API unreachable: {e}")

    outdir = Path(output_dir) if output_dir else Path(tempfile.mkdtemp(prefix="pdb_"))
    outdir.mkdir(parents=True, exist_ok=True)
    pdb_path = outdir / f"{pdb_id}.{ext}"
    pdb_path.write_bytes(resp.content)

    meta = _fetch_metadata(requests, pdb_id)
    meta["pdb_path"] = str(pdb_path)
    return meta


def _search_by_uniprot(requests, uniprot_id: str | None) -> dict:
    if not uniprot_id:
        raise ValueError("uniprot_id is required for query_type='search'")

    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers.uniprot_ids",
                "operator": "in",
                "value": [uniprot_id.upper()],
            },
        },
        "return_type": "entry",
        "request_options": {"return_all_hits": True},
    }

    try:
        resp = requests.post(
            _RCSB_SEARCH,
            headers={"Content-Type": "application/json"},
            data=json.dumps(query),
            timeout=30,
        )
        resp.raise_for_status()
    except requests.RequestException as e:
        raise RuntimeError(f"RCSB API unreachable: {e}")

    data = resp.json()
    pdb_ids = [r["identifier"] for r in data.get("result_set", [])]
    return {
        "uniprot_id": uniprot_id,
        "pdb_ids":    pdb_ids,
        "n_entries":  len(pdb_ids),
    }


def _fetch_metadata(requests, pdb_id: str | None) -> dict:
    if not pdb_id:
        return {}

    pdb_id = pdb_id.upper()
    graphql_query = """
    query($id: String!) {
      entry(entry_id: $id) {
        struct { title }
        refine { ls_d_res_high }
        exptl { method }
        polymer_entities {
          rcsb_polymer_entity_container_identifiers { uniprot_ids }
          entity_poly { pdbx_strand_id }
        }
        nonpolymer_entities {
          nonpolymer_comp { chem_comp { id name formula } }
        }
      }
    }
    """
    try:
        resp = requests.post(
            _RCSB_GRAPHQL,
            json={"query": graphql_query, "variables": {"id": pdb_id}},
            timeout=15,
        )
        resp.raise_for_status()
        data = resp.json().get("data", {}).get("entry") or {}
    except requests.RequestException as e:
        raise RuntimeError(f"RCSB API unreachable: {e}")

    chain_ids = []
    uniprot_ids = []
    for pe in data.get("polymer_entities") or []:
        ids = (pe.get("rcsb_polymer_entity_container_identifiers") or {})
        uniprot_ids.extend(ids.get("uniprot_ids") or [])
        strand = (pe.get("entity_poly") or {}).get("pdbx_strand_id") or ""
        chain_ids.extend([c.strip() for c in strand.split(",") if c.strip()])

    ligands = []
    for ne in data.get("nonpolymer_entities") or []:
        comp = ((ne.get("nonpolymer_comp") or {}).get("chem_comp") or {})
        if comp:
            ligands.append({
                "id":      comp.get("id"),
                "name":    comp.get("name"),
                "formula": comp.get("formula"),
            })

    refine = data.get("refine") or [{}]
    resolution = (refine[0] if refine else {}).get("ls_d_res_high")
    exptl = data.get("exptl") or [{}]
    method = (exptl[0] if exptl else {}).get("method")
    title = (data.get("struct") or {}).get("title")

    return {
        "pdb_id":          pdb_id,
        "title":           title,
        "resolution":      _safe_float(resolution),
        "experiment_type": method,
        "chain_ids":       chain_ids,
        "ligands":         ligands,
        "uniprot_ids":     list(set(uniprot_ids)),
    }


def _safe_float(val) -> float | None:
    try:
        return float(val)
    except (TypeError, ValueError):
        return None
