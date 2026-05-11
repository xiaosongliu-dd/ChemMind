# Skill: pdb_fetch

**Tier:** L1 — atomic tool skill
**Category:** data
**Tool:** RCSB PDB REST API — structure download and metadata retrieval
**GPU required:** No — REST API call; no compute

---

## What this skill does

`pdb_fetch` downloads protein structures and fetches metadata from the RCSB Protein Data Bank via its public REST and GraphQL APIs. It can resolve a UniProt accession to all matching PDB entries, download the full PDB/mmCIF file, and return key metadata (resolution, ligands, chain IDs, experiment type) in one call.

Key capabilities:
- Download PDB or mmCIF by four-letter PDB ID
- Resolve UniProt accession → list of deposited PDB IDs (SIFTS mapping)
- Retrieve structured metadata: resolution, experiment type, ligands, chains, UniProt IDs
- Three modes: `structure`, `search`, `metadata`

---

## When to call this skill

- Fetching the target structure before docking or pocket detection
- Looking up which PDB entries exist for a UniProt target
- Getting ligand information from a co-crystal structure
- Retrieving resolution to filter high-quality structures

Do NOT use if:
- Need predicted structures → use `esmfold` or `colabfold`
- Need bioactivity data → use `data/chembl`

---

## Inputs

| Field | Type | Required | Description |
|---|---|---|---|
| `pdb_id` | str | cond. | Four-letter PDB ID, e.g. `"2HYY"` |
| `uniprot_id` | str | cond. | UniProt accession, e.g. `"P00519"` (ABL1) |
| `query_type` | str | — | `"structure"` (default), `"search"`, or `"metadata"` |
| `download_format` | str | — | `"pdb"` (default) or `"cif"` |
| `output_dir` | str | — | Download directory; temp dir if None |

`pdb_id` required for `"structure"` and `"metadata"`; `uniprot_id` (or `pdb_id`) for `"search"`.

---

## Outputs

| Field | Type | Description |
|---|---|---|
| `pdb_path` | str \| None | Path to downloaded structure file (`structure` mode only) |
| `pdb_id` | str | PDB accession |
| `title` | str \| None | Structure title |
| `resolution` | float \| None | Crystallographic resolution (Å) |
| `experiment_type` | str \| None | `"X-RAY DIFFRACTION"`, `"CRYO-EM"`, etc. |
| `chain_ids` | list[str] | Polymer chain identifiers |
| `ligands` | list[dict] | `{id, name, formula}` for non-polymer entities |
| `uniprot_ids` | list[str] | Mapped UniProt accessions |
| `pdb_ids` | list[str] | (`search` mode) All PDB IDs for a UniProt ID |
| `runtime_s` | float | API call duration |

Raises `RuntimeError("RCSB API unreachable: ...")` on network failure.

---

## Implementation

`tools/data/pdb_fetch.py` — `run_pdb_fetch()`

`structure` mode: downloads `https://files.rcsb.org/download/{pdb_id}.pdb`.
`search` mode: POSTs to the RCSB Search API with a UniProt attribute query.
`metadata` mode: POSTs a GraphQL query to `https://data.rcsb.org/graphql`.

---

## Agent decision rules

- **Structure selection**: prefer `resolution < 2.5 Å` and `experiment_type == "X-RAY DIFFRACTION"` for docking prep; CryoEM at `< 3.5 Å` is acceptable
- **Ligand check**: inspect `ligands` to confirm whether a co-crystal ligand defines the binding site; if present, use its centre for docking box instead of running `pocket_detect`
- **UniProt → PDB workflow**: `search` mode returns multiple PDB IDs; pick the highest-resolution structure to download
- **Chain selection**: for multi-chain structures, check `chain_ids` and select the relevant chain for pocket detection
- **Upstream**: ChEMBL target metadata → UniProt ID → pdb_fetch search → best structure
- **Downstream**: `pocket_detect`, `docking/*`, `structure/colabfold` (as template)

---

## Install

```bash
# No install needed — pure requests
uv add requests  # if not already installed
```

---

## References

- Berman H.M. et al. *Nucleic Acids Res.* 2000, 28, 235 — "The Protein Data Bank"
- RCSB REST API: https://www.rcsb.org/docs/programmatic-access/web-services-overview
- RCSB Search API: https://search.rcsb.org/
- License: Open access (CC0 equivalent)
