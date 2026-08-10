# Skill: chembl

**Tier:** L1 — atomic tool skill
**Category:** data
**Tool:** ChEMBL REST API — bioactivity, target, and compound database
**GPU required:** No — REST API call; no compute

---

## What this skill does

`chembl` queries the European Bioinformatics Institute's ChEMBL database via its public REST API. ChEMBL contains >20 million bioactivity measurements for >2 million compounds against >15,000 targets. Use it to retrieve SAR data, identify known actives, or look up target metadata at the start of a drug design campaign.

Key capabilities:
- Fetch IC50/Ki/EC50 measurements for a target (by ChEMBL target ID)
- Return canonical SMILES + pChEMBL value for every active compound
- Compound lookup by SMILES (similarity search)
- Target metadata (gene name, organism, target type)

---

## When to call this skill

- Starting a new campaign — pull known actives to seed enumeration / generative models
- Looking up SAR context for a scaffold of interest
- Validating a target ID before fetching bioactivity
- Benchmarking predicted affinities against experimental data

Do NOT use if:
- Need 3D structures → use `data/pdb_fetch`
- Need retrosynthesis routes → use `data/askcos`

---

## Inputs

| Field | Type | Required | Description |
|---|---|---|---|
| `target_id` | str | cond. | ChEMBL target ID, e.g. `"CHEMBL217"` (ABL1) |
| `smiles` | str | cond. | Query SMILES for compound search |
| `query_type` | str | — | `"bioactivity"` (default), `"target"`, or `"compound"` |
| `standard_type` | str | — | Activity type filter for bioactivity, e.g. `"IC50"` (default) |
| `limit` | int | — | Max records to return (default 100, max 1000) |
| `api_url` | str | — | ChEMBL API base (default `"https://www.ebi.ac.uk/chembl/api/data"`) |

`target_id` required for `"bioactivity"` and `"target"`; `smiles` required for `"compound"`.

---

## Outputs

| Field | Type | Description |
|---|---|---|
| `records` | list[dict] | Query results |
| `n_records` | int | Number of records returned |
| `query_type` | str | Echo of input |
| `target_id` | str \| None | Echo of input |
| `runtime_s` | float | API call duration |

Bioactivity record: `{molecule_chembl_id, smiles, standard_value, standard_units, pchembl_value, assay_chembl_id, assay_description}`.
Target record: `{target_chembl_id, pref_name, organism, target_type, components}`.
Compound record: `{molecule_chembl_id, smiles, mw, alogp, hbd, hba, pref_name}`.

Raises `RuntimeError("ChEMBL API unreachable: ...")` on network failure.

---

## Implementation

`tools/data/chembl.py` — `run_chembl()`

Calls the ChEMBL REST endpoints (`/activity.json`, `/target/<id>.json`, `/molecule.json`) with lazy `import requests`.

---

## Agent decision rules

- **ABL1 ChEMBL ID**: `CHEMBL217`; look up others at https://www.ebi.ac.uk/chembl/target_report_card
- **pChEMBL gate**: `pchembl_value >= 6` means IC50/Ki ≤ 1 µM — use as the active/inactive split
- **standard_type**: prefer `"IC50"` for kinase targets; `"Ki"` for GPCR / protease
- **Limit**: 100 records is fine for seeding generation; use `limit=1000` for full SAR pull
- **After ChEMBL pull**: deduplicate by SMILES, apply `rdkit_props` to compute physicochemical descriptors, then filter with `ligand_filter` before feeding generative models
- **Upstream**: target name → ChEMBL target search → target ID
- **Downstream**: `enumeration/rdkit_enum`, `generation/*`, `filtering/ligand_filter`

---

## Install

```bash
# No install needed — pure requests
uv add requests  # if not already installed
```

---

## References

- Mendez D. et al. *Nucleic Acids Res.* 2019, 47, D930 — "ChEMBL: towards direct deposition of bioassay data"
- REST API docs: https://chembl.gitbook.io/chembl-interface-documentation/web-services/chembl-data-web-services
- License: CC BY-SA 3.0 (data), MIT (client)
