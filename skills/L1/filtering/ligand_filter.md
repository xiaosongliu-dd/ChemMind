# Skill: ligand_filter

**Tier:** L1 — atomic tool skill
**Category:** filtering
**Tool:** RDKit `FilterCatalog` (PAINS, BRENK, NIH, ZINC, CHEMBL alert sets) + custom SMARTS
**GPU required:** No — CPU only

---

## What this skill does

Flags problematic small molecules — pan-assay interference compounds (PAINS), reactive / unstable functionalities (BRENK), unwanted-substructure lists (NIH, ZINC) and CHEMBL medchem rules — using RDKit's built-in `FilterCatalog`. Optional user-supplied SMARTS patterns let the agent encode bespoke project alerts on top of the standard catalogs.

Use this skill as a **triage gate** between SMILES generation/enumeration and any expensive downstream tool (docking, BA prediction, FEP). Filtering out aggregators and Michael acceptors before docking saves substantial GPU time on hits that are guaranteed false positives.

---

## When to call this skill

- After `enumeration/rdkit_enum` or `enumeration/library_search` — triage the generated library before docking
- After de novo `generation/*` — many generative models produce reactive/strained structures
- Before `docking/*`, `binding_affinity/*`, or `md_fep/*` — any expensive downstream step
- When a project has bespoke structural alerts to enforce (`custom_smarts`)

---

## Inputs

| Field | Type | Required | Description |
|---|---|---|---|
| `smiles` | str \| list[str] | ✓ | Single SMILES or batch list |
| `filter_sets` | str \| list[str] | — | Preset (`"default"`, `"strict"`, `"all"`) or explicit list of catalog names; default `"default"` |
| `custom_smarts` | list[str] | — | Project-specific SMARTS patterns to flag (in addition to selected catalogs) |
| `return_canonical` | bool | — | Canonicalise SMILES on output (default True) |

**Presets**
- `"default"` → `["PAINS", "BRENK"]` — standard medchem triage
- `"strict"` → `["PAINS", "BRENK", "NIH", "ZINC"]` — stricter HTS-style filtering
- `"all"` → PAINS_A/B/C + BRENK + NIH + ZINC + all CHEMBL_* catalogs
- `"lilly"` → reserved for future `rd_filters` integration (raises `NotImplementedError`)

**Catalog names accepted in explicit list:** `PAINS`, `PAINS_A`, `PAINS_B`, `PAINS_C`, `BRENK`, `NIH`, `ZINC`, `CHEMBL`, `CHEMBL_Glaxo`, `CHEMBL_Dundee`, `CHEMBL_BMS`, `CHEMBL_LINT`, `CHEMBL_MLSMR`, `CHEMBL_SureChEMBL`.

---

## Outputs

### Single mode (input is `str`)

| Field | Type | Description |
|---|---|---|
| `passed` | bool | True if `n_alerts == 0` |
| `n_alerts` | int | Total alert hits across all selected catalogs + custom SMARTS |
| `failures` | list[dict] | One entry per alert hit (see schema below) |
| `smiles_canonical` | str | RDKit canonical SMILES |
| `filter_sets_used` | list[str] | Resolved catalog names actually applied |

### Batch mode (input is `list[str]`)

| Field | Type | Description |
|---|---|---|
| `passed_smiles` | list[str] | Canonical SMILES that passed all filters |
| `failed_smiles` | list[str] | Canonical SMILES that hit at least one alert |
| `n_passed` | int | `len(passed_smiles)` |
| `n_failed` | int | `len(failed_smiles)` |
| `results` | list[dict] | One single-mode dict per input, in order. Invalid SMILES yield `passed=False` plus an `error` key (no exception) |
| `filter_sets_used` | list[str] | Resolved catalog names actually applied |

### Failure entry schema

| Field | Type | Description |
|---|---|---|
| `catalog` | str | Catalog name (e.g. `"PAINS_A"`) or `"custom"` |
| `alert_name` | str | Short alert ID from RDKit |
| `description` | str | Human-readable alert description |
| `smarts` | str | SMARTS pattern that matched (may be empty for some catalogs) |

---

## Implementation

`tools/filtering/ligand_filter.py` — `run_ligand_filter()`

---

## Agent decision rules

- **Default for SAR triage**: `filter_sets="default"` (PAINS + BRENK) — fast, low false-positive
- **HTS-style screening**: `filter_sets="strict"` — adds NIH + ZINC reactive/aggregator lists
- **Property-based filtering** (MW, logP, Ro5, QED, SA): NOT this skill — call `admet/rdkit_props` or rely on `enumeration/rdkit_enum`'s built-in `filter_lipinski`
- **Standardisation / salt stripping**: NOT in this skill — input must be pre-standardised (future `tools/filtering/standardize.py` will cover this)
- **Custom alerts**: pass project-specific reactive/toxic SMARTS via `custom_smarts`
- **Upstream**: `enumeration/rdkit_enum`, `enumeration/library_search`, `generation/*`
- **Downstream**: `docking/gnina`, `docking/diffdock`, `binding_affinity/boltz2_affinity`, `md_fep/fep_openfe`

---

## Install

```bash
uv add rdkit  # already core dependency (>=2024.3)
```

---

## References

- Baell J. B. & Holloway G. A. *J. Med. Chem.* 2010, 53, 2719 — PAINS filters
- Brenk R. et al. *ChemMedChem* 2008, 3, 435 — unwanted functionality alerts
- RDKit FilterCatalog docs: https://www.rdkit.org/docs/source/rdkit.Chem.FilterCatalog.html
- Future work — Lilly MedChem rules via Pat Walters' `rd_filters`: https://github.com/PatWalters/rd_filters
- License: BSD (RDKit)
