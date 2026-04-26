# Skill: library_search

**Tier:** L1 — atomic tool skill
**Category:** enumeration
**Tool:** Enamine REAL / ZINC make-on-demand virtual library search
**GPU required:** No

---

## What this skill does

Searches make-on-demand virtual libraries for purchasable analogues similar to a query molecule:

- **Enamine REAL** — 48B+ synthetically accessible compounds, 2–4 week delivery
- **ZINC22** — 1.4B+ purchasable compounds, subset of Enamine and other vendors
- **Similarity search** — Tanimoto fingerprint search against the library
- **Substructure search** — retrieve all library members containing a scaffold/fragment

Use this when you want to expand a hit into a purchasable analogue series without synthesizing novel chemistry — compounds can be ordered immediately.

---

## When to call this skill

- Have a hit and want nearest purchasable neighbours for quick SAR
- Want a focused docking library anchored to a core substructure
- Need make-on-demand compounds for a fast experimental follow-up

---

## Inputs

| Field | Type | Required | Description |
|---|---|---|---|
| `smiles` | str | ✓ | Query molecule SMILES |
| `library` | str | — | `"enamine_real"` or `"zinc22"` (default `"enamine_real"`) |
| `mode` | str | — | `"similarity"` or `"substructure"` (default `"similarity"`) |
| `tanimoto_threshold` | float | — | Similarity cutoff 0–1 (default 0.7) |
| `max_results` | int | — | Max compounds to return (default 500) |
| `fingerprint` | str | — | `"ecfp4"` or `"fcfp4"` (default `"ecfp4"`) |

---

## Outputs

| Field | Type | Description |
|---|---|---|
| `smiles_list` | list[str] | Retrieved compound SMILES |
| `similarities` | list[float] | Tanimoto similarity to query (similarity mode) |
| `catalog_ids` | list[str] | Vendor catalog IDs |
| `n_returned` | int | Number of compounds returned |
| `library_searched` | str | Library name actually queried |

---

## Implementation

`tools/enumeration/library_search.py` — `run_library_search()`

---

## Agent decision rules

- **Default to Enamine REAL** for small molecule drug discovery — largest make-on-demand space
- **ZINC22** if you need free/open catalog IDs or cross-vendor coverage
- **Similarity threshold**: 0.7 for broad SAR coverage; 0.85+ for tight analogue series
- **Size**: keep `max_results` ≤ 500 before docking; filter further post-docking before ordering
- **Upstream**: `binding_affinity/` hit → `library_search` for purchasable analogues
- **Downstream**: `binding_affinity/gnina` or `autodock_gpu` → `admet/admetlab3` → order

---

## Install

```bash
uv add requests rdkit
# No API key required for Enamine public API or ZINC22
```

## References

- Enamine REAL database: https://enamine.net/compound-collections/real-compounds/real-database
- ZINC22: Tingle et al., J. Chem. Inf. Model. 2023
- Irwin et al., J. Med. Chem. 2020 (ZINC20)
- License: API access subject to Enamine / ZINC terms of service
