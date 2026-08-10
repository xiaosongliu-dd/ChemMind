# Skill: bindingdb

**Tier:** L1 — atomic tool skill
**Category:** data
**Tool:** BindingDB — experimental binding affinity database
**GPU required:** No

---

## What this skill does

Queries BindingDB (www.bindingdb.org) for experimentally measured binding affinities
(IC50, Ki, Kd, EC50) between small molecules and protein targets. Returns ligand SMILES,
affinity values in nM, assay types, and PubMed references. Use to:
- Seed a virtual screen with known actives (reference SMILES for similarity search)
- Validate docking predictions against experimental data
- Identify potency benchmarks for a target

---

## When to invoke

- Before Step 3 (library search) in virtual_screen: get reference SMILES for similarity query
- In lead_opt Step 3 (SAR context): retrieve known actives with IC50 data
- Whenever experimental potency data is needed for target X

---

## Signature

```python
bindingdb(
    target_name=<str | None>,       # e.g. "CDK2", "EGFR"
    smiles=<str | None>,            # query compound SMILES
    uniprot_id=<str | None>,        # e.g. "P24941" (CDK2)
    assay_type="IC50",              # "IC50", "Ki", "Kd", or "EC50"
    max_results=100,                # 1–500
    min_affinity_nm=None,           # filter: keep affinity ≤ X nM
    max_affinity_nm=None,           # filter: keep affinity ≥ X nM
)
→ {
    records: [
        {
            ligand_smiles: str | None,
            target_name:   str | None,
            affinity_nm:   float | None,    # in nM; None if not parseable
            affinity_type: str,             # IC50 / Ki / Kd / EC50
            uniprot_id:    str | None,
            pubmed_id:     str | None,
            ligand_name:   str | None,
        },
        ...    # sorted most-potent first (lowest affinity_nm first)
    ],
    n_records: int,
    query:     dict,
    runtime_s: float,
}
```

---

## Common UniProt IDs for drug targets

| Target | UniProt |
|---|---|
| CDK2 | P24941 |
| EGFR | P00533 |
| ABL1 | P00519 |
| BRAF | P15056 |
| JAK2 | O60674 |
| HIV-1 protease | P04585 |

---

## Usage pattern

```python
# Seed virtual screen with known CDK2 actives ≤ 100 nM
bindingdb(
    uniprot_id="P24941",
    assay_type="IC50",
    min_affinity_nm=100,
    max_results=50,
)
→ use records[0]["ligand_smiles"] as reference_smiles in library_search()
```

---

## References

- Gilson et al., *Nucleic Acids Res.* 2016 — BindingDB
- API docs: https://bindingdb.org/rwd/bind/chemsearch/marvin/SDFdownload.rwd
