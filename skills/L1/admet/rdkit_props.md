# Skill: rdkit_props

**Tier:** L1 — atomic tool skill
**Category:** admet
**Tool:** RDKit physicochemical descriptors + rule-based developability flags
**GPU required:** No — CPU only, pure RDKit

---

## What this skill does

Computes the canonical drug-likeness descriptor panel from RDKit on a single SMILES or a batch:

- **Physicochemical**: MW, logP (Crippen), HBD, HBA, TPSA, rotatable bonds, aromatic rings
- **Composite drug-likeness**: QED (Bickerton 2012), SA score (Ertl & Schuffenhauer 2009)
- **Rule-based pass/fail flags**: Lipinski Ro5, Veber, Egan
- Optional: Morgan-2 / Morgan-3 fingerprint bits

This is the bread-and-butter "is this molecule drug-like?" check. Fast, deterministic, no API calls. Use *before* docking / FEP to filter by physicochemistry, and use *after* generation to characterise produced libraries.

---

## When to call this skill

- After `enumeration/rdkit_enum` or `enumeration/library_search` — score libraries by QED + Ro5
- After `generation/*` — many generative models produce molecules with extreme MW or logP; rdkit_props flags them
- Before `docking/*` or `binding_affinity/*` — discard non-drug-like compounds before spending GPU time
- Whenever the agent needs to populate properties for a `visualization/ligand_viz` report

---

## Inputs

| Field | Type | Required | Description |
|---|---|---|---|
| `smiles` | str \| list[str] | ✓ | Single SMILES or batch list |
| `fingerprint` | str | — | `"morgan2"`, `"morgan3"`, or `None` (default `None` — skip FP) |
| `fp_bits` | int | — | Bit length when `fingerprint` is set (default 2048) |

---

## Outputs

### Single mode (input is `str`)

| Field | Type | Description |
|---|---|---|
| `smiles` | str | Echo of input SMILES |
| `mw` | float | Exact molecular weight (Da) |
| `logp` | float | Crippen logP |
| `hbd` | int | Hydrogen bond donors |
| `hba` | int | Hydrogen bond acceptors |
| `tpsa` | float | Topological polar surface area (Å²) |
| `rotatable_bonds` | int | Rotatable bond count |
| `aromatic_rings` | int | Aromatic ring count |
| `qed` | float | Quantitative drug-likeness ∈ [0, 1] |
| `sa_score` | float \| None | Synthetic accessibility ∈ [1, 10] (None if `Contrib/SA_Score` unavailable) |
| `lipinski_pass` | bool | MW ≤ 500, logP ≤ 5, HBD ≤ 5, HBA ≤ 10 |
| `veber_pass` | bool | RB ≤ 10, TPSA ≤ 140 |
| `egan_pass` | bool | −1 ≤ logP ≤ 6, TPSA ≤ 132 |
| `fingerprint` | list[int] | Bit array (only when `fingerprint` set) |

### Batch mode (input is `list[str]`)

| Field | Type | Description |
|---|---|---|
| `results` | list[dict] | One single-mode dict per input. Invalid SMILES yield `{"smiles": ..., "error": "Invalid SMILES"}` |
| `n_compounds` | int | Length of input list |

---

## Implementation

`tools/admet/rdkit_props.py` — `run_rdkit_props()`

---

## Agent decision rules

- **Default panel for triage**: don't supply a `fingerprint` (saves bytes); the rule-based flags + QED suffice for triage
- **Set `fingerprint="morgan2"`** only when downstream task is similarity / clustering / ML featurisation
- **`lipinski_pass`/`veber_pass`/`egan_pass`** are independent — log all three; agree-of-2 is a reasonable threshold for non-CNS targets
- **`qed > 0.5`** is a soft drug-likeness threshold; `qed > 0.7` for clinical-quality
- **`sa_score` returns None** when RDKit's Contrib SA_Score module isn't on the path — degrade gracefully, don't error
- **Property-based filter**: combine with `filtering/ligand_filter` (PAINS / BRENK) for full pre-docking gate
- **Upstream**: any tool producing SMILES — `rdkit_enum`, `library_search`, `generation/*`
- **Downstream**: `visualization/ligand_viz` (pass results as `properties`), `docking/*`, `binding_affinity/*`

---

## Install

```bash
uv add rdkit   # already a core dependency (>=2024.3)
```

---

## References

- Lipinski C.A. *Adv. Drug Deliv. Rev.* 2001, 46, 3 — Rule of 5
- Veber D.F. et al. *J. Med. Chem.* 2002, 45, 2615 — rotatable bonds + TPSA
- Egan W.J. et al. *J. Med. Chem.* 2000, 43, 3867 — logP × TPSA gate
- Bickerton G.R. et al. *Nature Chem.* 2012, 4, 90 — QED
- Ertl P. & Schuffenhauer A. *J. Cheminform.* 2009, 1, 8 — SA score
- License: BSD (RDKit)
