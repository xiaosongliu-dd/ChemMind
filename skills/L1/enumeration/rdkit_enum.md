# Skill: rdkit_enum

**Tier:** L1 — atomic tool skill
**Category:** enumeration
**Tool:** RDKit (BRICS, RECAP, R-group enumeration, scaffold decoration)
**GPU required:** No — CPU only

---

## What this skill does

Generates a combinatorial library of analogues from a known hit or scaffold using RDKit's built-in fragmentation and enumeration methods:

- **R-group enumeration** — fix a scaffold, vary substituents at defined attachment points
- **BRICS enumeration** — break molecule into BRICS fragments, recombine exhaustively
- **RECAP** — retrosynthetic combinatorial analysis procedure; generates synthetically accessible fragments
- **Scaffold decoration** — grow from a core by attaching user-supplied R-group SMILES lists

Use this skill immediately after identifying a hit to rapidly expand into an analogue series before scoring with docking or BA prediction.

---

## When to call this skill

- Have a hit molecule and want SAR exploration around a scaffold
- Need a focused library (100–50k compounds) for a docking campaign
- Upstream de novo generation produced a core; want systematic variation of side chains

---

## Inputs

| Field | Type | Required | Description |
|---|---|---|---|
| `smiles` | str | ✓ | Input hit/scaffold SMILES |
| `mode` | str | ✓ | `"rgroup"`, `"brics"`, or `"recap"` |
| `rgroup_smiles` | list[str] | — | R-group SMILES for `rgroup` mode; attachment via `[*]` wildcard |
| `attachment_points` | list[int] | — | Atom indices to vary in `rgroup` mode; auto-detected if None |
| `max_compounds` | int | — | Cap on output library size (default 10000) |
| `filter_lipinski` | bool | — | Apply Lipinski Ro5 filter to output (default True) |

---

## Outputs

| Field | Type | Description |
|---|---|---|
| `smiles_list` | list[str] | Enumerated compound SMILES |
| `n_generated` | int | Total generated before filtering |
| `n_returned` | int | After max_compounds / Lipinski filter |
| `scaffold_smiles` | str | MCS-derived scaffold used |

---

## Implementation

`tools/enumeration/rdkit_enum.py` — `run_rdkit_enum()`

---

## Agent decision rules

- **Mode selection**: use `brics` for scaffold-agnostic library expansion; use `rgroup` when you have a defined scaffold and an explicit substituent list from med chem
- **Size**: cap at 10k for docking campaigns; cap at 1k for FEP-destined libraries
- **Filter**: always apply Lipinski unless target is a natural product or macrocycle
- **Upstream**: `generation/` (hit) or `data/chembl` (known series) → `rdkit_enum`
- **Downstream**: `binding_affinity/gnina` or `binding_affinity/autodock_gpu` for scoring

---

## Install

```bash
uv add rdkit
```

## References

- RDKit BRICS: Degen et al., ChemMedChem 2008
- RDKit RECAP: Lewell et al., J. Chem. Inf. Comput. Sci. 1998
- GitHub: https://github.com/rdkit/rdkit
- License: BSD
