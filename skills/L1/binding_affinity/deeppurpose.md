# Skill: deeppurpose

**Tier:** L1 — atomic tool skill
**Category:** binding_affinity
**Tool:** DeepPurpose (DTI prediction)
**GPU required:** No — CPU inference, ~1 s/ligand

---

## What this skill does

DeepPurpose predicts drug-target binding affinity (pIC50 / pKd) using pretrained deep learning models (MPNN, Transformer, CNN) trained on BindingDB, DAVIS, and KIBA.

**Inputs are amino acid sequence + SMILES — no 3D structure, no PDB required.**

The model encodes the full protein sequence (not just the pocket) alongside a molecular graph of the ligand, then predicts affinity in a single forward pass. This makes it fast (CPU-only, ~1 s/ligand) but also pocket-agnostic: it cannot distinguish two proteins that differ only in their active site.

**Strength**: works for any target with a known sequence, even if no crystal structure exists.
**Weakness**: accuracy degrades for protein families outside the training distribution. BindingDB training data is dominated by kinases, proteases, and nuclear receptors. Use with caution for GPCRs, ion channels, or novel folds.

Use DeepPurpose for **fast relative rank ordering** (10–10k compounds) when you need to down-select before docking, or when no structure is available.

---

## When to call this skill

- Screening a series of analogues for relative potency ranking
- No GPU available (CPU-only environment)
- Need a fast pIC50 estimate before committing to docking

---

## Inputs

| Field | Type | Required | Description |
|---|---|---|---|
| `protein_sequence` | str | ✓ | Full amino acid sequence (single-letter, e.g. UniProt canonical sequence) |
| `ligand_smiles` | str or list[str] | ✓ | SMILES string or list for batch screening |
| `model` | str | — | Pretrained model ID (default `"MPNN_CNN_BindingDB_IC50"`) |

**Available models:**

| Model ID | Dataset | Output | Best for |
|---|---|---|---|
| `MPNN_CNN_BindingDB_IC50` | BindingDB | pIC50 | General small molecule (default) |
| `Transformer_CNN_BindingDB_Kd` | BindingDB | pKd | Equilibrium binding constant |
| `CNN_CNN_BindingDB_IC50` | BindingDB | pIC50 | Faster, lower accuracy |
| `MPNN_CNN_DAVIS` | DAVIS | pKd | Kinases specifically |

---

## Outputs

| Field | Type | Description |
|---|---|---|
| `predictions` | list[dict] | Per-ligand `{smiles, predicted_value, units}` sorted descending (higher = tighter) |
| `n_ligands` | int | Number of ligands scored |
| `units` | str | `"pIC50"` or `"pKd"` depending on model |
| `best_smiles` | str | Highest predicted affinity SMILES |
| `best_value` | float | Best predicted value in the batch |

---

## Implementation

`tools/binding_affinity/deeppurpose.py` — `run_deeppurpose()`

---

## Agent decision rules

- **Use when no structure is available**: if no PDB or AlphaFold model exists, this is the only fast BA option
- **Use for relative ranking, not absolute ΔG**: pIC50 outputs are only meaningful within a series against the same target
- **Do NOT use for novel protein families**: GPCRs, ion channels, transporters, or any family underrepresented in BindingDB will give unreliable results — use `boltz2_affinity` instead
- **Score threshold**: `predicted_value > 7.0` (pIC50 > 7 = IC50 < 100 nM) as initial pass filter
- **If structure is available**: prefer `docking/gnina` (pocket-aware) over DeepPurpose at similar throughput
- **Upstream**: `enumeration/rdkit_enum` or `enumeration/library_search` → `deeppurpose`
- **Downstream**: top 1% → `docking/gnina` for pose → `binding_affinity/fep_openfe` for gold standard

---

## Install

```bash
pip install DeepPurpose
```

## References

- Huang et al., Bioinformatics 2021 — "DeepPurpose: a deep learning library for drug-target interaction prediction"
- GitHub: https://github.com/kexinhuang12345/DeepPurpose
- License: BSD-3
