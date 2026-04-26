# Skill: boltz2_affinity

**Tier:** L1 — atomic tool skill
**Category:** binding_affinity
**Tool:** Boltz-2 (affinity-focused wrapper)
**GPU required:** Yes — single A100/H100; ~20 s per ligand via NIM

---

## What this skill does

Predicts binding affinity (ΔG kcal/mol) for one or a batch of ligands against a protein target using Boltz-2. Wraps `boltz2` but returns only affinity fields — no structure output. Supports **batch mode**: pass a list of SMILES to screen a series in one call, results sorted by affinity ascending (tighter binders first).

Call `boltz2` directly when you need the 3D complex structure. Call this skill when you only need the affinity number.

---

## When to call this skill

- Triage a hit series (5–50 analogues) for relative potency ranking
- Quick affinity estimate on a newly generated molecule before docking
- No known pocket — Boltz-2 predicts affinity without a pre-defined box

---

## Inputs

| Field | Type | Required | Description |
|---|---|---|---|
| `protein_fasta` | str | ✓ | Target amino acid sequence |
| `ligand_smiles` | str or list[str] | ✓ | Single SMILES or list for batch screening |
| `n_samples` | int | — | Diffusion samples per ligand (default 5) |
| `use_nim` | bool | — | Use NVIDIA NIM API (default True) |
| `nim_api_key` | str | — | Reads `NVIDIA_API_KEY` env if None |

---

## Outputs

**Single ligand:**

| Field | Type | Description |
|---|---|---|
| `affinity_kcal_mol` | float | Predicted ΔG — lower = tighter binder |
| `affinity_confidence` | float | Model confidence 0–1 |
| `iptm_score` | float | Interface pTM (complex quality gate) |

**Batch (list input):**

| Field | Type | Description |
|---|---|---|
| `predictions` | list[dict] | Per-ligand results sorted by `affinity_kcal_mol` ascending |
| `n_ligands` | int | Total ligands scored |
| `n_failed` | int | Ligands that errored |
| `best_smiles` | str | Tightest predicted binder |
| `best_affinity_kcal_mol` | float | Best ΔG in the batch |

---

## Implementation

`tools/binding_affinity/boltz2_affinity.py` — `run_boltz2_affinity()`

---

## Agent decision rules

- **Affinity triage**: keep `affinity_kcal_mol < -8.0`; escalate `< -10.0` directly to FEP
- **Quality gate**: discard if `iptm_score < 0.5` — complex poorly predicted, score unreliable
- **Batch size**: keep ≤ 50 ligands per call to manage cost and latency
- **vs. deeppurpose**: use `boltz2_affinity` when absolute ΔG matters; use `deeppurpose` when only rank order is needed and speed is critical
- **Upstream**: `enumeration/rdkit_enum` or `generation/` → `boltz2_affinity`
- **Downstream**: top hits → `docking/gnina` (pose) → `binding_affinity/fep_openfe` (ΔΔG)

---

## Install

```bash
export NVIDIA_API_KEY="nvapi-..."   # NIM mode — no local install
```

## References

- Passaro et al. 2025 — "Boltz-2: Towards Accurate and Efficient Binding Affinity Prediction" (bioRxiv)
- NVIDIA NIM: https://docs.nvidia.com/nim/bionemo/boltz2/latest/overview.html
- License: MIT
