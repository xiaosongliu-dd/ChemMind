# Skill: gnina

**Tier:** L1 — atomic tool skill
**Category:** docking
**Tool:** GNINA (CNN-scored molecular docking)
**GPU required:** Yes — CUDA GPU strongly recommended; CPU fallback available

---

## What this skill does

GNINA is a molecular docking program that extends AutoDock Vina/Smina with a convolutional neural network (CNN) scoring function trained on protein–ligand complexes. It outputs both a CNN affinity score and a Vina score for each pose.

Key strengths:
- **CNN rescoring** dramatically improves pose ranking over Vina alone
- **Flexible receptor** side-chain handling
- Supports covalent docking
- Best general-purpose docking tool for SBDD campaigns when GPU is available

Use GNINA as the **default docking backend** for hit triage and virtual screening campaigns up to ~100k compounds.

---

## When to call this skill

- Need pose + docking score for a single molecule or focused library
- Running an SBDD virtual screen (≤100k)
- Rescoring poses from another docking program
- Covalent docking with a defined warhead

---

## Inputs

| Field | Type | Required | Description |
|---|---|---|---|
| `receptor_pdb` | str | ✓ | Path to prepared receptor PDB (no ligand, no water) |
| `ligand_smiles` | str or list[str] | ✓ | Query ligand(s) as SMILES |
| `center_x/y/z` | float | ✓ | Docking box center (Å); from P2Rank output |
| `box_size` | float | — | Box edge length in Å (default 22.5) |
| `n_poses` | int | — | Poses per ligand (default 9) |
| `cnn_scoring` | str | — | `"rescore"` or `"refinement"` (default `"rescore"`) |
| `covalent_res` | str | — | Residue for covalent docking e.g. `"CYS145"` |
| `exhaustiveness` | int | — | Search exhaustiveness (default 8) |

---

## Outputs

| Field | Type | Description |
|---|---|---|
| `poses` | list[dict] | Per-ligand list of `{smiles, pose_sdf, cnn_affinity, cnn_pose_score, vina_score, rank}` |
| `best_pose_sdf` | str | SDF path for top-ranked pose per ligand |
| `cnn_affinity` | float | CNN-predicted affinity (higher = better) |
| `vina_score` | float | Vina score kcal/mol (lower = better) |
| `runtime_s` | float | Wall time |

---

## Implementation

`tools/docking/gnina.py` — `run_gnina()`

---

## Agent decision rules

- **Primary docking tool** for SBDD; prefer over AutoDock-GPU when library < 100k
- **Score threshold**: `cnn_affinity > 6.0` as initial pass filter
- **Covalent docking**: pass `covalent_res` when target has known reactive cysteine/serine/lysine
- **Upstream**: `pocket/p2rank` (box center) → `gnina`
- **Downstream**: top poses → `binding_affinity/fep_openfe` (ΔΔG) or `admet/admetlab3`

---

## Install

```bash
# Linux binary (recommended)
wget https://github.com/gnina/gnina/releases/latest/download/gnina -O gnina && chmod +x gnina

# Or via conda
conda install -c conda-forge gnina
```

## References

- McNutt et al., J. Cheminform. 2021 — "GNINA 1.0: molecular docking with deep learning"
- Ragoza et al., J. Chem. Inf. Model. 2017 — original CNN scoring
- GitHub: https://github.com/gnina/gnina
- License: Apache 2.0
