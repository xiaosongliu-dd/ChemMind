# Skill: diffdock

**Tier:** L1 — atomic tool skill
**Category:** docking
**Tool:** DiffDock (diffusion-based blind docking)
**GPU required:** Yes — A100/H100 recommended; ~30 s per ligand on A100

---

## What this skill does

DiffDock is a diffusion generative model for molecular docking. Instead of searching within a predefined box, it samples poses over the full protein surface — making it a **blind docking** method that requires no pocket specification.

Key strengths:
- **Blind docking** — no box center or pocket needed
- **Diffusion sampling** — generates a diverse ensemble of poses; better coverage of binding modes
- **Top-1 accuracy** competitive with or better than Vina/GNINA on PoseBusters benchmark
- Available via NVIDIA NIM API (no local GPU required)

Use DiffDock when:
- Pocket location is unknown and P2Rank confidence is low
- You want pose diversity for an ensemble docking analysis
- Allosteric or cryptic site investigation

---

## When to call this skill

- No pocket identified yet (blind screen)
- Ensemble of diverse poses needed for downstream MD or FEP
- Novel target with no known binder — avoid pocket bias

---

## Inputs

| Field | Type | Required | Description |
|---|---|---|---|
| `protein_pdb` | str | ✓ | Path to prepared receptor PDB |
| `ligand_smiles` | str | ✓ | Ligand SMILES |
| `n_samples` | int | — | Diffusion samples (default 10) |
| `inference_steps` | int | — | Diffusion steps (default 20) |
| `use_nim` | bool | — | Use NVIDIA NIM API (default True) |
| `nim_api_key` | str | — | Reads `NVIDIA_API_KEY` env if None |

---

## Outputs

| Field | Type | Description |
|---|---|---|
| `poses` | list[dict] | `{pose_pdb, confidence_score, rank}` sorted by confidence |
| `best_pose_pdb` | str | Path to top-ranked pose |
| `confidence_score` | float | Model confidence for top pose (higher = better) |
| `n_poses` | int | Number of poses returned |
| `runtime_s` | float | Wall time |

---

## Implementation

`tools/docking/diffdock.py` — `run_diffdock()`

---

## Agent decision rules

- **Use when pocket is unknown** — skip P2Rank when target is novel or flexible
- **Ensemble use**: run with `n_samples=10`, cluster poses by RMSD, pick representative from each cluster for FEP
- **Confidence threshold**: `confidence_score > 0` as rough pass; negative = poorly predicted
- **Do NOT use for large-library VS** (>10k) — too slow; use GNINA or AutoDock-GPU instead
- **Upstream**: `structure/boltz2` or `structure/colabfold` → `diffdock` (no pocket needed)
- **Downstream**: `binding_affinity/fep_openfe` or `admet/admetlab3`

---

## Install

```bash
export NVIDIA_API_KEY="nvapi-..."   # NIM mode — no local install

# Local
git clone https://github.com/gcorso/DiffDock
cd DiffDock && pip install -e .
```

## References

- Corso et al., ICLR 2023 — "DiffDock: Diffusion Steps, Twists, and Turns for Molecular Docking"
- NVIDIA NIM: https://docs.nvidia.com/nim/bionemo/diffdock/latest/overview.html
- GitHub: https://github.com/gcorso/DiffDock
- License: MIT
