# Skill: autodock_gpu

**Tier:** L1 — atomic tool skill
**Category:** docking
**Tool:** AutoDock-GPU (massively parallel GPU docking)
**GPU required:** Yes — designed for multi-GPU; scales linearly with GPU count

---

## What this skill does

AutoDock-GPU is a GPU-accelerated reimplementation of AutoDock4 that runs thousands of Lamarckian genetic algorithm (LGA) docking jobs in parallel. It is the fastest open-source docking tool for **large virtual screening** campaigns (100k–10M compounds).

Key strengths:
- **Throughput**: 100x faster than CPU AutoDock4; screens millions of compounds overnight on 4× A100
- **AutoDock4 scoring**: established empirical scoring function, widely validated
- **PDBQT format**: direct compatibility with AutoDock ecosystem tooling (ADFRsuite, MGLTools)
- Best choice when GNINA throughput is insufficient

Use AutoDock-GPU as the **high-throughput VS backend** for large library screening; use GNINA for accuracy on shortlisted hits.

---

## When to call this skill

- Library > 100k compounds
- Rapid first-pass screen before GNINA rescoring
- Multi-GPU cluster available

---

## Inputs

| Field | Type | Required | Description |
|---|---|---|---|
| `receptor_pdbqt` | str | ✓ | Path to prepared receptor PDBQT |
| `ligand_pdbqt_dir` | str | ✓ | Directory of ligand PDBQT files |
| `gpf_path` | str | ✓ | AutoGrid parameter file (.gpf) defining docking box |
| `n_runs` | int | — | LGA runs per ligand (default 20) |
| `heuristics` | bool | — | Use heuristic termination for speed (default True) |
| `npts` | tuple | — | Grid points (x,y,z); overrides gpf if provided |
| `center` | tuple | — | Box center (x,y,z) in Å |

---

## Outputs

| Field | Type | Description |
|---|---|---|
| `results` | list[dict] | Per-ligand `{ligand_id, best_energy, ki_estimate, pose_pdbqt}` |
| `best_energy` | float | Lowest docking energy kcal/mol (most negative = best) |
| `n_docked` | int | Ligands successfully docked |
| `failed` | list[str] | Ligand IDs that failed |
| `runtime_s` | float | Wall time |

---

## Implementation

`tools/docking/autodock_gpu.py` — `run_autodock_gpu()`

---

## Agent decision rules

- **Use for VS >100k** — below that, GNINA accuracy is worth the extra time
- **Two-stage screen**: AutoDock-GPU first pass → top 1% → GNINA rescore → top 0.1% → FEP
- **Score threshold**: `best_energy < -7.0 kcal/mol` as initial pass
- **Receptor prep**: must use ADFRsuite/AutoGrid4 to generate maps before docking
- **Upstream**: `pocket/p2rank` (box) → AutoGrid4 (maps) → `autodock_gpu`
- **Downstream**: GNINA rescore → `binding_affinity/fep_openfe`

---

## Install

```bash
# Pre-built GPU binary (Linux)
wget https://github.com/ccsb-scripps/AutoDock-GPU/releases/latest/download/autodock_gpu_128wi
chmod +x autodock_gpu_128wi

# Receptor & ligand prep
pip install meeko
conda install -c conda-forge adfrsuite
```

## References

- Santos-Martins et al., J. Chem. Theory Comput. 2021 — "Accelerating AutoDock4 with GPUs"
- Morris et al., J. Comput. Chem. 2009 — AutoDock4 scoring function
- GitHub: https://github.com/ccsb-scripps/AutoDock-GPU
- License: LGPL v2.1
