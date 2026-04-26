# Skill: vina

**Tier:** L1 — atomic tool skill
**Category:** docking
**Tool:** AutoDock Vina / Vina-GPU 2.0
**GPU required:** No (Vina); Yes for Vina-GPU 2.0 (optional speedup)

---

## What this skill does

AutoDock Vina is the most widely used open-source docking program, featuring an iterated local search global optimizer and an empirical scoring function. Vina-GPU 2.0 is a GPU-accelerated fork that achieves 10–30× speedup with identical results.

Key strengths:
- **Simplest setup** — takes PDB/PDBQT receptor directly, no map pre-generation
- **Well-validated** — decades of benchmarks, largest literature coverage
- **Vina-GPU 2.0** matches CPU Vina results while screening ~100k/day on a single GPU
- **Flexible ligand + flexible receptor** side-chain support via `--flex`
- Good CPU fallback when no GPU is available

Use Vina as the **baseline/reference docking tool** and when setting up a new target without AutoGrid infrastructure. Prefer GNINA for accuracy and AutoDock-GPU for large throughput.

---

## When to call this skill

- Quick baseline score on a new target with minimal setup
- CPU-only environment (no GPU available)
- Cross-validation: compare Vina score against GNINA CNN score
- Flexible receptor docking with known flexible residues

---

## Inputs

| Field | Type | Required | Description |
|---|---|---|---|
| `receptor_pdbqt` | str | ✓ | Receptor in PDBQT format |
| `ligand_smiles` | str or list[str] | ✓ | Ligand(s) as SMILES |
| `center_x/y/z` | float | ✓ | Docking box center (Å); from P2Rank |
| `box_size` | float | — | Box edge length in Å (default 22.5) |
| `exhaustiveness` | int | — | Search thoroughness (default 8; use 32 for accuracy) |
| `n_poses` | int | — | Poses per ligand (default 9) |
| `flex_residues_pdbqt` | str | — | Path to flexible residues PDBQT for induced-fit docking |
| `use_gpu` | bool | — | Use Vina-GPU 2.0 if available (default True) |

---

## Outputs

| Field | Type | Description |
|---|---|---|
| `poses` | list[dict] | Per-ligand `{smiles, pose_pdbqt, vina_score, rank}` |
| `best_pose_pdbqt` | str | Path to top pose |
| `vina_score` | float | Best Vina score kcal/mol (lower = better) |
| `runtime_s` | float | Wall time |

---

## Implementation

`tools/docking/vina.py` — `run_vina()`

---

## Agent decision rules

- **Baseline tool**: always run when benchmarking a new target; compare to GNINA CNN score
- **Score threshold**: `vina_score < -7.0 kcal/mol` as initial pass
- **Flexible docking**: pass `flex_residues_pdbqt` for known hinge/loop residues (e.g. kinase DFG loop)
- **Exhaustiveness**: use 8 for VS, 32 for lead optimization accuracy
- **Upstream**: `pocket/p2rank` (box center) → `vina`
- **Downstream**: GNINA rescore → `binding_affinity/fep_openfe`

---

## Install

```bash
pip install vina meeko          # CPU Vina via Python bindings

# Vina-GPU 2.0 (CUDA)
git clone https://github.com/DeltaGroupNJUPT/Vina-GPU-2.0
# follow build instructions for your CUDA version
```

## References

- Eberhardt et al., J. Chem. Inf. Model. 2021 — "AutoDock Vina 1.2.0"
- Trott & Olson, J. Comput. Chem. 2010 — original AutoDock Vina
- Ding et al., J. Chem. Inf. Model. 2023 — Vina-GPU 2.0
- GitHub: https://github.com/ccsb-scripps/AutoDock-Vina
- License: Apache 2.0
