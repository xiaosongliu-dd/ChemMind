# Skill: fep_openfe

**Tier:** L1 — atomic tool skill
**Category:** binding_affinity
**Tool:** OpenFE (relative binding free energy)
**GPU required:** Yes — multi-GPU recommended; ~4–12 h per edge on A100

---

## What this skill does

OpenFE runs relative binding free energy perturbation (RBFE) calculations to compute ΔΔG between pairs of ligands with uncertainty estimates. It uses the hybrid topology HREX protocol via OpenMM, with LoMap for automatic perturbation network planning.

This is the **gold standard binding affinity method** — use it only at the late lead optimization stage (5–20 compounds) after docking and ML-based triage. Results are quantitative ΔΔG in kcal/mol with sub-kcal uncertainty when converged.

---

## When to call this skill

- Final lead series (≤20 compounds) requiring quantitative ΔΔG
- Prioritizing which analogue to synthesize next
- Scaffold hop validation — confirm structural change doesn't lose potency

---

## Inputs

| Field | Type | Required | Description |
|---|---|---|---|
| `protein_pdb` | str | ✓ | Prepared receptor PDB (no ligand, no crystal waters) |
| `ligand_sdf_dir` | str | ✓ | Directory of ligand SDF files (one per file, with 3D coords) |
| `reference_ligand_sdf` | str | ✓ | SDF of reference/anchor ligand (known binder) |
| `n_replicas` | int | — | HREX replicas per λ window (default 3) |
| `n_lambda` | int | — | λ windows per edge (default 11) |
| `forcefield` | str | — | Small molecule FF: `"openff-2.1.0"` or `"gaff-2.11"` (default openff) |
| `solvent` | str | — | Water model: `"tip3p"` or `"tip4p-ew"` (default tip3p) |
| `work_dir` | str | — | Output directory (default `"fep_output"`) |

---

## Outputs

| Field | Type | Description |
|---|---|---|
| `edge_results` | list[dict] | Per-edge `{ligandA, ligandB, DDG_kcal_mol, uncertainty, converged}` |
| `n_edges` | int | Total perturbation edges run |
| `n_converged` | int | Edges with uncertainty < 0.5 kcal/mol |
| `work_dir` | str | Path to trajectories and result files |

---

## Implementation

`tools/binding_affinity/fep_openfe.py` — `run_fep_openfe()`

---

## Agent decision rules

- **Only use for ≤20 compounds** — cost is O(edges × n_lambda × n_replicas × simulation time)
- **Convergence gate**: discard edges where `uncertainty > 0.5 kcal/mol`; re-run with more replicas
- **Upstream**: `docking/gnina` top hits → 3D poses in SDF → `fep_openfe`
- **Downstream**: order synthesis of top-ranked compounds → experimental IC50

---

## Install

```bash
pip install openfe lomap2 openmmforcefields
```

## References

- Scheen et al., JCTC 2023 — "OpenFE: a robust, accessible framework for alchemical FEP"
- LoMap: Liu et al., J. Comput.-Aided Mol. Des. 2013
- GitHub: https://github.com/OpenFreeEnergy/openfe
- License: MIT
