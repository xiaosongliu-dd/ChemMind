# Skill: openmm

**Tier:** L1 — atomic tool skill
**Category:** md_fep
**Tool:** OpenMM — GPU-accelerated molecular dynamics simulation
**GPU required:** Yes — CUDA strongly recommended; CPU ~100× slower for production runs

---

## What this skill does

OpenMM runs classical molecular dynamics (MD) simulations for protein–ligand systems using AMBER14 or CHARMM36 forcefields with particle-mesh Ewald electrostatics. It executes a standard three-phase protocol: energy minimisation → NVT equilibration → NVT/NPT production run. Trajectory and energy data are written to disk for downstream analysis with `mdanalysis`.

Key capabilities:
- NVT (constant volume/temperature) and NPT (constant pressure/temperature) ensembles
- AMBER14 all-atom forcefield with TIP3P-FB or CHARMM36 + TIP3P water
- Langevin Middle integrator with configurable timestep (2 fs default)
- DCD trajectory + CSV energy reporter at user-defined intervals
- Automatic solvation (periodic box with 1 nm water padding)

---

## When to call this skill

- Equilibrating a docked protein–ligand complex before FEP
- Sampling conformational flexibility of a binding pocket
- Generating an ensemble of frames for MMGBSA re-scoring
- Studying protein stability or loop dynamics of a design candidate

Do NOT use if:
- Need relative binding free energies → use `binding_affinity/fep_openfe` (handles FEP setup)
- Need a quick structure prediction → use `structure/esmfold` or `structure/colabfold`
- No GPU available and system > 10,000 atoms → MD will be too slow

---

## Inputs

| Field | Type | Required | Description |
|---|---|---|---|
| `structure_pdb` | str | ✓ | Input PDB (protein ± ligand; must have all heavy atoms) |
| `forcefield` | str | — | OpenMM forcefield XML (default `"amber14-all.xml"`) |
| `water_model` | str | — | Water XML (default `"amber14/tip3pfb.xml"`) |
| `ensemble` | str | — | `"NPT"` (default) or `"NVT"` |
| `temperature_k` | float | — | Simulation temperature in Kelvin (default 300.0) |
| `pressure_bar` | float | — | Pressure for NPT barostat (default 1.0) |
| `n_steps` | int | — | Production steps (default 1,000,000 = 2 ns at 2 fs timestep) |
| `timestep_fs` | float | — | Integration timestep in femtoseconds (default 2.0) |
| `report_interval` | int | — | Steps between trajectory frames (default 10,000 = 20 ps) |
| `output_dir` | str | — | Results directory; temp dir if None |

---

## Outputs

| Field | Type | Description |
|---|---|---|
| `trajectory_path` | str | DCD trajectory file of production run |
| `energy_path` | str | CSV log: step, potential energy (kJ/mol), temperature (K), volume, density |
| `final_structure` | str | PDB of last production frame |
| `potential_energy_kcal_mean` | float \| None | Mean potential energy (kcal/mol) over production |
| `temperature_k_mean` | float \| None | Mean temperature (K) over production |
| `n_steps_run` | int | Echo of `n_steps` |
| `timestep_fs` | float | Echo of `timestep_fs` |
| `runtime_s` | float | Wall time in seconds |

Raises `RuntimeError("openmm not installed: ...")` if OpenMM is missing.

---

## Implementation

`tools/md_fep/openmm.py` — `run_openmm()`

Uses `openmm.app.Modeller` to add solvent (1 nm TIP3P padding), builds the force field system with PME electrostatics and HBonds constraints, runs Langevin equilibration then production. Energy is parsed from the StateDataReporter CSV and converted from kJ/mol to kcal/mol.

---

## Agent decision rules

- **Pre-processing**: input PDB must have all heavy atoms — run through `pdbfixer` first if needed
- **n_steps=1,000,000** (2 ns) is the minimum for equilibrated sampling; use 5,000,000+ (10 ns) for loop dynamics
- **NPT vs NVT**: use NPT for solvated systems (default); NVT for vacuum or pre-equilibrated boxes
- **report_interval**: 10,000 steps = 20 ps per frame = 100 frames for 2 ns trajectory — sufficient for RMSD/RMSF
- **Post-simulation**: always run `mdanalysis` on the trajectory to compute RMSD drift — if RMSD > 3 Å in the first 500 ps, the system has not equilibrated; extend equilibration
- **FEP prep**: use OpenMM to equilibrate the apo or holo structure before creating the FEP perturbation network in `fep_openfe`
- **Upstream**: `data/pdb_fetch` (structure), `structure/colabfold` (predicted structure)
- **Downstream**: `md_fep/mdanalysis`, `binding_affinity/fep_openfe`

---

## Install

```bash
conda install -c conda-forge openmm
# or:
pip install openmm
# GPU support requires CUDA ≥ 11.x; verify with:
python -m openmm.testInstallation
```

---

## References

- Eastman P. et al. *PLOS Comput. Biol.* 2017, 13, e1005659 — "OpenMM 7: Rapid development of high performance algorithms for molecular dynamics"
- GitHub: https://github.com/openmm/openmm
- Documentation: http://docs.openmm.org
- License: MIT / LGPL-3.0
