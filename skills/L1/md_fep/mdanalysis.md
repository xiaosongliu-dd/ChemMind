# Skill: mdanalysis

**Tier:** L1 — atomic tool skill
**Category:** md_fep
**Tool:** MDAnalysis — Python-based MD trajectory analysis
**GPU required:** No — CPU-only; fast for typical production trajectories

---

## What this skill does

MDAnalysis parses molecular dynamics trajectories in any standard format (DCD, XTC, TRR, DMS, etc.) and computes structural analysis metrics that quantify system stability, flexibility, and binding-pocket geometry. It accepts topology + trajectory pairs from OpenMM, GROMACS, NAMD, or AMBER outputs.

Key capabilities:
- **RMSD**: backbone heavy-atom root-mean-square deviation vs reference frame or external PDB
- **RMSF**: per-residue root-mean-square fluctuation (flexibility profile)
- **Contacts**: inter-residue heavy-atom contact frequency map
- **Pocket volume**: convex-hull estimate of binding-site volume over time

---

## When to call this skill

- After `md_fep/openmm` — analyse the DCD trajectory for equilibration and stability
- Characterising pocket flexibility before docking (use mean pocket volume for box sizing)
- Identifying flexible loops (high RMSF) for constraint in docking or FEP
- Verifying that a protein–ligand complex remains stable over the simulation

Do NOT use if:
- Need free energy estimates → use `binding_affinity/fep_openfe`
- Need to run an MD simulation (not analyse one) → use `md_fep/openmm`

---

## Inputs

| Field | Type | Required | Description |
|---|---|---|---|
| `topology` | str | ✓ | Topology file: PDB, PSF, GRO, or PRMTOP |
| `trajectory` | str | ✓ | Trajectory file: DCD, XTC, TRR, or multi-frame PDB |
| `analyses` | list[str] | — | Subset of `["rmsd","rmsf","contacts","pocket_volume"]` (default: `["rmsd","rmsf"]`) |
| `reference_pdb` | str | — | External PDB for RMSD alignment; uses first trajectory frame if None |
| `selection` | str | — | MDAnalysis selection string (default `"protein and backbone"`) |
| `output_dir` | str | — | Results directory; temp dir if None |

---

## Outputs

| Field | Type | Description |
|---|---|---|
| `rmsd` | dict \| absent | `{step: list[int], rmsd_angstrom: list[float]}` |
| `rmsf` | dict \| absent | `{resids: list[int], rmsf_angstrom: list[float]}` |
| `contacts` | dict \| absent | `{matrix: list[list[float]], resids: list[int]}` (normalised 0–1) |
| `pocket_volume` | dict \| absent | `{step: list[int], volume_angstrom3: list[float]}` (convex hull) |
| `n_frames` | int | Total frames in trajectory |
| `n_atoms` | int | Atoms in universe |
| `runtime_s` | float | Wall time in seconds |

Raises `RuntimeError("MDAnalysis not installed: ...")` if the package is missing.
Raises `ValueError` for unknown analysis names.

---

## Implementation

`tools/md_fep/mdanalysis.py` — `run_mdanalysis()`

Uses `MDAnalysis.analysis.rms.RMSD` and `.RMSF` for the standard metrics. Contact maps iterate over all trajectory frames with a 4.5 Å heavy-atom cutoff. Pocket volume uses SciPy `ConvexHull` on the selection coordinates per frame. RMSD is written to a CSV file in `output_dir`.

---

## Agent decision rules

- **RMSD interpretation**: production RMSD should plateau within 5–10 ns; drift > 3 Å or continuous increase signals instability — discard this simulation and troubleshoot setup
- **RMSF hot spots**: residues with RMSF > 2 Å in the binding pocket deserve attention in docking; consider restraining them in rigid-receptor docking
- **Contacts**: symmetric matrix — entry `[i][j]` is the fraction of frames where residues i and j have a heavy-atom contact within 4.5 Å. Persistent contacts (> 0.7) define the pharmacophore envelope
- **Pocket volume**: if volume fluctuates > 30%, the pocket is cryptic or allosteric — standard docking may miss true binding modes
- **Selection**: use `"protein and backbone"` for global stability; narrow to `"resid 48:86 and segid A"` for pocket-specific analysis
- **Upstream**: `md_fep/openmm` (produces DCD + topology), `data/pdb_fetch` (reference PDB for RMSD)
- **Downstream**: `docking/*` (pocket volume informs box size), `binding_affinity/fep_openfe`

---

## Install

```bash
pip install MDAnalysis
# For scipy convex hull (pocket_volume analysis):
pip install scipy
```

---

## References

- Michaud-Agrawal N. et al. *J. Comput. Chem.* 2011, 32, 2319 — "MDAnalysis: A toolkit for the analysis of molecular dynamics simulations"
- Gowers R.J. et al. *Proc. 15th Python Sci. Conf.* 2016 — "MDAnalysis: A Python Package for the Rapid Analysis of Molecular Dynamics Simulations"
- GitHub: https://github.com/MDAnalysis/mdanalysis
- License: GPL-2.0
