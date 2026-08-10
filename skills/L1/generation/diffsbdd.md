# Skill: diffsbdd

**Tier:** L1 — atomic tool skill
**Category:** generation
**Tool:** DiffSBDD — 3D structure-based molecule generation via score-based diffusion
**GPU required:** Yes — A100 recommended; ~5–30 min per 100 molecules

---

## What this skill does

DiffSBDD generates drug-like molecules **conditioned on a 3D protein binding pocket** using a score-based diffusion model that operates on 3D atom coordinates. Unlike SMILES-based generators, it samples molecules in 3D space simultaneously with the protein pocket, producing geometrically complementary scaffolds from the start. Output molecules include 3D coordinates in SDF format.

Key capabilities:
- Pocket-conditioned generation: molecules shaped to fit the target site
- Outputs SMILES + 3D SDF coordinates for each generated molecule
- Pocket specified by centre coordinates or residue IDs (Cα centroid computed)
- Validity filtering via RDKit (optional — falls back to raw output if RDKit unavailable)

---

## When to call this skill

- Structure-based hit generation: target has a known pocket and you want novel chemotypes
- No known actives / scaffolds — truly de novo generation
- Want 3D-aware generation rather than sequence-space sampling

Do NOT use if:
- You have a known scaffold and want analogues → use `enumeration/rdkit_enum`
- No pocket structure available → use `generation/molgpt` or `generation/reinvent4`
- Need antibody / peptide generation → use `biologics/rfantibody` or `biologics/rfdiffusion_pep`

---

## Inputs

| Field | Type | Required | Description |
|---|---|---|---|
| `protein_pdb` | str | ✓ | Target protein PDB (pocket atoms held fixed during generation) |
| `pocket_center` | tuple[float, float, float] | cond. | Pocket centre in Å (x, y, z); from `pocket_detect` output |
| `pocket_residues` | list[str] | cond. | Pocket lining residues `"<chain><resnum>"` — centroid computed automatically |
| `n_molecules` | int | — | Number of molecules to generate (default 100) |
| `diffsbdd_dir` | str | — | Path to DiffSBDD repository clone (default `"DiffSBDD"`) |
| `output_dir` | str | — | Results directory; temp dir if None |

Provide either `pocket_center` or `pocket_residues` (not both required; center takes priority).

---

## Outputs

| Field | Type | Description |
|---|---|---|
| `smiles` | list[str] | Valid generated SMILES (RDKit-validated) |
| `sdf_path` | str \| None | SDF file with 3D coordinates |
| `n_generated` | int | Number requested |
| `n_valid` | int | Number of valid SMILES |
| `validity_rate` | float | Fraction valid (0–1) |
| `pocket_center` | tuple \| None | Pocket centre used |
| `output_dir` | str | Output directory |

Raises `RuntimeError("DiffSBDD not found at ...")` if the repo clone is missing.
Raises `RuntimeError("DiffSBDD failed: ...")` on non-zero subprocess exit.

---

## Implementation

`tools/generation/diffsbdd.py` — `run_diffsbdd()`

Calls `test_conditioned_generation.py` from the DiffSBDD repo as a subprocess. Pocket centre is either passed directly or computed as the Cα centroid of the listed residues by parsing the PDB ATOM records. Output SDF is parsed by RDKit to extract canonical SMILES.

---

## Agent decision rules

- **Always use `pocket_detect` first**: pass `best_pocket["center_x/y/z"]` as `pocket_center` and `best_pocket["residues"]` as `pocket_residues`
- **n_molecules=100**: typical starting batch; ~60–80% validity rate expected
- **After generation**: apply `ligand_filter` (PAINS/BRENK) then `rdkit_props` (Lipinski/QED) before docking
- **Diversity**: if generated molecules cluster, run again with a different random seed or increase n_molecules
- **3D coordinates**: the SDF output can be used directly as a docking starting pose in `vina` or `gnina`
- **Upstream**: `pocket_detect` (pocket coordinates), `data/pdb_fetch` (target structure)
- **Downstream**: `ligand_filter`, `rdkit_props`, `docking/diffdock`, `docking/gnina`

---

## Install

```bash
git clone https://github.com/arneschneuing/DiffSBDD.git
cd DiffSBDD
conda env create -f environment.yml
conda activate DiffSBDD
# Download pre-trained checkpoint from repo README (~300 MB)
```

---

## References

- Schneuing A. et al. *ICLR* 2023 — "Structure-based Drug Design with Equivariant Diffusion Models"
- GitHub: https://github.com/arneschneuing/DiffSBDD
- License: MIT
