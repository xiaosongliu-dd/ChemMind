# Skill: colabfold

**Tier:** L1 — atomic tool skill
**Category:** structure
**Tool:** ColabFold — fast MSA-based structure prediction via AlphaFold2
**GPU required:** Yes — A100/V100 recommended; ~5–30 min per model × 5 models

---

## What this skill does

ColabFold accelerates AlphaFold2 by replacing the expensive Jackhmmer/HHblits MSA search with **MMseqs2 API lookup**, reducing compute time from hours to minutes while retaining AF2's accuracy. It supports monomer and **multimer** prediction, optional structural templates, and AMBER refinement.

Key capabilities:
- Monomer and multimer (protein complex) structure prediction
- MSA computed via fast MMseqs2 API (no local databases needed)
- 5 model ensemble with ranked outputs (rank_1 = best)
- Optional AMBER side-chain refinement
- Template search against PDB30 (optional)

---

## When to call this skill

- Need AF2-quality prediction with MSA
- Predicting protein–protein complex structure (multimer mode via `:` separator)
- Comparative modelling against known templates
- More accurate than ESMFold when homologous sequences exist

Do NOT use if:
- Need a fast single-sequence estimate → use `esmfold`
- Input is an antibody → use `igfold` or `abodybuilder3`
- Need protein–ligand complex → use `boltz2`

---

## Inputs

| Field | Type | Required | Description |
|---|---|---|---|
| `sequences` | str \| list[str] | ✓ | FASTA string(s); multimer: chains separated by `:` in one entry |
| `n_models` | int | — | Number of AF2 models to run (default 5) |
| `use_templates` | bool | — | Search PDB30 for templates (default True) |
| `use_amber` | bool | — | AMBER side-chain refinement (default False) |
| `colabfold_bin` | str | — | Path to `colabfold_batch` binary (default `"colabfold_batch"`) |
| `output_dir` | str | — | Results directory; temp dir if None |

---

## Outputs

| Field | Type | Description |
|---|---|---|
| `pdb_paths` | list[str] | All model PDBs, rank_1 first |
| `best_pdb_path` | str \| None | Highest-ranked model PDB |
| `mean_plddt` | float \| None | Mean pLDDT of best model |
| `ptm_scores` | list[float \| None] | Per-model pTM scores |
| `n_sequences` | int | Number of input sequences |
| `output_dir` | str | Path to colabfold output directory |

Raises `RuntimeError("colabfold_batch failed: ...")` on non-zero exit.

---

## Implementation

`tools/structure/colabfold.py` — `run_colabfold()`

Writes input sequences to a FASTA file, calls `colabfold_batch <fasta> <outdir>`, parses ranked PDBs and optional per-model JSON scores.

---

## Agent decision rules

- **Always prefer rank_1 model**: it has the highest pLDDT-ranked score across the ensemble
- **Multimer**: pass multiple chains as `sequence1:sequence2` in a single entry
- **n_models=1**: use only for quick scouts; 5 models give a consensus confidence estimate
- **Template toggle**: disable (`use_templates=False`) for de novo / engineered sequences with no known homologues
- **AMBER**: adds ~30% wall time; enables only when backbone accuracy matters for docking
- **pLDDT gate**: `mean_plddt > 80` is high confidence; 70–80 is acceptable; below 70 investigate
- **Upstream**: `data/pdb_fetch` (homologue structures), `data/chembl` (target IDs)
- **Downstream**: `pocket_detect`, `docking/*`, `boltz2` (affinity)

---

## Install

```bash
pip install "colabfold[alphafold]"
# Or via conda:
conda install -c conda-forge -c bioconda colabfold
# See https://github.com/sokrypton/ColabFold for full setup
```

---

## References

- Mirdita M. et al. *Nat. Methods* 2022, 19, 679 — "ColabFold: making protein folding accessible to all"
- GitHub: https://github.com/sokrypton/ColabFold
- AlphaFold2: Jumper J. et al. *Nature* 2021, 596, 583
- License: MIT
