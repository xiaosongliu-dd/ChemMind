# Skill: igfold

**Tier:** L1 — atomic tool skill
**Category:** biologics
**Tool:** IgFold (deep-learning antibody / nanobody structure prediction)
**GPU required:** Optional — CPU-runnable; GPU 5–10× faster

---

## What this skill does

IgFold is a fast, antibody-specific structure-prediction model from the Gray Lab. It folds VH+VL paired antibodies (or VHH single-domain nanobodies) in seconds without an MSA, with particular accuracy on the CDR loops — especially the notoriously variable CDR-H3.

Use this skill whenever you need an antibody 3D structure and don't already have a co-crystal: it's faster than ColabFold/AF2 for antibodies, doesn't need a multiple sequence alignment, and reports per-residue pLDDT + per-residue predicted RMSD (pRMSD) so the agent can gate on confidence.

---

## When to call this skill

- After `proteinmpnn` or `rfantibody` — fold the designed sequence to validate CDR geometry
- After querying ChEMBL / antibody databases for a candidate sequence
- Before docking an antibody as receptor — IgFold gives a clean Fv structure to feed into `docking/diffdock` or `binding_affinity/boltz2_affinity` (Ab–Ag mode)
- Whenever speed matters more than peak accuracy — IgFold ~30 s, ColabFold ~5–30 min

---

## Inputs

| Field | Type | Required | Description |
|---|---|---|---|
| `heavy_sequence` | str | ✓ | VH (or VHH for nanobody) amino acid sequence |
| `light_sequence` | str | — | VL sequence; omit for nanobody / single-domain |
| `n_models` | int | — | Ensemble size 1–4; more = slower but better pRMSD estimate (default 4) |
| `do_refine` | bool | — | Rosetta-based all-atom refinement (default False — requires PyRosetta) |
| `output_dir` | str | — | Where to write the PDB; uses temp dir if None |

---

## Outputs

| Field | Type | Description |
|---|---|---|
| `pdb_path` | str | Path to predicted structure PDB |
| `mean_plddt` | float \| None | Mean per-residue pLDDT confidence ∈ [0, 1] |
| `mean_prmsd` | float \| None | Mean per-residue predicted RMSD (Å); lower = higher confidence |
| `heavy_sequence` | str | Echo of input |
| `light_sequence` | str \| None | Echo of input |
| `mode` | str | `"antibody"` or `"nanobody"` (auto-detected from light chain presence) |
| `n_models` | int | Echo of input |

Raises `RuntimeError("igfold not installed")` if the `igfold` package can't be imported.

---

## Implementation

`tools/biologics/igfold.py` — `run_igfold()`

---

## Agent decision rules

- **Confidence gate**: discard structures with `mean_plddt < 0.85` or `mean_prmsd > 1.5 Å` — CDR-H3 below this is unreliable
- **Ensemble size**: default `n_models=4` for production; `n_models=1` only for ultra-fast triage when you already trust the sequence
- **Refinement**: skip `do_refine=True` unless you specifically need all-atom side-chain accuracy — adds 30–120 s and requires PyRosetta installed
- **Pair with abodybuilder3 for an ensemble**: IgFold and ABodyBuilder3 disagree on ~10% of CDR-H3 conformations; running both and taking the consensus boosts confidence
- **Upstream**: `biologics/proteinmpnn`, `biologics/rfantibody`, `data/chembl` (antibody DB)
- **Downstream**: `docking/*` (Ab as receptor), `binding_affinity/boltz2_affinity` (Ab–Ag complex), `biologics/abodybuilder3` (cross-check)

---

## Install

```bash
uv add igfold
# IgFold pulls torch + fair-esm transitively (~3 GB)
# Optional: pip install pyrosetta  # for do_refine=True
```

---

## References

- Ruffolo J.A. et al. *Nature Comm.* 2023, 14, 2389 — "Fast, accurate antibody structure prediction from deep learning on massive set of natural antibodies"
- GitHub: https://github.com/Graylab/IgFold
- License: JHU Software License (free for academic / non-commercial use)
