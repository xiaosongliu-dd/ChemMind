# Skill: pocket_detect

**Tier:** L1 — atomic tool skill
**Category:** pocket
**Tool:** P2Rank (primary) / fpocket (fallback) — druggable binding pocket detection
**GPU required:** No — CPU-only; P2Rank ~5–30 s per structure

---

## What this skill does

`pocket_detect` identifies and ranks druggable binding pockets on a target protein structure. P2Rank uses a random-forest model trained on protein surface features to score pockets; fpocket uses purely geometric Voronoi clustering. Both return ranked pocket centres, volumes, and the lining residues — the key inputs for docking box setup.

Key capabilities:
- Ranked pockets sorted by druggability score
- Pocket centre coordinates (`center_x/y/z`) ready for docking box
- Pocket volume for size filtering
- Residue list for hotspot specification
- P2Rank (ML) and fpocket (geometry) backends

---

## When to call this skill

- Need to identify the binding site before docking (no crystal ligand present)
- Characterising a newly folded structure from `esmfold` or `colabfold`
- Selecting hotspot residues for `rfantibody` or `bindcraft`
- Ranking multiple candidate pockets for multi-site analysis

Do NOT use if:
- Crystal structure with co-crystallised ligand already defines the pocket → use its coordinates directly
- Running AlphaFold predictions with known template pockets

---

## Inputs

| Field | Type | Required | Description |
|---|---|---|---|
| `protein_pdb` | str | ✓ | Path to protein structure PDB |
| `tool` | str | — | `"p2rank"` (default) or `"fpocket"` |
| `p2rank_bin` | str | — | P2Rank binary name / path (default `"prank"`) |
| `fpocket_bin` | str | — | fpocket binary name / path (default `"fpocket"`) |
| `output_dir` | str | — | Results directory; temp dir if None |

---

## Outputs

| Field | Type | Description |
|---|---|---|
| `pockets` | list[dict] | Pockets sorted by score descending |
| `n_pockets` | int | Number of pockets found |
| `best_pocket` | dict \| None | Top-ranked pocket |
| `tool` | str | Backend used (`"p2rank"` or `"fpocket"`) |
| `output_dir` | str | Path to tool output directory |

Each pocket dict: `{rank, score, center_x, center_y, center_z, volume, residues}`.

Raises `RuntimeError("P2Rank failed: ...")` on non-zero exit.

---

## Implementation

`tools/structure/pocket_detect.py` — `run_pocket_detect()`

P2Rank mode: calls `prank predict -f <pdb> -o <outdir>`, parses `*_predictions.csv`.
fpocket mode: calls `fpocket -f <pdb>`, parses `*_out/*_info.txt`.

---

## Agent decision rules

- **Always use P2Rank by default**: it outperforms fpocket on the DUD-E benchmark and is faster
- **Box setup**: use `best_pocket["center_x/y/z"]` as docking box centre; set box size = `sqrt(volume) * 1.5` (rough heuristic)
- **Hotspot residues**: pass `best_pocket["residues"]` to `rfantibody` or `bindcraft` as `hotspot_residues`
- **Multi-pocket screening**: if top pocket has `score < 0.4`, the target may be undruggable; consider allosteric sites (pocket rank 2–3)
- **Cross-validate**: run both P2Rank and fpocket; consensus pockets (appear in both top-3) are more reliable
- **Upstream**: `structure/esmfold`, `structure/colabfold`, `data/pdb_fetch`
- **Downstream**: `docking/vina`, `docking/gnina`, `docking/diffdock`, `biologics/rfantibody`

---

## Install

```bash
# P2Rank
wget https://github.com/rdk/p2rank/releases/latest/download/p2rank_<version>.tar.gz
tar xf p2rank_*.tar.gz && export PATH=$PATH:$(pwd)/p2rank_*/

# fpocket
sudo apt install fpocket  # Debian/Ubuntu
# or: brew install fpocket  # macOS
```

---

## References

- Krivák R. & Hoksza D. *J. Cheminform.* 2018, 10, 39 — "P2Rank: machine learning based tool for rapid and accurate prediction of ligand binding sites from protein structure"
- Le Guilloux V. et al. *BMC Bioinform.* 2009, 10, 168 — "Fpocket: an open source platform for ligand pocket detection"
- GitHub P2Rank: https://github.com/rdk/p2rank
- GitHub fpocket: https://github.com/Discngine/fpocket
- License: MIT (P2Rank), BSD (fpocket)
