# Skill: ligand_viz

**Tier:** L1 — atomic tool skill
**Category:** visualization
**Tool:** RDKit + pandas + (optional) mols2grid / plotly / umap-learn / matplotlib
**GPU required:** No — CPU only

---

## What this skill does

Renders a **self-contained HTML report** for a ligand set, combining 2D structure grid, sortable property table, scatter plot, and optional chemical-space UMAP — Datagrok / DataWarrior–style "browse a SAR table" output. The single HTML file has no external runtime dependencies (pulls jQuery / DataTables / plotly.js from CDN) and opens in any browser.

This is a **"publish results"** step at the end of an agentic workflow. Use it whenever the agent has produced a list of compounds with associated activity / docking / property data and the human needs to inspect them.

---

## When to call this skill

- After scoring a library with `docking/gnina`, `docking/diffdock`, `binding_affinity/boltz2_affinity`, `binding_affinity/deeppurpose` — render hits with their scores
- After enumerating with `enumeration/rdkit_enum` and filtering with `filtering/ligand_filter` — show what survived
- After ADMET prediction (`admet/rdkit_props`, `admet/admetlab3`) — visualise property distributions
- Whenever the user asks "show me the structures" / "what does the SAR look like"

---

## Inputs

| Field | Type | Required | Description |
|---|---|---|---|
| `smiles` | list[str] | ✓ | Compound SMILES (one per row) |
| `properties` | dict[str, list] | — | Column name → values aligned with `smiles`; e.g. `{"pIC50":[...], "logP":[...]}` |
| `plots` | list[str] | — | Subset of `["grid", "table", "scatter", "umap"]`; default `["grid", "table"]` |
| `scatter_x` | str | — | Property column for x-axis; defaults to first `properties` key |
| `scatter_y` | str | — | Property column for y-axis; defaults to second `properties` key |
| `color_by` | str | — | Property column to colour points (scatter / UMAP) |
| `title` | str | — | Report title shown in `<h1>`; default `"ChemMind ligand report"` |
| `output_dir` | str | — | Write report.html here; uses temp dir if None |

**Plot types**
- `"grid"` — 2D structure grid via `mols2grid` (interactive) if installed, else RDKit `Draw.MolsToGridImage` (static PNG)
- `"table"` — sortable / paginated DataTables HTML table
- `"scatter"` — Plotly scatter (interactive hover) if available, else matplotlib PNG
- `"umap"` — 2D UMAP of Morgan-2 fingerprints (Jaccard metric); requires `umap-learn`

---

## Outputs

| Field | Type | Description |
|---|---|---|
| `html_path` | str | Path to the rendered `report.html` |
| `output_dir` | str | Directory containing the report and any auxiliary PNGs |
| `n_compounds` | int | `len(smiles)` |
| `plots_generated` | list[str] | Plot types successfully rendered |
| `properties_indexed` | list[str] | Property column names included in the report |
| `skipped_plots` | list[dict] | `[{"plot": name, "reason": str}, ...]` for any plot that couldn't be drawn (missing optional dep, insufficient data, etc.) |

---

## Implementation

`tools/visualization/ligand_viz.py` — `run_ligand_viz()`

---

## Agent decision rules

- **Default end-of-workflow call**: `plots=["grid","table"]` covers 90% of "show me the hits" requests with zero optional deps
- **SAR plot**: add `"scatter"` with `scatter_x="logP"`, `scatter_y="pIC50"`, `color_by="QED"` when both axes are continuous and meaningful
- **Diversity check / cluster overview**: add `"umap"` only when the library is ≥20 compounds (smaller sets cluster trivially)
- **Don't over-spec**: if the user just asked "show structures", don't add UMAP/scatter — they noise up the report
- **Property alignment**: every list in `properties` must equal `len(smiles)`; the agent should compute / collect properties first (e.g. via `admet/rdkit_props`) then pass them in
- **Upstream**: any tool that produces SMILES + scores — `gnina`, `diffdock`, `boltz2_affinity`, `deeppurpose`, `rdkit_enum`, `ligand_filter`
- **Downstream**: human review (no further tool consumes this)

---

## Install

Core (already in repo):
```bash
uv add rdkit pandas
```

Optional — enable richer rendering:
```bash
uv add mols2grid plotly umap-learn matplotlib
```

The tool degrades gracefully — if an optional dep is missing, the affected plot lands in `skipped_plots` with a reason and the rest of the report still renders.

---

## References

- `mols2grid` — interactive HTML grid for RDKit molecules: https://github.com/cbouy/mols2grid
- `plotly` — interactive scientific charts: https://plotly.com/python/
- `umap-learn` — McInnes et al. 2018: https://umap-learn.readthedocs.io
- DataWarrior (inspiration) — Sander et al. *J. Chem. Inf. Model.* 2015, 55, 460
- License: BSD (RDKit), MIT (mols2grid), MIT (plotly), BSD (umap-learn)
