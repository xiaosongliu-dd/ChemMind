# Workflow: lead_opt

**Tier:** L2 — workflow skill
**Purpose:** Systematically optimise a known hit or lead compound — generate a focussed analogue library, filter for ADMET compliance, rank by predicted affinity, and select candidates for synthesis or FEP. Designed for iterative medicinal chemistry cycles.

---

## When to invoke this workflow

- A hit compound is known (SMILES provided) and needs improvement
- Goal is improving potency, selectivity, metabolic stability, or reducing toxicity flags
- SAR data exists (ChEMBL actives) and should be used to guide enumeration
- Keyword triggers: "optimise", "improve affinity", "scaffold decoration", "lead series", "SAR", "analogues", "R-group"

---

## Required inputs

| Input | Description |
|---|---|
| Hit SMILES | Starting compound for optimisation |
| Target PDB | For docking — resolved by `target_prep` if absent |
| Optimisation goal | e.g. "improve potency", "reduce hERG", "improve solubility" |
| SAR data (optional) | Known active SMILES from `chembl` for seed enumeration |

---

## L1 skill sequence

### Step 1 — Baseline characterisation
```
rdkit_props(smiles=hit_smiles)
→ {qed, sa_score, mw, logp, hbd, hba, tpsa, lipinski_pass, fingerprint}
```
Identify which physicochemical properties need improvement. This baseline guides Step 3.

### Step 2 — Target preparation (if not already done)
```
target_prep sub-workflow
→ {pdb_path, pocket.center_x/y/z, pocket.volume}
```

### Step 3 — SAR context (optional but recommended)
```
chembl(target_id=<CHEMBL_ID>, standard_type="IC50", limit=100)
→ {records: [{smiles, pchembl_value, standard_value}]}
```
Use the top-20 most potent actives (highest `pchembl_value`) as R-group enumeration seeds in Step 4.

### Step 4 — Analogue generation
```
rdkit_enum(
    smiles=hit_smiles,
    mode="brics",             # BRICS decomposition + R-group replacement
    rgroup_smiles=sar_actives,
    max_compounds=2000,
)
→ {smiles: list[str], n_products}
```
If goal is scaffold hop (not just decoration):
```
reinvent4(
    mode="sampling",
    scaffold_smarts=murcko_scaffold,
    n_steps=100,
    batch_size=64,
)
→ {smiles: list[str]}
```

### Step 5 — Physicochemical filtering
```
rdkit_props(smiles=analogue_smiles)
```
Apply target-goal-dependent filters:
- **Potency**: `mw ≤ 500`, `logp ≤ 5`, `qed ≥ 0.5`
- **Oral bioavailability**: add `tpsa ≤ 90`, `hbd ≤ 3`
- **CNS penetration**: add `mw ≤ 450`, `tpsa ≤ 90`, `logp 1–5`

### Step 6 — Structural alert filtering
```
ligand_filter(smiles=filtered_smiles, filter_sets="strict")
→ pass/fail per molecule
```
Use `"strict"` preset (PAINS + BRENK + NIH + ZINC) for lead opt (more stringent than hit gen).

### Step 7 — Docking (≤ 500 molecules)
```
gnina(
    receptor_pdb=pdb_path,
    ligand_smiles=passed_smiles,
    center_x=..., center_y=..., center_z=...,
    box_size=box_size,
    exhaustiveness=16,            # higher than VS — quality over throughput
)
→ {results sorted by docking_score_kcal_mol}
```
Take top 30.

### Step 8 — Affinity rescoring and Boltz-2 (top-30)
```
boltz2_affinity(protein_fasta=target_fasta, ligand_smiles=top30_smiles)
→ {predictions sorted by affinity_kcal_mol}
```
Take top 10.

### Step 9 — Full ADMET profiling (top-10)
```
admetlab3(smiles=top10_smiles)
```
Flag per optimisation goal:
- **hERG reduction goal**: keep `hERG_IC50_pred > 10 µM`
- **Metabolic stability**: keep `CLhep_pred < 10 mL/min/kg`
- **Solubility**: keep `LogS_pred > −4`

### Step 10 — MD equilibration of top-3 (optional, high-confidence leads)
```
openmm(structure_pdb=docked_pose_pdb, n_steps=500000, ensemble="NVT")
→ {trajectory_path, final_structure}

mdanalysis(topology=..., trajectory=trajectory_path, analyses=["rmsd", "rmsf"])
→ {rmsd, rmsf}
```
If backbone RMSD of docked ligand > 3 Å over 1 ns, the pose is likely unstable.

### Step 11 — Retrosynthesis (top-3)
```
askcos(smiles=top3_smiles, n_steps=3)
→ {routes}
```
Discard candidates with no buyable route unless the synthesis is strategically important.

### Step 12 — Visualisation
```
ligand_viz(
    smiles=top10_smiles,
    properties={"delta_qed": [...], "docking_score": [...], "boltz2_dG": [...]},
    plots=["grid", "table", "scatter", "umap"],
    color_by="docking_score",
    title="Lead opt — <hit SMILES truncated>",
)
```

---

## Decision branches

```
Optimisation goal is "reduce hERG"?
  → add hERG_IC50_pred filter in Step 9 before further shortlisting

Fewer than 20 analogues pass ADMET?
  → loosen one property gate (e.g. qed ≥ 0.4) or run reinvent4 sampling
    with higher temperature (1.2) to explore more chemical space

FEP warranted for top-3?
  → goal is accurate ΔΔG ranking → invoke fep_campaign sub-workflow on top-3

SAR data not available (no ChEMBL ID)?
  → skip Step 3; run BRICS enumeration on hit alone
```

---

## Expected outputs

Ranked table of ≤ 10 optimised candidates vs the parent hit:
```
Rank | SMILES | Δ(docking vs parent) | QED | SA | hERG flag | Synthesis route
```
Plus `ligand_viz` HTML report.

---

## Success criteria

- ≥ 1 candidate with predicted affinity ≥ 10-fold better than starting hit
- All final candidates pass strict PAINS/BRENK and target-specific ADMET criteria
- ≥ 1 candidate with a feasible ≤ 3-step synthesis route
- SA score improvement vs parent (lower is better)

---

## References

- BRICS enumeration: Degen et al., *ChemMedChem* 2008 — `rdkit_enum`
- REINVENT4: Loeffler et al., *J. Cheminform.* 2024 — `reinvent4`
- ADMETlab3 — `admetlab3`
- OpenFE / RBFE — `fep_openfe` (invoked via `fep_campaign`)
