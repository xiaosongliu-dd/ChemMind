# Workflow: fep_campaign

**Tier:** L2 — workflow skill
**Purpose:** Accurately rank the binding affinities of a lead series (3–10 compounds) using relative binding free energy perturbation (RBFE). Provides the most physically accurate ΔΔG estimates available from computation, suitable for rank-ordering before synthesis decisions.

---

## When to invoke this workflow

- Top compounds from `hit_gen_sbdd`, `virtual_screen`, or `lead_opt` need accurate ΔΔG ranking
- ≤ 10 compounds in the lead series (FEP scales poorly beyond this in a single campaign)
- Docking scores or Boltz-2 ΔG spread is narrow (< 2 kcal/mol) and you need better discrimination
- Keyword triggers: "FEP", "free energy", "RBFE", "rank affinity accurately", "perturbation"

---

## Required inputs

| Input | Description |
|---|---|
| Lead series SMILES | 3–10 structurally similar compounds (congeneric series preferred for RBFE) |
| Docked poses (PDB) | From `gnina` or `vina` — required as starting geometry for RBFE |
| Target PDB | Apo or holo protein structure |
| Reference compound | The best-characterised compound in the series (anchor for RBFE network) |

---

## L1 skill sequence

### Step 1 — Input validation
```
rdkit_props(smiles=lead_series_smiles)
```
Confirm: all compounds are congeneric (same core scaffold, ≤ 10 heavy-atom differences). If not congeneric, RBFE convergence will be poor — recommend `boltz2_affinity` instead.

Verify: docked poses exist for all compounds in the series. If missing:
```
gnina(receptor_pdb=target_pdb, ligand_smiles=missing_smiles, exhaustiveness=16)
→ {pose_pdb}
```

### Step 2 — Protein–ligand MD equilibration
For each compound (or just the reference compound for efficiency):
```
openmm(
    structure_pdb=docked_pose_pdb,    # protein + ligand complex PDB
    ensemble="NPT",
    temperature_k=300.0,
    n_steps=500000,                   # 1 ns equilibration
    timestep_fs=2.0,
    forcefield="amber14-all.xml",
)
→ {trajectory_path, final_structure, energy_path}
```

### Step 3 — Equilibration verification
```
mdanalysis(
    topology=docked_pose_pdb,
    trajectory=trajectory_path,
    analyses=["rmsd", "rmsf"],
    selection="protein and backbone",
)
→ {rmsd: {step, rmsd_angstrom}, rmsf: {resids, rmsf_angstrom}}
```
Gate: protein backbone RMSD plateau < 2 Å over last 500 ps. If not plateaued, extend simulation to 2 ns (n_steps=1,000,000) and re-check.

### Step 4 — RBFE calculation (OpenFE)
```
fep_openfe(
    protein_pdb=equilibrated_structure,
    ligand_sdf_dir=ligand_sdf_dir,           # directory of SDF files, one per ligand
    reference_ligand_sdf=ref_sdf,            # SDF of the reference/anchor compound
    n_replicas=3,
    n_lambda=11,
    work_dir=workdir,
)
→ {
    ddg_predictions: [{smiles_A, smiles_B, ddg_kcal_mol, uncertainty}],
    absolute_dg:     [{smiles, dg_kcal_mol, uncertainty}],
    perturbation_map: <Lomap network>
}
```
Convergence check: uncertainty < 0.5 kcal/mol per edge. Re-run longer if higher.

### Step 5 — Rank-ordering
Combine RBFE ΔΔG predictions with the reference compound's Boltz-2 ΔG to build absolute ΔG estimates:
```
final_dG[compound_i] = boltz2_dG[reference] + ddg_predictions[reference → compound_i]
```
Sort ascending (lower = tighter binder).

### Step 6 — ADMET gate on top-3
```
admetlab3(smiles=top3_smiles)
```
Final developability check before synthesis decision.

### Step 7 — Retrosynthesis for top-3
```
askcos(smiles=each_of_top3, n_steps=3)
→ {routes}
```

### Step 8 — Report
Text table (no small-molecule ligand_viz at this stage):
```
Rank | SMILES | RBFE ΔG (kcal/mol) | Uncertainty | Boltz-2 ΔG | QED | Synthesis route
```
Plus: perturbation map diagram (Lomap network topology from fep_openfe output).

---

## Decision branches

```
Series is non-congeneric (scaffold hops)?
  → RBFE convergence is unreliable; use boltz2_affinity instead for ranking
  → Report: "Non-congeneric series — using Boltz-2 for ranking; FEP not applied"

openmm equilibration RMSD doesn't plateau after 2 ns?
  → Check for steric clashes in the input pose (gnina may have produced a bad pose)
  → Re-dock with higher exhaustiveness (32) and re-equilibrate

fep_openfe uncertainty > 1.0 kcal/mol on key edges?
  → Extend simulation_time_ns to 10.0 for those edges
  → If still > 1.0, flag as unconverged; report best estimate with error bar

< 3 compounds in series?
  → FEP is not worth the cost; use boltz2_affinity for pairwise ranking
  → Report: "Series too small for FEP; using Boltz-2 for ranking"

> 10 compounds in series?
  → Run two FEP cycles: pre-screen with boltz2_affinity → pick top 8 → run FEP
```

---

## Expected outputs

Ranked table of lead series compounds:
```
Rank | SMILES | RBFE ΔG ± σ (kcal/mol) | Boltz-2 ΔG | QED | SA | ADMET flags | Synthesis
```
Plus perturbation network topology and convergence metrics per edge.

---

## Success criteria

- All RBFE edges have uncertainty < 0.5 kcal/mol
- Perturbation network is fully connected (no isolated nodes)
- Top-ranked compound also passes ADMET and has feasible synthesis
- Rank ordering is consistent between RBFE ΔG and Boltz-2 ΔG for top-3 (Spearman ρ > 0.7)
- Consistency check: if RBFE and Boltz-2 rankings disagree on #1 vs #2, note the discrepancy and recommend experimental validation

---

## References

- OpenFE: Ganguly et al., *J. Chem. Theory Comput.* 2023 — `fep_openfe`
- OpenMM: Eastman et al., *PLOS Comput. Biol.* 2017 — `openmm`
- MDAnalysis: Michaud-Agrawal et al., 2011 — `mdanalysis`
- LoMap perturbation maps: Liu et al., 2013
- Boltz-2 — `boltz2_affinity` (complementary cross-check)
