# Workflow: hit_gen_sbdd

**Tier:** L2 — workflow skill
**Purpose:** Generate novel, diverse small-molecule hit candidates conditioned on a 3D target binding pocket using structure-based diffusion. Covers the full pipeline from raw target information to a ranked, ADMET-filtered, docked hit list ready for experimental triage.

---

## When to invoke this workflow

- No known actives or scaffold for the target — truly de novo generation
- Want geometrically complementary molecules from the pocket shape
- Goal: 5–20 high-confidence candidates for experimental screening
- Keyword triggers: "generate hits", "de novo SBDD", "novel chemotypes", "structure-based drug design"

---

## Required inputs

| Input | Description |
|---|---|
| Target protein | PDB ID, UniProt ID, or sequence (resolved by `target_prep`) |
| Pocket specification | PDB ID with co-crystal ligand, or residue IDs for hotspot — resolved by `pocket_detect` if absent |
| n_designs | Optional — default 100 molecules from DiffSBDD |

---

## L1 skill sequence

### Step 1 — Target preparation (invoke `target_prep` sub-workflow)
```
pdb_fetch or esmfold/colabfold → pocket_detect
→ {pdb_path, pocket.center_x/y/z, pocket.residues, pocket.volume}
```
If pocket is already known (caller passes `pocket_center`), skip this step.

### Step 2 — De novo molecule generation
```
diffsbdd(
    protein_pdb=pdb_path,
    pocket_center=(center_x, center_y, center_z),
    pocket_residues=pocket.residues,
    n_molecules=100,
)
→ {smiles: list[str], sdf_path, n_valid, validity_rate}
```
Fallback: if DiffSBDD is unavailable or validity_rate < 0.5:
```
reinvent4(mode="sampling", n_steps=200, batch_size=64)
→ {smiles: list[str]}
```

### Step 3 — Physicochemical filtering
```
rdkit_props(smiles=generated_smiles)
→ {results: [{smiles, qed, sa_score, mw, logp, lipinski_pass, ...}]}
```
Keep molecules where:
- `lipinski_pass == True`
- `qed ≥ 0.4`
- `sa_score ≤ 5`
- `mw ≤ 600`

Typical retention: 40–70% of generated molecules.

### Step 4 — Structural alert filtering
```
ligand_filter(smiles=filtered_smiles, filter_sets="default")
→ {results: [{smiles, passed, failures}]}
```
Keep `passed == True` only. Typical retention: 80–95% of physicochemical-passed set.

### Step 5 — Docking (broad triage)
Dock all ADMET-passed molecules (up to 200; if > 200, sort by QED descending and take top 200):
```
gnina(
    protein_pdb=pdb_path,
    ligand_smiles=passed_smiles,       # batch mode
    center_x=center_x, center_y=center_y, center_z=center_z,
    box_size=box_size,                 # sqrt(pocket.volume) * 1.5
    exhaustiveness=8,
)
→ {results: [{smiles, docking_score_kcal_mol, pose_pdb}]}
```
Sort by `docking_score_kcal_mol` ascending (more negative = better). Take top 20.

Fallback if gnina unavailable:
```
vina(protein_pdb=..., ligand_smiles=..., center_x=..., ...)
```

### Step 6 — Affinity rescoring (top-20)
```
boltz2_affinity(
    protein_fasta=target_sequence,
    ligand_smiles=top20_smiles,        # batch mode
)
→ {predictions: [{smiles, affinity_kcal_mol, iptm_score}]}
```
Sort by `affinity_kcal_mol` ascending. Take top 10.

### Step 7 — ADMET profiling (top-10)
```
admetlab3(smiles=top10_smiles)
→ {predictions: [{Caco2, hERG, ..., BBB, T12, ...}]}
```
Flag molecules with predicted hERG inhibition (IC50 < 10 µM) or Caco2 < 10 nm/s.

### Step 8 — Retrosynthesis (top-1)
```
askcos(smiles=best_smiles, n_steps=3)
→ {routes: [{steps, buyable_leaves, overall_score}]}
```
Report whether the top candidate has a ≤ 3-step route with buyable precursors.

### Step 9 — Visualisation
```
ligand_viz(
    smiles=top10_smiles,
    properties={"docking_score": [...], "qed": [...], "pIC50_pred": [...]},
    plots=["grid", "table", "scatter"],
    scatter_x="mw", scatter_y="docking_score",
    title="SBDD hits — <target name>",
)
→ {html_path}
```

---

## Decision branches

```
validity_rate < 0.5 after DiffSBDD?
  → switch to reinvent4(sampling) + rdkit_props filter

Fewer than 10 molecules pass ADMET?
  → increase n_molecules to 300 and re-run DiffSBDD with different seed

gnina docking fails (binary not found)?
  → fall back to vina; note lower accuracy in output

boltz2_affinity unreachable (NIM API error)?
  → report gnina scores only; flag that Boltz-2 rescoring was skipped

askcos API not configured?
  → report SA scores as synthesisability proxy; note retrosynthesis skipped
```

---

## Expected outputs

Ranked table of ≤ 10 candidates:
```
Rank | SMILES | Docking (kcal/mol) | Boltz-2 ΔG | QED | SA | Lipinski | ADMET flags
```
Plus:
- `ligand_viz` HTML report
- ASKCOS route for top-1 candidate
- Confidence statement on each prediction

---

## Success criteria

- ≥ 5 candidates with docking score < −8 kcal/mol
- ≥ 3 candidates with Boltz-2 affinity < −9 kcal/mol and `iptm_score > 0.6`
- All final candidates pass Lipinski, QED ≥ 0.5, SA ≤ 5
- Top-1 candidate has a feasible retrosynthetic route (ASKCOS ≥ 1 route with buyable leaves)

---

## References

- DiffSBDD: Schneuing et al., *ICLR* 2023 — `diffsbdd`
- GNINA: McNutt et al., *J. Cheminform.* 2021 — `gnina`
- Boltz-2: MIT 2025 — `boltz2_affinity`
- ADMETlab3: Xiong et al., 2021 — `admetlab3`
