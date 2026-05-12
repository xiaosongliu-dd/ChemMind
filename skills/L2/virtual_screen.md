# Workflow: virtual_screen

**Tier:** L2 — workflow skill
**Purpose:** Screen commercially available compound libraries (Enamine REAL, ZINC22) against a target pocket to identify buyable, drug-like hits ready for experimental purchase and testing. Prioritises throughput — hundreds of thousands of compounds docked efficiently — over molecular novelty.

---

## When to invoke this workflow

- Want to find commercially available hits quickly (same-week delivery, no synthesis)
- Starting from a large chemical library (Enamine REAL, ZINC22, ChEMBL subset)
- Complement to `hit_gen_sbdd` — run both to get novel + buyable hit series
- Keyword triggers: "screen library", "VS campaign", "find buyable actives", "ZINC", "Enamine", "virtual screening"

---

## Required inputs

| Input | Description |
|---|---|
| Target protein | PDB ID, UniProt, or sequence |
| Library | `"enamine_real"`, `"zinc22"`, or a SMILES list |
| Query SMILES (similarity search) | Optional — if a reference compound is known |
| Property filters | Optional — default: Lipinski, QED ≥ 0.4, SA ≤ 5 |

---

## L1 skill sequence

### Step 1 — Target preparation
```
target_prep sub-workflow
→ {pdb_path, pocket.center_x/y/z, pocket.volume, pocket.residues}
```

### Step 2 — Library retrieval
```
library_search(
    library="enamine_real",            # or "zinc22"
    query_smiles=reference_smiles,     # similarity search; omit for random subset
    n_results=10000,
    similarity_threshold=0.4,
)
→ {smiles: list[str], catalog_ids: list[str], prices: list[float]}
```
If no reference compound: use `chembl(target_id=..., standard_type="IC50", limit=100)` to get known actives as similarity seeds.

### Step 3 — Rapid physicochemical pre-filter
```
rdkit_props(smiles=library_smiles)
→ {results: [{smiles, qed, sa_score, mw, logp, lipinski_pass}]}
```
Hard filters: `lipinski_pass == True`, `qed ≥ 0.4`, `sa_score ≤ 5`.
Expected retention: ~40–60% of library. Proceed with up to 5,000 molecules.

### Step 4 — Structural alert filter
```
ligand_filter(smiles=filtered_smiles, filter_sets="default")
```
Remove PAINS + BRENK flagged compounds. Expected retention: ~85–90%.

### Step 5 — High-throughput docking
For large sets (> 500 compounds), use AutoDock-GPU (fastest):
```
autodock_gpu(
    protein_pdb=pdb_path,
    ligand_smiles=passed_smiles,
    center_x=center_x, center_y=center_y, center_z=center_z,
    box_size=box_size,
    nrun=20,
)
→ {results: [{smiles, binding_energy, ki_nm}]}
```
Sort by `binding_energy` ascending. Take top 50.

Fallback (if AutoDock-GPU unavailable):
```
vina(protein_pdb=..., ligand_smiles=..., exhaustiveness=4)  # fast, lower accuracy
```

### Step 6 — Medium-throughput rescoring (top-50)
```
gnina(
    protein_pdb=pdb_path,
    ligand_smiles=top50_smiles,
    center_x=..., ...,
    exhaustiveness=16,              # higher accuracy
)
→ {results: [{smiles, docking_score_kcal_mol, cnn_score}]}
```
Sort by `docking_score_kcal_mol`. Take top 20.

### Step 7 — Affinity rescoring (top-20)
```
boltz2_affinity(protein_fasta=target_sequence, ligand_smiles=top20_smiles)
→ {predictions: [{smiles, affinity_kcal_mol, iptm_score}]}
```
Final ranking: sort by `affinity_kcal_mol`. Take top 10.

### Step 8 — ADMET profiling (top-10)
```
admetlab3(smiles=top10_smiles)
```
Flag hERG, Caco2, solubility issues.

### Step 9 — Retrosynthesis (top-1)
```
askcos(smiles=best_smiles, n_steps=2)
```
Note: library compounds are pre-synthesised — ASKCOS confirms feasibility if analogues are needed.

### Step 10 — Visualisation
```
ligand_viz(
    smiles=top10_smiles,
    properties={"docking_score": [...], "qed": [...], "catalog_id": [...]},
    plots=["grid", "table", "scatter"],
    title="VS hits — <target name>",
)
```

---

## Decision branches

```
Library > 5,000 after pre-filter?
  → use autodock_gpu for Step 5; vina is too slow for > 500

No reference SMILES for similarity search?
  → chembl search for known actives as seeds; or use structural fingerprint
    of a co-crystal ligand from pdb_fetch output

autodock_gpu unavailable?
  → vina with exhaustiveness=4 (fast mode); note lower accuracy

boltz2_affinity API error?
  → skip Step 7; report gnina scores only

< 5 compounds pass all filters?
  → loosen QED threshold to 0.35 or expand similarity search to 0.35
```

---

## Expected outputs

Ranked table of ≤ 10 buyable compounds:
```
Rank | SMILES | Catalog ID | Price | AutoDock-GPU (kcal/mol) | GNINA | Boltz-2 ΔG | QED | SA | ADMET flags
```
Plus `ligand_viz` HTML report.

---

## Success criteria

- ≥ 5 candidates with GNINA docking score < −8 kcal/mol
- All final candidates pass Lipinski and PAINS/BRENK filters
- ≥ 1 candidate available in the library for immediate purchase
- Boltz-2 `iptm_score > 0.55` for top-3 candidates

---

## References

- AutoDock-GPU: Santos-Martins et al., *J. Chem. Theory Comput.* 2021 — `autodock_gpu`
- GNINA: McNutt et al., *J. Cheminform.* 2021 — `gnina`
- Enamine REAL: https://enamine.net/compound-libraries/real-compounds
- ZINC22: Tingle et al., *J. Chem. Inf. Model.* 2023 — `library_search`
