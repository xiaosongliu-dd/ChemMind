# Workflow: target_prep

**Tier:** L2 — workflow skill
**Purpose:** Prepare a biological target for computational drug design — fetch or predict the 3D structure, identify druggable pockets, and produce the coordinates and metadata needed by all downstream docking and generation workflows.

---

## When to invoke this workflow

- Starting any drug design campaign when the pocket location is unknown
- Input is a PDB ID, UniProt accession, gene name, or bare amino acid sequence
- Implicitly required upstream of `hit_gen_sbdd`, `virtual_screen`, `lead_opt`, `peptide_design`, or any docking step

---

## Required inputs

| Input | Source | Notes |
|---|---|---|
| PDB ID or UniProt ID or AA sequence | User | One of the three required |
| Resolution cutoff | Optional | Default: < 2.5 Å (X-ray), < 3.5 Å (CryoEM) |
| Pocket hint (residues or co-crystal ligand) | Optional | If known, speeds up pocket step |

---

## L1 skill sequence

### Step 1 — Obtain a 3D structure

**If PDB ID known:**
```
pdb_fetch(pdb_id=<ID>, query_type="structure")
→ {pdb_path, resolution, experiment_type, ligands, chain_ids}
```
Check `ligands` — if a drug-like ligand is present (not solvent), its coordinates define the pocket; skip Step 3.

**If UniProt / gene name known, no PDB ID:**
```
pdb_fetch(uniprot_id=<accession>, query_type="search")
→ {pdb_ids: [...]}
```
Pick the lowest-resolution number entry (highest quality). Fetch it with `query_type="structure"`.

**If sequence only:**
```
esmfold(sequence=<aa_seq>, use_nim=True)
→ {pdb_path, mean_plddt}
```
Gate: `mean_plddt ≥ 70`. If below, switch to:
```
colabfold(sequences=<aa_seq>, n_models=5, use_templates=True)
→ {best_pdb_path, mean_plddt}
```

### Step 2 — Quality check

- X-ray / CryoEM: resolution ≤ 2.5 Å (X-ray) or ≤ 3.5 Å (CryoEM)
- Predicted: `mean_plddt ≥ 70` across the binding-site region
- If quality fails: try an alternate PDB entry or escalate to ColabFold

### Step 3 — Pocket detection

**No co-crystal ligand:**
```
pocket_detect(protein_pdb=<pdb_path>, tool="p2rank")
→ {best_pocket: {score, center_x, center_y, center_z, volume, residues}}
```
If `best_pocket.score < 0.4`, cross-check with fpocket:
```
pocket_detect(protein_pdb=<pdb_path>, tool="fpocket")
```
Use the pocket with the higher score that is spatially consistent between both tools.

**Co-crystal ligand present:** compute the ligand centroid from PDB HETATM records directly.

---

## Decision branches

```
Has PDB ID?
  yes → pdb_fetch structure
  no → Has UniProt ID?
         yes → pdb_fetch search → pick best PDB → fetch structure
         no → esmfold (pLDDT ≥ 70?) → yes: done
                                      → no: colabfold

Has co-crystal ligand in structure?
  yes → use ligand centroid as pocket centre (skip pocket_detect)
  no  → pocket_detect (P2Rank) → score < 0.4 → cross-check fpocket
```

---

## Outputs

```
pdb_path:    str            # absolute path to prepared PDB
source:      str            # "pdb_fetch" | "esmfold" | "colabfold"
resolution:  float | None   # Å; None for predicted structures
mean_plddt:  float | None   # predicted-structure confidence
pocket:
  center_x:  float          # Å
  center_y:  float          # Å
  center_z:  float          # Å
  volume:    float          # Å³
  residues:  list[str]      # ["A48", "A86", ...]
  score:     float          # P2Rank druggability score (0–1)
chain_ids:   list[str]
ligands:     list[dict]     # co-crystal ligands
```

---

## Success criteria

- `pdb_path` is a valid, non-empty PDB file
- `pocket.score ≥ 0.4` or ligand-defined pocket used
- `pocket.center_x/y/z` are finite non-zero floats
- Resolution / pLDDT quality gate passed

---

## Failure recovery

| Failure | Recovery |
|---|---|
| PDB not found | Try alternate PDB ID; search by UniProt |
| ESMFold pLDDT < 60 | Switch to ColabFold (`use_templates=True`) |
| P2Rank score < 0.3 | Run fpocket; check literature; try allosteric pockets (rank 2–3) |
| No pocket found | Target may be disordered; consider PPI inhibition or allosteric approach |

---

## Downstream handoff

- `pdb_path` → DiffSBDD, vina, gnina, diffdock, rfantibody, rfdiffusion_pep
- `pocket.center_x/y/z` → docking box centre
- `pocket.residues` → `hotspot_residues` for rfantibody, rfdiffusion_pep, bindcraft
- `pocket.volume` → box size estimate: `box_size ≈ sqrt(volume) * 1.5` Å

---

## References

- P2Rank: Krivák & Hoksza, *J. Cheminform.* 2018 — `pocket_detect`
- ESMFold: Lin et al., *Science* 2023 — `esmfold`
- ColabFold: Mirdita et al., *Nat. Methods* 2022 — `colabfold`
- RCSB PDB: Berman et al., *Nucleic Acids Res.* 2000 — `pdb_fetch`
