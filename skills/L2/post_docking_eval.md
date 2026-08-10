# Workflow: post_docking_eval

**Tier:** L2 — workflow skill
**Purpose:** Re-rank a docked compound set using multi-component consensus scoring
(Vina + GNINA + ProLIF contacts + Boltz-2 ΔG) to reduce false positives from any
single score. Run after the primary docking step and before ADMET filtering.

---

## When to invoke this workflow

- After Step 6 (docking) in `virtual_screen` for any library > 50 compounds
- After Step 7 (GNINA rescoring) when ProLIF contacts and Boltz-2 are also available
- Keyword triggers: "consensus ranking", "re-rank docked compounds", "post-docking",
  "interaction fingerprint scoring", "combine docking scores"

---

## Required inputs

| Input | Description |
|---|---|
| Receptor PDB | Protein structure (same as used for docking) |
| Docked SDF dir | Directory containing one SDF per compound (from vina/gnina output) |
| SMILES list | Compounds in same order as SDF files |
| Vina scores | list[float] — kcal/mol from Step 6 docking |
| GNINA CNN affinities | list[float] — from Step 7 rescoring (optional) |
| Boltz-2 ΔG | list[float] — from Step 8 affinity prediction (optional) |

---

## L1 skill sequence

### Step 1 — ProLIF interaction fingerprints

For each docked pose SDF in `ligand_sdf_dir/`:
```python
prolif(
    protein_pdb=receptor_pdb,
    ligand_sdf=f"{ligand_sdf_dir}/{compound_i}.sdf",
    output_dir=workdir,
)
→ {interactions, summary: {n_contacts}, fingerprint, n_poses}
```

Collect `summary["n_contacts"]` for each compound into `contact_counts: list[int]`.

**Parallelise:** run prolif on all SDFs concurrently if > 20 compounds.

**If no SDF dir (SMILES-only docking):** skip Step 1 and set `prolif_contact_counts=None`.

---

### Step 2 — Consensus ranking

```python
equiscore(
    smiles=compound_smiles,
    vina_scores=vina_scores,            # from Step 6
    gnina_affinities=gnina_affinities,  # from Step 7; None if not run
    boltz2_dg=boltz2_dg,               # from Step 8; None if not run
    prolif_contact_counts=contact_counts,  # from Step 1; None if no SDF
    weights={
        "vina":   0.35,
        "gnina":  0.35,
        "boltz2": 0.15,
        "prolif": 0.15,
    },
)
→ {
    ranked: [{smiles, vina_score, gnina_affinity, boltz2_dg,
              prolif_contacts, consensus, rank}],
    n_scored,
    weights_used,
}
```

Take `ranked[:20]` — these are the consensus top-20 for ADMET.

---

### Step 3 — Interaction profile report

For the top-10 compounds, emit a structured SAR summary:

```
Rank | SMILES | Vina | GNINA | Boltz-2 ΔG | Contacts | Key residues
```

Highlight compounds that engage all three of:
1. ≥ 1 H-bond (from prolif interactions)
2. Hydrophobic contact with gatekeeper residue
3. Fit within the box (prolif n_contacts > 3)

These are the highest-confidence hits.

---

## Decision branches

```
No SDF files available (SMILES-only docking)?
  → skip prolif step; run equiscore with vina + gnina only
  → weight redistribution is automatic

Only one score component (e.g. Vina only)?
  → equiscore with weights={"vina": 1.0}
  → equivalent to simple sort — still useful for pipeline uniformity

< 5 compounds pass consensus > 0.5?
  → lower docking exhaustiveness and re-dock with more poses
  → OR loosen QED/SA filters in virtual_screen Step 4

All top-10 have prolif_contacts < 3?
  → pocket may be too small or ligands not entering — review docking box
  → check docking box: box_size may need to be increased
```

---

## Expected outputs

Consensus-ranked table:
```
Rank | SMILES | Vina (kcal/mol) | GNINA CNN | Boltz-2 ΔG | ProLIF contacts | Consensus score
```

Plus per-residue interaction table for top-3.

---

## Success criteria

- ≥ 10 compounds scored with ≥ 2 active score components
- Top-1 compound has consensus score > 0.65
- Top-3 compounds each have ProLIF n_contacts ≥ 3

---

## References

- ProLIF: Bouysset & Fiorucci, *J. Cheminform.* 2021 — `prolif`
- Consensus scoring review: Clark, *J. Med. Chem.* 2005 — motivation for multi-component ranking
