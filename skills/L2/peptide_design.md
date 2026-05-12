# Workflow: peptide_design

**Tier:** L2 — workflow skill
**Purpose:** Design short peptide binders (8–50 residues) against a target protein epitope using RFdiffusion backbone generation, ProteinMPNN sequence design, ColabFold structure validation, and Boltz-2 complex scoring.

---

## When to invoke this workflow

- Targeting a protein–protein interface, surface groove, or shallow pocket
- Peptide therapeutic needed (< 50 residues; no full antibody)
- Goal: 5–10 peptide sequences with predicted nanomolar binding
- Keyword triggers: "design peptide binder", "stapled peptide", "cyclic peptide", "peptide drug", "macrocycle", "mini-protein binder"

---

## Required inputs

| Input | Description |
|---|---|
| Target PDB | Antigen / target structure; resolved by `pdb_fetch` if needed |
| Hotspot residues | Target surface residues for the peptide to engage; from `pocket_detect` if absent |
| Peptide length range | `(min_length, max_length)` in residues; default (15, 40) |
| Peptide type | `"linear"` (default), `"cyclic"`, or `"stapled"` — affects backbone topology |

---

## L1 skill sequence

### Step 1 — Target preparation
```
pdb_fetch(pdb_id=<target_pdb_id>, query_type="structure")
→ {pdb_path, chain_ids, ligands}
```
If hotspot residues are unknown:
```
pocket_detect(protein_pdb=pdb_path, tool="p2rank")
→ best_pocket.residues  # use as hotspot_residues
```

### Step 2 — Peptide backbone design (RFdiffusion)
```
rfdiffusion_pep(
    target_pdb=target_pdb_path,
    hotspot_residues=hotspot_residues,
    binder_length=(min_length, max_length),
    n_designs=50,
    use_nim=True,
)
→ {designs: [{pdb_path, mean_plddt, rank}], best_pdb_path}
```
Keep designs with `mean_plddt > 70`. Target 30–45 passing designs.

For a broader design space or higher hit rate, use BindCraft (end-to-end pipeline with AF2 validation):
```
bindcraft(
    target_pdb=target_pdb_path,
    hotspot_residues=hotspot_residues,
    binder_length=(min_length, max_length),
    n_designs=20,
    af2_filter_iptm=0.6,
)
→ {designs, n_passed, best_pdb_path, best_iptm}
```
Use BindCraft if a single-call, pre-validated output is preferred over step-by-step control.

### Step 3 — Sequence design on backbones (ProteinMPNN)
For each backbone with `mean_plddt > 70`:
```
proteinmpnn(
    pdb_path=backbone_pdb,
    chains_to_design=["B"],        # peptide chain (B by convention)
    n_sequences=8,
    sampling_temperature=0.1,
)
→ {sequences: [<peptide_seq>, ...], scores: [...]}
```
Total candidates: ~240 (30 backbones × 8 sequences).

### Step 4 — Structure prediction and validation
For top-50 sequences by ProteinMPNN score, fold the peptide–target complex:
```
colabfold(
    sequences=f"{target_seq}:{peptide_seq}",   # multimer prediction
    n_models=2,                                 # fast, 2 models sufficient for triage
    use_templates=True,
)
→ {best_pdb_path, mean_plddt, ptm_scores}
```
Filter: `mean_plddt > 75` and `ptm_scores[0] > 0.6`.

### Step 5 — Binding affinity scoring (top-20)
```
boltz2_affinity(
    protein_fasta=target_fasta,
    ligand_smiles=None,
    antibody_fasta=peptide_fasta,  # peptide treated as "antibody" in Boltz-2
)
→ {affinity_kcal_mol, iptm_score}
```
Keep `iptm_score > 0.6` and `affinity_kcal_mol < −9 kcal/mol`.

### Step 6 — Physicochemical analysis (top-10 peptides)
```
rdkit_props(smiles=[peptide_to_smiles(seq) for seq in top10_seqs])
```
Note: peptide SMILES are complex — use MW and logP as rough guides:
- Linear peptide MW < 5 kDa typical (~40 residues max)
- LogP: use `rdkit_props` on the amino acid composition as a proxy

### Step 7 — Synthesisability assessment
For each top-10 peptide, check:
- Standard amino acids only → solid-phase peptide synthesis (SPPS) compatible
- Non-standard residues (e.g. Aib, D-amino acids for stapled) → flag as custom synthesis
- Cys count: 0 or 2 (for disulfide-cyclic) → flag multiple Cys as oxidation risk

### Step 8 — Summary visualisation (text-based)
Report as table:
```
Rank | Sequence | Length | ipTM | ColabFold pLDDT | Boltz-2 ΔG | Synthesis complexity
```

---

## Decision branches

```
Peptide type is "cyclic"?
  → use rfdiffusion_pep with cyclic backbone constraints (via contig string)
  → proteinmpnn: add fixed positions at N/C termini for cyclisation

Peptide type is "stapled"?
  → design linear backbone first; add i,i+4 or i,i+7 staple constraint
  → requires custom Aib residue parameterisation — note limitation

BindCraft vs rfdiffusion_pep + proteinmpnn?
  → BindCraft: use when you want a validated single-step result and can wait
    4–12 h for AF2 in the loop
  → rfdiffusion_pep + proteinmpnn: use when you want fine-grained control
    or are running on a tight time budget

Hotspot residues define a PPI interface?
  → set binder_length range wider (25–50 residues) to cover more interface area

ColabFold fails for multimer?
  → fall back to boltz2 in protein-protein complex mode
```

---

## Expected outputs

Ranked table of ≤ 10 peptide candidates:
```
Rank | Sequence | Length (aa) | ipTM | Boltz-2 ΔG (kcal/mol) | Synthesis type | Notes
```

---

## Success criteria

- ≥ 5 candidates with `iptm_score > 0.6`
- ≥ 3 candidates with `boltz2_affinity < −9 kcal/mol`
- All candidates < 50 residues
- Synthesis feasibility assessed for all top-10

---

## References

- RFdiffusion: Watson et al., *Nature* 2023 — `rfdiffusion_pep`
- BindCraft: Pacesa et al., *bioRxiv* 2024 — `bindcraft`
- ProteinMPNN: Dauparas et al., *Science* 2022 — `proteinmpnn`
- ColabFold: Mirdita et al., *Nat. Methods* 2022 — `colabfold`
- Boltz-2 — `boltz2_affinity`
