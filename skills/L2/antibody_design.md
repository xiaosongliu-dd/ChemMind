# Workflow: antibody_design

**Tier:** L2 — workflow skill
**Purpose:** Design de novo antibodies (scFv or VHH nanobodies) against a user-specified antigen epitope using diffusion-based backbone generation followed by sequence design, structure prediction, complex scoring, and developability filtering.

---

## When to invoke this workflow

- Need a new antibody against a target antigen where no starting antibody exists
- Goal: ≥ 5 sequences passing developability filters, ready for mammalian expression
- Keyword triggers: "design antibody", "nanobody", "scFv", "VHH", "biologic", "antibody against"

---

## Required inputs

| Input | Description |
|---|---|
| Antigen PDB | Antigen structure; resolved by `pdb_fetch` if only PDB ID given |
| Epitope residues | List of `"<chain><resnum>"` strings defining the target epitope |
| Antibody type | `"scfv"` (default) or `"nanobody"` |
| n_designs | Optional — number of backbone designs (default 50) |

---

## L1 skill sequence

### Step 1 — Antigen structure
```
pdb_fetch(pdb_id=<antigen_pdb_id>, query_type="structure")
→ {pdb_path, chain_ids, ligands}
```
If epitope residues are not provided, run:
```
pocket_detect(protein_pdb=pdb_path, tool="p2rank")
→ best_pocket.residues  # use as hotspot_residues
```

### Step 2 — Antibody backbone design (RFantibody)
```
rfantibody(
    target_pdb=antigen_pdb_path,
    hotspot_residues=epitope_residues,
    cdr_lengths={"H3": (8, 15)},
    n_designs=50,
    use_nim=True,
)
→ {designs: [{pdb_path, mean_plddt, rank}], best_pdb_path}
```
Keep designs with `mean_plddt > 70`. Typical: 30–45 out of 50 pass.

### Step 3 — Sequence design on backbones (ProteinMPNN)
For each passing backbone design:
```
proteinmpnn(
    pdb_path=backbone_pdb,
    chains_to_design=["H", "L"],   # or ["H"] for nanobody
    n_sequences=8,
    sampling_temperature=0.1,
)
→ {sequences: [{"H": <vh_seq>, "L": <vl_seq>}, ...], scores: [...]}
```
Total candidates: ~240 (30 backbones × 8 sequences).

### Step 4 — Antibody structure prediction
For each unique sequence pair, predict structure for validation:
```
igfold(
    heavy_sequence=vh_seq,
    light_sequence=vl_seq,         # None for nanobody
    do_refine=True,
)
→ {pdb_path, plddt_score}
```
Filter: `plddt_score > 0.7`. Run `abodybuilder3` as cross-check on top-20 by pLDDT:
```
abodybuilder3(
    heavy_sequence=vh_seq,
    light_sequence=vl_seq,
    numbering_scheme="imgt",
)
→ {pdb_path, predicted_error}
```
Discard where `predicted_error > 1.5 Å`.

### Step 5 — Complex scoring (Boltz-2 Ab–Ag mode)
For the top-20 validated antibodies:
```
boltz2(
    protein_fasta=antigen_fasta,
    antibody_fasta=vh_vl_fasta,    # heavy:light or VHH
    n_samples=3,
    use_nim=True,
)
→ {iptm_score, pdb_path}
```
Keep candidates with `iptm_score > 0.6`. Typical retention: 5–10 candidates.

### Step 6 — Developability / ADMET filtering
```
admetlab3(smiles=None, sequences=top_ab_sequences)
```
Note: ADMETlab3 for biologics is limited — check manually or with BioPharma Finder:
- Hydrophobicity index (HI) < 0.5
- Net charge at pH 7 between −3 and +3
- No unpaired cysteines in VH/VL
- Predicted aggregation propensity (AP) flag

Report `admetlab3` ADMET predictions if the antibody is conjugated to a small-molecule warhead.

### Step 7 — CDR-H3 optimisation (optional, for top-3)
If top candidates need CDR-H3 diversity:
```
abdiffuser(
    antibody_framework_pdb=top_pdb,
    antigen_pdb=antigen_pdb_path,
    cdrs_to_design=["H3"],
    n_designs=20,
)
→ {designs: [{pdb_path, cdr_sequences, log_likelihood}]}
```

### Step 8 — Visualisation summary
Report as text table (no small-molecule ligand_viz):
```
Rank | VH sequence | VL sequence | ipTM | IgFold pLDDT | ABodyBuilder3 error | CDR-H3 length
```

---

## Decision branches

```
Antibody type is "nanobody"?
  → omit VL in rfantibody (use VHH mode)
  → igfold: light_sequence=None → triggers nanobody mode
  → abodybuilder3: light_sequence=None → NanoBodyBuilder2

Epitope residues unknown?
  → pocket_detect to get surface residues; alternatively use ConSurf
    conservation scores to identify functional surface patches

< 5 designs pass ipTM > 0.6?
  → increase n_designs to 100 in rfantibody
  → run abdiffuser on top-5 to diversify CDR-H3

ProteinMPNN unavailable (NIM error)?
  → fall back to a fixed-sequence design (VH germline framework + random CDR-H3)
    and note the limitation
```

---

## Expected outputs

Ranked table of ≤ 5 antibody candidates:
```
Rank | VH sequence | VL sequence | CDR-H3 | ipTM | IgFold pLDDT | AB3 error | Notes
```

---

## Success criteria

- ≥ 5 candidates with `iptm_score > 0.6`
- All candidates have `igfold.plddt_score > 0.7` and `abodybuilder3.predicted_error < 1.5 Å`
- CDR-H3 diversity: ≥ 3 distinct CDR-H3 sequences among top candidates
- No obvious developability red flags (unpaired Cys, extreme charge)

---

## References

- RFantibody: Bennett et al., 2025 — `rfantibody`
- ProteinMPNN: Dauparas et al., *Science* 2022 — `proteinmpnn`
- IgFold: Ruffolo et al., *Nat. Commun.* 2023 — `igfold`
- ABodyBuilder3: Abanades et al., *Nat. Commun. Biol.* 2023 — `abodybuilder3`
- AbDiffuser: Martinkus et al., *NeurIPS* 2023 — `abdiffuser`
- Boltz-2 — `boltz2`
