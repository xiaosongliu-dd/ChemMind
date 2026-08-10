# Skill: alphafold3

**Tier:** L1 — atomic tool skill
**Category:** structure
**Tool:** AlphaFold 3 (Google DeepMind, local weights)
**GPU required:** Yes — A100 (40 GB) minimum; ~15–60 min per complex depending on size

---

## What this skill does

AlphaFold 3 is Google DeepMind's unified structure prediction model that predicts the 3D structures of **any biomolecular complex** in a single diffusion-based forward pass. It accepts proteins, small molecules (as SMILES or CCD codes), DNA, and RNA simultaneously and outputs a co-complex structure with per-residue confidence scores.

Key capabilities beyond AF2/ColabFold:
- **Protein–small molecule co-complexes** — first major model to handle drug-like ligands end-to-end with MSA
- **Protein–DNA / protein–RNA** structures (e.g., transcription factors, CRISPR complexes)
- **Covalent modifications** — phosphorylation, glycosylation, disulfide bonds
- **Arbitrary assemblies** — heterodimers, oligomers, ligand + cofactor + protein + nucleic acid in one shot
- **Ranking score** — built-in quality metric for selecting the best model seed

---

## When to call this skill

- Need a protein–ligand co-complex structure with the best available accuracy
- Target has DNA/RNA binding partners (not supported by Boltz-2, ESMFold, or ColabFold)
- Modelling post-translational modifications on the target
- Generating diverse structural hypotheses for an ensemble docking campaign (run 5 seeds → 5 poses)
- Benchmarking: want the most accurate publicly available structure model

---

## When NOT to use this skill

| Situation | Better choice |
|-----------|--------------|
| Need fast screening of 100s of sequences | ESMFold (no GPU queue, ~2 s/protein) |
| Protein-only, need MSA-quality accuracy | ColabFold (same AF2 engine, 3× faster than AF3 for protein-only) |
| Need affinity ΔG simultaneously with structure | Boltz-2 (returns affinity in same call) |
| >50 residues of disordered loop regions | AF3 handles IDRs but confidence will be low — check pLDDT |

---

## Inputs

```python
run_alphafold3(
    protein_sequences = "MKTAYIAKQRQ...",        # str or list[str] — plain sequence, no FASTA header
    ligand_smiles     = "CC1=CC=CC=C1",          # optional SMILES; str or list[str]
    ligand_ccd_codes  = ["ATP", "MG"],           # optional CCD codes (alternative to SMILES)
    dna_sequences     = ["ATCGATCG"],            # optional 5'→3' DNA strand(s)
    rna_sequences     = ["AUCGAUCG"],            # optional RNA strand(s)
    seeds             = [1, 2, 3, 4, 5],        # explicit seeds; or use num_seeds=5
    num_seeds         = 5,                       # number of independent predictions
    model_dir         = "/path/to/af3/weights",  # or set $AF3_MODEL_DIR
    output_dir        = None,                    # auto tmpdir if not given
    af3_script        = "run_alphafold.py",      # or $AF3_SCRIPT
    timeout           = 3600,
)
```

**Chain ID assignment (automatic):**
Proteins get A, B, C…; ligands follow; DNA/RNA last. Order within each type follows input list order.

---

## Outputs

```python
{
    "top_cif_path":     "/tmp/alphafold3_xxx/chemmind_af3/chemmind_af3_model_0.cif",
    "all_cif_paths":    ["/tmp/.../model_0.cif", "/tmp/.../model_1.cif", ...],  # best first
    "ranking_score":    0.78,   # AF3 ranking score — higher is better; >0.7 is high confidence
    "ptm":              0.82,   # predicted TM-score (0–1); >0.8 = reliable global fold
    "iptm":             0.74,   # interface pTM — only for complexes; >0.7 = reliable interface
    "mean_plddt":       81.4,   # mean per-residue pLDDT (0–100); >80 = confident
    "has_ligand":       True,
    "has_nucleic_acid": False,
    "n_chains":         2,
    "output_dir":       "/tmp/alphafold3_xxx",
    "seeds_used":       [1, 2, 3, 4, 5],
}
```

Structures are returned in **CIF format** (mmCIF). Convert to PDB if needed using gemmi:
```bash
gemmi convert model.cif model.pdb
```

---

## Agent decision rules

### Quality gates
| Metric | Threshold | Action |
|--------|-----------|--------|
| `ranking_score` | ≥ 0.7 | Accept model for downstream use |
| `ranking_score` | 0.5–0.7 | Use with caution; check pLDDT per region |
| `ranking_score` | < 0.5 | Low confidence — try more seeds or different complex composition |
| `ptm` | ≥ 0.8 | Reliable global fold |
| `iptm` (complex) | ≥ 0.7 | Interface is likely correct |
| `mean_plddt` | ≥ 80 | High confidence overall |
| `mean_plddt` | 50–80 | Partial confidence — loop regions may be wrong |

### Choosing num_seeds
- **5 seeds** (default): standard; 5× slower than 1 but gives diversity and confidence spread
- **1 seed**: quick check or when compute is limited
- **10+ seeds**: ensemble-based downstream analyses (MD, FEP), statistical significance checks

### What to do after AF3
1. Visualise in PyMOL/ChimeraX — inspect ligand pose, pocket shape
2. Run `pocket_detect` to confirm/refine the binding site
3. If `iptm` < 0.6 for a complex: the interface geometry is unreliable — fall back to docking (GNINA) on the AF3 apo structure
4. Pass `top_cif_path` (converted to PDB) to docking tools or FEP

---

## Installation

```bash
# 1. Clone the AF3 repository
git clone https://github.com/google-deepmind/alphafold3
cd alphafold3 && pip install -e .

# 2. Request model weights (non-commercial only)
#    https://forms.gle/svvpY4u2jsHEwWYS6
#    Place downloaded weights at: /data/af3_model_weights/

# 3. Set environment variables
export AF3_MODEL_DIR=/data/af3_model_weights
export AF3_SCRIPT=/opt/alphafold3/run_alphafold.py

# 4. (Optional) gemmi for CIF→PDB conversion
pip install gemmi
```

**License:** AlphaFold 3 model weights are non-commercial only. Each user must independently request weights from Google DeepMind. The code itself is Apache 2.0. Do not use for clinical diagnostic purposes.

---

## References

- Abramson et al. "Accurate structure prediction of biomolecular interactions with AlphaFold 3." *Nature* 630, 493–500 (2024). https://doi.org/10.1038/s41586-024-07487-w
- AF3 GitHub: https://github.com/google-deepmind/alphafold3
- AlphaFold Server: https://alphafoldserver.com
