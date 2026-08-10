# Skill: protenix

**Tier:** L1 — atomic tool skill
**Category:** structure
**Tool:** Protenix v2 (ByteDance AI4Science)
**GPU required:** Yes — A100 40 GB minimum; ~17 s for 500 tokens, ~60 s for 1000 tokens
**License:** Apache 2.0 (code + weights) — fully open, no restrictions

---

## What this skill does

Protenix is ByteDance's open-source reimplementation of AlphaFold 3, with independent training and architectural improvements. Protenix-v2 (April 2026) **outperforms AF3** on antibody-antigen benchmarks and matches it on protein-ligand and protein-only tasks — while being fully open under Apache 2.0, meaning there is no non-commercial weight restriction.

Capabilities (same as AF3):
- Protein structure (single chain, multimer)
- Protein–small molecule co-complex (SMILES or CCD code)
- Protein–DNA / protein–RNA complexes
- Post-translational modifications (phosphorylation, glycosylation via CCD codes)
- Covalent bond specification between chains

**Why use Protenix instead of AF3:** No weight request required, Apache 2.0 license allows commercial use, v2 model is more accurate on Ab-Ag complexes, pip-installable.

---

## When to call this skill

- All use cases where AF3 would be called, but the user has not requested AF3 weights from Google DeepMind
- Commercial or production drug discovery workflows (Apache 2.0 removes licensing friction)
- Antibody-antigen complex prediction (Protenix-v2 substantially better than AF3 here)
- When you need the best openly licensed structure prediction for any biomolecule type

---

## When NOT to use this skill

| Situation | Better choice |
|-----------|--------------|
| Need affinity ΔG simultaneously | Boltz-2 (structure + affinity in one call) |
| Ultra-fast screening of many proteins | ESMFold (no MSA, ~2 s/protein) |
| Need MSA-quality protein-only, no ligand | ColabFold (3× faster than Protenix for protein-only) |
| Very large complexes (>4000 tokens) | May OOM — reduce n_sample to 1 first |

---

## Inputs

```python
run_protenix(
    protein_sequences = "MKTAYIAKQR...",        # str or list[str]
    ligand_smiles     = "CCO",                  # optional str or list[str]
    ligand_ccd_codes  = ["ATP", "MG"],          # optional list[str], prefixed as CCD_
    dna_sequences     = ["ATCGATCG"],           # optional 5'→3' single-stranded
    rna_sequences     = ["AUCGAUCG"],           # optional
    seed              = 101,                    # single random seed
    n_sample          = 5,                      # structural samples per seed
    model             = "protenix-v2",          # or "protenix_base_default_v1.0.0"
    output_dir        = None,                   # auto temp dir if None
    timeout           = 3600,
)
```

**Ligand input precedence:** SMILES > CCD code > SDF file. Use SMILES for novel ligands, CCD for standard cofactors (ATP, NAD, HEM, etc.).

---

## Outputs

```python
{
    "top_cif_path":     "/tmp/protenix_xxx/chemmind_ptx/seed-101_sample-0/model.cif",
    "all_cif_paths":    ["/tmp/.../model.cif", ...],  # best-first by ranking_score
    "ranking_score":    0.79,     # AF3-style ranking score; >0.7 = high confidence
    "ptm":              0.84,     # predicted TM-score (0–1)
    "iptm":             0.76,     # interface pTM (complexes); >0.7 = reliable interface
    "mean_plddt":       82.1,     # per-residue mean pLDDT (0–100)
    "has_ligand":       True,
    "has_nucleic_acid": False,
    "n_chains":         2,
    "output_dir":       "/tmp/protenix_xxx",
    "seed":             101,
    "n_sample":         5,
    "model":            "protenix-v2",
}
```

Structures are in **CIF format**. Convert to PDB with `gemmi convert model.cif model.pdb`.

---

## Agent decision rules

### Quality gates
| Metric | Threshold | Action |
|--------|-----------|--------|
| `ranking_score` | ≥ 0.7 | Accept for downstream use |
| `ranking_score` | 0.5–0.7 | Use with caution; inspect pLDDT per region |
| `ranking_score` | < 0.5 | Low confidence — try different complex composition or more samples |
| `ptm` | ≥ 0.8 | Reliable global fold |
| `iptm` (complex) | ≥ 0.7 | Reliable interface |
| `mean_plddt` | ≥ 80 | High confidence |
| `mean_plddt` | 50–80 | Partial — loop regions uncertain |

### Model selection
| Need | Model |
|------|-------|
| Best accuracy (default) | `protenix-v2` |
| Reproduce published benchmarks | `protenix_base_default_v1.0.0` |
| Recent PDB structures (post-2021) | `protenix_base_20250630_v1.0.0` |
| Low memory / fast check | `protenix-mini` (see docs) |

### Seed and sample strategy
- **n_sample=5, seed=101** (default): 5 independent structures, suitable for most tasks
- **n_sample=1**: fastest, for rapid structural hypothesis check
- **Multiple seeds** (run protenix multiple times with different seeds): for ensemble analyses or when ranking_score < 0.6

---

## Installation

```bash
pip install --upgrade protenix --index-url https://pypi.org/simple

# Verify
protenix pred --help

# Optional: pre-build CUDA kernels for 2–5× speedup
protenix build-kernels
```

**No weight request required.** Weights are downloaded automatically on first run.
**License:** Apache 2.0 — code and weights. No non-commercial restriction.

---

## References

- Zhang et al. "Protenix-v2: Broadening the Reach of Structure Prediction and Biomolecular Design." *bioRxiv* (2026). https://doi.org/10.64898/2026.04.10.717613
- Zhang et al. "Protenix-v1: Toward High-Accuracy Open-Source Biomolecular Structure Prediction." *bioRxiv* (2026). https://doi.org/10.64898/2026.02.05.703733
- GitHub: https://github.com/bytedance/protenix
- Web server (no GPU needed): https://protenix-server.com
