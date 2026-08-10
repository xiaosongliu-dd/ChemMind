# Skill: chai1

**Tier:** L1 — atomic tool skill
**Category:** structure
**Tool:** Chai-1 (Chai Discovery)
**GPU required:** Yes — A100/H100 recommended; A10/A30/RTX 4090 work for smaller complexes
**License:** Weights auto-download, no separate request required

---

## What this skill does

Chai-1 is a multi-modal foundation model from Chai Discovery that predicts structures of **any biomolecular complex** in a single pass. It provides a clean Python API (`run_inference`) with auto-downloaded weights and strong performance on protein–ligand complexes.

Unique capabilities vs other tools in ChemMind:
- **Python API** — no subprocess, cleaner integration and error handling
- **Automatic weight download** — `pip install chai_lab` then run directly
- **Glycosylation support** — modified residues encoded inline (e.g. `AAA(SEP)AAA` for phosphoserine)
- **Restraint-guided folding** — can provide inter-chain contact restraints and covalent bonds to guide the prediction (useful when partial experimental data is available)
- Handles proteins, small molecules (SMILES), DNA, RNA in one call

---

## When to call this skill

- Clean Python-only workflow (no subprocess management)
- Need glycoprotein structure (modified residues)
- Want to incorporate experimental distance restraints into the structure prediction
- Protein–ligand co-complex, protein–DNA/RNA — standard multi-chain task
- Rapid iteration: weights download once then cached, no setup beyond `pip install`

---

## When NOT to use this skill

| Situation | Better choice |
|-----------|--------------|
| Commercial production deployment | Protenix (Apache 2.0, no NC restriction) |
| Best antibody-antigen accuracy | Protenix-v2 (outperforms Chai-1 on Ab-Ag) |
| Need binding affinity ΔG | Boltz-2 (returns affinity in same call) |
| Fastest protein-only folding | ESMFold (~2 s, no GPU queue) |
| MSA-quality protein structure | ColabFold (pass --use-msa-server to chai1 or use ColabFold directly) |

---

## Inputs

```python
run_chai1(
    protein_sequences    = "MKTAYIAKQR...",   # str or list[str]
    ligand_smiles        = "CCO",             # optional str or list[str]
    dna_sequences        = ["ATCGATCG"],      # optional
    rna_sequences        = ["AUCGAUCG"],      # optional
    seeds                = [42, 137],         # list[int]; one prediction per seed
    num_trunk_recycles   = 3,                 # Evoformer recycles (3 = default)
    num_diffn_timesteps  = 200,               # diffusion steps (50 = fast, 200 = accurate)
    use_esm_embeddings   = True,              # ESM-2 protein LM augmentation
    device               = "cuda:0",
    output_dir           = None,
)
```

**Ligands as SMILES:** Pass any valid SMILES string; Chai-1 automatically featurizes it.
**Glycosylated residues:** Use the modified FASTA notation inline, e.g. `"MKTAA(HY3)K..."` where `HY3` is the CCD code of the modification.
**Multiple seeds = ensemble:** Each seed produces one structural prediction. Use 5 seeds for confident decisions, 1 for rapid hypothesis testing.

---

## Outputs

```python
{
    "top_cif_path":     "/tmp/chai1_xxx/seed_42/pred.model_idx_0.cif",
    "all_cif_paths":    ["/tmp/.../pred.model_idx_0.cif", ...],  # best-first
    "aggregate_scores": [0.81, 0.77, 0.74, ...],   # per-prediction (best first)
    "mean_plddt":       83.2,     # mean pLDDT of best prediction (0–100)
    "has_ligand":       True,
    "has_nucleic_acid": False,
    "n_chains":         2,
    "output_dir":       "/tmp/chai1_xxx",
    "seeds_used":       [42, 137],
}
```

Additional per-prediction data is in `scores.model_idx_{i}.npz` files in each seed output directory:
```python
import numpy as np
scores = np.load("seed_42/scores.model_idx_0.npz")
# keys: 'plddt', 'pae', 'pde', 'ptm', 'iptm'
ptm  = scores["ptm"].item()
iptm = scores["iptm"].item()
```

---

## Agent decision rules

### Quality gates
| Metric | Threshold | Action |
|--------|-----------|--------|
| `aggregate_score` | > 0.7 | Reliable prediction |
| `aggregate_score` | 0.5–0.7 | Use with care; inspect per-residue pLDDT |
| `aggregate_score` | < 0.5 | Uncertain; try more seeds or check complex composition |
| `mean_plddt` | ≥ 80 | High confidence overall |
| `mean_plddt` | 50–80 | Partial confidence; loops may be disordered |

### Speed vs accuracy tradeoff
| `num_diffn_timesteps` | Speed | Quality |
|-----------------------|-------|---------|
| 50 | ~3× faster | Screening mode |
| 200 (default) | standard | Publication quality |

### When to use multiple seeds
- Any decision that will advance to synthesis: use ≥ 3 seeds
- Ensemble mode: 5–10 seeds, take the consensus binding pocket across predictions
- Rapid hypothesis: 1 seed only

---

## Installation

```bash
pip install chai_lab   # latest stable
# OR latest dev:
pip install git+https://github.com/chaidiscovery/chai-lab.git

# Weights auto-downloaded on first run to <package_root>/downloads/
# Override with: export CHAI_DOWNLOADS_DIR=/your/path
```

Requires Linux, Python ≥ 3.10, CUDA GPU with bfloat16 support.

---

## References

- Chai-1 Technical Report (2024). https://www.biorxiv.org/content/10.1101/2024.10.10.615955
- GitHub: https://github.com/chaidiscovery/chai-lab
- Web server: https://lab.chaidiscovery.com
