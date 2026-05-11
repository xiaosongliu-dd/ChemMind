# Skill: esmfold

**Tier:** L1 — atomic tool skill
**Category:** structure
**Tool:** ESMFold — fast sequence-only protein structure prediction
**GPU required:** Optional — CPU works; NIM API removes local GPU requirement

---

## What this skill does

ESMFold predicts protein 3D structure directly from amino acid sequence using a language-model-based folding head (ESM-2). Unlike AlphaFold2 / ColabFold, ESMFold **requires no multiple sequence alignment** — folding happens in seconds per protein. Accuracy is slightly below AF2 on well-conserved proteins but competitive for novel sequences and dramatically faster.

Key capabilities:
- Single-sequence structure prediction (no MSA, no templates)
- Per-residue pLDDT confidence (B-factor column)
- pTM global fold quality score (NIM mode)
- Batch mode: list of sequences, errors captured per-item

---

## When to call this skill

- Quick structure needed without MSA computation
- Screening many novel sequences (fast throughput)
- No local GPU — use NIM API mode
- First-pass fold before committing to ColabFold

Do NOT use if:
- High accuracy required on well-characterised protein → use `colabfold`
- Predicting a protein–ligand or protein–protein complex → use `boltz2`

---

## Inputs

| Field | Type | Required | Description |
|---|---|---|---|
| `sequence` | str \| list[str] | ✓ | Amino acid sequence or list for batch |
| `use_nim` | bool | — | Use NVIDIA NIM API (default True) |
| `nim_api_key` | str | — | Reads `NVIDIA_API_KEY` env if None |
| `output_dir` | str | — | PDB output directory; temp dir if None |

---

## Outputs

| Field | Type | Description |
|---|---|---|
| `pdb_path` | str | Path to predicted PDB |
| `mean_plddt` | float \| None | Mean pLDDT (0–100); parsed from B-factor |
| `ptm_score` | float \| None | Predicted TM-score (NIM mode only) |
| `sequence` | str | Echo of input |
| `runtime_s` | float | Wall time in seconds |

Batch mode (`sequence: list[str]`): `{results: [single-mode dict, ...], n_sequences}`.

---

## Implementation

`tools/structure/esmfold.py` — `run_esmfold()`

NIM mode POSTs to `https://health.api.nvidia.com/v1/biology/nvidia/esmfold`; local mode calls `esm.pretrained.esmfold_v1().infer_pdb()`. Mean pLDDT is extracted from the B-factor column of the output PDB.

---

## Agent decision rules

- **Default to NIM mode**: local `esm` install requires ~3 GB of weights + GPU
- **pLDDT gate**: `mean_plddt > 70` = confident; 50–70 = uncertain; < 50 = unreliable → retry with `colabfold`
- **After folding**: pipe into `pocket_detect` for druggability, or directly to `docking/*`
- **Batch**: ESMFold runs one sequence at a time; for >20 sequences prefer `colabfold` (shares MSA)
- **Upstream**: sequence from `data/chembl` target lookup or `data/pdb_fetch` SIFTS

---

## Install

```bash
# NIM mode (recommended)
export NVIDIA_API_KEY=<your-key>
# No local install needed

# Local mode
pip install fair-esm  # downloads ~3 GB weights on first use
```

---

## References

- Lin Z. et al. *Science* 2023, 379, 1123 — "Evolutionary-scale prediction of atomic-level protein structure with a language model"
- GitHub: https://github.com/facebookresearch/esm
- NVIDIA NIM: https://build.nvidia.com/nvidia/esmfold
- License: MIT
