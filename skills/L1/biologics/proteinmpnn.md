# Skill: proteinmpnn

**Tier:** L1 — atomic tool skill
**Category:** biologics
**Tool:** ProteinMPNN (inverse folding / sequence design)
**GPU required:** No — CPU inference; NVIDIA NIM API available

---

## What this skill does

ProteinMPNN designs amino acid sequences for a given protein backbone using a message-passing neural network trained on PDB structures. Given any 3D backbone PDB, it outputs sequences predicted to fold into that backbone, ranked by log-probability score.

It is used in virtually every backbone-first protein design pipeline: RFdiffusion generates the backbone → ProteinMPNN assigns sequences → AlphaFold2/IgFold validates the fold.

Key capabilities:
- **Fixed-backbone sequence design** across one or multiple chains
- **Multi-chain**: design VH + VL simultaneously with cross-chain communication
- **Fixed residues**: lock catalytic/interface residues; design the rest
- **Soluble model**: variant biased toward soluble, expressible sequences

---

## When to call this skill

- After `rfdiffusion_pep` or `rfantibody` produce a backbone — always needed to get sequences
- Designing sequence variants of a known structure for expression or stability optimization
- Redesigning an interface while keeping the scaffold fixed

---

## Inputs

| Field | Type | Required | Description |
|---|---|---|---|
| `pdb_path` | str | ✓ | Backbone PDB to design sequences for |
| `chains_to_design` | list[str] | ✓ | Chain IDs to design, e.g. `["A"]` or `["H", "L"]` |
| `fixed_positions` | dict | — | `{chain: [residue_indices]}` — residues locked to original sequence |
| `n_sequences` | int | — | Sequences per design (default 8) |
| `sampling_temp` | float | — | Diversity: 0.1 = focused / confident, 1.0 = diverse (default 0.1) |
| `use_soluble_model` | bool | — | Soluble-biased variant (default False) |
| `use_nim` | bool | — | Use NVIDIA NIM API (default True) |
| `nim_api_key` | str | — | Reads `NVIDIA_API_KEY` env if None |

---

## Outputs

| Field | Type | Description |
|---|---|---|
| `sequences` | list[str] | Designed amino acid sequences |
| `scores` | list[float] | Global log-likelihood per sequence (higher = better predicted foldability) |
| `n_returned` | int | Sequences returned |
| `chains_designed` | list[str] | Which chains were designed |

---

## Implementation

`tools/biologics/proteinmpnn.py` — `run_proteinmpnn()`

---

## Agent decision rules

- **Always follow any backbone generation step** (`rfdiffusion_pep`, `rfantibody`) — backbone without sequence is unusable
- **Temperature**: 0.1 for CDR loops (conservative); 0.3–0.5 for scaffold regions (diversity)
- **Fixed positions**: always lock catalytic residues, disulfides, and known epitope contacts
- **n_sequences**: 8–50; then filter by AF2 ipTM or IgFold pLDDT
- **Upstream**: `rfdiffusion_pep` or `rfantibody` backbone → `proteinmpnn`
- **Downstream**: sequences → `igfold` or `abodybuilder3` (fold check) → `admet/admetlab3` (developability)

---

## Install

```bash
export NVIDIA_API_KEY="nvapi-..."   # NIM mode — no local install

# Local
git clone https://github.com/dauparas/ProteinMPNN
cd ProteinMPNN && pip install -e .
```

## References

- Dauparas et al., Science 2022 — "Robust deep learning-based protein sequence design using ProteinMPNN"
- NVIDIA NIM: https://docs.nvidia.com/nim/bionemo/proteinmpnn/latest/overview.html
- GitHub: https://github.com/dauparas/ProteinMPNN
- License: MIT
