# Skill: rfantibody

**Tier:** L1 — atomic tool skill
**Category:** biologics
**Tool:** RFantibody — RFdiffusion fine-tuned for antibody design
**License:** MIT (Baker Lab, University of Washington)
**GPU required:** Yes — A100 recommended; ~2–10 min per design batch

---

## What this skill does

RFantibody is a fine-tuned variant of RFdiffusion trained specifically on antibody loops. It performs **de novo design of human-like antibodies** (scFvs and VHH nanobodies) against a user-specified antigen epitope — no starting antibody structure needed.

Key capabilities:
- De novo scFv design (heavy + light chain variable fragments)
- VHH / nanobody design
- CDR-H3 loop generation (the primary binding determinant)
- Epitope-targeted design (specify which residues on the antigen to engage)
- Outputs backbone structures → feed to ProteinMPNN for sequence design

---

## When to call this skill

- Need a brand-new antibody against a target antigen
- No existing antibody scaffold available — truly de novo
- Want to target a specific epitope or interface region
- Generating a diverse library of CDR designs for experimental screening

Do NOT use if:
- You have an existing antibody and only need CDR optimisation → use `abdiffuser` or `proteinmpnn`
- You need full-atom structure → feed RFantibody backbone into `igfold` or `abodybuilder3`

---

## Inputs

| Field | Type | Required | Description |
|---|---|---|---|
| `target_pdb` | str | ✓ | Antigen / target PDB (held fixed during diffusion) |
| `hotspot_residues` | list[str] | ✓ | Epitope residues on target, format `"<chain><resnum>"`, e.g. `["B25", "B30", "B35"]` |
| `cdr_lengths` | dict[str, tuple[int,int]] | — | Per-CDR length range, e.g. `{"H3": (8, 15)}`. CDRs not listed use framework lengths (no design). Default: `{"H3": (8, 15)}` |
| `framework_pdb` | str | — | Optional antibody framework PDB to scaffold onto; if None, generates a full VH backbone de novo |
| `n_designs` | int | — | Number of antibody backbone designs (default 10) |
| `use_nim` | bool | — | Use NVIDIA NIM API endpoint (default True) |
| `nim_api_key` | str | — | Reads `NVIDIA_API_KEY` env if None |

---

## Outputs

| Field | Type | Description |
|---|---|---|
| `designs` | list[dict] | `{pdb_path, mean_plddt, rank}` sorted by pLDDT descending |
| `best_pdb_path` | str \| None | Top-ranked design's PDB path |
| `n_designs` | int | Number of designs returned |
| `hotspot_res` | list[str] | Echo of input |
| `cdr_lengths` | dict | Echo of resolved CDR lengths |

---

## Implementation

`tools/biologics/rfantibody.py` — `run_rfantibody()`

NIM mode hits the RFdiffusion endpoint (`https://health.api.nvidia.com/v1/biology/ipd/rfdiffusion`) with `model_runner="rf_antibody"`; local mode runs `RFdiffusion/scripts/run_inference.py` with antibody-specific contigs and the `models/rf_antibody.pt` checkpoint.

---

## Standard antibody design pipeline (recommended L2 sequence)

```
antigen PDB
    │
    ▼
rfantibody          ← this skill — generates backbone geometries
    │  50 backbones
    ▼
proteinmpnn         ← L1/biologics/proteinmpnn — designs sequences on backbones
    │  8 seqs / backbone = 400 candidates
    ▼
igfold / abodybuilder3   ← L1/biologics/igfold — fold + validate CDR-H3
    │
    ▼
boltz2 (Ab–Ag mode) ← L1/structure/boltz2 — score complex + affinity
    │
    ▼
admetlab3           ← L1/admet/admetlab3 — developability filters
```

This pipeline is encoded as the `antibody_design` L2 skill.

---

## Agent decision rules

- **Design count**: start with `n_designs=50`; ~5–10% will pass downstream filters
- **Epitope targeting**: always specify `epitope_residues` — blind design produces poor hit rates
- **scFv vs nanobody**: use `nanobody` for membrane/cryptic epitopes; `scfv` for most targets
- **Filtering**: after ProteinMPNN, filter by `boltz2.iptm_score > 0.6` before wet-lab
- **Diversity**: if designs cluster, increase `n_designs` or run with different random seeds

---

## Install

```bash
# Clone and install RFdiffusion (includes RFantibody fine-tuned weights)
git clone https://github.com/RosettaCommons/RFdiffusion.git
cd RFdiffusion
pip install -e .

# Download antibody-specific model weights (~2 GB)
bash scripts/download_models.sh RFantibody
```

---

## References

- Bennett et al. 2025 — "De novo design of human-like antibodies with RFdiffusion" (Baker Lab)
- GitHub: https://github.com/RosettaCommons/RFdiffusion
- License: MIT
- Chai-2 (June 2025) achieves ~50% hit rate using similar diffusion-based ab design approach
