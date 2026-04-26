# Skill: rfdiffusion_pep

**Tier:** L1 — atomic tool skill
**Category:** biologics
**Tool:** RFdiffusion (peptide / mini-protein binder backbone design)
**GPU required:** Yes — A100/H100 recommended; NVIDIA NIM API available

---

## What this skill does

RFdiffusion generates de novo protein backbone conformations geometrically complementary to a target epitope. In peptide/binder mode, it produces backbone-only PDB files — no sequences yet. Sequences are designed in the next step by ProteinMPNN.

Key capabilities:
- **Peptide binder design**: short peptides (10–30 aa) against a target hotspot
- **Mini-protein design**: small proteins (50–150 aa) for tighter, more stable binders
- **Motif scaffolding**: fix a functional motif (e.g. a binding helix), design scaffold around it
- **Partial diffusion**: refine an existing binder backbone without full redesign

---

## When to call this skill

- Starting a peptide or protein binder campaign from a target PDB
- Refining an existing binder backbone via partial diffusion
- Generating diverse binder geometries to explore binding mode space

---

## Inputs

| Field | Type | Required | Description |
|---|---|---|---|
| `target_pdb` | str | ✓ | Receptor PDB (held fixed during diffusion) |
| `hotspot_residues` | list[str] | ✓ | Epitope residues e.g. `["A25", "A30", "A35"]` (chain+resnum) |
| `binder_length` | tuple[int,int] | — | (min, max) residue count for designed binder (default `(10, 20)`) |
| `n_designs` | int | — | Number of backbone designs to generate (default 10) |
| `noise_scale` | float | — | 1.0 = full de novo; <1.0 = partial diffusion from existing structure |
| `partial_diffusion_pdb` | str | — | Starting backbone for partial diffusion |
| `partial_T` | int | — | Partial diffusion steps (used with `partial_diffusion_pdb`) |
| `use_nim` | bool | — | Use NVIDIA NIM API (default True) |
| `nim_api_key` | str | — | Reads `NVIDIA_API_KEY` env if None |

---

## Outputs

| Field | Type | Description |
|---|---|---|
| `designs` | list[dict] | `{pdb_path, mean_plddt, rank}` sorted by pLDDT descending |
| `best_pdb_path` | str | Top-ranked design PDB |
| `n_designs` | int | Number of designs generated |
| `hotspot_res` | list[str] | Hotspot residues used |

---

## Implementation

`tools/biologics/rfdiffusion_pep.py` — `run_rfdiffusion_pep()`

---

## Agent decision rules

- **Always follow with proteinmpnn**: RFdiffusion outputs backbone-only — no sequence
- **pLDDT gate**: discard designs with `mean_plddt < 0.7` before ProteinMPNN
- **n_designs**: 10 for quick triage; 50–100 for a real campaign (expect ~5–15% to pass AF2 filter)
- **Partial diffusion**: use `noise_scale=0.3` to refine a hit; keep geometry, change sequence compatibility
- **Upstream**: `pocket/p2rank` or `structure/boltz2` (hotspot definition) → `rfdiffusion_pep`
- **Downstream**: → `proteinmpnn` (sequence) → `igfold` (fold check) → `bindcraft` (AF2 ipTM filter)

---

## Install

```bash
export NVIDIA_API_KEY="nvapi-..."   # NIM mode — no local install

# Local
git clone https://github.com/RosettaCommons/RFdiffusion
cd RFdiffusion && pip install -e .
bash scripts/download_models.sh Base  # ~1.5 GB
```

## References

- Watson et al., Nature 2023 — "De novo design of protein structure and function with RFdiffusion"
- NVIDIA NIM: https://docs.nvidia.com/nim/bionemo/rfdiffusion/latest/overview.html
- GitHub: https://github.com/RosettaCommons/RFdiffusion
- License: BSD-3
