# Skill: abdiffuser

**Tier:** L1 — atomic tool skill
**Category:** biologics
**Tool:** AbDiffuser — joint sequence + structure CDR design via score-based diffusion
**GPU required:** Yes — A100/H100 recommended; ~5–30 min per design batch

---

## What this skill does

AbDiffuser performs **co-design** of antibody CDR loops: instead of fixing a backbone and designing sequence (ProteinMPNN-style) or fixing a sequence and folding it (IgFold-style), AbDiffuser samples sequence and 3D conformation *jointly* via score-based diffusion. This couples sequence to structure during sampling, often producing CDRs better matched to their predicted geometry.

Use this skill when you have a working antibody framework + antigen, and want to redesign one or more CDR loops (typically CDR-H3, sometimes H1/H2) for an existing target.

---

## When to call this skill

- Optimising CDRs of an existing antibody scaffold (you already have the framework)
- Redesigning CDR-H3 against a defined epitope where ProteinMPNN's sequence-on-fixed-backbone is too rigid
- When you want the design-time joint sequence/structure coupling that ProteinMPNN doesn't provide
- After `igfold` / `abodybuilder3` — fold the framework first, then redesign CDRs

Do NOT use if:
- No starting framework exists → use `rfantibody` first to generate one
- Designing the entire antibody from scratch → `rfantibody` is the right starting point
- Designing non-antibody binders → `rfdiffusion_pep` or `bindcraft`

---

## Inputs

| Field | Type | Required | Description |
|---|---|---|---|
| `antibody_framework_pdb` | str | ✓ | Antibody PDB with framework regions intact (CDRs to redesign should be marked / will be regenerated) |
| `antigen_pdb` | str | ✓ | Antigen / target PDB (held fixed during design) |
| `cdrs_to_design` | list[str] | ✓ | CDR loops to redesign, e.g. `["H1", "H2", "H3"]` or just `["H3"]` |
| `n_designs` | int | — | Number of designs to sample (default 20) |
| `n_diffusion_steps` | int | — | Reverse diffusion steps; 100 = default quality, 50 = faster (default 100) |
| `temperature` | float | — | Sampling diversity; lower = closer to training mean (default 0.5) |
| `abdiffuser_dir` | str | — | Path to AbDiffuser repository clone (default `"AbDiffuser"`) |
| `output_dir` | str | — | Results directory; uses temp dir if None |

---

## Outputs

| Field | Type | Description |
|---|---|---|
| `designs` | list[dict] | Designs sorted by `log_likelihood` descending |
| `best_pdb_path` | str \| None | Top-ranked design's PDB path |
| `n_designs` | int | Number of designs produced |
| `cdrs_designed` | list[str] | Echo of input |

Each design dict: `{pdb_path, cdr_sequences: {H1: ..., H2: ..., H3: ...}, log_likelihood, energy}`.

Raises `RuntimeError("AbDiffuser failed: ...")` if the subprocess returns non-zero.

---

## Implementation

`tools/biologics/abdiffuser.py` — `run_abdiffuser()`

---

## Agent decision rules

- **Always start with `cdrs_to_design=["H3"]`**: H3 contributes ~70% of binding specificity; redesigning all six CDRs at once explodes the search space and rarely improves over H3-only
- **Confidence gate**: `log_likelihood` > median is a soft pass; pair with downstream `igfold` re-folding to validate predicted CDR geometry
- **n_diffusion_steps**: 100 default; drop to 50 only for diversity scouts. Below 50, designs degrade noticeably
- **temperature**: 0.3–0.5 for production (close to training distribution); 0.7–1.0 for exploring exotic CDRs
- **Validate downstream**: AbDiffuser's energy is a proxy — feed top designs through `igfold` then `binding_affinity/boltz2_affinity` (Ab–Ag mode) for an orthogonal score
- **vs ProteinMPNN**: ProteinMPNN sequences-only, faster (~seconds per design); AbDiffuser is slower (minutes) but jointly samples geometry. Use ProteinMPNN for breadth, AbDiffuser for high-stakes leads
- **Upstream**: `biologics/rfantibody` (initial framework), `biologics/igfold` (fold check), antibody crystal structures
- **Downstream**: `biologics/igfold` / `abodybuilder3` (re-fold designs), `binding_affinity/boltz2_affinity`, experimental validation

---

## Install

```bash
git clone https://github.com/igashov/AbDiffuser.git
cd AbDiffuser
# See repo README for env setup and pretrained weights (~2 GB)
pip install -e .
```

The wrapper expects to find AbDiffuser's `sample.py` entry script under `abdiffuser_dir`.

---

## References

- Martinkus K. et al. *NeurIPS* 2023 — "AbDiffuser: full-atom generation of in-vitro functioning antibodies"
- GitHub: https://github.com/igashov/AbDiffuser
- License: MIT
