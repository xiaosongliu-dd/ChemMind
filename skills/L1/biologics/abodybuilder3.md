# Skill: abodybuilder3

**Tier:** L1 — atomic tool skill
**Category:** biologics
**Tool:** ABodyBuilder3 / ImmuneBuilder (antibody + nanobody structure prediction)
**GPU required:** Optional — fast on CPU; ~2–5× speedup with GPU

---

## What this skill does

ABodyBuilder3 (via the `ImmuneBuilder` package) predicts antibody and nanobody 3D structures from VH (+ VL) sequences using an ensemble of four IgFold-style networks fine-tuned with IMGT numbering and CDR annotation. Output PDBs include per-residue IMGT numbering and CDR loop labels in the B-factor column — useful for downstream analysis.

It's complementary to `igfold`: same task, similar speed, slightly different conformational biases. Running both and comparing gives a free-of-charge confidence cross-check on CDR-H3.

---

## When to call this skill

- After `proteinmpnn` or `rfantibody` — fold sequences with IMGT-annotated output
- When you need IMGT numbering preserved in the PDB (e.g. for downstream antibody-engineering tools that expect IMGT)
- As a cross-validation pass alongside `igfold` for high-stakes designs
- Whenever you want the fastest non-IgFold antibody folder

---

## Inputs

| Field | Type | Required | Description |
|---|---|---|---|
| `heavy_sequence` | str | ✓ | VH (or VHH for nanobody) amino acid sequence |
| `light_sequence` | str | — | VL sequence; omit for nanobody / single-domain |
| `numbering_scheme` | str | — | `"imgt"` (default), `"chothia"`, or `"kabat"` |
| `output_dir` | str | — | Where to write the PDB; uses temp dir if None |

---

## Outputs

| Field | Type | Description |
|---|---|---|
| `pdb_path` | str | Path to predicted structure PDB (with CDR labels in B-factor column) |
| `predicted_error` | float \| None | Model's self-estimated prediction error (Å); lower = more confident |
| `heavy_sequence` | str | Echo of input |
| `light_sequence` | str \| None | Echo of input |
| `mode` | str | `"antibody"` if light_sequence given, else `"nanobody"` |
| `numbering` | str | Echo of `numbering_scheme` |

Raises `RuntimeError("ImmuneBuilder not installed")` if `ImmuneBuilder` can't be imported.

---

## Implementation

`tools/biologics/abodybuilder3.py` — `run_abodybuilder3()`

---

## Agent decision rules

- **Confidence gate**: keep designs with `predicted_error < 1.5 Å`; CDR-H3 above this is unreliable
- **Numbering choice**: default `imgt` for modern ML pipelines; `chothia` only if downstream tools require it (older Rosetta scripts)
- **vs IgFold**: ABodyBuilder3 is slightly faster on CPU, IgFold has marginally better CDR-H3 accuracy on benchmark sets — preferences are dataset-dependent. Consensus of both is the safe bet
- **Nanobody mode dispatches automatically**: omit `light_sequence` and the tool routes to NanoBodyBuilder2 internally
- **Don't use as primary fold for non-antibody proteins**: ABodyBuilder3 is antibody-domain-specific; for general proteins use `structure/colabfold` or `structure/esmfold`
- **Upstream**: `biologics/proteinmpnn`, `biologics/rfantibody`, `data/chembl`
- **Downstream**: `docking/*`, `binding_affinity/boltz2_affinity` (Ab–Ag), `biologics/igfold` (cross-check)

---

## Install

```bash
uv add ImmuneBuilder
# Pulls torch + AbNumber for residue numbering
```

---

## References

- Abanades B. et al. *Nat. Commun. Biol.* 2023, 6, 575 — "ImmuneBuilder: Deep-Learning models for predicting the structures of immune proteins"
- GitHub: https://github.com/oxpig/ImmuneBuilder
- License: BSD-3 (Oxford Protein Informatics Group)
