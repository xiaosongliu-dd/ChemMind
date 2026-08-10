# Skill: molgpt

**Tier:** L1 — atomic tool skill
**Category:** generation
**Tool:** MolGPT — GPT-based conditional SMILES generation
**GPU required:** Yes for training; inference on CPU is feasible (~1 min / 100 molecules)

---

## What this skill does

MolGPT is a GPT language model trained on SMILES strings with optional **property conditioning** (QED, logP, SA score). At inference it autoregressively samples new SMILES tokens, optionally constrained to start from a scaffold fragment, and biased toward target property values via discrete condition tokens prepended to the prompt.

Key capabilities:
- Unconditional generation (sample from the learned SMILES distribution)
- Property-conditional generation: bias toward target QED / logP / SA
- Scaffold-constrained generation: fix a substructure and generate the rest
- Validity / uniqueness reporting

---

## When to call this skill

- Need diverse SMILES-space exploration without a 3D target pocket
- Want property-guided generation (high QED, specific logP window)
- Scaffold hopping: keep a pharmacophoric core, generate varied periphery
- Faster than REINVENT4 for quick diversity scouts (no RL loop)

Do NOT use if:
- Have a 3D pocket structure → use `generation/diffsbdd` for geometry-aware generation
- Need RL-guided multi-objective optimisation → use `generation/reinvent4`
- Need protein structure prediction → use `structure/esmfold` or `structure/colabfold`

---

## Inputs

| Field | Type | Required | Description |
|---|---|---|---|
| `n_molecules` | int | — | Number of molecules to generate (default 100) |
| `scaffold_smiles` | str | — | Scaffold SMILES fragment to fix (generation continues from this) |
| `temperature` | float | — | Sampling temperature (default 1.0); lower = closer to training mode |
| `target_qed` | float | — | Target QED (0–1); conditions the generation prompt |
| `target_logp` | float | — | Target logP |
| `target_sa` | float | — | Target SA score (1–10; lower = more synthesisable) |
| `molgpt_dir` | str | — | Path to MolGPT repository clone (default `"MolGPT"`) |
| `output_dir` | str | — | Results directory; temp dir if None |

---

## Outputs

| Field | Type | Description |
|---|---|---|
| `smiles` | list[str] | Unique valid SMILES |
| `n_generated` | int | Total tokens sampled |
| `n_valid` | int | Number passing RDKit parsing |
| `n_unique` | int | Unique valid SMILES |
| `validity_rate` | float | Fraction valid (0–1) |
| `uniqueness` | float | Unique / valid ratio |
| `scaffold` | str \| None | Echo of `scaffold_smiles` |
| `output_file` | str | Path to `generated.smi` |

Raises `RuntimeError("MolGPT not installed ...")` if the repo and `sample.py` are missing.

---

## Implementation

`tools/generation/molgpt.py` — `run_molgpt()`

Subprocess mode (default): calls `sample.py` from the MolGPT repo. Library mode: imports `model.GPT` directly (raises `NotImplementedError` for full weight setup — use subprocess mode). Generated SMILES are validated with RDKit and deduplicated.

---

## Agent decision rules

- **Temperature**: 0.7–1.0 for exploration; 0.5 for conservative sampling near training mean
- **Validity expectation**: ~85–95% validity at temperature=1.0 on a well-trained model
- **Property conditioning**: specify at most one or two properties — all three simultaneously pushes the model into a corner of property space that may be sparsely populated in training data
- **Scaffold**: provide as minimal valid SMILES (not SMARTS); the model continues after the last token
- **After generation**: apply `ligand_filter` (PAINS/BRENK) then `docking/*` on the top-ranked subset
- **Upstream**: training corpus from `data/chembl`; scaffold from existing actives
- **Downstream**: `rdkit_props`, `ligand_filter`, `docking/*`

---

## Install

```bash
git clone https://github.com/devalab/molgpt.git
cd molgpt
pip install -r requirements.txt
# Download pre-trained checkpoint (~100 MB) from repo releases
```

---

## References

- Bagal V. et al. *J. Chem. Inf. Model.* 2022, 62, 2064 — "MolGPT: Molecular Generation using a Transformer-Decoder Model"
- GitHub: https://github.com/devalab/molgpt
- License: MIT
