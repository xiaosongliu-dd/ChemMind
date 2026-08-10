# Skill: reinvent4

**Tier:** L1 — atomic tool skill
**Category:** generation
**Tool:** REINVENT 4 — reinforcement learning-guided molecular generation
**GPU required:** Yes for fast RL — A100 recommended; CPU feasible for sampling-only mode

---

## What this skill does

REINVENT 4 (Molecular AI, AstraZeneca) uses **reinforcement learning** to guide a pre-trained SMILES prior model toward molecules that score highly on a composite objective function. It supports scaffold-constrained generation (Scaffold Decorator mode) and multi-objective optimisation over arbitrary scoring components (docking, ADMET, similarity, custom SMARTS filters). The `sampling` mode draws directly from the prior without RL — useful for fast diverse generation.

Key capabilities:
- `sampling` mode: draw from prior, output scored SMILES in seconds
- `reinforcement` mode: RL loop with configurable scoring components (docking, properties, alerts)
- Scaffold SMARTS constraint: fix a core, generate variable attachments
- Diversity filter: prevents mode collapse via Murcko scaffold bucketing
- TOML configuration for full control over scoring components

---

## When to call this skill

- Multi-objective optimisation with a custom scoring function (docking + ADMET + synthesis)
- Scaffold-constrained lead optimisation (vary R-groups while keeping the core)
- Need RL-driven convergence toward a target property profile (vs random sampling)
- Want a well-validated, production-grade generative system

Do NOT use if:
- Just need fast unconditional generation → use `molgpt` (simpler, no RL overhead)
- Need 3D pocket-conditioned generation → use `diffsbdd`
- Need antibody or peptide generation → use `rfantibody` / `rfdiffusion_pep`

---

## Inputs

| Field | Type | Required | Description |
|---|---|---|---|
| `mode` | str | — | `"sampling"` (default) or `"reinforcement"` |
| `scaffold_smarts` | str | — | Scaffold SMARTS with attachment points (`[*]`) |
| `scoring_config` | dict | — | REINVENT4 scoring section dict (see docs for component library) |
| `n_steps` | int | — | RL steps or sampling batch count (default 100) |
| `batch_size` | int | — | SMILES per step (default 64); total = n_steps × batch_size for sampling |
| `model_file` | str | — | Path to prior `.prior` file (default `"priors/reinvent.prior"`) |
| `reinvent_bin` | str | — | REINVENT4 executable (default `"reinvent"`) |
| `output_dir` | str | — | Results directory; temp dir if None |

---

## Outputs

| Field | Type | Description |
|---|---|---|
| `smiles` | list[str] | Generated / optimised SMILES, sorted by score descending |
| `scores` | list[float] | Total score per SMILES (0–1) |
| `n_generated` | int | Number of SMILES returned |
| `n_steps_run` | int | Echo of `n_steps` |
| `mode` | str | Echo of mode |
| `output_dir` | str | Path containing CSV outputs and TensorBoard logs |

Raises `RuntimeError("REINVENT4 failed: ...")` on non-zero subprocess exit.

---

## Implementation

`tools/generation/reinvent4.py` — `run_reinvent4()`

Builds a TOML configuration, writes it to `output_dir/config.toml`, then calls the `reinvent` CLI. Output CSVs are parsed and results sorted by score descending. `tomli_w` is used for serialisation if available; falls back to a minimal custom TOML writer.

---

## Agent decision rules

- **sampling → reinforcement progression**: always scout with `sampling` first to check prior coverage, then run `reinforcement` to push toward the target
- **n_steps for RL**: 200–500 steps typical for convergence; 100 is a quick sanity check
- **Scoring config**: start with `"QED"` + `"SA_score"` + a custom alerts component; add docking only after property gates pass (docking scoring is ~10× slower per step)
- **Diversity filter**: leave on (default) to prevent mode collapse into one scaffold family
- **Score gate**: keep SMILES with `score > 0.6` for downstream docking; `score > 0.8` for high-priority leads
- **Upstream**: `data/chembl` (known actives as seed SMILES), `pocket_detect` (pocket for docking scorer)
- **Downstream**: `ligand_filter`, `rdkit_props`, `docking/gnina`, `binding_affinity/fep_openfe`

---

## Install

```bash
pip install reinvent  # PyPI package (REINVENT 4)
# or from source:
git clone https://github.com/MolecularAI/REINVENT4.git
cd REINVENT4 && pip install -e .
# Download priors (~200 MB) from repo releases
```

---

## References

- Loeffler H. et al. *J. Cheminform.* 2024, 16, 20 — "REINVENT 4: Modern AI–driven generative molecule design"
- GitHub: https://github.com/MolecularAI/REINVENT4
- License: Apache 2.0
