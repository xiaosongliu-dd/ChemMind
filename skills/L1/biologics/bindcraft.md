# Skill: bindcraft

**Tier:** L1 — atomic tool skill
**Category:** biologics
**Tool:** BindCraft — end-to-end protein binder design pipeline
**GPU required:** Yes — A100 strongly recommended; multi-hour wall time

---

## What this skill does

BindCraft is a turnkey pipeline that chains **RFdiffusion → ProteinMPNN → AlphaFold2** to design and validate de novo mini-protein binders against a user-specified target epitope. Unlike calling each step manually, BindCraft handles the design–score–filter loop end-to-end: it generates backbones, designs sequences on each, predicts the binder–target complex with AF2, and returns only designs passing user-specified ipTM / pTM thresholds.

Use this skill when you want a *single shot* at protein binder design with built-in quality gates, rather than orchestrating RFdiffusion + ProteinMPNN + AF2 yourself.

---

## When to call this skill

- Starting a fresh mini-protein binder campaign against a target PDB + epitope
- Want a turnkey "give me 20 high-confidence designs" pipeline rather than per-step orchestration
- Have time for a multi-hour run (typical: 50 designs → 4–12 h on a single A100)
- Comfortable with the BindCraft repo's opinionated defaults

Do NOT use if:
- You need fine-grained control of each step → call `rfdiffusion_pep` + `proteinmpnn` + `colabfold` individually
- Designing antibodies → use `rfantibody` or `abdiffuser`
- Designing peptide binders <30 aa → `rfdiffusion_pep` is more focused

---

## Inputs

| Field | Type | Required | Description |
|---|---|---|---|
| `target_pdb` | str | ✓ | Target protein PDB (held fixed during design) |
| `hotspot_residues` | list[str] | ✓ | Target epitope residues, e.g. `["A25", "A30", "A35"]` |
| `binder_length` | tuple[int,int] | — | (min, max) binder residue count (default `(50, 100)`) |
| `n_designs` | int | — | Backbone designs to generate (default 20) |
| `n_mpnn_sequences` | int | — | ProteinMPNN sequences per backbone (default 8) |
| `af2_filter_iptm` | float | — | Minimum ipTM to keep a design (default 0.55) |
| `af2_filter_ptm` | float | — | Minimum pTM to keep a design (default 0.55) |
| `bindcraft_dir` | str | — | Path to BindCraft repository clone (default `"BindCraft"`) |
| `output_dir` | str | — | Results directory; uses temp dir if None |

---

## Outputs

| Field | Type | Description |
|---|---|---|
| `designs` | list[dict] | Designs that passed both ipTM and pTM filters, sorted by ipTM descending |
| `all_designs` | list[dict] | All designs including those that failed filters |
| `n_total` | int | Total designs produced |
| `n_passed` | int | Number passing both filters |
| `best_pdb_path` | str \| None | Top-ranked passing design's PDB path |
| `best_iptm` | float \| None | ipTM of top design |
| `filter_iptm` | float | Echo of `af2_filter_iptm` |
| `filter_ptm` | float | Echo of `af2_filter_ptm` |

Each design dict: `{pdb_path, sequence, iptm, ptm, plddt, passed_filter}`.

Raises `RuntimeError("BindCraft failed: ...")` if the subprocess returns non-zero.

---

## Implementation

`tools/biologics/bindcraft.py` — `run_bindcraft()`

---

## Agent decision rules

- **Filter thresholds**: defaults `iptm=0.55, ptm=0.55` are loose. For wet-lab triage use `iptm=0.7, ptm=0.7` — much tighter, ~5–15% pass rate but real hits
- **n_designs**: 20 = quick scout; 100+ = serious campaign. Wall time scales linearly
- **Wall time budgeting**: ~10–25 min per (RFdiffusion + ProteinMPNN + AF2) cycle on A100. 100 designs = full overnight run
- **If `n_passed == 0`**: relax filters slightly first; if still zero, the epitope may be undruggable for this binder length — try `binder_length=(80, 150)` or different hotspots
- **Pair with `binding_affinity/fep_openfe`**: top BindCraft hits feed into FEP for affinity ranking before wet-lab
- **Upstream**: `pocket/p2rank` (epitope picking), target prep
- **Downstream**: experimental validation; `binding_affinity/boltz2_affinity` for additional scoring

---

## Install

```bash
git clone https://github.com/martinpacesa/BindCraft.git
cd BindCraft
# Follow repo README for AF2 weight download (~5 GB) and conda env setup
conda env create -f environment.yml
```

The wrapper expects to find BindCraft's `bindcraft.py` entry script under `bindcraft_dir`.

---

## References

- Pacesa M. et al. *bioRxiv* 2024 — "BindCraft: one-shot design of functional protein binders"
- GitHub: https://github.com/martinpacesa/BindCraft
- Underlying methods: RFdiffusion (Watson 2023), ProteinMPNN (Dauparas 2022), AlphaFold2 (Jumper 2021)
- License: MIT
