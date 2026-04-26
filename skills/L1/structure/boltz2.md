# Skill: boltz2

**Tier:** L1 — atomic tool skill
**Category:** structure
**Tool:** Boltz-2 (MIT open weights + NVIDIA NIM)
**GPU required:** Yes — single A100/H100; ~20 s per complex

---

## What this skill does

Boltz-2 is a structural biology foundation model from MIT CSAIL that produces **both** a predicted 3D complex structure **and** a binding affinity score in a single forward pass. It supports:

- Protein structure prediction (single and multi-chain)
- Protein–ligand complex prediction with ΔG binding affinity (kcal/mol)
- Protein–protein, protein–RNA, protein–DNA complexes
- Antibody–antigen complexes (notable CDR-H3 loop improvement over Boltz-1)

Use Boltz-2 as the **primary structure + early-triage affinity engine** in ChemMind. Only escalate to dedicated docking (GNINA, AutoDock-GPU) when you need large-library throughput (>10k compounds) or explicit pose ensemble diversity.

---

## When to call this skill

- Target has no PDB — need predicted structure
- Need a rapid affinity estimate for an analogue series
- Antibody–antigen complex modelling
- Checking structural plausibility of a generated molecule before expensive FEP

---

## Inputs

| Field | Type | Required | Description |
|---|---|---|---|
| `protein_fasta` | str | ✓ | Target sequence; use `/` between chains |
| `ligand_smiles` | str | — | Small molecule SMILES |
| `antibody_fasta` | str | — | Heavy+light chain FASTA (chain-break `/`) |
| `pocket_residues` | list[int] | — | Residue IDs to bias pocket; from P2Rank |
| `n_samples` | int | — | Diffusion samples (default 5) |
| `use_nim` | bool | — | Use NVIDIA NIM API (default True) |
| `nim_api_key` | str | — | Reads `NVIDIA_API_KEY` env var if None |

---

## Outputs

| Field | Type | Description |
|---|---|---|
| `pdb_path` | str | Path to best-ranked predicted structure |
| `all_pdb_paths` | list[str] | All n_samples structures |
| `affinity_kcal_mol` | float | Predicted ΔG — lower = tighter binder |
| `affinity_confidence` | float | Model confidence 0–1 |
| `ptm_score` | float | Predicted TM-score (structure quality) |
| `iptm_score` | float | Interface pTM (complex quality) |
| `runtime_s` | float | Wall time |

---

## Implementation

`tools/structure/boltz2.py` — `run_boltz2()`

---

## Agent decision rules

- **Affinity triage**: keep `affinity_kcal_mol < -8.0`; escalate `< -10.0` to FEP
- **Quality gate**: discard if `iptm_score < 0.5` — complex poorly predicted
- **Ab–Ag**: pass heavy+light as `antibody_fasta` with chain-break `/`
- **Blind docking**: fine for initial triage; pass `pocket_residues` from `pocket_detect` for speed in VS campaigns
- **Upstream**: `pocket_detect` (P2Rank) → `boltz2`
- **Downstream**: `gnina` (pose diversity) → `fep_openfe` (ΔΔG) → `admetlab3`

---

## Install

```bash
export NVIDIA_API_KEY="nvapi-..."   # NIM mode — no local install

# Local weights mode
uv add boltz                         # downloads ~5 GB weights on first run
```

## References

- Passaro et al. 2025 — "Boltz-2: Towards Accurate and Efficient Binding Affinity Prediction" (bioRxiv 2025.06.14)
- NVIDIA NIM: https://docs.nvidia.com/nim/bionemo/boltz2/latest/overview.html
- GitHub: https://github.com/jwohlwend/boltz
- License: MIT
