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

## Execution (Python wrapper)

```python
# dependencies = ["boltz>=0.4", "requests"]

import os, json, base64, time, tempfile, subprocess
from pathlib import Path

def run_boltz2(
    protein_fasta: str,
    ligand_smiles: str | None = None,
    antibody_fasta: str | None = None,
    pocket_residues: list | None = None,
    n_samples: int = 5,
    use_nim: bool = True,
    nim_api_key: str | None = None,
) -> dict:
    if use_nim:
        return _run_nim(protein_fasta, ligand_smiles, antibody_fasta,
                        n_samples, nim_api_key or os.environ.get("NVIDIA_API_KEY", ""))
    return _run_local(protein_fasta, ligand_smiles, n_samples)


def _run_nim(protein_fasta, ligand_smiles, antibody_fasta, n_samples, api_key):
    import requests
    NIM_URL = "https://health.api.nvidia.com/v1/biology/mit/boltz2"
    payload = {
        "sequences": [{"protein": {"id": "A", "sequence": protein_fasta}}],
        "diffusion_samples": n_samples,
    }
    if ligand_smiles:
        payload["sequences"].append({"ligand": {"id": "LIG", "smiles": ligand_smiles}})
    if antibody_fasta:
        payload["sequences"].append({"protein": {"id": "B", "sequence": antibody_fasta}})

    resp = requests.post(
        NIM_URL,
        headers={"Authorization": f"Bearer {api_key}", "Content-Type": "application/json"},
        json=payload, timeout=300,
    )
    resp.raise_for_status()
    data = resp.json()

    outdir = Path(tempfile.mkdtemp(prefix="boltz2_"))
    pdb_paths = []
    for i, s in enumerate(data.get("structures", [])):
        p = outdir / f"sample_{i}.pdb"
        p.write_bytes(base64.b64decode(s["pdb"]))
        pdb_paths.append(str(p))

    best = data["structures"][0] if data.get("structures") else {}
    return {
        "pdb_path":            pdb_paths[0] if pdb_paths else None,
        "all_pdb_paths":       pdb_paths,
        "affinity_kcal_mol":   best.get("affinity_kcal_mol"),
        "affinity_confidence": best.get("affinity_confidence"),
        "ptm_score":           best.get("ptm"),
        "iptm_score":          best.get("iptm"),
        "runtime_s":           data.get("runtime_s"),
    }


def _run_local(protein_fasta, ligand_smiles, n_samples):
    workdir = Path(tempfile.mkdtemp(prefix="boltz2_local_"))
    fasta_path = workdir / "input.fasta"
    fasta_path.write_text(f">A\n{protein_fasta}\n")
    extra = ["--ligand", str(workdir / "ligand.smi")] if ligand_smiles else []
    if ligand_smiles:
        (workdir / "ligand.smi").write_text(ligand_smiles)

    cmd = [
        "uvx", "--from", "boltz>=0.4", "boltz", "predict",
        str(fasta_path), "--out_dir", str(workdir / "out"),
        "--num_diffusion_samples", str(n_samples),
        "--device", "cuda",
    ] + extra
    subprocess.run(cmd, check=True, capture_output=True)

    summary_path = workdir / "out" / "summary.json"
    summary = json.loads(summary_path.read_text()) if summary_path.exists() else {}
    pdb_paths = sorted((workdir / "out").glob("*.pdb"))
    return {
        "pdb_path":            str(pdb_paths[0]) if pdb_paths else None,
        "all_pdb_paths":       [str(p) for p in pdb_paths],
        "affinity_kcal_mol":   summary.get("affinity_kcal_mol"),
        "affinity_confidence": summary.get("affinity_confidence"),
        "ptm_score":           summary.get("ptm"),
        "iptm_score":          summary.get("iptm"),
        "runtime_s":           summary.get("runtime_s"),
    }
```

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
