from __future__ import annotations

import base64
import os
import subprocess
import tempfile
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
        return _run_nim(
            protein_fasta, ligand_smiles, antibody_fasta,
            n_samples, nim_api_key or os.environ.get("NVIDIA_API_KEY", ""),
        )
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
        json=payload,
        timeout=300,
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
    import json

    workdir = Path(tempfile.mkdtemp(prefix="boltz2_local_"))
    fasta_path = workdir / "input.fasta"
    fasta_path.write_text(f">A\n{protein_fasta}\n")

    extra = []
    if ligand_smiles:
        lig_path = workdir / "ligand.smi"
        lig_path.write_text(ligand_smiles)
        extra = ["--ligand", str(lig_path)]

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
