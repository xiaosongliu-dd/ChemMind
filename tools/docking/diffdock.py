from __future__ import annotations

import os
import subprocess
import tempfile
from pathlib import Path


NIM_URL = "https://health.api.nvidia.com/v1/biology/mit/diffdock"


def run_diffdock(
    protein_pdb: str,
    ligand_smiles: str,
    n_samples: int = 10,
    inference_steps: int = 20,
    use_nim: bool = True,
    nim_api_key: str | None = None,
) -> dict:
    if use_nim:
        return _run_nim(
            protein_pdb, ligand_smiles, n_samples, inference_steps,
            nim_api_key or os.environ.get("NVIDIA_API_KEY", ""),
        )
    return _run_local(protein_pdb, ligand_smiles, n_samples, inference_steps)


def _run_nim(protein_pdb, ligand_smiles, n_samples, steps, api_key):
    import requests

    payload = {
        "protein":          Path(protein_pdb).read_text(),
        "ligand":           ligand_smiles,
        "ligand_file_type": "smi",
        "num_poses":        n_samples,
        "time_divisions":   steps,
    }
    resp = requests.post(
        NIM_URL,
        headers={"Authorization": f"Bearer {api_key}", "Content-Type": "application/json"},
        json=payload,
        timeout=300,
    )
    resp.raise_for_status()
    data = resp.json()

    outdir = Path(tempfile.mkdtemp(prefix="diffdock_"))
    poses = []
    for i, p in enumerate(data.get("ligand_positions", [])):
        path = outdir / f"pose_{i}.pdb"
        path.write_text(p)
        poses.append({
            "pose_pdb":         str(path),
            "confidence_score": data.get("position_confidence", [None])[i],
            "rank":             i,
        })
    poses.sort(key=lambda x: -(x["confidence_score"] or 0))
    best = poses[0] if poses else {}
    return {
        "poses":            poses,
        "best_pose_pdb":    best.get("pose_pdb"),
        "confidence_score": best.get("confidence_score"),
        "n_poses":          len(poses),
        "runtime_s":        data.get("runtime_s"),
    }


def _run_local(protein_pdb, ligand_smiles, n_samples, steps):
    workdir = Path(tempfile.mkdtemp(prefix="diffdock_local_"))
    lig_path = workdir / "ligand.smi"
    lig_path.write_text(ligand_smiles)
    cmd = [
        "python", "inference.py",
        "--protein_path",        protein_pdb,
        "--ligand",              str(lig_path),
        "--out_dir",             str(workdir / "out"),
        "--samples_per_complex", str(n_samples),
        "--inference_steps",     str(steps),
    ]
    subprocess.run(cmd, check=True, capture_output=True, cwd="DiffDock")
    pose_files = sorted((workdir / "out").glob("rank*.pdb"))
    poses = [{"pose_pdb": str(p), "confidence_score": None, "rank": i}
             for i, p in enumerate(pose_files)]
    best = poses[0] if poses else {}
    return {
        "poses":            poses,
        "best_pose_pdb":    best.get("pose_pdb"),
        "confidence_score": None,
        "n_poses":          len(poses),
        "runtime_s":        None,
    }
