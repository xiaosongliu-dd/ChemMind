from __future__ import annotations

import os
import subprocess
import tempfile
from pathlib import Path

NIM_URL = "https://health.api.nvidia.com/v1/biology/ipd/rfdiffusion"


def run_rfdiffusion_pep(
    target_pdb: str,
    hotspot_residues: list[str],
    binder_length: tuple[int, int] = (10, 20),
    n_designs: int = 10,
    noise_scale: float = 1.0,
    partial_diffusion_pdb: str | None = None,
    partial_T: int | None = None,
    use_nim: bool = True,
    nim_api_key: str | None = None,
) -> dict:
    """
    Design peptide / mini-protein binder backbones using RFdiffusion.

    Parameters
    ----------
    target_pdb        : receptor PDB (fixed during diffusion)
    hotspot_residues  : epitope residues to target, e.g. ["A25", "A30", "A35"]
                        format: "<chain><resnum>"
    binder_length     : (min, max) residue count for the designed binder
    n_designs         : number of backbone designs to generate
    noise_scale       : diffusion noise; 1.0 = full de novo, <1.0 = partial diffusion
    partial_diffusion_pdb : if provided, partially diffuse this starting backbone
    partial_T         : partial diffusion steps (used with partial_diffusion_pdb)
    """
    if use_nim:
        return _run_nim(
            target_pdb, hotspot_residues, binder_length, n_designs, noise_scale,
            nim_api_key or os.environ.get("NVIDIA_API_KEY", ""),
        )
    return _run_local(target_pdb, hotspot_residues, binder_length, n_designs,
                      noise_scale, partial_diffusion_pdb, partial_T)


def _run_nim(target_pdb, hotspot_residues, binder_length, n_designs, noise_scale, api_key):
    import requests

    pdb_str = Path(target_pdb).read_text()
    min_len, max_len = binder_length
    contig = f"{min_len}-{max_len}"

    payload = {
        "pdb_str":       pdb_str,
        "contigs":       [contig],
        "hotspot_res":   hotspot_residues,
        "num_designs":   n_designs,
        "noise_scale_ca":    noise_scale,
        "noise_scale_frame": noise_scale,
    }
    resp = requests.post(
        NIM_URL,
        headers={"Authorization": f"Bearer {api_key}", "Content-Type": "application/json"},
        json=payload,
        timeout=600,
    )
    resp.raise_for_status()
    data = resp.json()

    outdir = Path(tempfile.mkdtemp(prefix="rfdiff_pep_"))
    pdb_paths = []
    for i, pdb_str in enumerate(data.get("pdbs", [])):
        p = outdir / f"design_{i}.pdb"
        p.write_text(pdb_str)
        pdb_paths.append(str(p))

    plddts = data.get("plddt_array", [None] * len(pdb_paths))
    designs = [
        {"pdb_path": path, "mean_plddt": _mean(plddt), "rank": i}
        for i, (path, plddt) in enumerate(zip(pdb_paths, plddts))
    ]
    designs.sort(key=lambda x: -(x["mean_plddt"] or 0))
    return {
        "designs":       designs,
        "best_pdb_path": designs[0]["pdb_path"] if designs else None,
        "n_designs":     len(designs),
        "hotspot_res":   hotspot_residues,
        "binder_length": binder_length,
    }


def _run_local(target_pdb, hotspot_residues, binder_length, n_designs,
               noise_scale, partial_pdb, partial_T):
    workdir = Path(tempfile.mkdtemp(prefix="rfdiff_pep_local_"))
    out_dir = workdir / "designs"
    out_dir.mkdir()

    min_len, max_len = binder_length
    hotspot_str = ",".join(hotspot_residues)
    contig = f"[{min_len}-{max_len}/0 A1-{_count_residues(target_pdb)}]"

    cmd = [
        "python", "scripts/run_inference.py",
        f"inference.output_prefix={out_dir}/design",
        f"inference.input_pdb={target_pdb}",
        f"contigmap.contigs=['{contig}']",
        f"ppi.hotspot_res=[{hotspot_str}]",
        f"inference.num_designs={n_designs}",
        f"denoiser.noise_scale_ca={noise_scale}",
        f"denoiser.noise_scale_frame={noise_scale}",
    ]
    if partial_pdb and partial_T:
        cmd += [
            f"inference.input_pdb={partial_pdb}",
            f"diffuser.partial_T={partial_T}",
        ]

    subprocess.run(cmd, check=True, capture_output=True, cwd="RFdiffusion")
    pdb_paths = sorted(out_dir.glob("design*.pdb"))
    designs = [{"pdb_path": str(p), "mean_plddt": None, "rank": i}
               for i, p in enumerate(pdb_paths)]
    return {
        "designs":       designs,
        "best_pdb_path": designs[0]["pdb_path"] if designs else None,
        "n_designs":     len(designs),
        "hotspot_res":   hotspot_residues,
        "binder_length": binder_length,
    }


def _mean(values) -> float | None:
    if not values:
        return None
    try:
        vals = [float(v) for v in values if v is not None]
        return sum(vals) / len(vals) if vals else None
    except (TypeError, ValueError):
        return None


def _count_residues(pdb_path: str) -> int:
    residues = set()
    for line in Path(pdb_path).read_text().splitlines():
        if line.startswith("ATOM") and line[21] == "A":
            residues.add(int(line[22:26].strip()))
    return max(residues) if residues else 100
