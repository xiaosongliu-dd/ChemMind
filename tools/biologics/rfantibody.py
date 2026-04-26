from __future__ import annotations

import os
import subprocess
import tempfile
from pathlib import Path

# RFAntibody uses the RFdiffusion NIM endpoint with antibody-specific contigs
NIM_URL = "https://health.api.nvidia.com/v1/biology/ipd/rfdiffusion"


def run_rfantibody(
    target_pdb: str,
    hotspot_residues: list[str],
    cdr_lengths: dict[str, tuple[int, int]] | None = None,
    framework_pdb: str | None = None,
    n_designs: int = 10,
    use_nim: bool = True,
    nim_api_key: str | None = None,
) -> dict:
    """
    Design antibody CDR loops (especially CDR-H3) against a target epitope
    using RFdiffusion fine-tuned on antibody structures.

    Parameters
    ----------
    target_pdb       : antigen PDB; target is held fixed during diffusion
    hotspot_residues : epitope residues on the target, e.g. ["B25", "B30", "B35"]
    cdr_lengths      : per-CDR length range, e.g. {"H3": (8, 15), "H1": (5, 8)}
                       CDRs not listed use their framework lengths (no design)
    framework_pdb    : optional antibody framework PDB to scaffold onto;
                       if None, generates a full VH backbone de novo
    n_designs        : number of antibody backbone designs
    """
    if use_nim:
        return _run_nim(
            target_pdb, hotspot_residues, cdr_lengths, framework_pdb, n_designs,
            nim_api_key or os.environ.get("NVIDIA_API_KEY", ""),
        )
    return _run_local(target_pdb, hotspot_residues, cdr_lengths, framework_pdb, n_designs)


def _run_nim(target_pdb, hotspot_residues, cdr_lengths, framework_pdb, n_designs, api_key):
    import requests

    target_str = Path(target_pdb).read_text()

    # Build contig string for CDR design
    # Format: "[H_framework/CDR_H3:8-15/H_framework 0 target_chain:hotspots]"
    cdr_lengths = cdr_lengths or {"H3": (8, 15)}
    contig_parts = []
    for cdr, (lo, hi) in cdr_lengths.items():
        contig_parts.append(f"{lo}-{hi}")
    contig = "/".join(contig_parts) if contig_parts else "8-15"

    payload = {
        "pdb_str":      target_str,
        "contigs":      [contig],
        "hotspot_res":  hotspot_residues,
        "num_designs":  n_designs,
        "model_runner": "rf_antibody",  # antibody-specific diffusion model
    }
    if framework_pdb:
        payload["scaffold_pdb_str"] = Path(framework_pdb).read_text()

    resp = requests.post(
        NIM_URL,
        headers={"Authorization": f"Bearer {api_key}", "Content-Type": "application/json"},
        json=payload,
        timeout=600,
    )
    resp.raise_for_status()
    data = resp.json()

    outdir = Path(tempfile.mkdtemp(prefix="rfab_"))
    designs = []
    for i, pdb_str in enumerate(data.get("pdbs", [])):
        p = outdir / f"design_{i}.pdb"
        p.write_text(pdb_str)
        plddts = data.get("plddt_array", [[]])[i]
        designs.append({
            "pdb_path":   str(p),
            "mean_plddt": sum(plddts) / len(plddts) if plddts else None,
            "rank":       i,
        })

    designs.sort(key=lambda x: -(x["mean_plddt"] or 0))
    return {
        "designs":       designs,
        "best_pdb_path": designs[0]["pdb_path"] if designs else None,
        "n_designs":     len(designs),
        "hotspot_res":   hotspot_residues,
        "cdr_lengths":   cdr_lengths,
    }


def _run_local(target_pdb, hotspot_residues, cdr_lengths, framework_pdb, n_designs):
    workdir = Path(tempfile.mkdtemp(prefix="rfab_local_"))
    out_dir = workdir / "designs"
    out_dir.mkdir()

    cdr_lengths = cdr_lengths or {"H3": (8, 15)}
    hotspot_str = ",".join(hotspot_residues)

    # Build contig for CDR-H3 de novo with framework scaffold
    cdr_contig_parts = []
    for cdr, (lo, hi) in cdr_lengths.items():
        cdr_contig_parts.append(f"{lo}-{hi}")
    contig_str = "/".join(cdr_contig_parts)

    cmd = [
        "python", "scripts/run_inference.py",
        f"inference.output_prefix={out_dir}/rfab",
        f"inference.input_pdb={target_pdb}",
        f"contigmap.contigs=['{contig_str}/0 B1-200']",
        f"ppi.hotspot_res=[{hotspot_str}]",
        f"inference.num_designs={n_designs}",
        "inference.ckpt_override_path=models/rf_antibody.pt",
    ]
    if framework_pdb:
        cmd[4] = f"inference.input_pdb={framework_pdb}"

    subprocess.run(cmd, check=True, capture_output=True, cwd="RFdiffusion")
    pdb_paths = sorted(out_dir.glob("rfab*.pdb"))
    designs = [{"pdb_path": str(p), "mean_plddt": None, "rank": i}
               for i, p in enumerate(pdb_paths)]
    return {
        "designs":       designs,
        "best_pdb_path": designs[0]["pdb_path"] if designs else None,
        "n_designs":     len(designs),
        "hotspot_res":   hotspot_residues,
        "cdr_lengths":   cdr_lengths,
    }
