from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path


def run_abdiffuser(
    antibody_framework_pdb: str,
    antigen_pdb: str,
    cdrs_to_design: list[str],
    n_designs: int = 20,
    n_diffusion_steps: int = 100,
    temperature: float = 0.5,
    abdiffuser_dir: str = "AbDiffuser",
    output_dir: str | None = None,
) -> dict:
    """
    Co-design antibody CDR sequence + structure simultaneously using AbDiffuser.

    Unlike ProteinMPNN (sequence-only on a fixed backbone), AbDiffuser jointly
    samples CDR sequence and 3D conformation via score-based diffusion. This
    allows it to find sequences that are better coupled to their 3D structure.

    Parameters
    ----------
    antibody_framework_pdb : antibody PDB with framework regions (CDR residues
                              will be regenerated; mark them with chain B-factor)
    antigen_pdb            : antigen/target PDB (held fixed during design)
    cdrs_to_design         : list of CDR loops to redesign, e.g. ["H1", "H2", "H3"]
    n_designs              : number of CDR designs to sample
    n_diffusion_steps      : reverse diffusion steps (100 = default quality)
    temperature            : sampling diversity; lower = closer to training mean
    abdiffuser_dir         : path to AbDiffuser repository clone
    output_dir             : results directory; uses temp dir if None
    """
    outdir = Path(output_dir or tempfile.mkdtemp(prefix="abdiffuser_"))
    outdir.mkdir(parents=True, exist_ok=True)

    config = {
        "antibody_pdb":    str(Path(antibody_framework_pdb).resolve()),
        "antigen_pdb":     str(Path(antigen_pdb).resolve()),
        "cdrs_to_design":  cdrs_to_design,
        "n_samples":       n_designs,
        "n_steps":         n_diffusion_steps,
        "temperature":     temperature,
        "output_dir":      str(outdir),
    }
    config_path = outdir / "config.json"
    config_path.write_text(json.dumps(config, indent=2))

    cmd = [
        "python", "sample.py",
        "--config", str(config_path),
    ]
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=abdiffuser_dir,
    )
    if result.returncode != 0:
        raise RuntimeError(f"AbDiffuser failed:\n{result.stderr[-2000:]}")

    return _collect_results(outdir, cdrs_to_design)


def _collect_results(outdir: Path, cdrs_designed: list[str]) -> dict:
    designs = []
    for pdb_path in sorted(outdir.glob("sample_*.pdb")):
        score_path = pdb_path.with_suffix(".json")
        scores = {}
        if score_path.exists():
            try:
                scores = json.loads(score_path.read_text())
            except json.JSONDecodeError:
                pass
        cdr_seqs = {cdr: scores.get(f"seq_{cdr.lower()}") for cdr in cdrs_designed}
        designs.append({
            "pdb_path":    str(pdb_path),
            "cdr_sequences": cdr_seqs,
            "log_likelihood": scores.get("log_likelihood"),
            "energy":      scores.get("energy"),
        })

    designs.sort(key=lambda x: -(x["log_likelihood"] or float("-inf")))
    return {
        "designs":       designs,
        "best_pdb_path": designs[0]["pdb_path"] if designs else None,
        "n_designs":     len(designs),
        "cdrs_designed": cdrs_designed,
    }
