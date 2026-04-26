from __future__ import annotations

import tempfile
from pathlib import Path


def run_igfold(
    heavy_sequence: str,
    light_sequence: str | None = None,
    n_models: int = 4,
    do_refine: bool = False,
    output_dir: str | None = None,
) -> dict:
    """
    Predict 3D antibody or nanobody structure using IgFold.

    Parameters
    ----------
    heavy_sequence : VH (or VHH for nanobody) amino acid sequence
    light_sequence : VL sequence; omit for nanobody / single-domain antibody
    n_models       : ensemble size (1–4); more = slower but better pRMSD estimate
    do_refine      : Rosetta-based refinement (requires Rosetta); off by default
    output_dir     : path to write PDB; uses temp dir if None
    """
    try:
        from igfold import IgFoldRunner
    except ImportError:
        raise RuntimeError(
            "igfold not installed. Run: pip install igfold\n"
            "  First install: pip install torch torchvision"
        )

    outdir = Path(output_dir or tempfile.mkdtemp(prefix="igfold_"))
    outdir.mkdir(parents=True, exist_ok=True)
    pdb_path = str(outdir / "structure.pdb")

    sequences: dict[str, str] = {"H": heavy_sequence}
    if light_sequence:
        sequences["L"] = light_sequence

    runner = IgFoldRunner(num_models=n_models)
    pred = runner.fold(
        sequences=sequences,
        output_pdb=pdb_path,
        do_refine=do_refine,
    )

    mean_plddt = None
    mean_prmsd = None
    try:
        import torch
        mean_plddt = float(pred.plddt.mean())
        mean_prmsd = float(pred.prmsd.mean())
    except Exception:
        pass

    return {
        "pdb_path":       pdb_path,
        "mean_plddt":     mean_plddt,
        "mean_prmsd":     mean_prmsd,
        "heavy_sequence": heavy_sequence,
        "light_sequence": light_sequence,
        "mode":           "nanobody" if light_sequence is None else "antibody",
        "n_models":       n_models,
    }
