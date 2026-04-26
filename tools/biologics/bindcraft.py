from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path


def run_bindcraft(
    target_pdb: str,
    hotspot_residues: list[str],
    binder_length: tuple[int, int] = (50, 100),
    n_designs: int = 20,
    n_mpnn_sequences: int = 8,
    af2_filter_iptm: float = 0.55,
    af2_filter_ptm: float = 0.55,
    bindcraft_dir: str = "BindCraft",
    output_dir: str | None = None,
) -> dict:
    """
    End-to-end protein binder design using BindCraft:
    RFdiffusion → ProteinMPNN → ColabFold/AF2 filter → ranked outputs.

    BindCraft generates backbone designs with RFdiffusion, then sequences with
    ProteinMPNN, then validates/ranks each design by predicting the binder-target
    complex with AlphaFold2 and filtering by ipTM + pTM scores.

    Parameters
    ----------
    target_pdb         : target protein PDB (held fixed)
    hotspot_residues   : target epitope residues, e.g. ["A25", "A30"]
    binder_length      : (min, max) binder residue count
    n_designs          : RFdiffusion backbone designs to generate
    n_mpnn_sequences   : ProteinMPNN sequences per backbone
    af2_filter_iptm    : minimum ipTM to keep a design (complex quality)
    af2_filter_ptm     : minimum pTM to keep a design (binder quality)
    bindcraft_dir      : path to BindCraft repository clone
    output_dir         : output directory; uses temp dir if None
    """
    outdir = Path(output_dir or tempfile.mkdtemp(prefix="bindcraft_"))
    outdir.mkdir(parents=True, exist_ok=True)

    settings = {
        "target_pdb":        str(Path(target_pdb).resolve()),
        "hotspot_res":       hotspot_residues,
        "binder_length_min": binder_length[0],
        "binder_length_max": binder_length[1],
        "num_designs":       n_designs,
        "num_seqs":          n_mpnn_sequences,
        "af2_iptm_cutoff":   af2_filter_iptm,
        "af2_ptm_cutoff":    af2_filter_ptm,
        "output_dir":        str(outdir),
    }
    settings_path = outdir / "settings.json"
    settings_path.write_text(json.dumps(settings, indent=2))

    cmd = [
        "python", "bindcraft.py",
        "--settings", str(settings_path),
    ]
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=bindcraft_dir,
    )
    if result.returncode != 0:
        raise RuntimeError(f"BindCraft failed:\n{result.stderr[-2000:]}")

    return _collect_results(outdir, af2_filter_iptm, af2_filter_ptm)


def _collect_results(outdir: Path, iptm_cut: float, ptm_cut: float) -> dict:
    designs = []
    for pdb_path in sorted(outdir.glob("designs/*.pdb")):
        score_path = pdb_path.with_suffix(".json")
        scores = {}
        if score_path.exists():
            try:
                scores = json.loads(score_path.read_text())
            except json.JSONDecodeError:
                pass
        designs.append({
            "pdb_path":  str(pdb_path),
            "sequence":  scores.get("sequence"),
            "iptm":      scores.get("iptm"),
            "ptm":       scores.get("ptm"),
            "plddt":     scores.get("plddt"),
            "passed_filter": (
                (scores.get("iptm") or 0) >= iptm_cut and
                (scores.get("ptm")  or 0) >= ptm_cut
            ),
        })

    passed = [d for d in designs if d["passed_filter"]]
    passed.sort(key=lambda x: -(x["iptm"] or 0))

    return {
        "designs":        passed,
        "all_designs":    designs,
        "n_total":        len(designs),
        "n_passed":       len(passed),
        "best_pdb_path":  passed[0]["pdb_path"] if passed else None,
        "best_iptm":      passed[0]["iptm"] if passed else None,
        "filter_iptm":    iptm_cut,
        "filter_ptm":     ptm_cut,
    }
