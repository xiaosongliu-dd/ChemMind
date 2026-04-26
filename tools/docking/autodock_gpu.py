from __future__ import annotations

import math
import subprocess
import tempfile
from pathlib import Path


def run_autodock_gpu(
    receptor_pdbqt: str,
    ligand_pdbqt_dir: str,
    gpf_path: str,
    n_runs: int = 20,
    heuristics: bool = True,
    autodock_gpu_bin: str = "autodock_gpu_128wi",
) -> dict:
    workdir = Path(tempfile.mkdtemp(prefix="adgpu_"))
    ligand_files = sorted(Path(ligand_pdbqt_dir).glob("*.pdbqt"))
    if not ligand_files:
        raise FileNotFoundError(f"No PDBQT files in {ligand_pdbqt_dir}")

    flist_path = workdir / "ligands.flist"
    flist_path.write_text("\n".join(str(f) for f in ligand_files))

    fld_path = gpf_path.replace(".gpf", ".maps.fld")
    cmd = [
        autodock_gpu_bin,
        "--ffile",  fld_path,
        "--lfile",  str(flist_path),
        "--nrun",   str(n_runs),
        "--resnam", str(workdir / "results"),
    ]
    if heuristics:
        cmd += ["--heuristics", "1"]

    subprocess.run(cmd, check=True, capture_output=True)

    results = [r for r in (_parse_dlg(f) for f in sorted(workdir.glob("results*.dlg"))) if r]
    results.sort(key=lambda x: x["best_energy"])
    return {
        "results":     results,
        "best_energy": results[0]["best_energy"] if results else None,
        "n_docked":    len(results),
        "failed":      [],
        "runtime_s":   None,
    }


def _parse_dlg(dlg_path: Path) -> dict | None:
    energies = []
    for line in dlg_path.read_text().splitlines():
        if line.startswith("DOCKED: USER    Estimated Free Energy of Binding"):
            try:
                energies.append(float(line.split("=")[1].split()[0]))
            except (IndexError, ValueError):
                pass
    if not energies:
        return None
    best = min(energies)
    ki = math.exp(best * 1000 / (1.987 * 298.15))
    return {
        "ligand_id":   dlg_path.stem,
        "best_energy": best,
        "ki_estimate": ki,
        "pose_pdbqt":  str(dlg_path),
    }
