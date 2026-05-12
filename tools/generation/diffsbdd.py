from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path


def run_diffsbdd(
    protein_pdb: str,
    pocket_center: tuple[float, float, float] | None = None,
    pocket_residues: list[str] | None = None,
    n_molecules: int = 100,
    diffsbdd_dir: str = "DiffSBDD",
    output_dir: str | None = None,
) -> dict:
    """
    Generate drug-like molecules conditioned on a 3D protein binding pocket
    using DiffSBDD (score-based diffusion over atoms).

    Provide EITHER `pocket_center` (x, y, z Å) OR `pocket_residues` (list of
    "<chain><resnum>" strings, e.g. ["A48", "A86"]).  When residues are given
    the pocket centre is computed as their Cα centroid.

    Returns:
        {
            smiles:        list[str],            # generated SMILES (valid only)
            sdf_path:      str | None,           # output SDF with 3D coords
            n_generated:   int,
            n_valid:       int,
            validity_rate: float,
            pocket_center: tuple[float, float, float] | None,
        }

    Raises RuntimeError("DiffSBDD failed: ...") on non-zero subprocess exit.
    """
    outdir = Path(output_dir) if output_dir else Path(tempfile.mkdtemp(prefix="diffsbdd_"))
    outdir.mkdir(parents=True, exist_ok=True)

    center = pocket_center or _centroid_from_residues(protein_pdb, pocket_residues or [])

    _run_inference(protein_pdb, center, n_molecules, diffsbdd_dir, outdir)

    smiles, sdf_path = _collect_outputs(outdir)
    n_valid = len(smiles)

    return {
        "smiles":        smiles,
        "sdf_path":      sdf_path,
        "n_generated":   n_molecules,
        "n_valid":       n_valid,
        "validity_rate": round(n_valid / max(n_molecules, 1), 3),
        "pocket_center": center,
        "output_dir":    str(outdir),
    }


def _run_inference(
    protein_pdb: str,
    center: tuple | None,
    n_molecules: int,
    diffsbdd_dir: str,
    outdir: Path,
) -> None:
    script = Path(diffsbdd_dir) / "test_conditioned_generation.py"
    if not script.exists():
        raise RuntimeError(
            f"DiffSBDD not found at {diffsbdd_dir}. "
            "Clone from https://github.com/arneschneuing/DiffSBDD and set diffsbdd_dir."
        )

    cmd = [
        "python", str(script),
        "--protein",      protein_pdb,
        "--n_samples",    str(n_molecules),
        "--output_dir",   str(outdir),
    ]
    if center:
        cmd += ["--pocket_center", f"{center[0]},{center[1]},{center[2]}"]

    result = subprocess.run(cmd, capture_output=True, text=True, cwd=diffsbdd_dir)
    if result.returncode != 0:
        raise RuntimeError(f"DiffSBDD failed:\n{result.stderr[-2000:]}")


def _collect_outputs(outdir: Path) -> tuple[list[str], str | None]:
    sdf_files = sorted(outdir.glob("*.sdf"))
    sdf_path = str(sdf_files[0]) if sdf_files else None
    smiles: list[str] = []

    smiles_files = sorted(outdir.glob("*.smi")) + sorted(outdir.glob("*.smiles"))
    if smiles_files:
        for line in smiles_files[0].read_text().splitlines():
            smi = line.split()[0].strip()
            if smi:
                smiles.append(smi)
    elif sdf_path:
        smiles = _sdf_to_smiles(sdf_path)

    return smiles, sdf_path


def _sdf_to_smiles(sdf_path: str) -> list[str]:
    try:
        from rdkit import Chem
        suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
        return [Chem.MolToSmiles(m) for m in suppl if m is not None]
    except ImportError:
        return []


def _centroid_from_residues(
    protein_pdb: str,
    residue_ids: list[str],
) -> tuple[float, float, float] | None:
    """Compute Cα centroid for a list of '<chain><resnum>' residue IDs."""
    if not residue_ids:
        return None

    res_set = set()
    for r in residue_ids:
        chain = r[0]
        resnum = r[1:]
        res_set.add((chain, resnum))

    coords: list[tuple[float, float, float]] = []
    try:
        for line in Path(protein_pdb).read_text().splitlines():
            if not line.startswith("ATOM") and not line.startswith("HETATM"):
                continue
            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue
            chain = line[21]
            resnum = line[22:26].strip()
            if (chain, resnum) in res_set:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                coords.append((x, y, z))
    except (OSError, ValueError):
        return None

    if not coords:
        return None
    n = len(coords)
    return (
        round(sum(c[0] for c in coords) / n, 3),
        round(sum(c[1] for c in coords) / n, 3),
        round(sum(c[2] for c in coords) / n, 3),
    )
