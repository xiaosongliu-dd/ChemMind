from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Any


def run_prolif(
    protein_pdb: str,
    ligand_sdf: str,
    output_dir: str = ".",
    *,
    residue_radius: float = 6.0,
    count_contacts: bool = False,
) -> dict[str, Any]:
    """
    Compute a protein-ligand interaction fingerprint (IFP) using ProLIF.

    Parameters
    ----------
    protein_pdb : str
        Path to receptor PDB file.
    ligand_sdf : str
        Path to ligand SDF file (one or more poses, one per molecule block).
    output_dir : str
        Directory to write prolif_fingerprint.csv.
    residue_radius : float
        Radius (Å) around the ligand to select interacting protein residues.
    count_contacts : bool
        Return contact counts (True) vs binary presence/absence (False).

    Returns
    -------
    {
        interactions : dict[residue_id → dict[interaction_type → bool/int]]
        summary      : {n_hbond, n_hydrophobic, n_ionic, n_pistack, n_contacts}
        fingerprint  : list[int]   # binary IFP vector across all residue×interaction pairs
        n_poses      : int
        csv_path     : str | None
    }
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        import MDAnalysis  # noqa: F401
        import prolif  # noqa: F401
    except ImportError:
        return _prolif_subprocess(protein_pdb, ligand_sdf, output_dir)

    return _prolif_python(protein_pdb, ligand_sdf, output_dir, count_contacts)


def _prolif_python(
    protein_pdb: str,
    ligand_sdf: str,
    output_dir: Path,
    count_contacts: bool,
) -> dict[str, Any]:
    import MDAnalysis as mda
    import prolif as plf

    protein_u = mda.Universe(protein_pdb)
    ligand_u = mda.Universe(ligand_sdf)

    protein_mol = plf.Molecule.from_mda(protein_u.select_atoms("protein"))

    poses: list = []
    for _ts in ligand_u.trajectory:
        lig_mol = plf.Molecule.from_mda(ligand_u.atoms)
        poses.append(lig_mol)

    fp = plf.Fingerprint(
        interactions=[
            "HBDonor", "HBAcceptor", "Hydrophobic",
            "Ionic", "PiStacking", "EdgeToFace", "CationPi",
        ],
        count=count_contacts,
    )
    fp.run_from_iterable(poses, protein_mol)

    df = fp.to_dataframe()
    csv_path = output_dir / "prolif_fingerprint.csv"
    df.to_csv(csv_path)

    interactions: dict[str, dict[str, Any]] = {}
    summary = {
        "n_hbond": 0, "n_hydrophobic": 0,
        "n_ionic": 0, "n_pistack": 0, "n_contacts": 0,
    }
    fingerprint: list[int] = []

    if not df.empty:
        for col in df.columns:
            if not isinstance(col, tuple) or len(col) != 3:
                continue
            _lig, res, itype = col
            val = df[col].iloc[0]
            bit = int(bool(val))
            residue_key = str(res)
            if residue_key not in interactions:
                interactions[residue_key] = {}
            interactions[residue_key][itype] = int(val) if count_contacts else bool(val)
            fingerprint.append(bit)
            if bit:
                summary["n_contacts"] += 1
                itype_upper = itype.upper()
                if "HB" in itype_upper:
                    summary["n_hbond"] += 1
                elif "HYDROPHOBIC" in itype_upper:
                    summary["n_hydrophobic"] += 1
                elif "IONIC" in itype_upper:
                    summary["n_ionic"] += 1
                elif "PI" in itype_upper or "FACE" in itype_upper:
                    summary["n_pistack"] += 1

    return {
        "interactions": interactions,
        "summary":      summary,
        "fingerprint":  fingerprint,
        "n_poses":      len(poses),
        "csv_path":     str(csv_path),
    }


def _prolif_subprocess(
    protein_pdb: str,
    ligand_sdf: str,
    output_dir: Path,
) -> dict[str, Any]:
    """Fallback: run ProLIF in a child Python process."""
    csv_path = str(output_dir / "prolif_fingerprint.csv")
    script = (
        "import sys, json\n"
        "try:\n"
        "    import prolif as plf, MDAnalysis as mda\n"
        f"    u = mda.Universe({protein_pdb!r})\n"
        f"    lu = mda.Universe({ligand_sdf!r})\n"
        "    prot = plf.Molecule.from_mda(u.select_atoms('protein'))\n"
        "    poses = []\n"
        "    for ts in lu.trajectory:\n"
        "        poses.append(plf.Molecule.from_mda(lu.atoms))\n"
        "    fp = plf.Fingerprint()\n"
        "    fp.run_from_iterable(poses, prot)\n"
        "    df = fp.to_dataframe()\n"
        f"    df.to_csv({csv_path!r})\n"
        "    print(json.dumps({'n_poses': len(poses), 'n_contacts': int(df.values.sum())}))\n"
        "except Exception as e:\n"
        "    print(json.dumps({'error': str(e)}))\n"
    )
    proc = subprocess.run(
        ["python", "-c", script],
        capture_output=True, text=True, timeout=180,
    )
    if proc.returncode == 0 and proc.stdout.strip():
        import json
        data = json.loads(proc.stdout.strip())
        if "error" not in data:
            return {
                "interactions": {},
                "summary": {"n_contacts": data.get("n_contacts", 0),
                            "n_hbond": 0, "n_hydrophobic": 0,
                            "n_ionic": 0, "n_pistack": 0},
                "fingerprint":  [],
                "n_poses":      data.get("n_poses", 0),
                "csv_path":     csv_path,
            }

    raise RuntimeError(
        "ProLIF unavailable. Install with: pip install prolif MDAnalysis\n"
        f"stderr: {proc.stderr[:400]}"
    )
