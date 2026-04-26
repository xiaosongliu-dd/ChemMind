from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path


def run_vina(
    receptor_pdbqt: str,
    ligand_smiles: str | list,
    center_x: float,
    center_y: float,
    center_z: float,
    box_size: float = 22.5,
    exhaustiveness: int = 8,
    n_poses: int = 9,
    flex_residues_pdbqt: str | None = None,
    use_gpu: bool = True,
    vina_bin: str | None = None,
) -> dict:
    if isinstance(ligand_smiles, str):
        ligand_smiles = [ligand_smiles]

    if vina_bin is None:
        vina_bin = "vina_gpu" if use_gpu else "vina"

    workdir = Path(tempfile.mkdtemp(prefix="vina_"))
    all_poses = []

    for i, smi in enumerate(ligand_smiles):
        lig_pdbqt = _smiles_to_pdbqt(smi, workdir / f"lig_{i}.pdbqt")
        out_pdbqt = workdir / f"out_{i}.pdbqt"

        cmd = [
            vina_bin,
            "--receptor",       receptor_pdbqt,
            "--ligand",         str(lig_pdbqt),
            "--out",            str(out_pdbqt),
            "--center_x",       str(center_x),
            "--center_y",       str(center_y),
            "--center_z",       str(center_z),
            "--size_x",         str(box_size),
            "--size_y",         str(box_size),
            "--size_z",         str(box_size),
            "--exhaustiveness", str(exhaustiveness),
            "--num_modes",      str(n_poses),
        ]
        if flex_residues_pdbqt:
            cmd += ["--flex", flex_residues_pdbqt]

        subprocess.run(cmd, capture_output=True, text=True)
        poses = _parse_vina_pdbqt(str(out_pdbqt), smi)
        all_poses.append(poses)

    best = all_poses[0][0] if all_poses and all_poses[0] else {}
    return {
        "poses":           all_poses,
        "best_pose_pdbqt": best.get("pose_pdbqt"),
        "vina_score":      best.get("vina_score"),
        "runtime_s":       None,
    }


def _smiles_to_pdbqt(smiles: str, path: Path) -> Path:
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(mol)
    preparator = MoleculePreparation()
    preparator.prepare(mol)
    pdbqt_str = PDBQTWriterLegacy.write_string(preparator.setup)[0]
    path.write_text(pdbqt_str)
    return path


def _parse_vina_pdbqt(pdbqt_path: str, smiles: str) -> list[dict]:
    poses = []
    rank = 0
    current_score = None
    try:
        for line in open(pdbqt_path):
            if line.startswith("MODEL"):
                rank += 1
            elif "VINA RESULT" in line:
                try:
                    current_score = float(line.split()[3])
                except (IndexError, ValueError):
                    pass
            elif line.startswith("ENDMDL"):
                poses.append({
                    "smiles":     smiles,
                    "pose_pdbqt": pdbqt_path,
                    "vina_score": current_score,
                    "rank":       rank,
                })
    except FileNotFoundError:
        pass
    return poses
