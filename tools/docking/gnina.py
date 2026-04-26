from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path


def run_gnina(
    receptor_pdb: str,
    ligand_smiles: str | list,
    center_x: float,
    center_y: float,
    center_z: float,
    box_size: float = 22.5,
    n_poses: int = 9,
    cnn_scoring: str = "rescore",
    covalent_res: str | None = None,
    exhaustiveness: int = 8,
    gnina_bin: str = "gnina",
) -> dict:
    from rdkit import Chem
    from rdkit.Chem import AllChem

    if isinstance(ligand_smiles, str):
        ligand_smiles = [ligand_smiles]

    workdir = Path(tempfile.mkdtemp(prefix="gnina_"))
    all_poses = []

    for i, smi in enumerate(ligand_smiles):
        lig_sdf = _smiles_to_sdf(smi, workdir / f"lig_{i}.sdf")
        out_sdf = workdir / f"out_{i}.sdf"

        cmd = [
            gnina_bin,
            "--receptor",       receptor_pdb,
            "--ligand",         str(lig_sdf),
            "--out",            str(out_sdf),
            "--center_x",       str(center_x),
            "--center_y",       str(center_y),
            "--center_z",       str(center_z),
            "--size_x",         str(box_size),
            "--size_y",         str(box_size),
            "--size_z",         str(box_size),
            "--num_modes",      str(n_poses),
            "--cnn_scoring",    cnn_scoring,
            "--exhaustiveness", str(exhaustiveness),
        ]
        if covalent_res:
            cmd += ["--covalent_rec_atom", covalent_res]

        subprocess.run(cmd, capture_output=True, text=True)
        poses = _parse_gnina_output(str(out_sdf), smi)
        all_poses.append(poses)

    best = all_poses[0][0] if all_poses and all_poses[0] else {}
    return {
        "poses":         all_poses,
        "best_pose_sdf": best.get("pose_sdf"),
        "cnn_affinity":  best.get("cnn_affinity"),
        "vina_score":    best.get("vina_score"),
        "runtime_s":     None,
    }


def _smiles_to_sdf(smiles: str, path: Path) -> Path:
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(mol)
    w = Chem.SDWriter(str(path))
    w.write(mol)
    w.close()
    return path


def _parse_gnina_output(sdf_path: str, smiles: str) -> list[dict]:
    from rdkit import Chem

    poses = []
    try:
        suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
        for rank, mol in enumerate(suppl):
            if mol is None:
                continue
            props = mol.GetPropsAsDict()
            poses.append({
                "smiles":         smiles,
                "pose_sdf":       sdf_path,
                "rank":           rank,
                "cnn_affinity":   props.get("CNNaffinity"),
                "cnn_pose_score": props.get("CNNpose"),
                "vina_score":     props.get("minimizedAffinity"),
            })
    except Exception:
        pass
    return poses
