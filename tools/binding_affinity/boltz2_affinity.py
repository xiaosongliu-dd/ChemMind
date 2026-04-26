from __future__ import annotations

from tools.structure.boltz2 import run_boltz2


def run_boltz2_affinity(
    protein_fasta: str,
    ligand_smiles: str | list[str],
    n_samples: int = 5,
    use_nim: bool = True,
    nim_api_key: str | None = None,
) -> dict:
    """
    Predict binding affinity (ΔG kcal/mol) for one or more ligands against a target.

    Wraps run_boltz2() and returns only affinity-relevant fields.
    Structure outputs (pdb_path) are discarded — call boltz2 directly if you need the pose.

    Batch mode: pass ligand_smiles as a list to screen a series in one call.
    Results are sorted by affinity_kcal_mol ascending (tighter binders first).
    """
    if isinstance(ligand_smiles, str):
        return _predict_single(protein_fasta, ligand_smiles, n_samples, use_nim, nim_api_key)

    predictions = []
    for smi in ligand_smiles:
        try:
            result = _predict_single(protein_fasta, smi, n_samples, use_nim, nim_api_key)
            result["smiles"] = smi
            predictions.append(result)
        except Exception as e:
            predictions.append({
                "smiles":              smi,
                "affinity_kcal_mol":   None,
                "affinity_confidence": None,
                "iptm_score":          None,
                "error":               str(e),
            })

    predictions.sort(key=lambda x: (x["affinity_kcal_mol"] is None, x["affinity_kcal_mol"] or 0))
    return {
        "predictions": predictions,
        "n_ligands":   len(predictions),
        "n_failed":    sum(1 for p in predictions if p.get("error")),
        "best_smiles": predictions[0]["smiles"] if predictions else None,
        "best_affinity_kcal_mol": predictions[0]["affinity_kcal_mol"] if predictions else None,
    }


def _predict_single(
    protein_fasta: str,
    smiles: str,
    n_samples: int,
    use_nim: bool,
    nim_api_key: str | None,
) -> dict:
    result = run_boltz2(
        protein_fasta=protein_fasta,
        ligand_smiles=smiles,
        n_samples=n_samples,
        use_nim=use_nim,
        nim_api_key=nim_api_key,
    )
    return {
        "affinity_kcal_mol":   result["affinity_kcal_mol"],
        "affinity_confidence": result["affinity_confidence"],
        "iptm_score":          result["iptm_score"],
        "ptm_score":           result["ptm_score"],
        "runtime_s":           result["runtime_s"],
    }
