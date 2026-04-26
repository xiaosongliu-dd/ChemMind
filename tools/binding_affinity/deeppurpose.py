from __future__ import annotations

# Output units per model — used to annotate predictions
_MODEL_UNITS = {
    "MPNN_CNN_BindingDB_IC50":       "pIC50",
    "Transformer_CNN_BindingDB_Kd":  "pKd",
    "CNN_CNN_BindingDB_IC50":        "pIC50",
    "MPNN_CNN_DAVIS":                "pKd",
    "Transformer_CNN_BindingDB_IC50":"pIC50",
}


def run_deeppurpose(
    protein_sequence: str,
    ligand_smiles: str | list[str],
    model: str = "MPNN_CNN_BindingDB_IC50",
) -> dict:
    """
    Predict binding affinity (pIC50 / pKd) from amino acid sequence + SMILES.

    Input notes
    -----------
    protein_sequence : full amino acid sequence as a string — no PDB or structure needed.
                       The model encodes the entire sequence, not just the pocket.
                       Accuracy is best for targets similar to BindingDB training data
                       (kinases, proteases, nuclear receptors). Degrades for novel folds.
    ligand_smiles    : single SMILES or a list for batch screening.

    Output units are pIC50 or pKd depending on model (higher = tighter binder).
    All predictions use the field name 'predicted_value'; see 'units' field for interpretation.
    """
    try:
        from DeepPurpose import DTI as models
        from DeepPurpose import utils
    except ImportError:
        raise RuntimeError("DeepPurpose not installed. Run: pip install DeepPurpose")

    if isinstance(ligand_smiles, str):
        ligand_smiles = [ligand_smiles]

    net   = models.model_pretrained(model=model)
    units = _MODEL_UNITS.get(model, "predicted_value")

    # data_process_repurpose_virtual_screening is the correct inference API.
    # data_process() is for training and requires labels — do not use it here.
    X_pred = utils.data_process_repurpose_virtual_screening(
        ligand_smiles,
        [protein_sequence] * len(ligand_smiles),
        net.drug_encoding,
        net.target_encoding,
    )
    y_pred = net.predict(X_pred)

    predictions = [
        {
            "smiles":          smi,
            "predicted_value": float(score),
            "units":           units,
        }
        for smi, score in zip(ligand_smiles, y_pred)
    ]
    predictions.sort(key=lambda x: -x["predicted_value"])

    return {
        "predictions": predictions,
        "n_ligands":   len(predictions),
        "model":       model,
        "units":       units,
        "best_smiles": predictions[0]["smiles"] if predictions else None,
        "best_value":  predictions[0]["predicted_value"] if predictions else None,
    }
