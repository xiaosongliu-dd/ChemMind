from __future__ import annotations

from types import ModuleType
from unittest.mock import MagicMock, patch

import pytest

from tools.binding_affinity.deeppurpose import _MODEL_UNITS, run_deeppurpose


def _mock_deeppurpose(predictions: list[float]):
    """Build a fake DeepPurpose module stack."""
    net = MagicMock()
    net.drug_encoding  = "MPNN"
    net.target_encoding = "CNN"
    net.predict.return_value = predictions

    dti_mod = MagicMock()
    dti_mod.model_pretrained.return_value = net

    utils_mod = MagicMock()
    utils_mod.data_process_repurpose_virtual_screening.return_value = MagicMock()

    return dti_mod, utils_mod, net


# ── import guard ──────────────────────────────────────────────────────────────

def test_raises_if_deeppurpose_not_installed(imatinib, abl1_fasta):
    with patch.dict("sys.modules", {"DeepPurpose": None, "DeepPurpose.DTI": None,
                                     "DeepPurpose.utils": None}):
        with pytest.raises(RuntimeError, match="DeepPurpose not installed"):
            run_deeppurpose(abl1_fasta, imatinib)


# ── single ligand ─────────────────────────────────────────────────────────────

def test_single_ligand_returns_one_prediction(imatinib, abl1_fasta):
    dti, utils, net = _mock_deeppurpose([8.3])
    with patch.dict("sys.modules", {"DeepPurpose": MagicMock(), "DeepPurpose.DTI": dti,
                                     "DeepPurpose.utils": utils}):
        result = run_deeppurpose(abl1_fasta, imatinib)

    assert len(result["predictions"]) == 1
    assert result["predictions"][0]["predicted_value"] == pytest.approx(8.3)


def test_single_ligand_best_smiles_matches(imatinib, abl1_fasta):
    dti, utils, net = _mock_deeppurpose([8.3])
    with patch.dict("sys.modules", {"DeepPurpose": MagicMock(), "DeepPurpose.DTI": dti,
                                     "DeepPurpose.utils": utils}):
        result = run_deeppurpose(abl1_fasta, imatinib)

    assert result["best_smiles"] == imatinib
    assert result["best_value"] == pytest.approx(8.3)


# ── batch mode ────────────────────────────────────────────────────────────────

def test_batch_repeats_target_for_each_ligand(all_smiles, abl1_fasta):
    scores = [7.0 + i * 0.1 for i in range(len(all_smiles))]
    dti, utils, net = _mock_deeppurpose(scores)
    with patch.dict("sys.modules", {"DeepPurpose": MagicMock(), "DeepPurpose.DTI": dti,
                                     "DeepPurpose.utils": utils}):
        run_deeppurpose(abl1_fasta, all_smiles)

    call_args = utils.data_process_repurpose_virtual_screening.call_args[0]
    target_list = call_args[1]
    assert len(target_list) == len(all_smiles)
    assert all(t == abl1_fasta for t in target_list)


def test_batch_n_ligands_matches_input(all_smiles, abl1_fasta):
    scores = [7.0] * len(all_smiles)
    dti, utils, net = _mock_deeppurpose(scores)
    with patch.dict("sys.modules", {"DeepPurpose": MagicMock(), "DeepPurpose.DTI": dti,
                                     "DeepPurpose.utils": utils}):
        result = run_deeppurpose(abl1_fasta, all_smiles)

    assert result["n_ligands"] == len(all_smiles)


# ── sorting ───────────────────────────────────────────────────────────────────

def test_predictions_sorted_descending(all_smiles, abl1_fasta):
    scores = [6.0, 9.5, 7.2, 8.8, 5.1, 10.0][:len(all_smiles)]
    dti, utils, net = _mock_deeppurpose(scores)
    with patch.dict("sys.modules", {"DeepPurpose": MagicMock(), "DeepPurpose.DTI": dti,
                                     "DeepPurpose.utils": utils}):
        result = run_deeppurpose(abl1_fasta, all_smiles)

    values = [p["predicted_value"] for p in result["predictions"]]
    assert values == sorted(values, reverse=True)


def test_best_value_is_max(all_smiles, abl1_fasta):
    scores = [6.0, 9.5, 7.2, 8.8, 5.1, 10.0][:len(all_smiles)]
    dti, utils, net = _mock_deeppurpose(scores)
    with patch.dict("sys.modules", {"DeepPurpose": MagicMock(), "DeepPurpose.DTI": dti,
                                     "DeepPurpose.utils": utils}):
        result = run_deeppurpose(abl1_fasta, all_smiles)

    assert result["best_value"] == pytest.approx(max(scores))


# ── output schema ─────────────────────────────────────────────────────────────

def test_output_has_required_keys(imatinib, abl1_fasta):
    dti, utils, net = _mock_deeppurpose([8.0])
    with patch.dict("sys.modules", {"DeepPurpose": MagicMock(), "DeepPurpose.DTI": dti,
                                     "DeepPurpose.utils": utils}):
        result = run_deeppurpose(abl1_fasta, imatinib)

    for key in ("predictions", "n_ligands", "model", "units", "best_smiles", "best_value"):
        assert key in result


def test_prediction_entry_has_required_keys(imatinib, abl1_fasta):
    dti, utils, net = _mock_deeppurpose([8.0])
    with patch.dict("sys.modules", {"DeepPurpose": MagicMock(), "DeepPurpose.DTI": dti,
                                     "DeepPurpose.utils": utils}):
        result = run_deeppurpose(abl1_fasta, imatinib)

    pred = result["predictions"][0]
    for key in ("smiles", "predicted_value", "units"):
        assert key in pred


# ── model units ───────────────────────────────────────────────────────────────

def test_ic50_model_reports_pic50_units(imatinib, abl1_fasta):
    dti, utils, net = _mock_deeppurpose([8.0])
    with patch.dict("sys.modules", {"DeepPurpose": MagicMock(), "DeepPurpose.DTI": dti,
                                     "DeepPurpose.utils": utils}):
        result = run_deeppurpose(abl1_fasta, imatinib, model="MPNN_CNN_BindingDB_IC50")

    assert result["units"] == "pIC50"
    assert result["predictions"][0]["units"] == "pIC50"


def test_kd_model_reports_pkd_units(imatinib, abl1_fasta):
    dti, utils, net = _mock_deeppurpose([8.0])
    with patch.dict("sys.modules", {"DeepPurpose": MagicMock(), "DeepPurpose.DTI": dti,
                                     "DeepPurpose.utils": utils}):
        result = run_deeppurpose(abl1_fasta, imatinib,
                                 model="Transformer_CNN_BindingDB_Kd")

    assert result["units"] == "pKd"
