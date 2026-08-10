from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from tools.biologics.abodybuilder3 import run_abodybuilder3


_VH = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAR"
_VL = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPT"


def _mock_immune_builder(error_val=1.2):
    ab = MagicMock()
    ab.error = error_val
    ab.save = MagicMock()

    predictor = MagicMock()
    predictor.predict.return_value = ab

    abodybuilder2 = MagicMock(return_value=predictor)
    nanobodybuilder2 = MagicMock(return_value=predictor)

    immune_mock = MagicMock()
    immune_mock.ABodyBuilder2 = abodybuilder2
    immune_mock.NanoBodyBuilder2 = nanobodybuilder2

    return abodybuilder2, nanobodybuilder2, predictor, ab, immune_mock


# ── import guard ──────────────────────────────────────────────────────────────

def test_raises_if_immune_builder_not_installed():
    with patch.dict("sys.modules", {"ImmuneBuilder": None}):
        with pytest.raises(RuntimeError, match="ImmuneBuilder not installed"):
            run_abodybuilder3(_VH)


# ── antibody mode ─────────────────────────────────────────────────────────────

def test_antibody_mode_when_light_sequence_given():
    abb2, nbb2, predictor, ab, immune_mock = _mock_immune_builder()
    with patch.dict("sys.modules", {"ImmuneBuilder": immune_mock}):
        result = run_abodybuilder3(_VH, _VL)
    assert result["mode"] == "antibody"


def test_antibody_uses_abodybuilder2():
    abb2, nbb2, predictor, ab, immune_mock = _mock_immune_builder()
    with patch.dict("sys.modules", {"ImmuneBuilder": immune_mock}):
        run_abodybuilder3(_VH, _VL)
    abb2.assert_called_once()
    nbb2.assert_not_called()


def test_antibody_passes_both_chains():
    abb2, nbb2, predictor, ab, immune_mock = _mock_immune_builder()
    with patch.dict("sys.modules", {"ImmuneBuilder": immune_mock}):
        run_abodybuilder3(_VH, _VL)
    call_args = predictor.predict.call_args[0][0]
    assert call_args["H"] == _VH
    assert call_args["L"] == _VL


# ── nanobody mode ─────────────────────────────────────────────────────────────

def test_nanobody_mode_when_no_light_sequence():
    abb2, nbb2, predictor, ab, immune_mock = _mock_immune_builder()
    with patch.dict("sys.modules", {"ImmuneBuilder": immune_mock}):
        result = run_abodybuilder3(_VH)
    assert result["mode"] == "nanobody"


def test_nanobody_uses_nanobodybuilder2():
    abb2, nbb2, predictor, ab, immune_mock = _mock_immune_builder()
    with patch.dict("sys.modules", {"ImmuneBuilder": immune_mock}):
        run_abodybuilder3(_VH)
    nbb2.assert_called_once()
    abb2.assert_not_called()


def test_nanobody_passes_only_h_chain():
    abb2, nbb2, predictor, ab, immune_mock = _mock_immune_builder()
    with patch.dict("sys.modules", {"ImmuneBuilder": immune_mock}):
        run_abodybuilder3(_VH)
    call_args = predictor.predict.call_args[0][0]
    assert call_args == {"H": _VH}


# ── numbering scheme ──────────────────────────────────────────────────────────

def test_numbering_scheme_passed_to_builder():
    abb2, nbb2, predictor, ab, immune_mock = _mock_immune_builder()
    with patch.dict("sys.modules", {"ImmuneBuilder": immune_mock}):
        run_abodybuilder3(_VH, _VL, numbering_scheme="chothia")
    abb2.assert_called_once_with(numbering_scheme="chothia")


def test_default_numbering_is_imgt():
    abb2, nbb2, predictor, ab, immune_mock = _mock_immune_builder()
    with patch.dict("sys.modules", {"ImmuneBuilder": immune_mock}):
        result = run_abodybuilder3(_VH, _VL)
    assert result["numbering"] == "imgt"


# ── PDB is saved ──────────────────────────────────────────────────────────────

def test_save_called_on_prediction():
    abb2, nbb2, predictor, ab, immune_mock = _mock_immune_builder()
    with patch.dict("sys.modules", {"ImmuneBuilder": immune_mock}):
        run_abodybuilder3(_VH, _VL)
    ab.save.assert_called_once()


# ── output schema ─────────────────────────────────────────────────────────────

def test_output_has_required_keys():
    abb2, nbb2, predictor, ab, immune_mock = _mock_immune_builder()
    with patch.dict("sys.modules", {"ImmuneBuilder": immune_mock}):
        result = run_abodybuilder3(_VH, _VL)
    for key in ("pdb_path", "predicted_error", "heavy_sequence",
                 "light_sequence", "mode", "numbering"):
        assert key in result


def test_heavy_sequence_in_output():
    abb2, nbb2, predictor, ab, immune_mock = _mock_immune_builder()
    with patch.dict("sys.modules", {"ImmuneBuilder": immune_mock}):
        result = run_abodybuilder3(_VH, _VL)
    assert result["heavy_sequence"] == _VH
    assert result["light_sequence"] == _VL


def test_nanobody_light_sequence_none():
    abb2, nbb2, predictor, ab, immune_mock = _mock_immune_builder()
    with patch.dict("sys.modules", {"ImmuneBuilder": immune_mock}):
        result = run_abodybuilder3(_VH)
    assert result["light_sequence"] is None


def test_predicted_error_extracted():
    abb2, nbb2, predictor, ab, immune_mock = _mock_immune_builder(error_val=0.95)
    with patch.dict("sys.modules", {"ImmuneBuilder": immune_mock}):
        result = run_abodybuilder3(_VH, _VL)
    assert result["predicted_error"] == pytest.approx(0.95)
