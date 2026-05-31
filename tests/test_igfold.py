from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from tools.biologics.igfold import run_igfold


_VH = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAR"
_VL = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPT"


def _mock_runner(mean_plddt=0.88, mean_prmsd=0.42):
    pred = MagicMock()
    pred.plddt.mean.return_value = mean_plddt   # plain float; igfold.py calls float() on it
    pred.prmsd.mean.return_value = mean_prmsd

    runner = MagicMock()
    runner.fold.return_value = pred

    runner_cls = MagicMock(return_value=runner)

    igfold_mock = MagicMock()
    igfold_mock.IgFoldRunner = runner_cls

    return runner_cls, runner, pred, igfold_mock


# ── import guard ──────────────────────────────────────────────────────────────

def test_raises_if_igfold_not_installed():
    with patch.dict("sys.modules", {"igfold": None}):
        with pytest.raises(RuntimeError, match="igfold not installed"):
            run_igfold(_VH)


# ── antibody mode ─────────────────────────────────────────────────────────────

def test_antibody_mode_when_light_sequence_given():
    runner_cls, runner, pred, igfold_mock = _mock_runner()
    with patch.dict("sys.modules", {"igfold": igfold_mock}):
        result = run_igfold(_VH, _VL)
    assert result["mode"] == "antibody"


def test_antibody_passes_both_chains():
    runner_cls, runner, pred, igfold_mock = _mock_runner()
    with patch.dict("sys.modules", {"igfold": igfold_mock}):
        run_igfold(_VH, _VL)
    call_kwargs = runner.fold.call_args.kwargs
    assert "H" in call_kwargs["sequences"]
    assert "L" in call_kwargs["sequences"]


def test_antibody_sequences_match_input():
    runner_cls, runner, pred, igfold_mock = _mock_runner()
    with patch.dict("sys.modules", {"igfold": igfold_mock}):
        run_igfold(_VH, _VL)
    seqs = runner.fold.call_args.kwargs["sequences"]
    assert seqs["H"] == _VH
    assert seqs["L"] == _VL


# ── nanobody mode ─────────────────────────────────────────────────────────────

def test_nanobody_mode_when_no_light_sequence():
    runner_cls, runner, pred, igfold_mock = _mock_runner()
    with patch.dict("sys.modules", {"igfold": igfold_mock}):
        result = run_igfold(_VH)
    assert result["mode"] == "nanobody"


def test_nanobody_passes_only_h_chain():
    runner_cls, runner, pred, igfold_mock = _mock_runner()
    with patch.dict("sys.modules", {"igfold": igfold_mock}):
        run_igfold(_VH)
    seqs = runner.fold.call_args.kwargs["sequences"]
    assert "H" in seqs
    assert "L" not in seqs


# ── n_models ──────────────────────────────────────────────────────────────────

def test_n_models_passed_to_runner():
    runner_cls, runner, pred, igfold_mock = _mock_runner()
    with patch.dict("sys.modules", {"igfold": igfold_mock}):
        run_igfold(_VH, n_models=2)
    runner_cls.assert_called_once_with(num_models=2)


# ── do_refine ─────────────────────────────────────────────────────────────────

def test_do_refine_false_by_default():
    runner_cls, runner, pred, igfold_mock = _mock_runner()
    with patch.dict("sys.modules", {"igfold": igfold_mock}):
        run_igfold(_VH, _VL)
    call_kwargs = runner.fold.call_args.kwargs
    assert call_kwargs["do_refine"] is False


def test_do_refine_propagated():
    runner_cls, runner, pred, igfold_mock = _mock_runner()
    with patch.dict("sys.modules", {"igfold": igfold_mock}):
        run_igfold(_VH, _VL, do_refine=True)
    call_kwargs = runner.fold.call_args.kwargs
    assert call_kwargs["do_refine"] is True


# ── output schema ─────────────────────────────────────────────────────────────

def test_output_has_required_keys():
    runner_cls, runner, pred, igfold_mock = _mock_runner()
    with patch.dict("sys.modules", {"igfold": igfold_mock}):
        result = run_igfold(_VH, _VL)
    for key in ("pdb_path", "mean_plddt", "mean_prmsd", "heavy_sequence",
                 "light_sequence", "mode", "n_models"):
        assert key in result


def test_heavy_sequence_in_output():
    runner_cls, runner, pred, igfold_mock = _mock_runner()
    with patch.dict("sys.modules", {"igfold": igfold_mock}):
        result = run_igfold(_VH, _VL)
    assert result["heavy_sequence"] == _VH
    assert result["light_sequence"] == _VL


def test_nanobody_light_sequence_none_in_output():
    runner_cls, runner, pred, igfold_mock = _mock_runner()
    with patch.dict("sys.modules", {"igfold": igfold_mock}):
        result = run_igfold(_VH)
    assert result["light_sequence"] is None


def test_n_models_in_output():
    runner_cls, runner, pred, igfold_mock = _mock_runner()
    with patch.dict("sys.modules", {"igfold": igfold_mock}):
        result = run_igfold(_VH, n_models=3)
    assert result["n_models"] == 3
