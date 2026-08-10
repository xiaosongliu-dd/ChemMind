from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from tools.structure.chai1 import run_chai1, _build_fasta


_PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRP"
_LIGAND  = "CCO"
_DNA     = "ATCGATCG"
_RNA     = "AUCGAUCG"


def _make_candidates(cif_paths: list[Path], scores: list[float]) -> MagicMock:
    """Build a fake StructureCandidates object matching Chai-1's interface."""
    candidates = MagicMock()
    candidates.cif_paths = cif_paths

    ranking_data = []
    for s in scores:
        rd = MagicMock()
        rd.aggregate_score = MagicMock()
        rd.aggregate_score.item.return_value = s
        ranking_data.append(rd)
    candidates.ranking_data = ranking_data

    # plddt mock: supports plddt[0].mean().item() → 83.0 (no torch required)
    mean_item = MagicMock()
    mean_item.item.return_value = 83.0
    plddt_row = MagicMock()
    plddt_row.mean.return_value = mean_item
    plddt_mock = MagicMock()
    plddt_mock.__getitem__ = MagicMock(return_value=plddt_row)
    candidates.plddt = plddt_mock

    # sorted() returns a copy ordered best-first
    sorted_c = MagicMock()
    sorted_c.cif_paths = cif_paths  # assume already best-first in mock
    sorted_c.plddt = plddt_mock
    candidates.sorted.return_value = sorted_c

    return candidates


def _patched_inference(tmp_path: Path, scores: list[float] = None):
    """Return a fake run_inference that creates CIF files and returns candidates."""
    if scores is None:
        scores = [0.81, 0.74]

    def fake_inference(fasta_file, output_dir, **kwargs):
        cif_paths = []
        for i, s in enumerate(scores):
            cif = output_dir / f"pred.model_idx_{i}.cif"
            cif.write_text(f"# fake CIF {i}\n")
            cif_paths.append(cif)
        return _make_candidates(cif_paths, scores)

    return fake_inference


@contextmanager
def _chai_patch(fake_fn):
    """Patch chai_lab.chai1.run_inference via sys.modules (imported inside function body)."""
    chai1_mod = MagicMock()
    chai1_mod.run_inference = fake_fn
    with patch.dict("sys.modules", {"chai_lab": MagicMock(), "chai_lab.chai1": chai1_mod}):
        yield


# ── import guard ──────────────────────────────────────────────────────────────

def test_raises_if_chai_lab_not_installed():
    with patch.dict("sys.modules", {"chai_lab": None, "chai_lab.chai1": None}):
        with pytest.raises(RuntimeError, match="chai_lab not installed"):
            run_chai1(_PROTEIN)


# ── basic run ─────────────────────────────────────────────────────────────────

def test_returns_top_cif(tmp_path):
    with _chai_patch(_patched_inference(tmp_path)):
        result = run_chai1(_PROTEIN, output_dir=str(tmp_path))
    assert result["top_cif_path"] is not None
    assert result["top_cif_path"].endswith(".cif")


def test_all_cif_paths_non_empty(tmp_path):
    with _chai_patch(_patched_inference(tmp_path)):
        result = run_chai1(_PROTEIN, output_dir=str(tmp_path))
    assert len(result["all_cif_paths"]) > 0


def test_top_cif_is_first_in_all_cifs(tmp_path):
    with _chai_patch(_patched_inference(tmp_path)):
        result = run_chai1(_PROTEIN, output_dir=str(tmp_path))
    assert result["all_cif_paths"][0] == result["top_cif_path"]


# ── flags ─────────────────────────────────────────────────────────────────────

def test_has_ligand_false_for_protein_only(tmp_path):
    with _chai_patch(_patched_inference(tmp_path)):
        result = run_chai1(_PROTEIN, output_dir=str(tmp_path))
    assert result["has_ligand"] is False
    assert result["has_nucleic_acid"] is False


def test_has_ligand_true_with_smiles(tmp_path):
    with _chai_patch(_patched_inference(tmp_path)):
        result = run_chai1(_PROTEIN, ligand_smiles=_LIGAND, output_dir=str(tmp_path))
    assert result["has_ligand"] is True


def test_has_nucleic_acid_dna(tmp_path):
    with _chai_patch(_patched_inference(tmp_path)):
        result = run_chai1(_PROTEIN, dna_sequences=[_DNA], output_dir=str(tmp_path))
    assert result["has_nucleic_acid"] is True


def test_has_nucleic_acid_rna(tmp_path):
    with _chai_patch(_patched_inference(tmp_path)):
        result = run_chai1(_PROTEIN, rna_sequences=[_RNA], output_dir=str(tmp_path))
    assert result["has_nucleic_acid"] is True


def test_n_chains_protein_only(tmp_path):
    with _chai_patch(_patched_inference(tmp_path)):
        result = run_chai1(_PROTEIN, output_dir=str(tmp_path))
    assert result["n_chains"] == 1


def test_n_chains_complex(tmp_path):
    with _chai_patch(_patched_inference(tmp_path)):
        result = run_chai1(_PROTEIN, ligand_smiles=_LIGAND, dna_sequences=[_DNA],
                           output_dir=str(tmp_path))
    assert result["n_chains"] == 3


# ── scores ────────────────────────────────────────────────────────────────────

def test_aggregate_scores_non_empty(tmp_path):
    with _chai_patch(_patched_inference(tmp_path)):
        result = run_chai1(_PROTEIN, output_dir=str(tmp_path))
    assert len(result["aggregate_scores"]) > 0


def test_aggregate_scores_best_first(tmp_path):
    with _chai_patch(_patched_inference(tmp_path, scores=[0.65, 0.81])):
        result = run_chai1(_PROTEIN, output_dir=str(tmp_path))
    scores = result["aggregate_scores"]
    assert scores == sorted(scores, reverse=True)


def test_mean_plddt_extracted(tmp_path):
    with _chai_patch(_patched_inference(tmp_path)):
        result = run_chai1(_PROTEIN, output_dir=str(tmp_path))
    assert result["mean_plddt"] == pytest.approx(83.0, abs=1.0)


# ── seeds ─────────────────────────────────────────────────────────────────────

def test_default_seed_is_42(tmp_path):
    with _chai_patch(_patched_inference(tmp_path)):
        result = run_chai1(_PROTEIN, output_dir=str(tmp_path))
    assert result["seeds_used"] == [42]


def test_custom_seeds_used(tmp_path):
    call_seeds: list = []
    def capturing_inference(fasta_file, output_dir, seed, **kwargs):
        call_seeds.append(seed)
        cif = output_dir / "pred.model_idx_0.cif"
        cif.write_text("# fake")
        return _make_candidates([cif], [0.75])
    with _chai_patch(capturing_inference):
        result = run_chai1(_PROTEIN, seeds=[10, 20, 30], output_dir=str(tmp_path))
    assert result["seeds_used"] == [10, 20, 30]
    assert call_seeds == [10, 20, 30]


def test_multiple_seeds_gives_more_candidates(tmp_path):
    with _chai_patch(_patched_inference(tmp_path, scores=[0.80, 0.72])):
        result = run_chai1(_PROTEIN, seeds=[42, 99], output_dir=str(tmp_path))
    # 2 seeds × 2 predictions each = 4 total
    assert len(result["all_cif_paths"]) == 4


# ── FASTA construction ────────────────────────────────────────────────────────

def test_fasta_contains_protein_header():
    fasta = _build_fasta(["MKTAY"], [], [], [])
    assert ">protein|name=prot_0" in fasta
    assert "MKTAY" in fasta


def test_fasta_contains_ligand_smiles():
    fasta = _build_fasta(["MKTAY"], ["CCO"], [], [])
    assert ">ligand|name=lig_0" in fasta
    assert "CCO" in fasta


def test_fasta_contains_dna_header():
    fasta = _build_fasta(["MKTAY"], [], ["ATCG"], [])
    assert ">dna|name=dna_0" in fasta


def test_fasta_contains_rna_header():
    fasta = _build_fasta(["MKTAY"], [], [], ["AUCG"])
    assert ">rna|name=rna_0" in fasta


def test_fasta_multiple_proteins():
    fasta = _build_fasta(["MKTAY", "AKLMN"], [], [], [])
    assert ">protein|name=prot_0" in fasta
    assert ">protein|name=prot_1" in fasta


def test_fasta_written_to_disk(tmp_path):
    with _chai_patch(_patched_inference(tmp_path)):
        run_chai1(_PROTEIN[:10], output_dir=str(tmp_path))
    fasta_files = list(tmp_path.glob("*.fasta"))
    assert fasta_files, "no FASTA file written"


# ── output schema ─────────────────────────────────────────────────────────────

def test_output_has_required_keys(tmp_path):
    with _chai_patch(_patched_inference(tmp_path)):
        result = run_chai1(_PROTEIN, output_dir=str(tmp_path))
    for key in ("top_cif_path", "all_cif_paths", "aggregate_scores", "mean_plddt",
                "has_ligand", "has_nucleic_acid", "n_chains", "output_dir", "seeds_used"):
        assert key in result, f"missing key: {key}"


def test_output_dir_in_result(tmp_path):
    with _chai_patch(_patched_inference(tmp_path)):
        result = run_chai1(_PROTEIN, output_dir=str(tmp_path))
    assert result["output_dir"] == str(tmp_path)
