from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from tools.structure.alphafold3 import run_alphafold3


_PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL"
_LIGAND  = "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1"
_DNA     = "ATCGATCGATCG"
_RNA     = "AUCGAUCGAUCG"


def _make_af3_output(tmp_path: Path, job_name: str = "chemmind_af3",
                     n_models: int = 2, best_idx: int = 0) -> Path:
    """Create a minimal AF3 output directory for testing."""
    job_dir = tmp_path / job_name
    job_dir.mkdir(parents=True)
    scores = [0.82, 0.71, 0.65][:n_models]
    for i in range(n_models):
        cif = job_dir / f"{job_name}_model_{i}.cif"
        cif.write_text(f"# fake CIF {i}\ndata_model_{i}\n")
        conf = {
            "ranking_score": scores[i],
            "ptm":           0.85 - i * 0.05,
            "iptm":          0.78 - i * 0.04,
            "mean_plddt":    83.0 - i * 2.0,
        }
        (job_dir / f"{job_name}_summary_confidences_{i}.json").write_text(
            json.dumps(conf)
        )
    return job_dir


def _patched_run(tmp_path: Path, job_name: str = "chemmind_af3", n_models: int = 2):
    """Return a subprocess.run patch that creates fake AF3 output."""
    outdir_holder: list[Path] = []

    def fake_run(cmd, **kwargs):
        # Reconstruct output_dir from the --output_dir= argument
        for arg in cmd:
            if arg.startswith("--output_dir="):
                outdir = Path(arg.split("=", 1)[1])
                _make_af3_output(outdir, job_name, n_models)
                break
        return MagicMock(returncode=0)

    return fake_run


# ── no model_dir → RuntimeError ───────────────────────────────────────────────

def test_raises_if_no_model_dir(tmp_path, monkeypatch):
    monkeypatch.delenv("AF3_MODEL_DIR", raising=False)
    with pytest.raises(RuntimeError, match="model directory"):
        run_alphafold3(_PROTEIN, output_dir=str(tmp_path))


# ── basic protein-only run ────────────────────────────────────────────────────

def test_protein_only_returns_top_cif(tmp_path, monkeypatch):
    monkeypatch.setenv("AF3_MODEL_DIR", "/fake/weights")
    with patch("subprocess.run", side_effect=_patched_run(tmp_path)):
        result = run_alphafold3(_PROTEIN, output_dir=str(tmp_path))
    assert result["top_cif_path"] is not None
    assert result["top_cif_path"].endswith(".cif")


def test_protein_only_has_ligand_false(tmp_path, monkeypatch):
    monkeypatch.setenv("AF3_MODEL_DIR", "/fake/weights")
    with patch("subprocess.run", side_effect=_patched_run(tmp_path)):
        result = run_alphafold3(_PROTEIN, output_dir=str(tmp_path))
    assert result["has_ligand"] is False
    assert result["has_nucleic_acid"] is False


def test_n_chains_protein_only(tmp_path, monkeypatch):
    monkeypatch.setenv("AF3_MODEL_DIR", "/fake/weights")
    with patch("subprocess.run", side_effect=_patched_run(tmp_path)):
        result = run_alphafold3(_PROTEIN, output_dir=str(tmp_path))
    assert result["n_chains"] == 1


# ── protein–ligand complex ────────────────────────────────────────────────────

def test_protein_ligand_has_ligand_true(tmp_path, monkeypatch):
    monkeypatch.setenv("AF3_MODEL_DIR", "/fake/weights")
    with patch("subprocess.run", side_effect=_patched_run(tmp_path)):
        result = run_alphafold3(_PROTEIN, ligand_smiles=_LIGAND, output_dir=str(tmp_path))
    assert result["has_ligand"] is True


def test_protein_ligand_n_chains(tmp_path, monkeypatch):
    monkeypatch.setenv("AF3_MODEL_DIR", "/fake/weights")
    with patch("subprocess.run", side_effect=_patched_run(tmp_path)):
        result = run_alphafold3(_PROTEIN, ligand_smiles=_LIGAND, output_dir=str(tmp_path))
    assert result["n_chains"] == 2


def test_protein_ccd_ligand(tmp_path, monkeypatch):
    monkeypatch.setenv("AF3_MODEL_DIR", "/fake/weights")
    with patch("subprocess.run", side_effect=_patched_run(tmp_path)):
        result = run_alphafold3(_PROTEIN, ligand_ccd_codes=["ATP"], output_dir=str(tmp_path))
    assert result["has_ligand"] is True
    assert result["n_chains"] == 2


# ── nucleic acids ─────────────────────────────────────────────────────────────

def test_protein_dna_has_nucleic_acid(tmp_path, monkeypatch):
    monkeypatch.setenv("AF3_MODEL_DIR", "/fake/weights")
    with patch("subprocess.run", side_effect=_patched_run(tmp_path)):
        result = run_alphafold3(_PROTEIN, dna_sequences=[_DNA], output_dir=str(tmp_path))
    assert result["has_nucleic_acid"] is True
    assert result["n_chains"] == 2


def test_protein_rna_has_nucleic_acid(tmp_path, monkeypatch):
    monkeypatch.setenv("AF3_MODEL_DIR", "/fake/weights")
    with patch("subprocess.run", side_effect=_patched_run(tmp_path)):
        result = run_alphafold3(_PROTEIN, rna_sequences=[_RNA], output_dir=str(tmp_path))
    assert result["has_nucleic_acid"] is True


# ── multi-chain proteins ──────────────────────────────────────────────────────

def test_multi_chain_protein(tmp_path, monkeypatch):
    monkeypatch.setenv("AF3_MODEL_DIR", "/fake/weights")
    with patch("subprocess.run", side_effect=_patched_run(tmp_path)):
        result = run_alphafold3([_PROTEIN, _PROTEIN[:50]], output_dir=str(tmp_path))
    assert result["n_chains"] == 2


# ── confidence scores ─────────────────────────────────────────────────────────

def test_ranking_score_is_best(tmp_path, monkeypatch):
    monkeypatch.setenv("AF3_MODEL_DIR", "/fake/weights")
    with patch("subprocess.run", side_effect=_patched_run(tmp_path, n_models=2)):
        result = run_alphafold3(_PROTEIN, output_dir=str(tmp_path))
    assert result["ranking_score"] == pytest.approx(0.82)


def test_ptm_extracted(tmp_path, monkeypatch):
    monkeypatch.setenv("AF3_MODEL_DIR", "/fake/weights")
    with patch("subprocess.run", side_effect=_patched_run(tmp_path)):
        result = run_alphafold3(_PROTEIN, output_dir=str(tmp_path))
    assert result["ptm"] == pytest.approx(0.85)


def test_iptm_extracted(tmp_path, monkeypatch):
    monkeypatch.setenv("AF3_MODEL_DIR", "/fake/weights")
    with patch("subprocess.run", side_effect=_patched_run(tmp_path)):
        result = run_alphafold3(_PROTEIN, output_dir=str(tmp_path))
    assert result["iptm"] == pytest.approx(0.78)


def test_mean_plddt_extracted(tmp_path, monkeypatch):
    monkeypatch.setenv("AF3_MODEL_DIR", "/fake/weights")
    with patch("subprocess.run", side_effect=_patched_run(tmp_path)):
        result = run_alphafold3(_PROTEIN, output_dir=str(tmp_path))
    assert result["mean_plddt"] == pytest.approx(83.0)


def test_best_model_selected_by_ranking_score(tmp_path, monkeypatch):
    """Best model is the one with ranking_score=0.82, not the first by filename."""
    monkeypatch.setenv("AF3_MODEL_DIR", "/fake/weights")
    with patch("subprocess.run", side_effect=_patched_run(tmp_path, n_models=2)):
        result = run_alphafold3(_PROTEIN, output_dir=str(tmp_path))
    assert "model_0" in result["top_cif_path"]  # model_0 has 0.82 > model_1's 0.71


def test_all_cif_paths_sorted_best_first(tmp_path, monkeypatch):
    monkeypatch.setenv("AF3_MODEL_DIR", "/fake/weights")
    with patch("subprocess.run", side_effect=_patched_run(tmp_path, n_models=2)):
        result = run_alphafold3(_PROTEIN, output_dir=str(tmp_path))
    assert len(result["all_cif_paths"]) == 2
    assert result["all_cif_paths"][0] == result["top_cif_path"]


# ── seeds ─────────────────────────────────────────────────────────────────────

def test_custom_seeds_used(tmp_path, monkeypatch):
    monkeypatch.setenv("AF3_MODEL_DIR", "/fake/weights")
    custom = [42, 99, 7]
    with patch("subprocess.run", side_effect=_patched_run(tmp_path)):
        result = run_alphafold3(_PROTEIN, seeds=custom, output_dir=str(tmp_path))
    assert result["seeds_used"] == custom


def test_default_num_seeds(tmp_path, monkeypatch):
    monkeypatch.setenv("AF3_MODEL_DIR", "/fake/weights")
    with patch("subprocess.run", side_effect=_patched_run(tmp_path)):
        result = run_alphafold3(_PROTEIN, output_dir=str(tmp_path))
    assert result["seeds_used"] == [1, 2, 3, 4, 5]


# ── subprocess call arguments ─────────────────────────────────────────────────

def test_subprocess_called_with_json_path(tmp_path, monkeypatch):
    monkeypatch.setenv("AF3_MODEL_DIR", "/fake/weights")
    calls: list = []

    def fake_run(cmd, **kwargs):
        calls.append(cmd)
        for arg in cmd:
            if arg.startswith("--output_dir="):
                outdir = Path(arg.split("=", 1)[1])
                _make_af3_output(outdir)
        return MagicMock(returncode=0)

    with patch("subprocess.run", side_effect=fake_run):
        run_alphafold3(_PROTEIN, output_dir=str(tmp_path))

    assert any("--json_path=" in a for a in calls[0])


def test_subprocess_called_with_model_dir(tmp_path, monkeypatch):
    monkeypatch.setenv("AF3_MODEL_DIR", "/fake/weights")
    calls: list = []

    def fake_run(cmd, **kwargs):
        calls.append(cmd)
        for arg in cmd:
            if arg.startswith("--output_dir="):
                outdir = Path(arg.split("=", 1)[1])
                _make_af3_output(outdir)
        return MagicMock(returncode=0)

    with patch("subprocess.run", side_effect=fake_run):
        run_alphafold3(_PROTEIN, model_dir="/my/weights", output_dir=str(tmp_path))

    assert "--model_dir=/my/weights" in calls[0]


def test_subprocess_failure_raises_runtime_error(tmp_path, monkeypatch):
    monkeypatch.setenv("AF3_MODEL_DIR", "/fake/weights")
    with patch("subprocess.run", return_value=MagicMock(returncode=1, stderr="OOM error")):
        with pytest.raises(RuntimeError, match="AlphaFold 3 failed"):
            run_alphafold3(_PROTEIN, output_dir=str(tmp_path))


# ── output schema ─────────────────────────────────────────────────────────────

def test_output_has_required_keys(tmp_path, monkeypatch):
    monkeypatch.setenv("AF3_MODEL_DIR", "/fake/weights")
    with patch("subprocess.run", side_effect=_patched_run(tmp_path)):
        result = run_alphafold3(_PROTEIN, output_dir=str(tmp_path))
    for key in ("top_cif_path", "all_cif_paths", "ranking_score", "ptm",
                "iptm", "mean_plddt", "has_ligand", "has_nucleic_acid",
                "n_chains", "output_dir", "seeds_used"):
        assert key in result, f"missing key: {key}"


def test_output_dir_in_result(tmp_path, monkeypatch):
    monkeypatch.setenv("AF3_MODEL_DIR", "/fake/weights")
    with patch("subprocess.run", side_effect=_patched_run(tmp_path)):
        result = run_alphafold3(_PROTEIN, output_dir=str(tmp_path))
    assert result["output_dir"] == str(tmp_path)


# ── input JSON construction ───────────────────────────────────────────────────

def test_input_json_written_to_disk(tmp_path, monkeypatch):
    monkeypatch.setenv("AF3_MODEL_DIR", "/fake/weights")
    with patch("subprocess.run", side_effect=_patched_run(tmp_path)):
        run_alphafold3(_PROTEIN, output_dir=str(tmp_path))
    json_files = list(tmp_path.glob("*.json"))
    assert json_files, "no input JSON written"


def test_input_json_contains_protein_sequence(tmp_path, monkeypatch):
    monkeypatch.setenv("AF3_MODEL_DIR", "/fake/weights")
    with patch("subprocess.run", side_effect=_patched_run(tmp_path)):
        run_alphafold3(_PROTEIN[:20], output_dir=str(tmp_path))
    json_files = list(tmp_path.glob("*.json"))
    payload = json.loads(json_files[0].read_text())
    sequences = payload["sequences"]
    protein_seqs = [s["protein"]["sequence"] for s in sequences if "protein" in s]
    assert _PROTEIN[:20] in protein_seqs


def test_input_json_contains_ligand_smiles(tmp_path, monkeypatch):
    monkeypatch.setenv("AF3_MODEL_DIR", "/fake/weights")
    with patch("subprocess.run", side_effect=_patched_run(tmp_path)):
        run_alphafold3(_PROTEIN[:20], ligand_smiles="CCO", output_dir=str(tmp_path))
    json_files = list(tmp_path.glob("*.json"))
    payload = json.loads(json_files[0].read_text())
    ligand_smiles = [s["ligand"]["smiles"] for s in payload["sequences"] if "ligand" in s and "smiles" in s.get("ligand", {})]
    assert "CCO" in ligand_smiles


def test_input_json_dialect_alphafold3(tmp_path, monkeypatch):
    monkeypatch.setenv("AF3_MODEL_DIR", "/fake/weights")
    with patch("subprocess.run", side_effect=_patched_run(tmp_path)):
        run_alphafold3(_PROTEIN[:20], output_dir=str(tmp_path))
    json_files = list(tmp_path.glob("*.json"))
    payload = json.loads(json_files[0].read_text())
    assert payload["dialect"] == "alphafold3"
    assert payload["version"] == 2
