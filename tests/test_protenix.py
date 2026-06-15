from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from tools.structure.protenix import run_protenix


_PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAI"
_LIGAND  = "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1"
_DNA     = "ATCGATCGATCG"
_RNA     = "AUCGAUCGAUCG"


def _make_protenix_output(outdir: Path, n_samples: int = 2) -> None:
    """Create a minimal Protenix output structure for testing."""
    for i in range(n_samples):
        sample_dir = outdir / f"seed-101_sample-{i}"
        sample_dir.mkdir(parents=True, exist_ok=True)
        (sample_dir / "model.cif").write_text(f"# fake CIF sample {i}\n")
        scores = {
            "ranking_score": 0.82 - i * 0.08,
            "ptm":           0.85 - i * 0.05,
            "iptm":          0.77 - i * 0.04,
            "mean_plddt":    83.5 - i * 2.0,
        }
        (sample_dir / "summary_confidences.json").write_text(json.dumps(scores))


def _fake_run(cmd, **kwargs):
    """Intercept subprocess.run and create fake Protenix output."""
    outdir = None
    for arg in cmd:
        if arg in ("-o", "--output_dir"):
            idx = cmd.index(arg)
            outdir = Path(cmd[idx + 1])
            break
    # The -o flag: 'protenix pred -i json -o outdir -s seed ...'
    if outdir is None:
        # Parse '-o' positional
        for i, a in enumerate(cmd):
            if a == "-o" and i + 1 < len(cmd):
                outdir = Path(cmd[i + 1])
                break
    if outdir:
        _make_protenix_output(outdir)
    return MagicMock(returncode=0)


# ── basic run ─────────────────────────────────────────────────────────────────

def test_returns_top_cif(tmp_path):
    with patch("subprocess.run", side_effect=_fake_run):
        result = run_protenix(_PROTEIN, output_dir=str(tmp_path))
    assert result["top_cif_path"] is not None
    assert result["top_cif_path"].endswith(".cif")


def test_all_cif_paths_non_empty(tmp_path):
    with patch("subprocess.run", side_effect=_fake_run):
        result = run_protenix(_PROTEIN, output_dir=str(tmp_path))
    assert len(result["all_cif_paths"]) > 0


def test_has_ligand_false_for_protein_only(tmp_path):
    with patch("subprocess.run", side_effect=_fake_run):
        result = run_protenix(_PROTEIN, output_dir=str(tmp_path))
    assert result["has_ligand"] is False
    assert result["has_nucleic_acid"] is False


# ── ligand / DNA / RNA ────────────────────────────────────────────────────────

def test_has_ligand_true_with_smiles(tmp_path):
    with patch("subprocess.run", side_effect=_fake_run):
        result = run_protenix(_PROTEIN, ligand_smiles=_LIGAND, output_dir=str(tmp_path))
    assert result["has_ligand"] is True
    assert result["n_chains"] == 2


def test_has_ligand_true_with_ccd(tmp_path):
    with patch("subprocess.run", side_effect=_fake_run):
        result = run_protenix(_PROTEIN, ligand_ccd_codes=["ATP"], output_dir=str(tmp_path))
    assert result["has_ligand"] is True


def test_dna_sets_nucleic_acid_flag(tmp_path):
    with patch("subprocess.run", side_effect=_fake_run):
        result = run_protenix(_PROTEIN, dna_sequences=[_DNA], output_dir=str(tmp_path))
    assert result["has_nucleic_acid"] is True


def test_rna_sets_nucleic_acid_flag(tmp_path):
    with patch("subprocess.run", side_effect=_fake_run):
        result = run_protenix(_PROTEIN, rna_sequences=[_RNA], output_dir=str(tmp_path))
    assert result["has_nucleic_acid"] is True


def test_n_chains_protein_only(tmp_path):
    with patch("subprocess.run", side_effect=_fake_run):
        result = run_protenix(_PROTEIN, output_dir=str(tmp_path))
    assert result["n_chains"] == 1


def test_n_chains_complex(tmp_path):
    with patch("subprocess.run", side_effect=_fake_run):
        result = run_protenix(_PROTEIN, ligand_smiles="CCO", dna_sequences=[_DNA],
                              output_dir=str(tmp_path))
    assert result["n_chains"] == 3


# ── confidence scores ─────────────────────────────────────────────────────────

def test_ranking_score_is_best(tmp_path):
    with patch("subprocess.run", side_effect=_fake_run):
        result = run_protenix(_PROTEIN, output_dir=str(tmp_path))
    assert result["ranking_score"] == pytest.approx(0.82)


def test_ptm_extracted(tmp_path):
    with patch("subprocess.run", side_effect=_fake_run):
        result = run_protenix(_PROTEIN, output_dir=str(tmp_path))
    assert result["ptm"] == pytest.approx(0.85)


def test_iptm_extracted(tmp_path):
    with patch("subprocess.run", side_effect=_fake_run):
        result = run_protenix(_PROTEIN, output_dir=str(tmp_path))
    assert result["iptm"] == pytest.approx(0.77)


def test_mean_plddt_extracted(tmp_path):
    with patch("subprocess.run", side_effect=_fake_run):
        result = run_protenix(_PROTEIN, output_dir=str(tmp_path))
    assert result["mean_plddt"] == pytest.approx(83.5)


def test_all_cifs_sorted_best_first(tmp_path):
    with patch("subprocess.run", side_effect=_fake_run):
        result = run_protenix(_PROTEIN, output_dir=str(tmp_path))
    assert result["all_cif_paths"][0] == result["top_cif_path"]


# ── subprocess arguments ──────────────────────────────────────────────────────

def test_subprocess_uses_protenix_pred(tmp_path):
    calls: list = []
    def capturing_run(cmd, **kwargs):
        calls.append(cmd)
        _make_protenix_output(tmp_path)
        return MagicMock(returncode=0)
    with patch("subprocess.run", side_effect=capturing_run):
        run_protenix(_PROTEIN, output_dir=str(tmp_path))
    assert calls[0][0] == "protenix"
    assert calls[0][1] == "pred"


def test_subprocess_passes_seed(tmp_path):
    calls: list = []
    def capturing_run(cmd, **kwargs):
        calls.append(cmd)
        _make_protenix_output(tmp_path)
        return MagicMock(returncode=0)
    with patch("subprocess.run", side_effect=capturing_run):
        run_protenix(_PROTEIN, seed=777, output_dir=str(tmp_path))
    assert "777" in calls[0]


def test_subprocess_passes_model(tmp_path):
    calls: list = []
    def capturing_run(cmd, **kwargs):
        calls.append(cmd)
        _make_protenix_output(tmp_path)
        return MagicMock(returncode=0)
    with patch("subprocess.run", side_effect=capturing_run):
        run_protenix(_PROTEIN, model="protenix_base_default_v1.0.0", output_dir=str(tmp_path))
    assert "protenix_base_default_v1.0.0" in calls[0]


def test_subprocess_failure_raises_runtime_error(tmp_path):
    with patch("subprocess.run", return_value=MagicMock(returncode=1, stderr="CUDA OOM")):
        with pytest.raises(RuntimeError, match="Protenix failed"):
            run_protenix(_PROTEIN, output_dir=str(tmp_path))


# ── input JSON construction ───────────────────────────────────────────────────

def test_input_json_written(tmp_path):
    with patch("subprocess.run", side_effect=_fake_run):
        run_protenix(_PROTEIN, output_dir=str(tmp_path))
    json_files = list(tmp_path.glob("*.json"))
    assert json_files, "no input JSON written"


def test_input_json_is_list(tmp_path):
    with patch("subprocess.run", side_effect=_fake_run):
        run_protenix(_PROTEIN[:20], output_dir=str(tmp_path))
    payload = json.loads(next(tmp_path.glob("*.json")).read_text())
    assert isinstance(payload, list)


def test_input_json_protein_chain_key(tmp_path):
    with patch("subprocess.run", side_effect=_fake_run):
        run_protenix(_PROTEIN[:20], output_dir=str(tmp_path))
    payload = json.loads(next(tmp_path.glob("*.json")).read_text())
    seqs = payload[0]["sequences"]
    assert any("proteinChain" in s for s in seqs)


def test_input_json_ligand_smiles_raw(tmp_path):
    with patch("subprocess.run", side_effect=_fake_run):
        run_protenix(_PROTEIN[:20], ligand_smiles="CCO", output_dir=str(tmp_path))
    payload = json.loads(next(tmp_path.glob("*.json")).read_text())
    lig_entries = [s["ligand"]["ligand"] for s in payload[0]["sequences"] if "ligand" in s]
    assert "CCO" in lig_entries


def test_input_json_ccd_prefixed(tmp_path):
    with patch("subprocess.run", side_effect=_fake_run):
        run_protenix(_PROTEIN[:20], ligand_ccd_codes=["ATP"], output_dir=str(tmp_path))
    payload = json.loads(next(tmp_path.glob("*.json")).read_text())
    lig_entries = [s["ligand"]["ligand"] for s in payload[0]["sequences"] if "ligand" in s]
    assert "CCD_ATP" in lig_entries


def test_input_json_dna_key(tmp_path):
    with patch("subprocess.run", side_effect=_fake_run):
        run_protenix(_PROTEIN[:20], dna_sequences=[_DNA], output_dir=str(tmp_path))
    payload = json.loads(next(tmp_path.glob("*.json")).read_text())
    seqs = payload[0]["sequences"]
    assert any("dnaSequence" in s for s in seqs)


# ── output schema ─────────────────────────────────────────────────────────────

def test_output_has_required_keys(tmp_path):
    with patch("subprocess.run", side_effect=_fake_run):
        result = run_protenix(_PROTEIN, output_dir=str(tmp_path))
    for key in ("top_cif_path", "all_cif_paths", "ranking_score", "ptm",
                "iptm", "mean_plddt", "has_ligand", "has_nucleic_acid",
                "n_chains", "output_dir", "seed", "n_sample", "model"):
        assert key in result, f"missing key: {key}"


def test_seed_in_result(tmp_path):
    with patch("subprocess.run", side_effect=_fake_run):
        result = run_protenix(_PROTEIN, seed=42, output_dir=str(tmp_path))
    assert result["seed"] == 42


def test_model_in_result(tmp_path):
    with patch("subprocess.run", side_effect=_fake_run):
        result = run_protenix(_PROTEIN, output_dir=str(tmp_path))
    assert result["model"] == "protenix-v2"
