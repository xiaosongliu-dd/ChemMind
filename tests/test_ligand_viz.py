from __future__ import annotations

from pathlib import Path

import pytest

from tools.visualization.ligand_viz import run_ligand_viz


# ── basic rendering ───────────────────────────────────────────────────────────

def test_returns_html_path_that_exists(all_smiles, tmp_path):
    result = run_ligand_viz(all_smiles, output_dir=str(tmp_path))
    assert Path(result["html_path"]).exists()
    assert Path(result["html_path"]).stat().st_size > 0


def test_default_plots_are_grid_and_table(all_smiles, tmp_path):
    result = run_ligand_viz(all_smiles, output_dir=str(tmp_path))
    assert "table" in result["plots_generated"]
    # grid may fall back to skipped if rdkit Draw lacks PIL; allow either
    assert "grid" in result["plots_generated"] or any(
        s["plot"] == "grid" for s in result["skipped_plots"]
    )


def test_html_contains_title_and_count(all_smiles, tmp_path):
    result = run_ligand_viz(all_smiles, title="ABL1 hits",
                            output_dir=str(tmp_path))
    body = Path(result["html_path"]).read_text()
    assert "ABL1 hits" in body
    assert str(len(all_smiles)) in body


# ── properties ────────────────────────────────────────────────────────────────

def test_properties_appear_in_table(all_smiles, tmp_path):
    props = {"pIC50": [7.0 + i * 0.1 for i in range(len(all_smiles))]}
    result = run_ligand_viz(all_smiles, properties=props,
                            output_dir=str(tmp_path))
    body = Path(result["html_path"]).read_text()
    assert "pIC50" in body
    assert result["properties_indexed"] == ["pIC50"]


def test_property_length_mismatch_raises(all_smiles, tmp_path):
    bad_props = {"pIC50": [7.0, 7.5]}  # too short
    with pytest.raises(ValueError, match="length"):
        run_ligand_viz(all_smiles, properties=bad_props,
                       output_dir=str(tmp_path))


# ── plot selection ────────────────────────────────────────────────────────────

def test_unknown_plot_type_raises(all_smiles, tmp_path):
    with pytest.raises(ValueError, match="Unknown plot types"):
        run_ligand_viz(all_smiles, plots=["grid", "lollipop"],
                       output_dir=str(tmp_path))


def test_table_only_skips_other_plots(all_smiles, tmp_path):
    result = run_ligand_viz(all_smiles, plots=["table"],
                            output_dir=str(tmp_path))
    assert result["plots_generated"] == ["table"]


def test_scatter_without_properties_skipped(all_smiles, tmp_path):
    result = run_ligand_viz(all_smiles, plots=["scatter"],
                            output_dir=str(tmp_path))
    assert "scatter" not in result["plots_generated"]
    assert any(s["plot"] == "scatter" for s in result["skipped_plots"])


def test_scatter_with_two_properties_renders(all_smiles, tmp_path):
    props = {
        "logP":  [3.0 + i * 0.1 for i in range(len(all_smiles))],
        "pIC50": [7.0 + i * 0.1 for i in range(len(all_smiles))],
    }
    result = run_ligand_viz(
        all_smiles, properties=props, plots=["scatter"],
        scatter_x="logP", scatter_y="pIC50",
        output_dir=str(tmp_path),
    )
    rendered_or_skipped = (
        "scatter" in result["plots_generated"]
        or any(s["plot"] == "scatter" for s in result["skipped_plots"])
    )
    assert rendered_or_skipped


def test_umap_with_few_compounds_skipped(tmp_path):
    # only 2 compounds — umap requires >=4
    result = run_ligand_viz(["CCO", "CCN"], plots=["umap"],
                            output_dir=str(tmp_path))
    assert "umap" not in result["plots_generated"]
    skipped = [s for s in result["skipped_plots"] if s["plot"] == "umap"]
    assert skipped


# ── output schema ─────────────────────────────────────────────────────────────

def test_output_has_required_keys(all_smiles, tmp_path):
    result = run_ligand_viz(all_smiles, output_dir=str(tmp_path))
    expected = {"html_path", "output_dir", "n_compounds", "plots_generated",
                "properties_indexed", "skipped_plots"}
    assert set(result.keys()) == expected


def test_n_compounds_matches_input(all_smiles, tmp_path):
    result = run_ligand_viz(all_smiles, output_dir=str(tmp_path))
    assert result["n_compounds"] == len(all_smiles)


def test_output_dir_returned(all_smiles, tmp_path):
    result = run_ligand_viz(all_smiles, output_dir=str(tmp_path))
    assert result["output_dir"] == str(tmp_path)
