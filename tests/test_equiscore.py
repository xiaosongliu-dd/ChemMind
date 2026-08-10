from __future__ import annotations

import pytest

from tools.analysis.equiscore import _norm, run_equiscore

SMILES = [
    "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cnnc2)n1",  # imatinib
    "CC#Cc1cc(-c2nc(Nc3ccc(N4CCC(N5CCN(C)CC5)CC4)cc3F)c(=O)[nH]2)ccn1",
    "COc1cc2c(Nc3ccc(F)c(Cl)c3)ncnc2cc1OCCCN1CCOCC1",
    "Cc1cc(Nc2ncc3cc(-c4ccccc4)c(=O)n(C)c3n2)ccc1F",
]

VINA   = [-9.2, -8.5, -7.8, -10.1]   # kcal/mol, most negative = best
GNINA  = [0.82, 0.75, 0.68, 0.91]     # [0,1], higher = best
BOLTZ  = [-11.3, -10.2, -9.1, -12.4]  # kcal/mol, most negative = best
PROLIF = [8, 6, 5, 10]                # contacts, higher = best


# ── normalisation ─────────────────────────────────────────────────────────────

def test_norm_higher_is_better_max_gives_one():
    result = _norm([1.0, 5.0, 3.0], higher_is_better=True)
    assert result[1] == pytest.approx(1.0)


def test_norm_lower_is_better_min_gives_one():
    result = _norm([-10.0, -5.0, -8.0], higher_is_better=False)
    assert result[0] == pytest.approx(1.0)


def test_norm_all_equal_returns_half():
    result = _norm([3.0, 3.0, 3.0], higher_is_better=True)
    assert all(v == pytest.approx(0.5) for v in result)


def test_norm_none_returns_empty():
    assert _norm(None, higher_is_better=True) == []


# ── full pipeline ─────────────────────────────────────────────────────────────

def test_all_components_returns_n_records():
    result = run_equiscore(
        SMILES, vina_scores=VINA, gnina_affinities=GNINA,
        boltz2_dg=BOLTZ, prolif_contact_counts=PROLIF,
    )
    assert result["n_scored"] == len(SMILES)
    assert len(result["ranked"]) == len(SMILES)


def test_ranked_sorted_descending():
    result = run_equiscore(SMILES, vina_scores=VINA, gnina_affinities=GNINA)
    scores = [r["consensus"] for r in result["ranked"]]
    assert scores == sorted(scores, reverse=True)


def test_ranks_are_sequential():
    result = run_equiscore(SMILES, vina_scores=VINA)
    ranks = [r["rank"] for r in result["ranked"]]
    assert sorted(ranks) == list(range(1, len(SMILES) + 1))


def test_best_vina_gets_rank_one_when_vina_only():
    # index 3 has the most negative Vina score → should rank first
    result = run_equiscore(SMILES, vina_scores=VINA)
    assert result["ranked"][0]["smiles"] == SMILES[3]


def test_consensus_in_zero_one_range():
    result = run_equiscore(
        SMILES, vina_scores=VINA, gnina_affinities=GNINA,
        boltz2_dg=BOLTZ, prolif_contact_counts=PROLIF,
    )
    for r in result["ranked"]:
        assert 0.0 <= r["consensus"] <= 1.0


# ── output schema ─────────────────────────────────────────────────────────────

def test_output_has_required_keys():
    result = run_equiscore(SMILES, vina_scores=VINA)
    assert "ranked" in result
    assert "n_scored" in result
    assert "weights_used" in result


def test_record_has_required_keys():
    result = run_equiscore(SMILES, vina_scores=VINA)
    rec = result["ranked"][0]
    for key in ("smiles", "vina_score", "gnina_affinity", "boltz2_dg",
                "prolif_contacts", "consensus", "rank"):
        assert key in rec


def test_missing_components_are_none():
    result = run_equiscore(SMILES, vina_scores=VINA)
    rec = result["ranked"][0]
    assert rec["gnina_affinity"] is None
    assert rec["boltz2_dg"] is None
    assert rec["prolif_contacts"] is None


# ── weight handling ───────────────────────────────────────────────────────────

def test_weights_sum_to_one_after_redistribution():
    result = run_equiscore(SMILES, vina_scores=VINA)  # only vina
    w = result["weights_used"]
    total = sum(w.values())
    assert total == pytest.approx(1.0)


def test_custom_weights_accepted():
    result = run_equiscore(
        SMILES, vina_scores=VINA, gnina_affinities=GNINA,
        weights={"vina": 0.6, "gnina": 0.4, "boltz2": 0.0, "prolif": 0.0},
    )
    assert result["n_scored"] == len(SMILES)


# ── error handling ────────────────────────────────────────────────────────────

def test_mismatched_score_length_raises():
    with pytest.raises(ValueError, match="vina_scores"):
        run_equiscore(SMILES, vina_scores=VINA[:2])


def test_empty_smiles_returns_empty():
    result = run_equiscore([])
    assert result["ranked"] == []
    assert result["n_scored"] == 0


def test_no_scores_returns_zero_consensus():
    result = run_equiscore(SMILES)  # no score lists at all
    for r in result["ranked"]:
        assert r["consensus"] == pytest.approx(0.0)
