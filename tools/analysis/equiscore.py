from __future__ import annotations

from typing import Any


def run_equiscore(
    smiles: list[str],
    vina_scores: list[float] | None = None,
    gnina_affinities: list[float] | None = None,
    boltz2_dg: list[float] | None = None,
    prolif_contact_counts: list[int] | None = None,
    weights: dict[str, float] | None = None,
) -> dict[str, Any]:
    """
    Consensus docking scorer: combine multiple scoring components into a
    single normalised rank for each compound.

    All score lists must match len(smiles) if provided.

    Score conventions (higher raw consensus = better):
        vina_scores            : kcal/mol, more negative = better  (flipped)
        gnina_affinities       : [0, 1] CNN affinity, higher = better
        boltz2_dg              : kcal/mol, more negative = better  (flipped)
        prolif_contact_counts  : integer contact count, higher = better

    Default weights (MolClaw-aligned, redistributed for missing components):
        vina=0.35, gnina=0.35, boltz2=0.15, prolif=0.15

    Returns
    -------
    {
        ranked   : list[{smiles, vina_score, gnina_affinity, boltz2_dg,
                         prolif_contacts, consensus, rank}]
        n_scored : int
        weights_used : dict   # actual weights after redistribution
    }
    """
    if not smiles:
        return {"ranked": [], "n_scored": 0, "weights_used": {}}

    n = len(smiles)
    _check_len("vina_scores",           vina_scores,           n)
    _check_len("gnina_affinities",      gnina_affinities,      n)
    _check_len("boltz2_dg",             boltz2_dg,             n)
    _check_len("prolif_contact_counts", prolif_contact_counts, n)

    default_w = {"vina": 0.35, "gnina": 0.35, "boltz2": 0.15, "prolif": 0.15}
    raw_weights = {**default_w, **(weights or {})}

    # Zero out weights for missing components, then renormalise
    present = {
        "vina":   vina_scores is not None,
        "gnina":  gnina_affinities is not None,
        "boltz2": boltz2_dg is not None,
        "prolif": prolif_contact_counts is not None,
    }
    active_total = sum(raw_weights[k] for k, v in present.items() if v)
    if active_total == 0:
        # No scores at all — return unranked
        ranked = [
            {"smiles": s, "vina_score": None, "gnina_affinity": None,
             "boltz2_dg": None, "prolif_contacts": None, "consensus": 0.0, "rank": i + 1}
            for i, s in enumerate(smiles)
        ]
        return {"ranked": ranked, "n_scored": n, "weights_used": {}}

    weights_used = {
        k: raw_weights[k] / active_total if present[k] else 0.0
        for k in raw_weights
    }

    # Normalise each component to [0, 1]
    vina_n  = _norm(vina_scores,           higher_is_better=False)
    gnina_n = _norm(gnina_affinities,      higher_is_better=True)
    boltz_n = _norm(boltz2_dg,             higher_is_better=False)
    plif_n  = _norm(prolif_contact_counts, higher_is_better=True)

    records = []
    for i, smi in enumerate(smiles):
        consensus = (
            weights_used["vina"]   * (vina_n[i]  if vina_n  else 0.0) +
            weights_used["gnina"]  * (gnina_n[i] if gnina_n else 0.0) +
            weights_used["boltz2"] * (boltz_n[i] if boltz_n else 0.0) +
            weights_used["prolif"] * (plif_n[i]  if plif_n  else 0.0)
        )
        records.append({
            "smiles":          smi,
            "vina_score":      vina_scores[i]           if vina_scores           else None,
            "gnina_affinity":  gnina_affinities[i]      if gnina_affinities      else None,
            "boltz2_dg":       boltz2_dg[i]             if boltz2_dg             else None,
            "prolif_contacts": prolif_contact_counts[i] if prolif_contact_counts else None,
            "consensus":       round(consensus, 4),
        })

    records.sort(key=lambda r: r["consensus"], reverse=True)
    for i, r in enumerate(records):
        r["rank"] = i + 1

    return {
        "ranked":       records,
        "n_scored":     n,
        "weights_used": weights_used,
    }


def _norm(values: list | None, *, higher_is_better: bool) -> list[float]:
    """Min-max normalise to [0, 1]; flip if lower-is-better."""
    if values is None:
        return []
    floats = [float(v) for v in values]
    lo, hi = min(floats), max(floats)
    if hi == lo:
        return [0.5] * len(floats)
    normed = [(v - lo) / (hi - lo) for v in floats]
    if not higher_is_better:
        normed = [1.0 - v for v in normed]
    return normed


def _check_len(name: str, values: list | None, expected: int) -> None:
    if values is not None and len(values) != expected:
        raise ValueError(
            f"{name} has {len(values)} entries but smiles has {expected}."
        )
