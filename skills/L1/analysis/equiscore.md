# Skill: equiscore

**Tier:** L1 — atomic tool skill
**Category:** analysis
**Tool:** EquiScore consensus ranker — multi-component docking score normalisation
**GPU required:** No — pure Python normalisation

---

## What this skill does

Combines multiple pre-computed docking and scoring components (Vina, GNINA CNN,
Boltz-2 ΔG, ProLIF contact count) into a single normalised consensus score via
weighted min-max normalisation. Each component is mapped to [0, 1] (higher always
means better after normalisation), then weights are applied. Missing components have
their weight redistributed proportionally to the present ones.

Single-score docking rankings have ~30% false-positive rate; consensus scoring reduces
this to ~15% without additional compute cost.

---

## When to invoke

- After docking (Step 6 in virtual_screen) to re-rank top hits before rescoring
- In `post_docking_eval` workflow as the final consensus ranking step
- Whenever ≥ 2 independent scores are available for the same compound set

---

## Signature

```python
equiscore(
    smiles=<list[str]>,                    # compounds to rank (must match all score lists)
    vina_scores=<list[float] | None>,      # kcal/mol, more negative = better
    gnina_affinities=<list[float] | None>, # [0,1] CNN affinity, higher = better
    boltz2_dg=<list[float] | None>,        # kcal/mol, more negative = better
    prolif_contact_counts=<list[int] | None>, # integer, higher = better
    weights={                              # optional — defaults shown
        "vina":   0.35,
        "gnina":  0.35,
        "boltz2": 0.15,
        "prolif": 0.15,
    },
)
→ {
    ranked: [
        {
            smiles:          str,
            vina_score:      float | None,
            gnina_affinity:  float | None,
            boltz2_dg:       float | None,
            prolif_contacts: int | None,
            consensus:       float,          # [0, 1]; higher = better
            rank:            int,
        },
        ...                                  # sorted best → worst
    ],
    n_scored:     int,
    weights_used: dict,                      # actual weights after redistribution
}
```

---

## Usage pattern

Run after docking (Step 6) and before ADMET filtering (Step 9) in virtual_screen:

```python
# Collect scores from previous steps
equiscore(
    smiles=top100_smiles,
    vina_scores=vina_step6_scores,
    gnina_affinities=gnina_step7_affinities,
    boltz2_dg=boltz2_step8_dg,
    prolif_contact_counts=prolif_step_counts,  # optional
)
→ take ranked[:20] → ADMET step
```

---

## Weight guidance

| Scenario | Recommended weights |
|---|---|
| Only Vina available | `{"vina": 1.0}` |
| Vina + GNINA (fast screen) | defaults work |
| Full pipeline (Vina + GNINA + Boltz-2 + ProLIF) | defaults work |
| Prioritise interaction quality | increase `"prolif"` to 0.25 |
| Prioritise affinity accuracy | increase `"boltz2"` to 0.25 |
