# Skill: prolif

**Tier:** L1 — atomic tool skill
**Category:** analysis
**Tool:** ProLIF — Protein-Ligand Interaction Fingerprinting
**GPU required:** No — CPU-only; typically < 10 s per pose

---

## What this skill does

Computes a binary interaction fingerprint (IFP) between a docked ligand pose and its
protein receptor. Returns per-residue interaction types (H-bond donor/acceptor,
hydrophobic, ionic, π-stacking, edge-to-face, cation-π), a flat binary vector, and a
contact summary. Use after any docking step to understand *why* a compound scores well
and to enable SAR reasoning from structure.

---

## When to invoke

- After docking: characterise the interaction profile of a top-scoring compound
- SAR analysis: compare IFPs of actives vs inactives to identify pharmacophore contacts
- Selectivity: compare fingerprints against on-target vs off-target pocket
- Clustering: group docked poses by IFP similarity to remove redundant hits

---

## Signature

```python
prolif(
    protein_pdb=<str>,        # path to receptor PDB
    ligand_sdf=<str>,         # path to SDF file (one or more poses)
    output_dir=<str>,         # write prolif_fingerprint.csv here
    residue_radius=6.0,       # Å radius for residue selection
    count_contacts=False,     # True → return counts; False → binary
)
→ {
    interactions: {
        "<ChainResNum>": {
            "HBDonor":     bool,
            "HBAcceptor":  bool,
            "Hydrophobic": bool,
            "Ionic":       bool,
            "PiStacking":  bool,
            "EdgeToFace":  bool,
            "CationPi":    bool,
        },
        ...
    },
    summary: {
        n_hbond:       int,
        n_hydrophobic: int,
        n_ionic:       int,
        n_pistack:     int,
        n_contacts:    int,   # total non-zero contacts
    },
    fingerprint: list[int],  # flat binary IFP vector
    n_poses:     int,
    csv_path:    str,        # prolif_fingerprint.csv
}
```

---

## Interpretation guide

| Contact type | Drug-like significance |
|---|---|
| HBDonor / HBAcceptor | Strong directional interaction; essential for selectivity |
| Hydrophobic | Fills lipophilic sub-pockets; drives affinity |
| Ionic | Salt bridge; pH-dependent; contributes ~2–5 kcal/mol |
| PiStacking / EdgeToFace | Aromatic contacts; π-cation from Lys/Arg |
| CationPi | Basic amine over aromatic ring; potency boost |

**Quality gates:**
- `n_contacts == 0` → ligand is outside the pocket; check docking box center
- `n_hbond == 0` and target has a catalytic residue (Asp, Glu) → verify protonation state

---

## Dependencies

```bash
pip install prolif MDAnalysis
```

ProLIF requires MDAnalysis ≥ 2.7. Both are optional; the tool falls back to a
subprocess call if not importable in the current environment.

---

## References

- Bouysset C & Fiorucci S, *J. Cheminform.* 2021 — `prolif`
- MDAnalysis: Michaud-Agrawal et al., *J. Comput. Chem.* 2011
