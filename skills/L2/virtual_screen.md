# Workflow: virtual_screen

**Tier:** L2 — workflow skill
**Purpose:** Screen commercially available compound libraries (Enamine REAL, ZINC22) or a
user-provided SMILES list against a target pocket to find buyable, drug-like hits. A cascade
of cheap-to-expensive filters (properties → alerts → docking → CNN rescoring → affinity →
ADMET) maximises throughput while controlling compute cost.

---

## When to invoke this workflow

- Find commercially available hits quickly (same-week delivery, no synthesis required)
- Dock and rank an existing SMILES list against a known target
- Complement to `hit_gen_sbdd` — run both to get novel + buyable series in parallel
- Keyword triggers: "screen library", "VS campaign", "find buyable actives", "dock these
  compounds", "ZINC", "Enamine REAL", "virtual screening", "ranked compound list"

---

## Required inputs

| Input | Description |
|---|---|
| Target protein | PDB ID, PDB/CIF file path, or amino-acid sequence |
| Compound source | `"enamine_real"`, `"zinc22"`, or an explicit SMILES list |
| Query SMILES | Optional — reference ligand for similarity-guided retrieval |
| Property filters | Optional — defaults: Lipinski, QED ≥ 0.4, SA ≤ 5 |
| Binding site hint | Optional — known residue list or X/Y/Z coords; auto-detected otherwise |

---

## L1 skill sequence

### Step 1 — Acquire target structure

**1a. If a PDB ID or file is provided:**
```
pdb_fetch(pdb_id=pdb_id, include_ligands=True)
→ {pdb_path: str, sequence: str, ligands: [{smiles, name, center_x, center_y, center_z}]}
```
If a co-crystal ligand is present, store its geometric center — use it as the binding-site
seed in Step 2 instead of running pocket detection.

**1b. If only an amino-acid sequence is provided:**
```
# Default: Protenix (Apache 2.0, best openly-licensed accuracy)
protenix(
    protein_sequences=sequence,
    n_sample=3,
    output_dir=workdir,
)
→ {top_cif_path: str, mean_plddt: float, ranking_score: float}

# Alternative for protein-only / speed priority (< 2 s):
esmfold(sequence=sequence, output_dir=workdir)
→ {pdb_path: str, mean_plddt: float}
```
Convert CIF → PDB if needed: `gemmi convert model.cif model.pdb`

**Quality gate:** if `mean_plddt < 50` or `ranking_score < 0.4`, halt and request an
experimental structure — low-confidence models give unreliable docking boxes.

---

### Step 2 — Detect binding pocket

```
pocket_detect(
    protein_pdb=pdb_path,
    tool="p2rank",                  # fall back to "fpocket" if p2rank unavailable
    output_dir=workdir,
)
→ {
    best_pocket: {center_x, center_y, center_z, volume, score, residues},
    pockets: list[dict],
    n_pockets: int,
}
```

Derive docking box:
```python
import math
center_x = best_pocket["center_x"]
center_y = best_pocket["center_y"]
center_z = best_pocket["center_z"]
box_size  = max(22.5, math.cbrt(best_pocket["volume"]) + 4)
```
Override with co-crystal ligand center from Step 1a if available.

**If n_pockets == 0:** the protein may be disordered — try `fpocket` as fallback; if still
no pocket, ask user to provide binding site coordinates manually.

---

### Step 3 — Assemble compound library

**Option A — commercial library (similarity search):**
```
library_search(
    library="enamine_real",         # or "zinc22"
    query_smiles=reference_smiles,  # omit for diverse random subset
    n_results=10000,
    similarity_threshold=0.4,
)
→ {smiles: list[str], catalog_ids: list[str], prices: list[float]}
```
If no reference SMILES: mine known actives first —
```
chembl(target_id=uniprot_or_chembl_id, standard_type="IC50", limit=50)
→ {compounds: [{smiles, pchembl_value}]}
```
use the highest-potency compound as the similarity query.

**Option B — user SMILES list:** proceed directly with the provided list; skip Step 3.

---

### Step 4 — Physicochemical pre-filter

```
rdkit_props(smiles=library_smiles)
→ {results: [{smiles, qed, sa_score, mw, logp, hbd, hba, lipinski_pass, veber_pass}]}
```

Apply ALL hard filters:
- `lipinski_pass == True` (MW ≤ 500, logP ≤ 5, HBD ≤ 5, HBA ≤ 10)
- `qed ≥ 0.4`
- `sa_score ≤ 5`

Expected retention: 40–60 % of input. Cap at **5 000 compounds** before proceeding to docking.

---

### Step 5 — Structural alert filter

```
ligand_filter(
    smiles=step4_passed_smiles,
    filter_sets="default",          # PAINS + BRENK
)
→ {passed_smiles: list[str], failed_smiles: list[str], n_passed: int, n_failed: int}
```
Expected retention: 85–90 % of Step 4 output.

---

### Step 6 — High-throughput docking (scale-dependent)

Choose route based on compound count after Step 5:

**Route A — ≤ 200 compounds: GNINA (simplest, CNN-scored, accepts PDB + SMILES directly)**
```
gnina(
    receptor_pdb=pdb_path,
    ligand_smiles=passed_smiles,    # list accepted
    center_x=center_x,
    center_y=center_y,
    center_z=center_z,
    box_size=box_size,
    n_poses=3,
    cnn_scoring="rescore",
    exhaustiveness=4,               # fast screening mode
)
→ {poses: list, best_pose_sdf: str, cnn_affinity: float, vina_score: float, runtime_s: float}
```
Sort by `vina_score` ascending (most negative = best). Take top 50.

**Route B — 201–5 000 compounds: Vina (CPU/GPU batch, accepts SMILES)**
```
vina(
    receptor_pdbqt=receptor_pdbqt,  # convert PDB → PDBQT: obabel receptor.pdb -O receptor.pdbqt
    ligand_smiles=passed_smiles,
    center_x=center_x,
    center_y=center_y,
    center_z=center_z,
    box_size=box_size,
    exhaustiveness=4,
    n_poses=3,
    use_gpu=True,
)
→ {poses: list, best_pose_pdbqt: str, vina_score: float, runtime_s: float}
```
Sort by `vina_score` ascending. Take top 100.

**Route C — > 5 000 compounds: AutoDock-GPU (GPU batch, highest throughput)**
Requires PDBQT preparation for all inputs and a grid parameter file (GPF):
- Receptor: `prepare_receptor4.py -r receptor.pdb -o receptor.pdbqt` (MGLTools)
- Ligands: `mk_prepare_ligand.py -i compound.sdf -o compound.pdbqt` (Meeko), one file per compound, all in `ligand_dir/`
- GPF: `prepare_gpf4.py -r receptor.pdbqt -l representative_ligand.pdbqt -x {center_x},{center_y},{center_z} -s {box_size}`
```
autodock_gpu(
    receptor_pdbqt=receptor_pdbqt,
    ligand_pdbqt_dir=ligand_dir,
    gpf_path=gpf_path,
    n_runs=20,
    heuristics=True,
)
→ {results: [{smiles, binding_energy, ki_nm}], best_energy: float, n_docked: int, failed: list}
```
Sort by `binding_energy` ascending. Take top 100.

---

### Step 7 — Accurate CNN rescoring (top 50–100 from Step 6)

```
gnina(
    receptor_pdb=pdb_path,
    ligand_smiles=top_smiles,
    center_x=center_x,
    center_y=center_y,
    center_z=center_z,
    box_size=box_size,
    n_poses=9,
    cnn_scoring="rescore",
    exhaustiveness=16,              # high-accuracy mode
)
→ {poses: list, cnn_affinity: float, vina_score: float, runtime_s: float}
```
Sort by `cnn_affinity` descending (higher = more likely active). Take top 20.

---

### Step 8 — Binding affinity estimation (top 20)

```
boltz2_affinity(
    protein_fasta=target_sequence,
    ligand_smiles=top20_smiles,
)
→ {predictions: [{smiles, affinity_kcal_mol, iptm_score}]}
```
Accept predictions where `iptm_score > 0.55` (interface confidence).
Sort by `affinity_kcal_mol` ascending (most negative = tightest). Take top 10.

---

### Step 9 — ADMET profiling (top 10)

```
admetlab3(smiles=top10_smiles)
```
Flag and deprioritize (do not eliminate) compounds with:
- hERG risk ≥ "Medium"
- Caco-2 permeability < 5 nm/s
- Aqueous solubility < 10 μg/mL
- Pgp substrate = "Yes" combined with low permeability

Compounds with ≥ 2 serious ADMET flags are moved to a secondary list.

---

### Step 10 — Visualisation and report

```
ligand_viz(
    smiles=top10_smiles,
    properties={
        "docking_score_kcal_mol": docking_scores,
        "gnina_cnn_affinity":     cnn_affinities,
        "boltz2_dg_kcal_mol":     boltz2_dgs,
        "qed":                    qed_values,
        "admet_flags":            flag_counts,
        "catalog_id":             catalog_ids,
        "price_usd":              prices,
    },
    plots=["grid", "table", "scatter"],
    title="VS hits — <target name>",
    output_dir=workdir,
)
→ {html_path: str}
```

---

## Decision branches

```
Target is sequence only, ligand in library:
  → protenix (default — Apache 2.0, no license friction)
  → chai1 if glycoprotein target or restraints available
  → esmfold if protein-only and speed matters (< 2 s)

mean_plddt < 50 or ranking_score < 0.4:
  → Halt. Request experimental PDB or AlphaFold DB entry.

No co-crystal ligand + pocket_detect returns 0 pockets:
  → Re-run pocket_detect with tool="fpocket"
  → If still empty: ask user for binding site coordinates

No reference SMILES for library search:
  → chembl(target_id, standard_type="IC50", limit=50) for seeds
  → Fallback: co-crystal ligand from pdb_fetch as seed

Compound count after Step 5:
  ≤ 200     → Route A (GNINA)
  201–5000  → Route B (Vina)
  > 5000    → Route C (AutoDock-GPU); requires PDBQT prep — warn user

< 5 compounds pass Steps 4–5:
  → Loosen QED to 0.35, similarity_threshold to 0.35
  → Consider zinc22 if using enamine_real, or vice-versa

autodock_gpu unavailable or no GPU:
  → Fall back to vina (Route B) regardless of library size

boltz2_affinity fails or times out:
  → deeppurpose(smiles=top20_smiles, target_fasta=target_sequence) as fallback
  → If also unavailable: rank by gnina cnn_affinity only

≥ 8 / 10 final compounds have hERG risk ≥ Medium:
  → Filter to hERG-negative, report ≤ 5 hits
  → Note liability; recommend patch-clamp counter-screen before purchase
```

---

## Expected outputs

Ranked table of ≤ 10 buyable hits:
```
Rank | SMILES | Catalog ID | Price (USD) | Vina (kcal/mol) | GNINA CNN | Boltz-2 ΔG | QED | SA | ADMET flags
```
Plus `ligand_viz` HTML report at `workdir/vs_report.html`.

---

## Success criteria

- ≥ 5 final candidates with GNINA `cnn_affinity > 0.7` (0–1 scale)
- All final candidates pass Lipinski and PAINS/BRENK filters
- ≥ 1 candidate available for immediate purchase (Enamine/ZINC in-stock)
- Boltz-2 `iptm_score > 0.55` for top-3 (if Step 8 was run)
- Run time: Steps 1–7 complete in < 2 h on a single A100 GPU for ≤ 5 000 compounds

---

## Tool summary

| Step | Tool | Scale | Input |
|---|---|---|---|
| 1a | `pdb_fetch` | any | PDB ID |
| 1b | `protenix` / `esmfold` | any | sequence |
| 2 | `pocket_detect` | any | PDB path |
| 3 | `library_search` / `chembl` | library | target ID / ref SMILES |
| 4 | `rdkit_props` | batch | SMILES list |
| 5 | `ligand_filter` | batch | SMILES list |
| 6A | `gnina` | ≤ 200 | PDB + SMILES |
| 6B | `vina` | 201–5 000 | PDBQT + SMILES |
| 6C | `autodock_gpu` | > 5 000 | PDBQT dir + GPF |
| 7 | `gnina` | top 100 | PDB + SMILES |
| 8 | `boltz2_affinity` | top 20 | FASTA + SMILES |
| 9 | `admetlab3` | top 10 | SMILES |
| 10 | `ligand_viz` | top 10 | SMILES + properties |

---

## References

- AutoDock-GPU: Santos-Martins et al., *J. Chem. Theory Comput.* 2021 — `autodock_gpu`
- GNINA: McNutt et al., *J. Cheminform.* 2021 — `gnina`
- Vina-GPU 2.0: Ding et al., *J. Chem. Inf. Model.* 2023 — `vina`
- Boltz-2: Wohlwend et al., *bioRxiv* 2025 — `boltz2_affinity`
- Protenix v2: Zhang et al., *bioRxiv* 2026 — `protenix`
- P2Rank: Krivák & Hoksza, *J. Cheminform.* 2018 — `pocket_detect`
- ADMETlab 3.0: Xiong et al., *Nucleic Acids Res.* 2024 — `admetlab3`
- Enamine REAL: https://enamine.net/compound-libraries/real-compounds
- ZINC22: Tingle et al., *J. Chem. Inf. Model.* 2023 — `library_search`
