# Skill: rfantibody

**Tier:** L1 — atomic tool skill
**Category:** biologics
**Tool:** RFantibody — RFdiffusion fine-tuned for antibody design
**License:** MIT (Baker Lab, University of Washington)
**GPU required:** Yes — A100 recommended; ~2–10 min per design batch

---

## What this skill does

RFantibody is a fine-tuned variant of RFdiffusion trained specifically on antibody loops. It performs **de novo design of human-like antibodies** (scFvs and VHH nanobodies) against a user-specified antigen epitope — no starting antibody structure needed.

Key capabilities:
- De novo scFv design (heavy + light chain variable fragments)
- VHH / nanobody design
- CDR-H3 loop generation (the primary binding determinant)
- Epitope-targeted design (specify which residues on the antigen to engage)
- Outputs backbone structures → feed to ProteinMPNN for sequence design

---

## When to call this skill

- Need a brand-new antibody against a target antigen
- No existing antibody scaffold available — truly de novo
- Want to target a specific epitope or interface region
- Generating a diverse library of CDR designs for experimental screening

Do NOT use if:
- You have an existing antibody and only need CDR optimisation → use `abdiffuser` or `proteinmpnn`
- You need full-atom structure → feed RFantibody backbone into `igfold` or `abodybuilder3`

---

## Inputs

| Field | Type | Required | Description |
|---|---|---|---|
| `antigen_pdb` | str | ✓ | Path to antigen PDB file |
| `epitope_residues` | list[int] | ✓ | Residue numbers on antigen to target |
| `antibody_type` | str | — | `"scfv"` (default) or `"nanobody"` |
| `n_designs` | int | — | Number of backbones to generate (default 50) |
| `partial_diffusion` | bool | — | Start from existing framework (default False) |
| `framework_pdb` | str | — | Required if `partial_diffusion=True` |
| `hotspot_weights` | dict | — | Per-residue epitope weights (advanced) |
| `output_dir` | str | — | Where to write PDB outputs |

---

## Outputs

| Field | Type | Description |
|---|---|---|
| `backbone_pdbs` | list[str] | Paths to designed antibody backbone PDBs |
| `n_generated` | int | Number of successful designs |
| `epitope_contacts` | list[dict] | Per-design contact residue lists |
| `runtime_s` | float | Wall time |

---

## Execution (Python wrapper)

```python
# dependencies = ["rfdiffusion>=1.1", "torch>=2.0"]

import subprocess, time, tempfile, json
from pathlib import Path

def run_rfantibody(
    antigen_pdb: str,
    epitope_residues: list[int],
    antibody_type: str = "scfv",
    n_designs: int = 50,
    partial_diffusion: bool = False,
    framework_pdb: str | None = None,
    output_dir: str | None = None,
) -> dict:

    t0 = time.perf_counter()
    outdir = Path(output_dir or tempfile.mkdtemp(prefix="rfantibody_"))
    outdir.mkdir(parents=True, exist_ok=True)

    # Build hotspot string from epitope residues, e.g. "A10,A15,A22"
    hotspot_str = ",".join(f"A{r}" for r in epitope_residues)

    # RFantibody is invoked via the rfdiffusion CLI
    # see: github.com/RosettaCommons/RFdiffusion/tree/main/examples/antibody
    cmd = [
        "python", "-m", "rfdiffusion.run_inference",
        f"inference.input_pdb={antigen_pdb}",
        f"inference.num_designs={n_designs}",
        f"inference.output_prefix={outdir}/design",
        f"potentials.guiding_potentials=['type:antibody_contacts,weight:1.0,epitope:{hotspot_str}']",
        "contigmap.contigs=[A1-200]",
        f"antibody.antibody_type={antibody_type}",
    ]

    if partial_diffusion and framework_pdb:
        cmd += [
            f"diffuser.partial_T=20",
            f"inference.input_pdb={framework_pdb}",
        ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"RFantibody failed:\n{result.stderr}")

    pdbs = sorted(outdir.glob("design_*.pdb"))
    elapsed = round(time.perf_counter() - t0, 1)

    # Parse per-design contact summary if written by CLI
    contacts_file = outdir / "contacts.json"
    contacts = json.loads(contacts_file.read_text()) if contacts_file.exists() else []

    return {
        "backbone_pdbs":    [str(p) for p in pdbs],
        "n_generated":      len(pdbs),
        "epitope_contacts": contacts,
        "runtime_s":        elapsed,
        "output_dir":       str(outdir),
    }
```

---

## Standard antibody design pipeline (recommended L2 sequence)

```
antigen PDB
    │
    ▼
rfantibody          ← this skill — generates backbone geometries
    │  50 backbones
    ▼
proteinmpnn         ← L1/biologics/proteinmpnn — designs sequences on backbones
    │  8 seqs / backbone = 400 candidates
    ▼
igfold / abodybuilder3   ← L1/biologics/igfold — fold + validate CDR-H3
    │
    ▼
boltz2 (Ab–Ag mode) ← L1/structure/boltz2 — score complex + affinity
    │
    ▼
admetlab3           ← L1/admet/admetlab3 — developability filters
```

This pipeline is encoded as the `antibody_design` L2 skill.

---

## Agent decision rules

- **Design count**: start with `n_designs=50`; ~5–10% will pass downstream filters
- **Epitope targeting**: always specify `epitope_residues` — blind design produces poor hit rates
- **scFv vs nanobody**: use `nanobody` for membrane/cryptic epitopes; `scfv` for most targets
- **Filtering**: after ProteinMPNN, filter by `boltz2.iptm_score > 0.6` before wet-lab
- **Diversity**: if designs cluster, increase `n_designs` or run with different random seeds

---

## Install

```bash
# Clone and install RFdiffusion (includes RFantibody fine-tuned weights)
git clone https://github.com/RosettaCommons/RFdiffusion.git
cd RFdiffusion
pip install -e .

# Download antibody-specific model weights (~2 GB)
bash scripts/download_models.sh RFantibody
```

---

## References

- Bennett et al. 2025 — "De novo design of human-like antibodies with RFdiffusion" (Baker Lab)
- GitHub: https://github.com/RosettaCommons/RFdiffusion
- License: MIT
- Chai-2 (June 2025) achieves ~50% hit rate using similar diffusion-based ab design approach
