from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path


_DEFAULT_MODEL = "protenix-v2"


def run_protenix(
    protein_sequences: str | list[str],
    ligand_smiles: str | list[str] | None = None,
    ligand_ccd_codes: list[str] | None = None,
    dna_sequences: list[str] | None = None,
    rna_sequences: list[str] | None = None,
    seed: int = 101,
    n_sample: int = 5,
    model: str = _DEFAULT_MODEL,
    output_dir: str | None = None,
    timeout: int = 3600,
) -> dict:
    """
    Predict biomolecular structure with Protenix (ByteDance AF3 reimplementation).

    Protenix is fully open-source (Apache 2.0) with no restrictions on commercial use —
    unlike AF3 and Chai-1 whose weights have non-commercial licenses. Protenix-v2
    outperforms AF3 on antibody-antigen benchmarks and is pip-installable.

    Install: pip install protenix

    Parameters
    ----------
    protein_sequences : str or list[str]
        Amino-acid sequence(s). Each becomes a separate protein chain.
    ligand_smiles : str or list[str], optional
        Small-molecule ligand(s) as SMILES strings.
    ligand_ccd_codes : list[str], optional
        PDB CCD codes (e.g. ["ATP", "MG"]). Passed as "CCD_<code>" to Protenix.
    dna_sequences : list[str], optional
        Single-stranded DNA sequences (A/T/G/C/N).
    rna_sequences : list[str], optional
        RNA sequences (A/U/G/C/N).
    seed : int
        Random seed (default 101). Controls diffusion stochasticity.
    n_sample : int
        Number of structural samples to generate per seed (default 5).
    model : str
        Protenix model name. Options: "protenix-v2" (default, best accuracy),
        "protenix_base_default_v1.0.0", "protenix_base_20250630_v1.0.0".
    output_dir : str, optional
        Directory for Protenix output. Auto-generated temp dir if not given.
    timeout : int
        Subprocess timeout in seconds (default 3600).

    Returns
    -------
    {
        "top_cif_path":     str,          # best-ranked structure (CIF)
        "all_cif_paths":    list[str],    # all predictions, best first
        "ranking_score":    float | None, # AF3-style ranking score
        "ptm":              float | None,
        "iptm":             float | None,
        "mean_plddt":       float | None,
        "has_ligand":       bool,
        "has_nucleic_acid": bool,
        "n_chains":         int,
        "output_dir":       str,
        "seed":             int,
        "n_sample":         int,
        "model":            str,
    }
    """
    if isinstance(protein_sequences, str):
        protein_sequences = [protein_sequences]
    if isinstance(ligand_smiles, str):
        ligand_smiles = [ligand_smiles]

    outdir = Path(output_dir) if output_dir else Path(tempfile.mkdtemp(prefix="protenix_"))
    outdir.mkdir(parents=True, exist_ok=True)

    job_name = "chemmind_ptx"
    payload = _build_input_json(
        job_name=job_name,
        protein_sequences=protein_sequences,
        ligand_smiles=ligand_smiles or [],
        ligand_ccd_codes=ligand_ccd_codes or [],
        dna_sequences=dna_sequences or [],
        rna_sequences=rna_sequences or [],
    )
    json_path = outdir / f"{job_name}.json"
    json_path.write_text(json.dumps(payload, indent=2))

    _run_subprocess(json_path, outdir, seed, n_sample, model, timeout)

    top_cif, all_cifs, ranking_score, ptm, iptm, plddt_mean = _parse_results(outdir)

    has_ligand       = bool(ligand_smiles or ligand_ccd_codes)
    has_nucleic_acid = bool(dna_sequences or rna_sequences)
    n_chains = (
        len(protein_sequences)
        + len(ligand_smiles or [])
        + len(ligand_ccd_codes or [])
        + len(dna_sequences or [])
        + len(rna_sequences or [])
    )

    return {
        "top_cif_path":     top_cif,
        "all_cif_paths":    all_cifs,
        "ranking_score":    ranking_score,
        "ptm":              ptm,
        "iptm":             iptm,
        "mean_plddt":       plddt_mean,
        "has_ligand":       has_ligand,
        "has_nucleic_acid": has_nucleic_acid,
        "n_chains":         n_chains,
        "output_dir":       str(outdir),
        "seed":             seed,
        "n_sample":         n_sample,
        "model":            model,
    }


# ── Input construction ────────────────────────────────────────────────────────

def _build_input_json(
    job_name: str,
    protein_sequences: list[str],
    ligand_smiles: list[str],
    ligand_ccd_codes: list[str],
    dna_sequences: list[str],
    rna_sequences: list[str],
) -> list[dict]:
    """Build Protenix input JSON (top-level is a list, even for single jobs)."""
    sequences = []

    for seq in protein_sequences:
        sequences.append({
            "proteinChain": {
                "sequence":      seq.strip().upper(),
                "count":         1,
                "modifications": [],
            }
        })

    for smi in ligand_smiles:
        # Protenix accepts raw SMILES directly in the "ligand" field
        sequences.append({"ligand": {"ligand": smi, "count": 1}})

    for ccd in ligand_ccd_codes:
        # CCD codes must be prefixed with "CCD_"
        sequences.append({"ligand": {"ligand": f"CCD_{ccd}", "count": 1}})

    for dna in dna_sequences:
        sequences.append({
            "dnaSequence": {
                "sequence": dna.strip().upper(),
                "count":    1,
            }
        })

    for rna in rna_sequences:
        sequences.append({
            "rnaSequence": {
                "sequence": rna.strip().upper(),
                "count":    1,
            }
        })

    return [{"name": job_name, "sequences": sequences, "covalent_bonds": []}]


# ── Subprocess runner ─────────────────────────────────────────────────────────

def _run_subprocess(
    json_path: Path,
    output_dir: Path,
    seed: int,
    n_sample: int,
    model: str,
    timeout: int,
) -> None:
    cmd = [
        "protenix", "pred",
        "-i", str(json_path),
        "-o", str(output_dir),
        "-s", str(seed),
        "-n", model,
        f"sample_diffusion.N_sample={n_sample}",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    if result.returncode != 0:
        raise RuntimeError(
            f"Protenix failed (exit {result.returncode}):\n"
            f"{result.stderr[-3000:]}"
        )


# ── Output parsing ────────────────────────────────────────────────────────────

def _parse_results(
    outdir: Path,
) -> tuple[str | None, list[str], float | None, float | None, float | None, float | None]:
    """
    Locate CIF and confidence JSON files produced by Protenix.

    Protenix (AF3-style) writes per-sample directories:
        {outdir}/{job_name}/seed-{seed}_sample-{i}/model.cif
        {outdir}/{job_name}/seed-{seed}_sample-{i}/summary_confidences.json
    or flat CIFs directly in outdir for single-sample runs.
    """
    cif_files = sorted(outdir.rglob("*.cif"))
    if not cif_files:
        return None, [], None, None, None, None

    ranked: list[tuple[float, Path, dict]] = []
    for cif in cif_files:
        conf = _load_confidence(cif)
        score = conf.get("ranking_score", 0.0)
        ranked.append((score, cif, conf))

    ranked.sort(key=lambda x: -x[0])
    best_score, best_cif, best_conf = ranked[0]

    return (
        str(best_cif),
        [str(r[1]) for r in ranked],
        best_score or None,
        best_conf.get("ptm"),
        best_conf.get("iptm"),
        best_conf.get("mean_plddt"),
    )


def _load_confidence(cif_path: Path) -> dict:
    """Look for a summary_confidences.json adjacent to a CIF file."""
    candidates = [
        cif_path.parent / "summary_confidences.json",
        cif_path.with_suffix(".json"),
        cif_path.with_name(cif_path.stem + "_summary_confidences.json"),
    ]
    for p in candidates:
        if p.exists():
            try:
                return json.loads(p.read_text())
            except json.JSONDecodeError:
                pass
    return {}
