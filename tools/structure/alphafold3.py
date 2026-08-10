from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path


def run_alphafold3(
    protein_sequences: str | list[str],
    ligand_smiles: str | list[str] | None = None,
    ligand_ccd_codes: list[str] | None = None,
    dna_sequences: list[str] | None = None,
    rna_sequences: list[str] | None = None,
    seeds: list[int] | None = None,
    num_seeds: int = 5,
    model_dir: str | None = None,
    output_dir: str | None = None,
    af3_script: str = "run_alphafold.py",
    timeout: int = 3600,
) -> dict:
    """
    Predict biomolecular structure with AlphaFold 3 (local installation).

    Supports protein-only, protein-ligand co-complex, protein-DNA/RNA, and
    arbitrary multi-chain assemblies in a single call.

    Parameters
    ----------
    protein_sequences : str or list[str]
        One or more amino-acid sequences (plain string, no FASTA header).
        Multiple sequences are treated as separate chains in the complex.
    ligand_smiles : str or list[str], optional
        Small-molecule ligand(s) as SMILES. Each becomes a separate ligand chain.
    ligand_ccd_codes : list[str], optional
        CCD/PDB ligand codes (e.g. ["ATP", "MG"]) — alternative to SMILES.
    dna_sequences : list[str], optional
        DNA strand sequences (5'→3', nucleotide single-letter codes).
    rna_sequences : list[str], optional
        RNA strand sequences.
    seeds : list[int], optional
        Explicit seed list. Overrides num_seeds.
    num_seeds : int
        Number of random seeds to run (default 5 → 5 independent predictions).
    model_dir : str, optional
        Path to downloaded AF3 model weights directory.
        Falls back to $AF3_MODEL_DIR environment variable.
    output_dir : str, optional
        Directory for AF3 output files. Uses a temp directory if not given.
    af3_script : str
        Path to AlphaFold 3 run_alphafold.py entry-point script.
        Falls back to $AF3_SCRIPT environment variable.
    timeout : int
        Subprocess timeout in seconds (default 3600 = 1 h).

    Returns
    -------
    {
        "top_cif_path":     str,          # best-ranked structure (CIF format)
        "all_cif_paths":    list[str],    # all predicted structures, best-first
        "ranking_score":    float,        # AF3 ranking score for best model
        "ptm":              float,        # predicted TM-score (0–1)
        "iptm":             float | None, # interface pTM (complexes only)
        "mean_plddt":       float,        # mean per-residue pLDDT (0–100)
        "has_ligand":       bool,
        "has_nucleic_acid": bool,
        "n_chains":         int,
        "output_dir":       str,
        "seeds_used":       list[int],
    }
    """
    resolved_model_dir = model_dir or os.environ.get("AF3_MODEL_DIR")
    if not resolved_model_dir:
        raise RuntimeError(
            "AlphaFold 3 model directory not specified. Provide model_dir= or "
            "set the AF3_MODEL_DIR environment variable. "
            "Request model weights at https://github.com/google-deepmind/alphafold3"
        )

    resolved_script = af3_script if af3_script != "run_alphafold.py" else \
        os.environ.get("AF3_SCRIPT", "run_alphafold.py")

    if isinstance(protein_sequences, str):
        protein_sequences = [protein_sequences]
    if isinstance(ligand_smiles, str):
        ligand_smiles = [ligand_smiles]

    seeds_used = seeds if seeds is not None else list(range(1, num_seeds + 1))

    outdir = Path(output_dir) if output_dir else Path(tempfile.mkdtemp(prefix="alphafold3_"))
    outdir.mkdir(parents=True, exist_ok=True)

    job_name = "chemmind_af3"
    input_payload = _build_input_json(
        job_name=job_name,
        protein_sequences=protein_sequences,
        ligand_smiles=ligand_smiles or [],
        ligand_ccd_codes=ligand_ccd_codes or [],
        dna_sequences=dna_sequences or [],
        rna_sequences=rna_sequences or [],
        seeds=seeds_used,
    )
    json_path = outdir / f"{job_name}.json"
    json_path.write_text(json.dumps(input_payload, indent=2))

    _run_subprocess(resolved_script, json_path, resolved_model_dir, outdir, timeout)

    top_cif, all_cifs, ranking_score, ptm, iptm, mean_plddt = _parse_results(outdir, job_name)

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
        "mean_plddt":       mean_plddt,
        "has_ligand":       has_ligand,
        "has_nucleic_acid": has_nucleic_acid,
        "n_chains":         n_chains,
        "output_dir":       str(outdir),
        "seeds_used":       seeds_used,
    }


# ── Input construction ────────────────────────────────────────────────────────

def _build_input_json(
    job_name: str,
    protein_sequences: list[str],
    ligand_smiles: list[str],
    ligand_ccd_codes: list[str],
    dna_sequences: list[str],
    rna_sequences: list[str],
    seeds: list[int],
) -> dict:
    """Build the AF3 v2 JSON input format."""
    chain_ids = _chain_id_generator()
    sequences = []

    for seq in protein_sequences:
        sequences.append({
            "protein": {
                "id":           next(chain_ids),
                "sequence":     seq.strip().upper(),
                "modifications": [],
            }
        })

    for smi in ligand_smiles:
        sequences.append({
            "ligand": {
                "id":     next(chain_ids),
                "smiles": smi,
            }
        })

    for ccd in ligand_ccd_codes:
        sequences.append({
            "ligand": {
                "id":       next(chain_ids),
                "ccdCodes": [ccd],
            }
        })

    for dna in dna_sequences:
        sequences.append({
            "dna": {
                "id":       next(chain_ids),
                "sequence": dna.strip().upper(),
            }
        })

    for rna in rna_sequences:
        sequences.append({
            "rna": {
                "id":       next(chain_ids),
                "sequence": rna.strip().upper(),
            }
        })

    return {
        "name":        job_name,
        "modelSeeds":  seeds,
        "sequences":   sequences,
        "dialect":     "alphafold3",
        "version":     2,
    }


def _chain_id_generator():
    """Yield A, B, C … Z, AA, AB … indefinitely."""
    import string
    letters = string.ascii_uppercase
    for c in letters:
        yield c
    for c1 in letters:
        for c2 in letters:
            yield c1 + c2


# ── Subprocess runner ─────────────────────────────────────────────────────────

def _run_subprocess(
    af3_script: str,
    json_path: Path,
    model_dir: str,
    output_dir: Path,
    timeout: int,
) -> None:
    cmd = [
        "python", af3_script,
        f"--json_path={json_path}",
        f"--model_dir={model_dir}",
        f"--output_dir={output_dir}",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    if result.returncode != 0:
        raise RuntimeError(
            f"AlphaFold 3 failed (exit {result.returncode}):\n"
            f"{result.stderr[-3000:]}"
        )


# ── Output parsing ────────────────────────────────────────────────────────────

def _parse_results(
    outdir: Path, job_name: str
) -> tuple[str | None, list[str], float | None, float | None, float | None, float | None]:
    """
    Locate AF3 output files and select the best model by ranking_score.

    AF3 writes files as:
        {outdir}/{job_name}/{job_name}_model_{i}.cif
        {outdir}/{job_name}/{job_name}_summary_confidences_{i}.json
    """
    job_dir = outdir / job_name
    if not job_dir.exists():
        job_dir = outdir  # some AF3 versions write directly to outdir

    cif_files = sorted(job_dir.glob(f"{job_name}_model_*.cif"))
    if not cif_files:
        cif_files = sorted(job_dir.glob("*.cif"))

    if not cif_files:
        return None, [], None, None, None, None

    # Load per-seed confidence summaries and rank by ranking_score
    ranked: list[tuple[float, Path, dict]] = []
    for cif in cif_files:
        idx = _index_from_stem(cif.stem)
        conf_path = job_dir / f"{job_name}_summary_confidences_{idx}.json"
        if not conf_path.exists():
            conf_path = cif.with_name(cif.stem.replace("_model_", "_summary_confidences_") + ".json")

        conf: dict = {}
        if conf_path.exists():
            try:
                conf = json.loads(conf_path.read_text())
            except json.JSONDecodeError:
                pass

        score = conf.get("ranking_score", 0.0)
        ranked.append((score, cif, conf))

    ranked.sort(key=lambda x: -x[0])
    best_score, best_cif, best_conf = ranked[0]

    all_cifs = [str(r[1]) for r in ranked]
    ptm       = best_conf.get("ptm")
    iptm      = best_conf.get("iptm")
    plddt     = best_conf.get("mean_plddt")

    return str(best_cif), all_cifs, best_score or None, ptm, iptm, plddt


def _index_from_stem(stem: str) -> str:
    """Extract the numeric index from a filename stem like 'job_model_3'."""
    parts = stem.rsplit("_", 1)
    return parts[-1] if len(parts) == 2 else "0"
