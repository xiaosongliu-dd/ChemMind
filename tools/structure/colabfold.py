from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path


def run_colabfold(
    sequences: str | list[str],
    n_models: int = 5,
    use_templates: bool = True,
    use_amber: bool = False,
    colabfold_bin: str = "colabfold_batch",
    output_dir: str | None = None,
) -> dict:
    """
    Predict protein structure(s) via ColabFold (AlphaFold2 + MMseqs2 MSA).

    `sequences` can be a single FASTA string or a list of FASTA strings.
    Multiple sequences in a single entry are treated as a multimer (joined by ':').

    Returns:
        {pdb_paths, best_pdb_path, mean_plddt, ptm_scores, n_sequences}

    Raises RuntimeError if colabfold_batch returns non-zero exit code.
    """
    if isinstance(sequences, str):
        sequences = [sequences]

    outdir = Path(output_dir) if output_dir else Path(tempfile.mkdtemp(prefix="colabfold_"))
    outdir.mkdir(parents=True, exist_ok=True)

    fasta_path = _write_fasta(sequences, outdir)
    _run_colabfold_batch(fasta_path, outdir, n_models, use_templates, use_amber, colabfold_bin)

    pdb_paths, ptm_scores = _collect_results(outdir)
    best_pdb = pdb_paths[0] if pdb_paths else None
    mean_plddt = _mean_plddt_from_pdb(best_pdb) if best_pdb else None

    return {
        "pdb_paths":    pdb_paths,
        "best_pdb_path": best_pdb,
        "mean_plddt":   mean_plddt,
        "ptm_scores":   ptm_scores,
        "n_sequences":  len(sequences),
        "output_dir":   str(outdir),
    }


def _write_fasta(sequences: list[str], outdir: Path) -> Path:
    fasta_path = outdir / "input.fasta"
    lines = []
    for i, seq in enumerate(sequences):
        seq = seq.strip()
        if seq.startswith(">"):
            lines.append(seq)
        else:
            lines.append(f">seq_{i}")
            lines.append(seq)
    fasta_path.write_text("\n".join(lines) + "\n")
    return fasta_path


def _run_colabfold_batch(
    fasta_path: Path,
    outdir: Path,
    n_models: int,
    use_templates: bool,
    use_amber: bool,
    colabfold_bin: str,
) -> None:
    cmd = [
        colabfold_bin,
        str(fasta_path),
        str(outdir),
        "--num-models", str(n_models),
    ]
    if not use_templates:
        cmd.append("--no-templates")
    if use_amber:
        cmd.append("--amber")

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"colabfold_batch failed:\n{result.stderr[-2000:]}")


def _collect_results(outdir: Path) -> tuple[list[str], list[float]]:
    """Return PDB paths sorted by rank (rank_1 first) and their pTM scores."""
    pdbs = sorted(outdir.glob("*rank_*.pdb"))
    if not pdbs:
        pdbs = sorted(outdir.glob("*.pdb"))

    pdb_paths = [str(p) for p in pdbs]
    ptm_scores = []
    for p in pdbs:
        score_file = p.with_suffix(".json")
        if score_file.exists():
            import json
            try:
                data = json.loads(score_file.read_text())
                ptm_scores.append(data.get("ptm") or data.get("pTM"))
            except Exception:
                ptm_scores.append(None)
        else:
            ptm_scores.append(None)

    return pdb_paths, ptm_scores


def _mean_plddt_from_pdb(pdb_path: str) -> float | None:
    plddts = []
    try:
        for line in Path(pdb_path).read_text().splitlines():
            if line.startswith(("ATOM", "HETATM")):
                try:
                    plddts.append(float(line[60:66]))
                except (ValueError, IndexError):
                    pass
    except OSError:
        return None
    return round(sum(plddts) / len(plddts), 2) if plddts else None
