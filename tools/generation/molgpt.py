from __future__ import annotations

import tempfile
from pathlib import Path


def run_molgpt(
    n_molecules: int = 100,
    scaffold_smiles: str | None = None,
    temperature: float = 1.0,
    target_qed: float | None = None,
    target_logp: float | None = None,
    target_sa: float | None = None,
    molgpt_dir: str = "MolGPT",
    output_dir: str | None = None,
) -> dict:
    """
    Generate molecules with a GPT language model conditioned on a scaffold
    SMILES fragment and / or property targets.

    Property conditioning (target_qed / target_logp / target_sa) discretises
    the target into a vocabulary token that is prepended to the generation
    prompt, biasing sampling toward the desired property range.

    Returns:
        {
            smiles:        list[str],   # unique valid SMILES
            n_generated:   int,         # total sampled (before dedup/validity)
            n_valid:       int,
            n_unique:      int,
            validity_rate: float,
            uniqueness:    float,
            scaffold:      str | None,
        }

    Raises RuntimeError("MolGPT failed: ...") on import / subprocess failure.
    """
    outdir = Path(output_dir) if output_dir else Path(tempfile.mkdtemp(prefix="molgpt_"))
    outdir.mkdir(parents=True, exist_ok=True)

    smiles = _generate(
        molgpt_dir=molgpt_dir,
        n_molecules=n_molecules,
        scaffold=scaffold_smiles,
        temperature=temperature,
        target_qed=target_qed,
        target_logp=target_logp,
        target_sa=target_sa,
        outdir=outdir,
    )

    valid = _validate_smiles(smiles)
    unique = list(dict.fromkeys(valid))

    out_path = outdir / "generated.smi"
    out_path.write_text("\n".join(unique))

    return {
        "smiles":        unique,
        "n_generated":   n_molecules,
        "n_valid":       len(valid),
        "n_unique":      len(unique),
        "validity_rate": round(len(valid) / max(n_molecules, 1), 3),
        "uniqueness":    round(len(unique) / max(len(valid), 1), 3),
        "scaffold":      scaffold_smiles,
        "output_file":   str(out_path),
    }


def _generate(
    molgpt_dir: str,
    n_molecules: int,
    scaffold: str | None,
    temperature: float,
    target_qed: float | None,
    target_logp: float | None,
    target_sa: float | None,
    outdir: Path,
) -> list[str]:
    molgpt_path = Path(molgpt_dir)
    sample_script = molgpt_path / "sample.py"

    if sample_script.exists():
        return _run_subprocess(
            sample_script, molgpt_path, n_molecules, scaffold,
            temperature, target_qed, target_logp, target_sa, outdir,
        )

    return _run_library(n_molecules, scaffold, temperature, molgpt_path)


def _run_subprocess(
    script: Path,
    molgpt_dir: Path,
    n_molecules: int,
    scaffold: str | None,
    temperature: float,
    target_qed: float | None,
    target_logp: float | None,
    target_sa: float | None,
    outdir: Path,
) -> list[str]:
    import subprocess

    cmd = [
        "python", str(script),
        "--num_samples", str(n_molecules),
        "--temperature", str(temperature),
        "--output_dir",  str(outdir),
    ]
    if scaffold:
        cmd += ["--scaffold", scaffold]
    if target_qed is not None:
        cmd += ["--qed", str(target_qed)]
    if target_logp is not None:
        cmd += ["--logp", str(target_logp)]
    if target_sa is not None:
        cmd += ["--sa", str(target_sa)]

    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(molgpt_dir))
    if result.returncode != 0:
        raise RuntimeError(f"MolGPT failed:\n{result.stderr[-2000:]}")

    smi_files = sorted(outdir.glob("*.smi"))
    if smi_files:
        return [l.split()[0] for l in smi_files[0].read_text().splitlines() if l.strip()]
    return []


def _run_library(
    n_molecules: int,
    scaffold: str | None,
    temperature: float,
    molgpt_dir: Path,
) -> list[str]:
    """Attempt to load MolGPT as a Python module (molgpt_dir must be on sys.path)."""
    import sys
    if str(molgpt_dir) not in sys.path:
        sys.path.insert(0, str(molgpt_dir))

    try:
        from model import GPT, GPTConfig  # type: ignore
        from utils import SmilesDataset  # type: ignore
    except ImportError as e:
        raise RuntimeError(
            f"MolGPT not installed and sample.py not found at {molgpt_dir}: {e}. "
            "Clone from https://github.com/devalab/molgpt."
        )

    raise NotImplementedError(
        "MolGPT library-mode generation requires model weights and dataset config. "
        "Use the subprocess mode by pointing molgpt_dir to a cloned MolGPT repo."
    )


def _validate_smiles(smiles: list[str]) -> list[str]:
    try:
        from rdkit import Chem
        return [s for s in smiles if Chem.MolFromSmiles(s) is not None]
    except ImportError:
        return smiles
