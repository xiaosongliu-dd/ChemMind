from __future__ import annotations

import os
import tempfile
import time
from pathlib import Path


def run_esmfold(
    sequence: str | list[str],
    use_nim: bool = True,
    nim_api_key: str | None = None,
    output_dir: str | None = None,
) -> dict:
    """
    Predict protein 3D structure from sequence with ESMFold.

    Single mode (`sequence: str`):
        {pdb_path, mean_plddt, ptm_score, sequence, runtime_s}

    Batch mode (`sequence: list[str]`):
        {results: [single-mode dict, ...], n_sequences}

    NIM mode hits the NVIDIA ESMFold endpoint; local mode runs ESMFold
    via the `esm` Python package (requires GPU and the esm package installed).
    """
    if isinstance(sequence, str):
        return _fold_one(sequence, use_nim, nim_api_key, output_dir)

    results = []
    for seq in sequence:
        try:
            r = _fold_one(seq, use_nim, nim_api_key, output_dir)
        except Exception as e:
            r = {"sequence": seq, "error": str(e), "pdb_path": None, "mean_plddt": None}
        results.append(r)
    return {"results": results, "n_sequences": len(sequence)}


def _fold_one(
    sequence: str,
    use_nim: bool,
    nim_api_key: str | None,
    output_dir: str | None,
) -> dict:
    if use_nim:
        return _run_nim(sequence, nim_api_key, output_dir)
    return _run_local(sequence, output_dir)


def _run_nim(sequence: str, api_key: str | None, output_dir: str | None) -> dict:
    import requests

    key = api_key or os.environ.get("NVIDIA_API_KEY", "")
    NIM_URL = "https://health.api.nvidia.com/v1/biology/nvidia/esmfold"

    t0 = time.perf_counter()
    try:
        resp = requests.post(
            NIM_URL,
            headers={"Authorization": f"Bearer {key}", "Content-Type": "application/json"},
            json={"sequence": sequence},
            timeout=120,
        )
        resp.raise_for_status()
    except requests.RequestException as e:
        raise RuntimeError(f"ESMFold NIM API unreachable: {e}")

    data = resp.json()
    runtime = round(time.perf_counter() - t0, 1)

    pdb_str = data.get("pdbs", [data.get("pdb", "")])[0]
    outdir = Path(output_dir) if output_dir else Path(tempfile.mkdtemp(prefix="esmfold_"))
    outdir.mkdir(parents=True, exist_ok=True)
    pdb_path = outdir / "esmfold.pdb"
    pdb_path.write_text(pdb_str)

    mean_plddt = _extract_mean_plddt(pdb_str)

    return {
        "pdb_path":    str(pdb_path),
        "mean_plddt":  mean_plddt,
        "ptm_score":   data.get("ptm"),
        "sequence":    sequence,
        "runtime_s":   runtime,
    }


def _run_local(sequence: str, output_dir: str | None) -> dict:
    try:
        import torch
        import esm
    except ImportError as e:
        raise RuntimeError(
            f"esm package not installed (required for local ESMFold): {e}. "
            "Install with: pip install fair-esm"
        )

    t0 = time.perf_counter()
    model = esm.pretrained.esmfold_v1()
    model = model.eval()
    if torch.cuda.is_available():
        model = model.cuda()

    with torch.no_grad():
        output = model.infer_pdb(sequence)

    runtime = round(time.perf_counter() - t0, 1)

    outdir = Path(output_dir) if output_dir else Path(tempfile.mkdtemp(prefix="esmfold_"))
    outdir.mkdir(parents=True, exist_ok=True)
    pdb_path = outdir / "esmfold.pdb"
    pdb_path.write_text(output)

    mean_plddt = _extract_mean_plddt(output)

    return {
        "pdb_path":   str(pdb_path),
        "mean_plddt": mean_plddt,
        "ptm_score":  None,
        "sequence":   sequence,
        "runtime_s":  runtime,
    }


def _extract_mean_plddt(pdb_str: str) -> float | None:
    """Parse mean pLDDT from ATOM B-factor column."""
    plddts = []
    for line in pdb_str.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            try:
                plddts.append(float(line[60:66]))
            except (ValueError, IndexError):
                pass
    return round(sum(plddts) / len(plddts), 2) if plddts else None
