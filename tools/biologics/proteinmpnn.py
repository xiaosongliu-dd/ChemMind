from __future__ import annotations

import base64
import os
import subprocess
import tempfile
from pathlib import Path

NIM_URL = "https://health.api.nvidia.com/v1/biology/ipd/proteinmpnn"


def run_proteinmpnn(
    pdb_path: str,
    chains_to_design: list[str],
    fixed_positions: dict[str, list[int]] | None = None,
    n_sequences: int = 8,
    sampling_temp: float = 0.1,
    use_soluble_model: bool = False,
    use_nim: bool = True,
    nim_api_key: str | None = None,
) -> dict:
    """
    Design amino acid sequences on a fixed protein backbone using ProteinMPNN.

    Parameters
    ----------
    pdb_path          : receptor / scaffold PDB
    chains_to_design  : chain IDs to design, e.g. ["A"] or ["H", "L"]
    fixed_positions   : {chain_id: [residue_indices]} for residues that must
                        keep their original sequence (e.g. active site residues)
    n_sequences       : number of designed sequences per chain
    sampling_temp     : sequence diversity temperature; 0.1 = low diversity
                        (high confidence), 1.0 = high diversity
    use_soluble_model : use soluble-only ProteinMPNN variant (IPD)
    """
    if use_nim:
        return _run_nim(
            pdb_path, chains_to_design, fixed_positions,
            n_sequences, sampling_temp, use_soluble_model,
            nim_api_key or os.environ.get("NVIDIA_API_KEY", ""),
        )
    return _run_local(pdb_path, chains_to_design, fixed_positions, n_sequences, sampling_temp)


def _run_nim(pdb_path, chains_to_design, fixed_positions, n_sequences, temp, soluble, api_key):
    import requests

    pdb_str = Path(pdb_path).read_text()
    payload = {
        "pdb_str":           pdb_str,
        "chains_to_design":  ",".join(chains_to_design),
        "num_seq_per_target": n_sequences,
        "sampling_temp":     str(temp),
        "use_soluble_model": soluble,
    }
    if fixed_positions:
        payload["fixed_positions_jsonl"] = _encode_fixed_positions(fixed_positions)

    resp = requests.post(
        NIM_URL,
        headers={"Authorization": f"Bearer {api_key}", "Content-Type": "application/json"},
        json=payload,
        timeout=300,
    )
    resp.raise_for_status()
    data = resp.json()

    sequences = data.get("sequences", [])
    return {
        "sequences":     [s["sequence"] for s in sequences],
        "scores":        [s.get("global_score") for s in sequences],
        "n_returned":    len(sequences),
        "chains_designed": chains_to_design,
        "sampling_temp": temp,
    }


def _run_local(pdb_path, chains_to_design, fixed_positions, n_sequences, temp):
    workdir = Path(tempfile.mkdtemp(prefix="mpnn_"))
    out_dir = workdir / "outputs"
    out_dir.mkdir()

    chain_str = " ".join(chains_to_design)
    cmd = [
        "python", "protein_mpnn_run.py",
        "--pdb_path",            pdb_path,
        "--pdb_path_chains",     chain_str,
        "--out_folder",          str(out_dir),
        "--num_seq_per_target",  str(n_sequences),
        "--sampling_temp",       str(temp),
        "--seed",                "37",
        "--batch_size",          "1",
    ]
    if fixed_positions:
        fixed_path = workdir / "fixed.json"
        import json
        fixed_path.write_text(json.dumps(fixed_positions))
        cmd += ["--fixed_positions_jsonl", str(fixed_path)]

    subprocess.run(cmd, check=True, capture_output=True, cwd="ProteinMPNN")

    seqs, scores = _parse_mpnn_fasta(out_dir)
    return {
        "sequences":      seqs,
        "scores":         scores,
        "n_returned":     len(seqs),
        "chains_designed": chains_to_design,
        "sampling_temp":  temp,
    }


def _parse_mpnn_fasta(out_dir: Path) -> tuple[list[str], list[float | None]]:
    seqs, scores = [], []
    for fasta in out_dir.glob("*.fa"):
        current_score = None
        for line in fasta.read_text().splitlines():
            if line.startswith(">"):
                try:
                    current_score = float(line.split("score=")[1].split(",")[0])
                except (IndexError, ValueError):
                    pass
            elif line.strip():
                seqs.append(line.strip())
                scores.append(current_score)
    return seqs, scores


def _encode_fixed_positions(fixed: dict) -> str:
    import json
    return json.dumps(fixed)
