from __future__ import annotations

import csv
import subprocess
import tempfile
from pathlib import Path


def run_pocket_detect(
    protein_pdb: str,
    tool: str = "p2rank",
    p2rank_bin: str = "prank",
    fpocket_bin: str = "fpocket",
    output_dir: str | None = None,
) -> dict:
    """
    Detect and rank druggable binding pockets on a target protein structure.

    `tool`: "p2rank" (default, ML-based) or "fpocket" (geometry-based).

    Returns:
        {
            pockets: list[{rank, score, center_x, center_y, center_z, volume, residues}],
            n_pockets: int,
            best_pocket: dict | None,
            tool: str,
        }

    Raises RuntimeError if the external binary returns non-zero exit code.
    """
    tool = tool.lower()
    if tool not in ("p2rank", "fpocket"):
        raise ValueError(f"Unknown tool: {tool!r}. Choose 'p2rank' or 'fpocket'.")

    outdir = Path(output_dir) if output_dir else Path(tempfile.mkdtemp(prefix="pocket_"))
    outdir.mkdir(parents=True, exist_ok=True)

    if tool == "p2rank":
        pockets = _run_p2rank(protein_pdb, p2rank_bin, outdir)
    else:
        pockets = _run_fpocket(protein_pdb, fpocket_bin, outdir)

    return {
        "pockets":     pockets,
        "n_pockets":   len(pockets),
        "best_pocket": pockets[0] if pockets else None,
        "tool":        tool,
        "output_dir":  str(outdir),
    }


def _run_p2rank(protein_pdb: str, p2rank_bin: str, outdir: Path) -> list[dict]:
    cmd = [p2rank_bin, "predict", "-f", protein_pdb, "-o", str(outdir)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"P2Rank failed:\n{result.stderr[-2000:]}")

    pred_files = list(outdir.rglob("*_predictions.csv"))
    if not pred_files:
        return []

    pockets = []
    with open(pred_files[0]) as f:
        reader = csv.DictReader(f, skipinitialspace=True)
        for rank, row in enumerate(reader, start=1):
            try:
                pocket = {
                    "rank":     rank,
                    "score":    _safe_float(row.get("score") or row.get("probability")),
                    "center_x": _safe_float(row.get("center_x")),
                    "center_y": _safe_float(row.get("center_y")),
                    "center_z": _safe_float(row.get("center_z")),
                    "volume":   _safe_float(row.get("surf_atom_ids") or row.get("volume")),
                    "residues": _parse_residues(row.get("residue_ids", "")),
                }
                pockets.append(pocket)
            except Exception:
                continue

    pockets.sort(key=lambda p: -(p["score"] or 0))
    for i, p in enumerate(pockets, start=1):
        p["rank"] = i
    return pockets


def _run_fpocket(protein_pdb: str, fpocket_bin: str, outdir: Path) -> list[dict]:
    import shutil
    pdb_path = Path(protein_pdb)
    dest = outdir / pdb_path.name
    if not dest.exists():
        shutil.copy(str(pdb_path), str(dest))

    cmd = [fpocket_bin, "-f", str(dest)]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(outdir))
    if result.returncode != 0:
        raise RuntimeError(f"fpocket failed:\n{result.stderr[-2000:]}")

    stem = pdb_path.stem
    info_files = list(outdir.rglob(f"{stem}_out/{stem}_info.txt"))
    if not info_files:
        return []

    return _parse_fpocket_info(info_files[0])


def _parse_fpocket_info(info_file: Path) -> list[dict]:
    text = info_file.read_text()
    pockets = []
    current: dict | None = None

    for line in text.splitlines():
        line = line.strip()
        if line.startswith("Pocket"):
            if current:
                pockets.append(current)
            rank = int(line.split()[1]) if len(line.split()) > 1 else len(pockets) + 1
            current = {"rank": rank, "score": None, "center_x": None, "center_y": None,
                       "center_z": None, "volume": None, "residues": []}
        elif current and ":" in line:
            key, _, val = line.partition(":")
            val = val.strip()
            key = key.strip().lower()
            if "score" in key or "druggability" in key:
                current["score"] = _safe_float(val)
            elif "volume" in key:
                current["volume"] = _safe_float(val)
            elif "x_barycenter" in key or "center x" in key:
                current["center_x"] = _safe_float(val)
            elif "y_barycenter" in key or "center y" in key:
                current["center_y"] = _safe_float(val)
            elif "z_barycenter" in key or "center z" in key:
                current["center_z"] = _safe_float(val)

    if current:
        pockets.append(current)

    pockets.sort(key=lambda p: -(p["score"] or 0))
    for i, p in enumerate(pockets, start=1):
        p["rank"] = i
    return pockets


def _safe_float(val) -> float | None:
    try:
        return float(str(val).strip())
    except (TypeError, ValueError):
        return None


def _parse_residues(residue_str: str) -> list[str]:
    return [r.strip() for r in residue_str.split(",") if r.strip()]
