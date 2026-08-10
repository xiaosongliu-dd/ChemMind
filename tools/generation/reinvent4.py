from __future__ import annotations

import subprocess
import tempfile
import tomllib
from pathlib import Path


def run_reinvent4(
    mode: str = "sampling",
    scaffold_smarts: str | None = None,
    scoring_config: dict | None = None,
    n_steps: int = 100,
    batch_size: int = 64,
    model_file: str | None = None,
    reinvent_bin: str = "reinvent",
    output_dir: str | None = None,
) -> dict:
    """
    Generate and optimise molecules with REINVENT 4 (Molecular AI, AstraZeneca).

    `mode`:
        "sampling"       — sample from a prior model (fastest, no RL)
        "reinforcement"  — RL-guided optimisation with a scoring function

    `scoring_config`: dict passed verbatim into the TOML [scoring] section.
    Typical scoring components include component dicts with `type` and `params`
    keys; see the REINVENT4 documentation for the full component library.

    Returns:
        {
            smiles:      list[str],    # generated / optimised SMILES
            scores:      list[float],  # per-SMILES total score (0–1)
            n_generated: int,
            n_steps_run: int,
            mode:        str,
            output_dir:  str,
        }

    Raises RuntimeError("REINVENT4 failed: ...") on non-zero subprocess exit.
    """
    outdir = Path(output_dir) if output_dir else Path(tempfile.mkdtemp(prefix="reinvent4_"))
    outdir.mkdir(parents=True, exist_ok=True)

    config = _build_config(
        mode=mode,
        scaffold_smarts=scaffold_smarts,
        scoring_config=scoring_config or {},
        n_steps=n_steps,
        batch_size=batch_size,
        model_file=model_file,
        outdir=outdir,
    )
    config_path = outdir / "config.toml"
    _write_toml(config, config_path)

    _run_reinvent(reinvent_bin, config_path, outdir)

    smiles, scores = _parse_results(outdir)
    return {
        "smiles":      smiles,
        "scores":      scores,
        "n_generated": len(smiles),
        "n_steps_run": n_steps,
        "mode":        mode,
        "output_dir":  str(outdir),
    }


def _build_config(
    mode: str,
    scaffold_smarts: str | None,
    scoring_config: dict,
    n_steps: int,
    batch_size: int,
    model_file: str | None,
    outdir: Path,
) -> dict:
    if mode not in ("sampling", "reinforcement"):
        raise ValueError(f"Unknown mode: {mode!r}. Choose 'sampling' or 'reinforcement'.")

    base = {
        "run_type":  mode if mode != "reinforcement" else "reinforcement_learning",
        "device":    "cuda",
        "tb_logdir": str(outdir / "tb_logs"),
        "json_out_config": str(outdir / "effective_config.json"),
    }

    if mode == "sampling":
        base["parameters"] = {
            "model_file":  model_file or "priors/reinvent.prior",
            "output_file": str(outdir / "sampled.csv"),
            "num_smiles":  n_steps * batch_size,
            "batch_size":  batch_size,
            "with_prior_nll": True,
        }
    else:
        base["parameters"] = {
            "prior":        model_file or "priors/reinvent.prior",
            "agent":        model_file or "priors/reinvent.prior",
            "summary_csv_prefix": str(outdir / "results"),
            "batch_size":   batch_size,
            "n_steps":      n_steps,
        }
        if scoring_config:
            base["scoring"] = scoring_config
        else:
            base["scoring"] = {
                "type":  "custom_product",
                "params": {"transform": True},
                "component": [
                    {
                        "custom_alerts": {
                            "endpoint": [{"name": "Alerts", "weight": 1.0,
                                          "params": {"smarts": []}}]
                        }
                    }
                ],
            }

    if scaffold_smarts:
        base["diversity_filter"] = {
            "type": "IdenticalMurckoScaffold",
            "bucket_size": 25,
            "minscore": 0.4,
            "minsimilarity": 0.4,
        }
        base["inception"] = {
            "memory_size": 20,
            "sample_size": 5,
            "smiles": [scaffold_smarts],
        }

    return base


def _write_toml(config: dict, path: Path) -> None:
    try:
        import tomli_w
        path.write_bytes(tomli_w.dumps(config).encode())
    except ImportError:
        import json
        json_path = path.with_suffix(".json")
        json_path.write_text(json.dumps(config, indent=2))
        path.write_text(_dict_to_toml(config))


def _dict_to_toml(d: dict, prefix: str = "") -> str:
    """Minimal TOML serialiser (no arrays-of-tables, scalars only in nested)."""
    lines = []
    nested = {}
    for k, v in d.items():
        key = f"{prefix}.{k}" if prefix else k
        if isinstance(v, dict):
            nested[key] = v
        elif isinstance(v, bool):
            lines.append(f"{k} = {'true' if v else 'false'}")
        elif isinstance(v, str):
            lines.append(f'{k} = "{v}"')
        elif isinstance(v, (int, float)):
            lines.append(f"{k} = {v}")
        elif isinstance(v, list):
            items = ", ".join(f'"{x}"' if isinstance(x, str) else str(x) for x in v)
            lines.append(f"{k} = [{items}]")
    result = "\n".join(lines)
    for key, sub in nested.items():
        result += f"\n\n[{key}]\n" + _dict_to_toml(sub)
    return result


def _run_reinvent(reinvent_bin: str, config_path: Path, outdir: Path) -> None:
    cmd = [reinvent_bin, "-f", str(config_path)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"REINVENT4 failed (exit {result.returncode}):\n{result.stderr[-2000:]}"
        )


def _parse_results(outdir: Path) -> tuple[list[str], list[float]]:
    """Parse CSV output from REINVENT4 (sampled.csv or results*.csv)."""
    csv_files = sorted(outdir.glob("sampled.csv")) + sorted(outdir.glob("results*.csv"))
    if not csv_files:
        return [], []

    import csv as csv_mod
    smiles: list[str] = []
    scores: list[float] = []

    for csv_path in csv_files[:1]:
        with open(csv_path) as f:
            reader = csv_mod.DictReader(f)
            for row in reader:
                smi = row.get("SMILES") or row.get("Smiles") or row.get("smiles") or ""
                score = row.get("Score") or row.get("score") or row.get("TotalScore") or "0"
                if smi.strip():
                    smiles.append(smi.strip())
                    try:
                        scores.append(float(score))
                    except (ValueError, TypeError):
                        scores.append(0.0)

    paired = sorted(zip(scores, smiles), reverse=True)
    if paired:
        scores, smiles = zip(*paired)
        return list(smiles), list(scores)
    return [], []
