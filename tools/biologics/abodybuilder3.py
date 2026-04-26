from __future__ import annotations

import tempfile
from pathlib import Path


def run_abodybuilder3(
    heavy_sequence: str,
    light_sequence: str | None = None,
    numbering_scheme: str = "imgt",
    output_dir: str | None = None,
) -> dict:
    """
    Predict antibody or nanobody 3D structure using ABodyBuilder3 (ImmuneBuilder).

    Parameters
    ----------
    heavy_sequence   : VH (or VHH for nanobody) amino acid sequence
    light_sequence   : VL sequence; omit for nanobody / VHH single-domain
    numbering_scheme : antibody numbering: "imgt" (default), "chothia", or "kabat"
    output_dir       : write PDB here; uses temp dir if None

    Notes
    -----
    ABodyBuilder3 / ImmuneBuilder uses an ensemble of 4 IGFold-like networks
    fine-tuned on antibody structures with IMGT numbering + CDR annotation.
    It produces IMGT-annotated CDR loop residue labels in the B-factor column.
    """
    try:
        from ImmuneBuilder import ABodyBuilder2, NanoBodyBuilder2
    except ImportError:
        raise RuntimeError(
            "ImmuneBuilder not installed. Run: pip install ImmuneBuilder"
        )

    outdir = Path(output_dir or tempfile.mkdtemp(prefix="abodybuilder3_"))
    outdir.mkdir(parents=True, exist_ok=True)
    pdb_path = str(outdir / "antibody.pdb")

    if light_sequence is None:
        predictor = NanoBodyBuilder2(numbering_scheme=numbering_scheme)
        ab = predictor.predict({"H": heavy_sequence})
        mode = "nanobody"
    else:
        predictor = ABodyBuilder2(numbering_scheme=numbering_scheme)
        ab = predictor.predict({"H": heavy_sequence, "L": light_sequence})
        mode = "antibody"

    ab.save(pdb_path)

    error = None
    try:
        error = float(ab.error)
    except Exception:
        pass

    return {
        "pdb_path":        pdb_path,
        "predicted_error": error,
        "heavy_sequence":  heavy_sequence,
        "light_sequence":  light_sequence,
        "mode":            mode,
        "numbering":       numbering_scheme,
    }
