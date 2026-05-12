from __future__ import annotations

import tempfile
import time
from pathlib import Path


def run_mdanalysis(
    topology: str,
    trajectory: str,
    analyses: list[str] | None = None,
    reference_pdb: str | None = None,
    selection: str = "protein and backbone",
    output_dir: str | None = None,
) -> dict:
    """
    Analyse a molecular dynamics trajectory with MDAnalysis.

    `analyses` selects which quantities to compute:
        "rmsd"          — backbone RMSD vs reference frame (or reference_pdb)
        "rmsf"          — per-residue backbone RMSF
        "contacts"      — heavy-atom contact map averaged over trajectory
        "pocket_volume" — convex-hull volume of the selection over time

    Default analyses when None: ["rmsd", "rmsf"].

    Returns:
        {
            rmsd?:          {step: list[int], rmsd_angstrom: list[float]},
            rmsf?:          {resids: list[int], rmsf_angstrom: list[float]},
            contacts?:      {matrix: list[list[float]], resids: list[int]},
            pocket_volume?: {step: list[int], volume_angstrom3: list[float]},
            n_frames:       int,
            n_atoms:        int,
            runtime_s:      float,
        }

    Raises RuntimeError("MDAnalysis not installed: ...") if the package is missing.
    """
    try:
        import MDAnalysis as mda
        import numpy as np
    except ImportError as e:
        raise RuntimeError(
            f"MDAnalysis not installed: {e}. "
            "Install with: pip install MDAnalysis"
        )

    if analyses is None:
        analyses = ["rmsd", "rmsf"]

    analyses = [a.lower() for a in analyses]
    unknown = set(analyses) - {"rmsd", "rmsf", "contacts", "pocket_volume"}
    if unknown:
        raise ValueError(f"Unknown analyses: {unknown}. "
                         "Choose from 'rmsd', 'rmsf', 'contacts', 'pocket_volume'.")

    outdir = Path(output_dir) if output_dir else Path(tempfile.mkdtemp(prefix="mdanalysis_"))
    outdir.mkdir(parents=True, exist_ok=True)

    t0 = time.perf_counter()
    u = mda.Universe(topology, trajectory)
    sel = u.select_atoms(selection)
    n_frames = len(u.trajectory)
    n_atoms  = len(u.atoms)

    result: dict = {"n_frames": n_frames, "n_atoms": n_atoms}

    if "rmsd" in analyses:
        result["rmsd"] = _compute_rmsd(u, sel, reference_pdb, outdir)

    if "rmsf" in analyses:
        result["rmsf"] = _compute_rmsf(u, sel)

    if "contacts" in analyses:
        result["contacts"] = _compute_contacts(u, sel)

    if "pocket_volume" in analyses:
        result["pocket_volume"] = _compute_pocket_volume(u, sel)

    result["runtime_s"] = round(time.perf_counter() - t0, 1)
    return result


def _compute_rmsd(universe, sel, reference_pdb: str | None, outdir: Path) -> dict:
    import MDAnalysis as mda
    import numpy as np
    from MDAnalysis.analysis import rms

    if reference_pdb:
        ref = mda.Universe(reference_pdb)
        ref_sel = ref.select_atoms(sel.segids[0] if sel.segids else "protein and backbone")
    else:
        ref = universe
        ref_sel = sel

    R = rms.RMSD(sel, ref_sel, select="backbone", groupselections=None)
    R.run()

    data = R.results.rmsd.T
    steps = [int(s) for s in data[0]]
    rmsd_vals = [round(float(r), 4) for r in data[2]]

    csv_path = outdir / "rmsd.csv"
    csv_path.write_text("step,rmsd_angstrom\n" +
                        "\n".join(f"{s},{r}" for s, r in zip(steps, rmsd_vals)))
    return {"step": steps, "rmsd_angstrom": rmsd_vals}


def _compute_rmsf(universe, sel) -> dict:
    import numpy as np
    from MDAnalysis.analysis import rms

    R = rms.RMSF(sel)
    R.run()

    resids = [int(r) for r in sel.resids]
    rmsf_vals = [round(float(v), 4) for v in R.results.rmsf]

    per_res_resids: list[int] = []
    per_res_rmsf: list[float] = []
    seen = {}
    for resid, val in zip(resids, rmsf_vals):
        if resid not in seen:
            seen[resid] = []
        seen[resid].append(val)
    for resid, vals in seen.items():
        per_res_resids.append(resid)
        per_res_rmsf.append(round(sum(vals) / len(vals), 4))

    return {"resids": per_res_resids, "rmsf_angstrom": per_res_rmsf}


def _compute_contacts(universe, sel) -> dict:
    import numpy as np

    residues = list({a.resid for a in sel})
    residues.sort()
    n = len(residues)
    res_idx = {r: i for i, r in enumerate(residues)}

    contact_matrix = [[0.0] * n for _ in range(n)]
    n_frames = len(universe.trajectory)
    cutoff = 4.5  # Å

    for ts in universe.trajectory:
        positions = sel.positions
        atom_resids = [a.resid for a in sel]
        for ii, (pos_i, res_i) in enumerate(zip(positions, atom_resids)):
            for jj, (pos_j, res_j) in enumerate(zip(positions, atom_resids)):
                if jj <= ii or abs(res_idx[res_i] - res_idx[res_j]) < 2:
                    continue
                dist = float(np.linalg.norm(pos_i - pos_j))
                if dist < cutoff:
                    ri, rj = res_idx[res_i], res_idx[res_j]
                    contact_matrix[ri][rj] += 1.0
                    contact_matrix[rj][ri] += 1.0

    matrix = [[round(v / max(n_frames, 1), 4) for v in row] for row in contact_matrix]
    return {"matrix": matrix, "resids": residues}


def _compute_pocket_volume(universe, sel) -> dict:
    import numpy as np

    steps, volumes = [], []
    for ts in universe.trajectory:
        positions = sel.positions
        if len(positions) < 4:
            volumes.append(0.0)
        else:
            try:
                from scipy.spatial import ConvexHull
                hull = ConvexHull(positions)
                volumes.append(round(float(hull.volume), 2))
            except Exception:
                volumes.append(0.0)
        steps.append(int(ts.frame))

    return {"step": steps, "volume_angstrom3": volumes}
