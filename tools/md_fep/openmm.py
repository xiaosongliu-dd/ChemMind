from __future__ import annotations

import tempfile
import time
from pathlib import Path


_KCAL_PER_KJ = 1 / 4.184


def run_openmm(
    structure_pdb: str,
    forcefield: str = "amber14-all.xml",
    water_model: str = "amber14/tip3pfb.xml",
    ensemble: str = "NPT",
    temperature_k: float = 300.0,
    pressure_bar: float = 1.0,
    n_steps: int = 1_000_000,
    timestep_fs: float = 2.0,
    report_interval: int = 10_000,
    output_dir: str | None = None,
) -> dict:
    """
    Run a molecular dynamics simulation with OpenMM.

    Performs the standard three-phase protocol:
      1. Energy minimisation (until convergence or 10,000 steps max)
      2. NVT equilibration for n_steps // 10 steps
      3. Production run (NVT or NPT) for n_steps steps

    `ensemble`: "NPT" (constant pressure/temperature) or "NVT" (constant volume/temperature).

    Returns:
        {
            trajectory_path:    str,          # DCD trajectory of production run
            energy_path:        str,          # CSV of energy / temperature vs step
            final_structure:    str,          # PDB of last frame
            potential_energy_kcal_mean: float | None,
            temperature_k_mean: float | None,
            n_steps_run:        int,
            timestep_fs:        float,
            runtime_s:          float,
        }

    Raises RuntimeError("OpenMM not installed: ...") if openmm is missing.
    """
    try:
        import openmm
        import openmm.app as app
        import openmm.unit as unit
    except ImportError as e:
        raise RuntimeError(
            f"openmm not installed: {e}. "
            "Install with: conda install -c conda-forge openmm"
        )

    outdir = Path(output_dir) if output_dir else Path(tempfile.mkdtemp(prefix="openmm_"))
    outdir.mkdir(parents=True, exist_ok=True)

    t0 = time.perf_counter()

    pdb = app.PDBFile(structure_pdb)
    ff  = app.ForceField(forcefield, water_model)

    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.addSolvent(ff, model="tip3p", padding=1.0 * unit.nanometers)

    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometers,
        constraints=app.HBonds,
    )

    integrator = openmm.LangevinMiddleIntegrator(
        temperature_k * unit.kelvin,
        1.0 / unit.picoseconds,
        timestep_fs * unit.femtoseconds,
    )

    if ensemble.upper() == "NPT":
        system.addForce(openmm.MonteCarloBarostat(
            pressure_bar * unit.bar,
            temperature_k * unit.kelvin,
        ))

    simulation = app.Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)

    simulation.minimizeEnergy(maxIterations=10_000)

    eq_steps = max(n_steps // 10, 1_000)
    simulation.context.setVelocitiesToTemperature(temperature_k * unit.kelvin)
    simulation.step(eq_steps)

    traj_path   = outdir / "trajectory.dcd"
    energy_path = outdir / "energy.csv"
    final_pdb   = outdir / "final_frame.pdb"

    simulation.reporters.append(app.DCDReporter(str(traj_path), report_interval))
    simulation.reporters.append(app.StateDataReporter(
        str(energy_path),
        report_interval,
        step=True,
        potentialEnergy=True,
        temperature=True,
        volume=True,
        density=True,
    ))

    simulation.step(n_steps)

    positions = simulation.context.getState(getPositions=True).getPositions()
    with open(str(final_pdb), "w") as f:
        app.PDBFile.writeFile(simulation.topology, positions, f)

    pe_mean, temp_mean = _parse_energy_csv(energy_path)

    return {
        "trajectory_path":           str(traj_path),
        "energy_path":               str(energy_path),
        "final_structure":           str(final_pdb),
        "potential_energy_kcal_mean": pe_mean,
        "temperature_k_mean":        temp_mean,
        "n_steps_run":               n_steps,
        "timestep_fs":               timestep_fs,
        "runtime_s":                 round(time.perf_counter() - t0, 1),
    }


def _parse_energy_csv(energy_path: Path) -> tuple[float | None, float | None]:
    import csv as csv_mod
    pe_vals, temp_vals = [], []
    try:
        with open(energy_path) as f:
            reader = csv_mod.DictReader(f)
            for row in reader:
                for key, container in (
                    ("Potential Energy (kJ/mole)", pe_vals),
                    ("Temperature (K)", temp_vals),
                ):
                    val = row.get(key)
                    if val is not None:
                        try:
                            container.append(float(val))
                        except ValueError:
                            pass
    except (OSError, StopIteration):
        pass

    pe_mean = round(sum(pe_vals) / len(pe_vals) * _KCAL_PER_KJ, 2) if pe_vals else None
    temp_mean = round(sum(temp_vals) / len(temp_vals), 2) if temp_vals else None
    return pe_mean, temp_mean
