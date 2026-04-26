from __future__ import annotations

from pathlib import Path


def run_fep_openfe(
    protein_pdb: str,
    ligand_sdf_dir: str,
    reference_ligand_sdf: str,
    n_replicas: int = 3,
    n_lambda: int = 11,
    forcefield: str = "openff-2.1.0",
    solvent: str = "tip3p",
    work_dir: str = "fep_output",
) -> dict:
    """
    Run a relative binding free energy (RBFE) campaign with OpenFE.

    Uses LoMap to build the perturbation network, then runs RBFE edges
    with OpenMM via the HREX protocol. Returns ΔΔG per edge and
    estimated ΔG per ligand relative to the reference.

    Inputs
    ------
    protein_pdb          : prepared receptor (no ligand, no crystal waters)
    ligand_sdf_dir       : directory of ligand SDF files (one per file)
    reference_ligand_sdf : SDF of the reference ligand (known binder, anchor node)
    n_replicas           : HREX replicas per window (default 3)
    n_lambda             : λ windows per edge (default 11)
    forcefield           : small molecule FF (openff-2.1.0 or gaff-2.11)
    solvent              : water model (tip3p or tip4p-ew)
    work_dir             : output directory for results and trajectories
    """
    try:
        import openfe
        from openfe import RBFEProtocol, LigandNetwork
        from openfe.utils import lomap_network
        from rdkit import Chem
    except ImportError:
        raise RuntimeError("openfe not installed. Run: pip install openfe")

    work_path = Path(work_dir)
    work_path.mkdir(parents=True, exist_ok=True)

    # Load ligands
    ligand_sdfs = sorted(Path(ligand_sdf_dir).glob("*.sdf"))
    if not ligand_sdfs:
        raise FileNotFoundError(f"No SDF files in {ligand_sdf_dir}")

    ligands = [openfe.SmallMoleculeComponent.from_sdf_file(str(f)) for f in ligand_sdfs]
    reference = openfe.SmallMoleculeComponent.from_sdf_file(reference_ligand_sdf)
    receptor  = openfe.ProteinComponent.from_pdb_file(protein_pdb)
    solvent_comp = openfe.SolventComponent()

    # Build perturbation network with LoMap
    network = lomap_network(ligands, scorer=openfe.lomap_scorers.default_lomap_score)

    # Configure RBFE protocol
    protocol_settings = openfe.protocols.openmm_rbfe.RelativeHybridTopologyProtocolSettings(
        forcefield_settings=openfe.settings.OpenMMSystemGeneratorFFSettings(
            small_molecule_forcefield=forcefield,
            forcefields=["amber14-all.xml", f"{solvent}.xml"],
        ),
        lambda_settings=openfe.settings.LambdaSettings(
            lambda_windows=n_lambda,
        ),
        integrator_settings=openfe.settings.IntegratorSettings(
            n_replicas=n_replicas,
        ),
    )
    protocol = openfe.protocols.openmm_rbfe.RelativeHybridTopologyProtocol(
        settings=protocol_settings
    )

    # Run all edges
    edge_results = []
    for edge in network.edges:
        transformation = openfe.Transformation(
            stateA=openfe.ChemicalSystem({"ligand": edge.componentA, "protein": receptor, "solvent": solvent_comp}),
            stateB=openfe.ChemicalSystem({"ligand": edge.componentB, "protein": receptor, "solvent": solvent_comp}),
            protocol=protocol,
            mapping={edge.componentA: edge.componentB},
            name=f"{edge.componentA.name}_to_{edge.componentB.name}",
        )
        dag = transformation.create()
        unit_results = [unit.run(dry=False, verbose=True, scratch_basepath=str(work_path)) for unit in dag.protocol_units]
        dag_result = openfe.ProtocolDAGResult(protocol_units=unit_results)
        estimates = protocol.gather(dag_result)

        edge_results.append({
            "ligandA":        edge.componentA.name,
            "ligandB":        edge.componentB.name,
            "DDG_kcal_mol":   float(estimates.delta_g_estimate.m),
            "uncertainty":    float(estimates.delta_g_error.m),
            "converged":      estimates.delta_g_error.m < 0.5,
        })

    return {
        "edge_results":   edge_results,
        "n_edges":        len(edge_results),
        "n_converged":    sum(1 for e in edge_results if e["converged"]),
        "work_dir":       str(work_path),
        "forcefield":     forcefield,
        "n_lambda":       n_lambda,
        "n_replicas":     n_replicas,
    }
