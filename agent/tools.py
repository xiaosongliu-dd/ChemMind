"""
ChemMind — agent/tools.py
==========================
ToolRegistry: maps tool names to Python callables and exposes
Anthropic-compatible tool specs for native tool_use API calls.

Adding a new tool
-----------------
1. Implement `tools/{category}/mytool.py::run_mytool(**kwargs) -> dict`
2. Register it below with @registry.register(...)
3. Write the L1 skill markdown in skills/L1/{category}/mytool.md
"""

from __future__ import annotations
import importlib
from typing import Any, Callable


class ToolRegistry:
    def __init__(self):
        self._tools: dict[str, dict] = {}
        self._register_all()

    def register(
        self,
        name: str,
        description: str,
        parameters: dict,
        fn: Callable,
    ) -> None:
        self._tools[name] = {
            "description": description,
            "parameters":  parameters,
            "fn":          fn,
        }

    def call(self, name: str, inputs: dict[str, Any]) -> dict:
        if name not in self._tools:
            raise ValueError(f"Unknown tool: {name}. Available: {list(self._tools)}")
        return self._tools[name]["fn"](**inputs)

    def anthropic_tool_specs(self) -> list[dict]:
        """Return tool definitions in Anthropic tool_use format."""
        return [
            {
                "name":         name,
                "description":  spec["description"],
                "input_schema": {"type": "object", "properties": spec["parameters"]},
            }
            for name, spec in self._tools.items()
        ]

    def _register_all(self) -> None:
        """Register all ChemMind L1 tools."""
        self._reg("boltz2",         "tools.structure.boltz2",       "run_boltz2",
            "Predict protein/ligand/antibody-antigen structure and binding affinity with Boltz-2")
        self._reg("esmfold",        "tools.structure.esmfold",      "run_esmfold",
            "Fast protein structure prediction with ESMFold (no MSA required)")
        self._reg("colabfold",      "tools.structure.colabfold",    "run_colabfold",
            "MSA-based structure prediction via ColabFold / AlphaFold2")
        self._reg("pocket_detect",  "tools.structure.pocket_detect","run_pocket_detect",
            "Detect and rank druggable binding pockets with P2Rank and fpocket")
        self._reg("rfantibody",     "tools.biologics.rfantibody",   "run_rfantibody",
            "De novo antibody (scFv/VHH) backbone design with RFdiffusion-antibody")
        self._reg("abdiffuser",     "tools.biologics.abdiffuser",   "run_abdiffuser",
            "Full-atom antibody CDR co-design (sequence + structure) with AbDiffuser")
        self._reg("igfold",         "tools.biologics.igfold",       "run_igfold",
            "Antibody structure prediction and CDR-H3 loop modelling with IgFold")
        self._reg("abodybuilder3",  "tools.biologics.abodybuilder3","run_abodybuilder3",
            "Fast antibody and nanobody structure prediction with ABodyBuilder3")
        self._reg("proteinmpnn",    "tools.biologics.proteinmpnn",  "run_proteinmpnn",
            "Design amino acid sequences on a fixed protein backbone with ProteinMPNN")
        self._reg("bindcraft",      "tools.biologics.bindcraft",    "run_bindcraft",
            "Design peptide mini-binders and small protein binders with BindCraft")
        self._reg("rfdiffusion_pep","tools.biologics.rfdiffusion_pep","run_rfdiffusion_pep",
            "Design peptide binders with motif scaffolding via RFdiffusion")
        self._reg("diffsbdd",       "tools.generation.diffsbdd",    "run_diffsbdd",
            "3D structure-based molecule generation conditioned on a binding pocket")
        self._reg("reinvent4",      "tools.generation.reinvent4",   "run_reinvent4",
            "RL-based SMILES generation and lead optimisation with REINVENT 4")
        self._reg("molgpt",         "tools.generation.molgpt",      "run_molgpt",
            "GPT-based conditional SMILES generation with MolGPT")
        self._reg("gnina",          "tools.docking.gnina",          "run_gnina",
            "CNN-scored molecular docking for pose prediction and rescoring with GNINA")
        self._reg("autodock_gpu",   "tools.docking.autodock_gpu",   "run_autodock_gpu",
            "High-throughput GPU-accelerated virtual screening with AutoDock-GPU")
        self._reg("diffdock",       "tools.docking.diffdock",       "run_diffdock",
            "Diffusion-based blind docking and pose prediction with DiffDock")
        self._reg("openmm",         "tools.md_fep.openmm",          "run_openmm",
            "Molecular dynamics simulation with OpenMM (CUDA, NVT/NPT/metadynamics)")
        self._reg("fep_openfe",     "tools.md_fep.fep_openfe",      "run_fep_openfe",
            "Relative free energy perturbation with OpenFE and Lomap perturbation graph")
        self._reg("mdanalysis",     "tools.md_fep.mdanalysis",      "run_mdanalysis",
            "MD trajectory analysis: RMSD, RMSF, contact maps, pocket volume")
        self._reg("rdkit_props",    "tools.admet.rdkit_props",      "run_rdkit_props",
            "Compute QED, SA score, Lipinski, PAINS filters and Morgan fingerprints")
        self._reg("admetlab3",      "tools.admet.admetlab3",        "run_admetlab3",
            "Predict 70+ ADMET endpoints via ADMETlab3 API (batch mode)")
        self._reg("deeppurpose",    "tools.admet.deeppurpose",      "run_deeppurpose",
            "Drug-target interaction prediction and pIC50 estimation with DeepPurpose")
        self._reg("askcos",         "tools.data.askcos",            "run_askcos",
            "Retrosynthetic route planning and buyability scoring with ASKCOS")
        self._reg("chembl",         "tools.data.chembl",            "run_chembl",
            "Query ChEMBL for bioactivity data, SAR analysis and target information")
        self._reg("pdb_fetch",      "tools.data.pdb_fetch",         "run_pdb_fetch",
            "Fetch PDB structures, SIFTS annotations and ligand data from RCSB")

    def _reg(self, name: str, module: str, fn_name: str, description: str) -> None:
        """Lazy-import helper — tool modules are only loaded when first called."""
        def _caller(**kwargs):
            mod = importlib.import_module(module)
            return getattr(mod, fn_name)(**kwargs)
        self.register(name, description, {}, _caller)
