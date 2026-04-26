from __future__ import annotations


def run_rdkit_enum(
    smiles: str,
    mode: str = "brics",
    rgroup_smiles: list | None = None,
    attachment_points: list | None = None,
    max_compounds: int = 10_000,
    filter_lipinski: bool = True,
) -> dict:
    from rdkit.Chem import MolFromSmiles, MolToSmiles

    mol = MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    if mode == "brics":
        results = _brics_enum(mol, max_compounds)
    elif mode == "recap":
        results = _recap_enum(mol, max_compounds)
    elif mode == "rgroup":
        if not rgroup_smiles:
            raise ValueError("rgroup_smiles required for rgroup mode")
        results = _rgroup_enum(mol, rgroup_smiles, attachment_points, max_compounds)
    else:
        raise ValueError(f"Unknown mode: {mode}")

    n_generated = len(results)
    if filter_lipinski:
        results = [s for s in results if _passes_ro5(s)]

    return {
        "smiles_list":    results[:max_compounds],
        "n_generated":    n_generated,
        "n_returned":     min(len(results), max_compounds),
        "scaffold_smiles": _get_scaffold(smiles),
    }


def _brics_enum(mol, max_compounds: int) -> list[str]:
    from rdkit.Chem import BRICS, MolFromSmiles, MolToSmiles

    frags = list(BRICS.BRICSDecompose(mol))
    enum = list(BRICS.BRICSBuild([MolFromSmiles(f) for f in frags if MolFromSmiles(f)]))
    return [MolToSmiles(m) for m in enum if m][:max_compounds]


def _recap_enum(mol, max_compounds: int) -> list[str]:
    from rdkit.Chem import BRICS, MolFromSmiles, MolToSmiles
    from rdkit.Chem import Recap

    tree = Recap.RecapDecompose(mol)
    leaves = tree.GetLeaves()
    frags = [MolFromSmiles(s) for s in leaves if MolFromSmiles(s)]
    enum = list(BRICS.BRICSBuild(frags))
    return [MolToSmiles(m) for m in enum if m][:max_compounds]


def _rgroup_enum(mol, rgroup_smiles: list, attachment_points, max_compounds: int) -> list[str]:
    from rdkit.Chem import AllChem, MolFromSmiles, MolToSmiles

    results = []
    for smi in rgroup_smiles:
        rg = MolFromSmiles(smi)
        if rg is None:
            continue
        try:
            combos = AllChem.ReplaceSubstructs(mol, MolFromSmiles("[*]"), rg, replaceAll=False)
            for c in combos:
                results.append(MolToSmiles(c))
                if len(results) >= max_compounds:
                    return results
        except Exception:
            continue
    return results


def _passes_ro5(smiles: str) -> bool:
    from rdkit.Chem import MolFromSmiles
    from rdkit.Chem import rdMolDescriptors

    mol = MolFromSmiles(smiles)
    if mol is None:
        return False
    mw   = rdMolDescriptors.CalcExactMolWt(mol)
    hbd  = rdMolDescriptors.CalcNumHBD(mol)
    hba  = rdMolDescriptors.CalcNumHBA(mol)
    logp = rdMolDescriptors.CalcCrippenDescriptors(mol)[0]
    return mw <= 500 and hbd <= 5 and hba <= 10 and logp <= 5


def _get_scaffold(smiles: str) -> str:
    from rdkit.Chem import MolFromSmiles, MolToSmiles
    from rdkit.Chem.Scaffolds import MurckoScaffold

    mol = MolFromSmiles(smiles)
    return MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol)) if mol else smiles
