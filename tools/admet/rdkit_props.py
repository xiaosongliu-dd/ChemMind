from __future__ import annotations


def run_rdkit_props(
    smiles: str | list[str],
    fingerprint: str | None = None,
    fp_bits: int = 2048,
) -> dict:
    """
    Compute RDKit physicochemical descriptors and rule-based developability
    flags (Lipinski, Veber, Egan) for a SMILES — single or batch.

    Single mode (`smiles: str`):
        {smiles, mw, logp, hbd, hba, tpsa, rotatable_bonds, aromatic_rings,
         qed, sa_score, lipinski_pass, veber_pass, egan_pass, fingerprint?}

    Batch mode (`smiles: list[str]`):
        {results: [single-mode dict, ...], n_compounds}
        invalid SMILES yield {"smiles": <input>, "error": "Invalid SMILES"}
    """
    if isinstance(smiles, str):
        return _props_one(smiles, fingerprint, fp_bits, raise_invalid=True)

    results = [_props_one(s, fingerprint, fp_bits, raise_invalid=False)
               for s in smiles]
    return {
        "results":     results,
        "n_compounds": len(smiles),
    }


def _props_one(smiles: str, fingerprint, fp_bits, raise_invalid: bool) -> dict:
    from rdkit import Chem
    from rdkit.Chem import Crippen, QED, rdMolDescriptors

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        if raise_invalid:
            raise ValueError(f"Invalid SMILES: {smiles}")
        return {"smiles": smiles, "error": "Invalid SMILES"}

    mw   = rdMolDescriptors.CalcExactMolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd  = rdMolDescriptors.CalcNumHBD(mol)
    hba  = rdMolDescriptors.CalcNumHBA(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    rb   = rdMolDescriptors.CalcNumRotatableBonds(mol)

    props = {
        "smiles":          smiles,
        "mw":              mw,
        "logp":            logp,
        "hbd":             hbd,
        "hba":             hba,
        "tpsa":            tpsa,
        "rotatable_bonds": rb,
        "aromatic_rings":  rdMolDescriptors.CalcNumAromaticRings(mol),
        "qed":             QED.qed(mol),
        "sa_score":        _sa_score(mol),
        "lipinski_pass":   mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10,
        "veber_pass":      rb <= 10 and tpsa <= 140,
        "egan_pass":       -1 <= logp <= 6 and tpsa <= 132,
    }

    if fingerprint:
        from rdkit.Chem import AllChem
        radius = 3 if fingerprint == "morgan3" else 2
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=fp_bits)
        props["fingerprint"]      = list(fp)
        props["fingerprint_type"] = fingerprint
        props["fingerprint_bits"] = fp_bits

    return props


def _sa_score(mol) -> float | None:
    """SA score lives in rdkit/Contrib/SA_Score and isn't on sys.path by
    default. Add it lazily; return None if anything fails."""
    try:
        import os
        import sys
        from rdkit.Chem import RDConfig
        sa_dir = os.path.join(RDConfig.RDContribDir, "SA_Score")
        if sa_dir not in sys.path:
            sys.path.append(sa_dir)
        import sascorer
        return float(sascorer.calculateScore(mol))
    except Exception:
        return None
