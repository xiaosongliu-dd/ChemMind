from __future__ import annotations


_PRESETS = {
    "default": ["PAINS", "BRENK"],
    "strict":  ["PAINS", "BRENK", "NIH", "ZINC"],
    "all":     [
        "PAINS_A", "PAINS_B", "PAINS_C",
        "BRENK", "NIH", "ZINC",
        "CHEMBL_Glaxo", "CHEMBL_Dundee", "CHEMBL_BMS",
        "CHEMBL_LINT", "CHEMBL_MLSMR", "CHEMBL_SureChEMBL",
    ],
}


def run_ligand_filter(
    smiles: str | list[str],
    filter_sets: str | list[str] = "default",
    custom_smarts: list[str] | None = None,
    return_canonical: bool = True,
) -> dict:
    """
    Flag SMILES against PAINS / BRENK / reactive-group structural alerts using
    RDKit's FilterCatalog, plus optional user-supplied SMARTS.

    Single mode (`smiles: str`):
        {passed, n_alerts, failures, smiles_canonical, filter_sets_used}

    Batch mode (`smiles: list[str]`, even length 1):
        {passed_smiles, failed_smiles, n_passed, n_failed, results, filter_sets_used}

    Property-based filters (Ro5, QED, MW) are intentionally out of scope —
    use tools.admet.rdkit_props or tools.enumeration.rdkit_enum's _passes_ro5.
    """
    resolved = _resolve_filter_sets(filter_sets)
    catalogs = _build_catalogs(resolved)
    custom   = _compile_custom_smarts(custom_smarts or [])

    if isinstance(smiles, str):
        return _filter_one(smiles, catalogs, custom, resolved, return_canonical,
                           raise_on_invalid=True)
    return _filter_batch(smiles, catalogs, custom, resolved, return_canonical)


# ── filter-set resolution ─────────────────────────────────────────────────────

def _resolve_filter_sets(filter_sets: str | list[str]) -> list[str]:
    if isinstance(filter_sets, str):
        if filter_sets == "lilly":
            raise NotImplementedError(
                "Lilly MedChem rules planned via rd_filters; not yet wired"
            )
        if filter_sets in _PRESETS:
            return list(_PRESETS[filter_sets])
        raise ValueError(f"Unknown filter set: {filter_sets}")

    resolved: list[str] = []
    for name in filter_sets:
        if name == "lilly":
            raise NotImplementedError(
                "Lilly MedChem rules planned via rd_filters; not yet wired"
            )
        if name in _PRESETS:
            resolved.extend(_PRESETS[name])
        else:
            resolved.append(name)
    return resolved


# ── catalog construction ──────────────────────────────────────────────────────

def _build_catalogs(names: list[str]) -> list[tuple[str, object]]:
    from rdkit.Chem import FilterCatalog as FC

    catalogs_enum = FC.FilterCatalogParams.FilterCatalogs
    out = []
    for name in names:
        enum_val = getattr(catalogs_enum, name, None)
        if enum_val is None:
            raise ValueError(f"Unknown filter set: {name}")
        params = FC.FilterCatalogParams()
        params.AddCatalog(enum_val)
        out.append((name, FC.FilterCatalog(params)))
    return out


def _compile_custom_smarts(patterns: list[str]) -> list[tuple[str, object]]:
    from rdkit import Chem

    compiled = []
    for i, smarts in enumerate(patterns):
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            raise ValueError(f"Invalid SMARTS pattern: {smarts}")
        compiled.append((smarts, patt))
    return compiled


# ── single-molecule scan ──────────────────────────────────────────────────────

def _filter_one(
    smiles: str,
    catalogs: list[tuple[str, object]],
    custom: list[tuple[str, object]],
    resolved: list[str],
    return_canonical: bool,
    raise_on_invalid: bool,
) -> dict:
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        if raise_on_invalid:
            raise ValueError(f"Invalid SMILES: {smiles}")
        return {
            "passed":           False,
            "n_alerts":         0,
            "failures":         [],
            "smiles_canonical": smiles,
            "filter_sets_used": list(resolved),
            "error":            "Invalid SMILES",
        }

    canonical = Chem.MolToSmiles(mol) if return_canonical else smiles
    failures: list[dict] = []

    for cat_name, catalog in catalogs:
        for entry in catalog.GetMatches(mol):
            failures.append(_entry_to_failure(cat_name, entry))

    for i, (smarts, patt) in enumerate(custom):
        if mol.HasSubstructMatch(patt):
            failures.append({
                "catalog":     "custom",
                "alert_name":  f"custom_{i}",
                "description": f"User-supplied SMARTS: {smarts}",
                "smarts":      smarts,
            })

    return {
        "passed":           len(failures) == 0,
        "n_alerts":         len(failures),
        "failures":         failures,
        "smiles_canonical": canonical,
        "filter_sets_used": list(resolved),
    }


def _entry_to_failure(catalog_name: str, entry) -> dict:
    try:
        alert_name = entry.GetDescription()
    except Exception:
        alert_name = ""
    description = alert_name
    smarts = ""
    try:
        smarts = entry.GetSmarts() or ""
    except Exception:
        pass
    return {
        "catalog":     catalog_name,
        "alert_name":  alert_name,
        "description": description,
        "smarts":      smarts,
    }


# ── batch scan ────────────────────────────────────────────────────────────────

def _filter_batch(
    smiles_list: list[str],
    catalogs: list[tuple[str, object]],
    custom: list[tuple[str, object]],
    resolved: list[str],
    return_canonical: bool,
) -> dict:
    results = [
        _filter_one(s, catalogs, custom, resolved, return_canonical,
                    raise_on_invalid=False)
        for s in smiles_list
    ]
    passed = [r["smiles_canonical"] for r in results if r["passed"]]
    failed = [r["smiles_canonical"] for r in results if not r["passed"]]
    return {
        "passed_smiles":    passed,
        "failed_smiles":    failed,
        "n_passed":         len(passed),
        "n_failed":         len(failed),
        "results":          results,
        "filter_sets_used": list(resolved),
    }
