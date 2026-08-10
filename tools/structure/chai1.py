from __future__ import annotations

import tempfile
from pathlib import Path


def run_chai1(
    protein_sequences: str | list[str],
    ligand_smiles: str | list[str] | None = None,
    dna_sequences: list[str] | None = None,
    rna_sequences: list[str] | None = None,
    seeds: list[int] | None = None,
    num_trunk_recycles: int = 3,
    num_diffn_timesteps: int = 200,
    use_esm_embeddings: bool = True,
    device: str = "cuda:0",
    output_dir: str | None = None,
) -> dict:
    """
    Predict biomolecular structure with Chai-1 (Chai Discovery).

    Chai-1 is a multi-modal foundation model that supports proteins, small
    molecules (SMILES), DNA, RNA, and glycosylations in a single pass.
    Model weights are auto-downloaded on first run — no separate request required.

    Install: pip install chai_lab

    Parameters
    ----------
    protein_sequences : str or list[str]
        Amino-acid sequence(s). Each is a separate protein chain.
    ligand_smiles : str or list[str], optional
        Small-molecule ligand(s) as SMILES. Each becomes a ligand entity.
    dna_sequences : list[str], optional
        Single-stranded DNA sequences.
    rna_sequences : list[str], optional
        RNA sequences.
    seeds : list[int], optional
        Random seeds — one prediction per seed (default: [42]).
        More seeds → more diverse structural candidates.
    num_trunk_recycles : int
        Number of Evoformer recycling iterations (default 3).
    num_diffn_timesteps : int
        Diffusion steps (default 200; reduce to 50 for speed at cost of accuracy).
    use_esm_embeddings : bool
        Whether to augment with ESM-2 protein language model embeddings (default True).
    device : str
        PyTorch device string, e.g. "cuda:0", "cuda:1", "cpu".
    output_dir : str, optional
        Directory for Chai-1 output files. Auto temp dir if not given.

    Returns
    -------
    {
        "top_cif_path":     str,          # highest aggregate-score structure (CIF)
        "all_cif_paths":    list[str],    # all predictions, sorted best-first
        "aggregate_scores": list[float],  # per-prediction aggregate scores (best first)
        "mean_plddt":       float | None, # mean pLDDT of best prediction
        "has_ligand":       bool,
        "has_nucleic_acid": bool,
        "n_chains":         int,
        "output_dir":       str,
        "seeds_used":       list[int],
    }
    """
    try:
        from chai_lab.chai1 import run_inference
    except ImportError:
        raise RuntimeError(
            "chai_lab not installed. Run: pip install chai_lab"
        )

    if isinstance(protein_sequences, str):
        protein_sequences = [protein_sequences]
    if isinstance(ligand_smiles, str):
        ligand_smiles = [ligand_smiles]

    seeds_used = seeds if seeds is not None else [42]

    outdir = Path(output_dir) if output_dir else Path(tempfile.mkdtemp(prefix="chai1_"))
    outdir.mkdir(parents=True, exist_ok=True)

    fasta_path = outdir / "input.fasta"
    fasta_path.write_text(
        _build_fasta(
            protein_sequences=protein_sequences,
            ligand_smiles=ligand_smiles or [],
            dna_sequences=dna_sequences or [],
            rna_sequences=rna_sequences or [],
        )
    )

    all_candidates = []
    for seed in seeds_used:
        seed_outdir = outdir / f"seed_{seed}"
        seed_outdir.mkdir(exist_ok=True)
        # Chai-1 requires an empty output directory per run
        candidates = run_inference(
            fasta_file=fasta_path,
            output_dir=seed_outdir,
            num_trunk_recycles=num_trunk_recycles,
            num_diffn_timesteps=num_diffn_timesteps,
            seed=seed,
            device=device,
            use_esm_embeddings=use_esm_embeddings,
        )
        all_candidates.append(candidates)

    top_cif, all_cifs, agg_scores, mean_plddt = _collect_results(all_candidates)

    has_ligand       = bool(ligand_smiles)
    has_nucleic_acid = bool(dna_sequences or rna_sequences)
    n_chains = (
        len(protein_sequences)
        + len(ligand_smiles or [])
        + len(dna_sequences or [])
        + len(rna_sequences or [])
    )

    return {
        "top_cif_path":     top_cif,
        "all_cif_paths":    all_cifs,
        "aggregate_scores": agg_scores,
        "mean_plddt":       mean_plddt,
        "has_ligand":       has_ligand,
        "has_nucleic_acid": has_nucleic_acid,
        "n_chains":         n_chains,
        "output_dir":       str(outdir),
        "seeds_used":       seeds_used,
    }


# ── FASTA construction ────────────────────────────────────────────────────────

def _build_fasta(
    protein_sequences: list[str],
    ligand_smiles: list[str],
    dna_sequences: list[str],
    rna_sequences: list[str],
) -> str:
    """
    Build a Chai-1 FASTA string.

    Format per chain:
        >protein|name=prot_0       sequence of amino acids
        >ligand|name=lig_0         SMILES string (one-liner)
        >dna|name=dna_0            nucleotide sequence
        >rna|name=rna_0            nucleotide sequence
    """
    lines: list[str] = []

    for i, seq in enumerate(protein_sequences):
        lines.append(f">protein|name=prot_{i}")
        lines.append(seq.strip().upper())

    for i, smi in enumerate(ligand_smiles):
        lines.append(f">ligand|name=lig_{i}")
        lines.append(smi.strip())

    for i, dna in enumerate(dna_sequences):
        lines.append(f">dna|name=dna_{i}")
        lines.append(dna.strip().upper())

    for i, rna in enumerate(rna_sequences):
        lines.append(f">rna|name=rna_{i}")
        lines.append(rna.strip().upper())

    return "\n".join(lines) + "\n"


# ── Result collection ─────────────────────────────────────────────────────────

def _collect_results(
    all_candidates: list,
) -> tuple[str | None, list[str], list[float], float | None]:
    """
    Flatten all per-seed StructureCandidates, sort by aggregate_score, return top.
    """
    combined: list[tuple[float, Path]] = []
    for candidates in all_candidates:
        for cif, rd in zip(candidates.cif_paths, candidates.ranking_data):
            score = float(rd.aggregate_score.item())
            combined.append((score, cif))

    if not combined:
        return None, [], [], None

    combined.sort(key=lambda x: -x[0])  # descending by score
    all_cifs = [str(p) for _, p in combined]
    agg_scores = [s for s, _ in combined]

    # Best prediction pLDDT: use the plddt tensor from the first seed's best candidate
    top_plddt: float | None = None
    for candidates in all_candidates:
        sorted_c = candidates.sorted()
        if sorted_c.cif_paths and hasattr(sorted_c, "plddt"):
            try:
                top_plddt = float(sorted_c.plddt[0].mean().item())
            except Exception:
                pass
        break

    return all_cifs[0] if all_cifs else None, all_cifs, agg_scores, top_plddt
