"""ChemMind CLI — entry point for the `chemmind` command."""
from __future__ import annotations

import sys
from pathlib import Path

import click


@click.group()
@click.version_option(version="0.1.0", prog_name="chemmind")
def main() -> None:
    """ChemMind — LLM drug design agent with hierarchical skill architecture."""


@main.command()
@click.argument("task")
@click.option("--model", default="claude-sonnet-4-6", show_default=True,
              help="Anthropic model ID.")
@click.option("--max-steps", default=40, show_default=True,
              help="Maximum ReAct steps before forced stop.")
@click.option("--workdir", default="./chemmind_output", show_default=True,
              help="Working directory for tool file outputs.")
@click.option("--quiet", is_flag=True, default=False,
              help="Suppress step-by-step progress output.")
def run(task: str, model: str, max_steps: int, workdir: str, quiet: bool) -> None:
    """Run TASK through the ChemMind drug design agent."""
    from agent.core import ChemMindAgent, MAX_STEPS

    Path(workdir).mkdir(parents=True, exist_ok=True)
    agent = ChemMindAgent(model=model, verbose=not quiet)
    result = agent.run(task)

    click.echo(result.answer)
    click.echo(
        f"\n[{len(result.steps)} steps | {result.total_time_s:.1f}s | "
        f"skills: {', '.join(result.skills_used)}]",
        err=True,
    )


@main.command()
@click.option("--skills-dir", default=None,
              help="Path to skills/ directory (defaults to repo skills/).")
@click.option("--json", "as_json", is_flag=True, default=False,
              help="Emit results as JSON.")
def validate(skills_dir: str | None, as_json: bool) -> None:
    """Validate all L2 skill files for correct tool parameter names."""
    import json as _json

    from tools.eval.l2_validator import L2Validator, SKILLS_DIR

    target = Path(skills_dir) if skills_dir else SKILLS_DIR
    validator = L2Validator()
    results = validator.validate_all(target)

    errors   = [e for r in results for e in r.errors]
    warnings = [w for r in results for w in r.warnings]

    if as_json:
        click.echo(_json.dumps({
            "n_files":   len(results),
            "n_errors":  len(errors),
            "n_warnings": len(warnings),
            "errors": [
                {"file": e.skill_file.name, "step": e.step_hint, "message": e.message}
                for e in errors
            ],
        }, indent=2))
    else:
        for w in warnings:
            click.echo(f"WARN  {w.skill_file.name}: {w.message}", err=True)
        for e in errors:
            click.echo(f"ERROR {e.skill_file.name} [{e.step_hint}]: {e.message}", err=True)

        n = len(results)
        if errors:
            click.echo(f"\n✗ {len(errors)} error(s) across {n} skill file(s).", err=True)
            sys.exit(1)
        click.echo(f"✓ {n} skill file(s) — 0 errors, {len(warnings)} warning(s).")


@main.command()
@click.option("--query", required=True, help="Target or compound SMILES query.")
@click.option("--db", default="chembl",
              type=click.Choice(["chembl", "pubmed", "bindingdb"]),
              help="Database to query.")
@click.option("--limit", default=20, show_default=True)
def search(query: str, db: str, limit: int) -> None:
    """Search a bioactivity or literature database."""
    if db == "chembl":
        from tools.data.chembl import run_chembl
        result = run_chembl(target_id=query, standard_type="IC50", limit=limit)
        click.echo(f"ChEMBL: {result['n_records']} records for {query}")
        for r in result["records"][:5]:
            click.echo(f"  {r['smiles']}  pChEMBL={r['pchembl_value']}")
    elif db == "pubmed":
        from tools.data.pubmed import run_pubmed
        result = run_pubmed(query=query, max_results=limit)
        click.echo(f"PubMed: {result['n_results']} articles (of {result['total_found']} total)")
        for a in result["articles"][:3]:
            click.echo(f"  [{a['year']}] {a['title'][:80]}")
    elif db == "bindingdb":
        from tools.data.bindingdb import run_bindingdb
        result = run_bindingdb(target_name=query, max_results=limit)
        click.echo(f"BindingDB: {result['n_records']} records for '{query}'")
        for r in result["records"][:5]:
            click.echo(f"  {(r['ligand_smiles'] or '')[:40]}  {r['affinity_nm']} nM")
