"""
ChemMind L2 Skill Validator
===========================
Static analysis of L2 workflow markdowns: checks that every tool name in a
fenced code block is registered in agent/tools.py, and every keyword argument
matches the actual Python function signature.

Usage (CLI):
    python tools/eval/l2_validator.py               # validate all L2 skills
    python tools/eval/l2_validator.py path/to/skill.md  # single file

Returns exit code 1 if any errors are found (suitable for CI).
"""
from __future__ import annotations

import ast
import importlib
import inspect
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path


# ── Paths ─────────────────────────────────────────────────────────────────────

_REPO_ROOT = Path(__file__).parent.parent.parent
SKILLS_DIR = _REPO_ROOT / "skills" / "L2"
_TOOLS_PY  = _REPO_ROOT / "agent" / "tools.py"


# ── Data structures ───────────────────────────────────────────────────────────

@dataclass
class ValidationError:
    severity:   str    # "error" | "warning"
    skill_file: Path
    step_hint:  str    # e.g. "tool `gnina`"
    message:    str


@dataclass
class ValidationResult:
    skill_file: Path
    errors:   list[ValidationError] = field(default_factory=list)
    warnings: list[ValidationError] = field(default_factory=list)

    @property
    def passed(self) -> bool:
        return len(self.errors) == 0


REQUIRED_SECTIONS = [
    "When to invoke",
    "Required inputs",
    "L1 skill sequence",
    "Decision branches",
    "Expected outputs",
    "Success criteria",
]


# ── Registry introspection ────────────────────────────────────────────────────

def _parse_registry(tools_py: Path) -> dict[str, tuple[str, str]]:
    """
    AST-walk agent/tools.py and extract every self._reg(name, module, fn, ...)
    call, returning {tool_name: (module_path, function_name)}.
    """
    tree = ast.parse(tools_py.read_text())
    result: dict[str, tuple[str, str]] = {}
    for node in ast.walk(tree):
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == "_reg"
            and len(node.args) >= 3
        ):
            try:
                name   = ast.literal_eval(node.args[0])
                module = ast.literal_eval(node.args[1])
                fn     = ast.literal_eval(node.args[2])
                result[name] = (module, fn)
            except (ValueError, TypeError):
                pass
    return result


def _build_params(registry: dict[str, tuple[str, str]]) -> dict[str, set[str]]:
    """
    Import each tool module and introspect the function signature to get the
    set of valid parameter names.  Returns empty set for tools whose module
    cannot be imported (e.g. optional GPU dependencies) — param validation is
    skipped for those tools.
    """
    params: dict[str, set[str]] = {}
    for name, (mod_path, fn_name) in registry.items():
        try:
            mod = importlib.import_module(mod_path)
            fn  = getattr(mod, fn_name)
            sig = inspect.signature(fn)
            params[name] = {
                k
                for k, p in sig.parameters.items()
                if p.kind
                not in (
                    inspect.Parameter.VAR_POSITIONAL,
                    inspect.Parameter.VAR_KEYWORD,
                )
            }
        except Exception:
            params[name] = set()  # skip param check for this tool
    return params


# ── Markdown helpers ──────────────────────────────────────────────────────────

_CODE_BLOCK = re.compile(r"```[^\n]*\n(.*?)```", re.DOTALL)
_KWARG      = re.compile(r"\b([a-z_][a-z0-9_]*)\s*=(?!=)")   # word= but not word==


def _extract_code_blocks(text: str) -> list[str]:
    return _CODE_BLOCK.findall(text)


def _extract_call(text: str, open_pos: int) -> str:
    """
    Starting from the '(' at open_pos, scan forward tracking depth until the
    matching ')' is found.  Returns the full substring including both parens.
    """
    depth = 0
    for i in range(open_pos, len(text)):
        ch = text[i]
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth == 0:
                return text[open_pos : i + 1]
    return text[open_pos:]  # unbalanced — return rest


def _kwargs_in_call(call_body: str) -> list[str]:
    """
    Extract keyword argument names from a function call string that starts
    with '('.  Strips the opening paren before searching.
    """
    inner = call_body[1:] if call_body.startswith("(") else call_body
    return _KWARG.findall(inner)


# ── Core validation ───────────────────────────────────────────────────────────

def validate_skill(
    md_path: Path,
    tool_params: dict[str, set[str]],
) -> ValidationResult:
    result = ValidationResult(skill_file=md_path)
    text   = md_path.read_text()
    known  = set(tool_params.keys())

    # 1. Required-section completeness (warnings only)
    for section in REQUIRED_SECTIONS:
        if section not in text:
            result.warnings.append(ValidationError(
                "warning", md_path, "structure",
                f"Missing recommended section: '{section}'"
            ))

    # 2. Validate tool calls inside fenced code blocks
    for block in _extract_code_blocks(text):
        # Process longer names first so e.g. "boltz2_affinity" is matched
        # before the shorter "boltz2" on the same text.
        for tool_name in sorted(known, key=len, reverse=True):
            pattern = re.compile(r"\b" + re.escape(tool_name) + r"\s*\(")
            for m in pattern.finditer(block):
                # m.end()-1 is the position of the '('
                call_body = _extract_call(block, m.end() - 1)
                kwargs    = _kwargs_in_call(call_body)
                params    = tool_params.get(tool_name, set())

                if not params:
                    # Module not importable — can't validate params for this tool
                    continue

                for kw in kwargs:
                    if kw not in params:
                        result.errors.append(ValidationError(
                            "error",
                            md_path,
                            f"tool `{tool_name}`",
                            f"unknown parameter `{kw}=` "
                            f"(valid: {sorted(params)})",
                        ))

    return result


def validate_all(
    skills_dir: Path,
    tool_params: dict[str, set[str]],
) -> list[ValidationResult]:
    return [
        validate_skill(f, tool_params)
        for f in sorted(skills_dir.rglob("*.md"))
    ]


# ── Public class (used by tests) ──────────────────────────────────────────────

class L2Validator:
    """
    Convenience wrapper that builds the tool-params map once and exposes
    validate_skill / validate_all methods.
    """

    def __init__(self, tools_py: Path | None = None) -> None:
        tools_py = tools_py or _TOOLS_PY
        registry       = _parse_registry(tools_py)
        self._params   = _build_params(registry)
        self._registry = registry

    @property
    def registered_tools(self) -> set[str]:
        return set(self._registry.keys())

    def validate_skill(self, md_path: Path) -> ValidationResult:
        return validate_skill(md_path, self._params)

    def validate_all(self, skills_dir: Path = SKILLS_DIR) -> list[ValidationResult]:
        return validate_all(skills_dir, self._params)


# ── CLI ───────────────────────────────────────────────────────────────────────

def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(
        description="Validate ChemMind L2 skill markdowns against the tool registry."
    )
    parser.add_argument(
        "path",
        nargs="?",
        default=str(SKILLS_DIR),
        help="Path to a single .md file or the skills/L2/ directory (default).",
    )
    args = parser.parse_args()

    target   = Path(args.path)
    registry = _parse_registry(_TOOLS_PY)
    params   = _build_params(registry)

    if target.is_file():
        results = [validate_skill(target, params)]
    else:
        results = validate_all(target, params)

    total_errors   = sum(len(r.errors)   for r in results)
    total_warnings = sum(len(r.warnings) for r in results)

    for r in results:
        has_issues = r.errors or r.warnings
        prefix = "  ✗" if r.errors else ("  ⚠" if r.warnings else "  ✓")
        print(f"{prefix}  {r.skill_file.name}")
        for e in r.errors:
            print(f"       ERROR   [{e.step_hint}] {e.message}")
        for w in r.warnings:
            print(f"       WARN    [{w.step_hint}] {w.message}")

    print(
        f"\n{len(results)} file(s) — "
        f"{total_errors} error(s), {total_warnings} warning(s)"
    )
    sys.exit(1 if total_errors else 0)


if __name__ == "__main__":
    sys.path.insert(0, str(_REPO_ROOT))
    main()
