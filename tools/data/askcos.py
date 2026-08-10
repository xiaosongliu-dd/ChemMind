from __future__ import annotations

import os
import time


def run_askcos(
    smiles: str,
    n_steps: int = 3,
    max_branching: int = 25,
    api_url: str | None = None,
    api_key: str | None = None,
    timeout: int = 120,
) -> dict:
    """
    Plan retrosynthetic routes for a target molecule using ASKCOS.

    `api_url` defaults to the ASKCOS_API_URL environment variable.
    ASKCOS is typically self-hosted; the public MIT endpoint requires
    registration. See https://askcos.mit.edu for access.

    Returns:
        {
            routes: list[{steps, buyable_leaves, overall_score}],
            n_routes: int,
            target_smiles: str,
            runtime_s: float,
        }

    Raises RuntimeError("ASKCOS API unreachable: ...") on network failure.
    Raises ValueError if no api_url is configured.
    """
    try:
        import requests
    except ImportError as e:
        raise RuntimeError(f"requests not installed (required for ASKCOS): {e}")

    url = api_url or os.environ.get("ASKCOS_API_URL", "")
    if not url:
        raise ValueError(
            "No ASKCOS API URL configured. Set ASKCOS_API_URL env var or pass api_url=."
        )
    url = url.rstrip("/")

    headers = {"Content-Type": "application/json"}
    if api_key:
        headers["Authorization"] = f"Bearer {api_key}"

    payload = {
        "smiles":        smiles,
        "max_depth":     n_steps,
        "max_branching": max_branching,
    }

    t0 = time.perf_counter()
    try:
        resp = requests.post(
            f"{url}/api/tree-builder/",
            json=payload,
            headers=headers,
            timeout=timeout,
        )
        resp.raise_for_status()
    except requests.RequestException as e:
        raise RuntimeError(f"ASKCOS API unreachable: {e}")

    data = resp.json()
    runtime = round(time.perf_counter() - t0, 1)

    raw_routes = data.get("result") or data.get("routes") or data.get("trees") or []
    routes = [_parse_route(r) for r in raw_routes]
    routes.sort(key=lambda r: -(r["overall_score"] or 0))

    return {
        "routes":        routes,
        "n_routes":      len(routes),
        "target_smiles": smiles,
        "runtime_s":     runtime,
    }


def _parse_route(route: dict) -> dict:
    """Flatten an ASKCOS tree node into a summary dict."""
    steps = _collect_steps(route)
    buyable = [s for s in steps if s.get("is_buyable")]
    score = route.get("score") or route.get("overall_score")

    return {
        "steps":          steps,
        "n_steps":        len(steps),
        "buyable_leaves": buyable,
        "n_buyable":      len(buyable),
        "overall_score":  _safe_float(score),
    }


def _collect_steps(node: dict, depth: int = 0) -> list[dict]:
    """Recursively collect reaction steps from a tree node."""
    steps = []
    smiles = node.get("smiles") or node.get("smi") or ""
    is_buyable = node.get("is_chemical") and node.get("purchase_sources")
    steps.append({
        "smiles":        smiles,
        "depth":         depth,
        "is_buyable":    bool(is_buyable),
        "reaction_smarts": node.get("smarts") or node.get("reaction_smarts"),
        "template_score":  _safe_float(node.get("template_score")),
    })
    for child in node.get("children") or []:
        steps.extend(_collect_steps(child, depth + 1))
    return steps


def _safe_float(val) -> float | None:
    try:
        return float(val)
    except (TypeError, ValueError):
        return None
