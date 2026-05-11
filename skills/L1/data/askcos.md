# Skill: askcos

**Tier:** L1 — atomic tool skill
**Category:** data
**Tool:** ASKCOS — AI-powered retrosynthetic route planning
**GPU required:** No — REST API call to self-hosted or cloud ASKCOS instance

---

## What this skill does

ASKCOS (Automated System for Knowledge-based Continuous Optimization of Synthesis) predicts retrosynthetic routes to a target molecule using Monte Carlo Tree Search guided by learned reaction templates. Given a SMILES, it returns ranked multi-step routes with estimated buyability of leaf building blocks.

Key capabilities:
- Multi-step retrosynthetic tree search (1–5 steps typical)
- Buyability check against commercial building-block catalogs
- Per-reaction template confidence scores
- Multiple ranked route alternatives

---

## When to call this skill

- Assessing synthetic accessibility of a design candidate before committing to wet-lab
- Finding commercially available building blocks for a novel scaffold
- Comparing synthetic complexity across lead series
- Triage step before FEP or experimental validation

Do NOT use if:
- Need a simple SA score proxy → use `rdkit_props` (`sa_score` field)
- Target molecule is a known drug — the route likely exists in ChEMBL / Reaxys already

---

## Inputs

| Field | Type | Required | Description |
|---|---|---|---|
| `smiles` | str | ✓ | Target molecule SMILES |
| `n_steps` | int | — | Maximum retrosynthetic depth (default 3) |
| `max_branching` | int | — | Max reaction branches per step (default 25) |
| `api_url` | str | — | ASKCOS API base URL; reads `ASKCOS_API_URL` env if None |
| `api_key` | str | — | Bearer token for authenticated endpoints |
| `timeout` | int | — | HTTP timeout in seconds (default 120) |

`api_url` is required; set `ASKCOS_API_URL` environment variable or pass explicitly.

---

## Outputs

| Field | Type | Description |
|---|---|---|
| `routes` | list[dict] | Ranked routes, best score first |
| `n_routes` | int | Number of routes returned |
| `target_smiles` | str | Echo of input |
| `runtime_s` | float | API call duration |

Each route: `{steps: list[dict], n_steps, buyable_leaves: list[dict], n_buyable, overall_score}`.
Each step: `{smiles, depth, is_buyable, reaction_smarts, template_score}`.

Raises `RuntimeError("ASKCOS API unreachable: ...")` on network failure.
Raises `ValueError` if no `api_url` is configured.

---

## Implementation

`tools/data/askcos.py` — `run_askcos()`

POSTs to `{api_url}/api/tree-builder/` with SMILES and search parameters. Parses the returned route tree recursively to collect steps and buyable leaves.

---

## Agent decision rules

- **SA score first**: `rdkit_props.sa_score > 4` is a strong signal to check ASKCOS — don't run ASKCOS on every compound
- **n_steps = 3**: default covers most lead-like molecules; increase to 5 for complex natural-product-like scaffolds
- **Buyability gate**: `n_buyable == n_steps` means all precursors are commercially available — this is the ideal scenario
- **Route ranking**: `overall_score` is a composite; prefer routes with `template_score > 0.5` at each step
- **API setup**: ASKCOS requires a running server (MIT public instance or self-hosted). Set `ASKCOS_API_URL=https://askcos.mit.edu` for the MIT public endpoint (requires free registration)
- **Upstream**: `generation/*` or `enumeration/rdkit_enum` → `askcos` for synthetic feasibility triage
- **Downstream**: experimental synthesis, `binding_affinity/fep_openfe`

---

## Install

```bash
# Self-hosted ASKCOS (Docker):
git clone https://github.com/ASKCOS/ASKCOS.git
cd ASKCOS && docker-compose up
# then set ASKCOS_API_URL=http://localhost:9100

# Or use MIT public endpoint (requires registration):
# https://askcos.mit.edu
```

---

## References

- Coley C.W. et al. *ACS Cent. Sci.* 2019, 5, 1236 — "ASKCOS: a library of algorithms for organic synthesis planning"
- GitHub: https://github.com/ASKCOS/ASKCOS
- License: MIT
