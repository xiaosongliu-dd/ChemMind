# Skill: admetlab3

**Tier:** L1 — atomic tool skill
**Category:** admet
**Tool:** ADMETlab3 web service (~70 ADMET endpoint predictions)
**GPU required:** No — REST API call

---

## What this skill does

Calls the ADMETlab3 web service (Scientific Computing & Big Data Drug Discovery Lab, Central South University) to predict ~70 ADMET endpoints from SMILES, including absorption, distribution, metabolism, excretion, and toxicity panels — Caco-2, HIA, BBB, P-gp inhibition, CYP isoforms, hERG, AMES, hepatotoxicity, LD50, etc.

This is the heavyweight cousin of `rdkit_props`: instead of cheap rule-based flags, you get neural-network endpoint predictions trained on curated experimental datasets. Use *after* basic Ro5/QED triage when you need high-confidence ADMET signals before committing to a wet-lab campaign.

---

## When to call this skill

- After narrowing a hit list (typically <100 compounds) where wet-lab validation cost justifies the API call
- When property-of-interest is a specific ADMET endpoint not covered by RDKit (e.g. hERG, BBB, CYP)
- For lead-optimisation triage — pick analogues with the best toxicity profile to advance
- For developability triage on biologics-derived small molecules where ADMET is binary go/no-go

---

## Inputs

| Field | Type | Required | Description |
|---|---|---|---|
| `smiles` | str \| list[str] | ✓ | Single SMILES or batch list |
| `api_url` | str | — | ADMETlab3 endpoint (default placeholder — verify before production use) |
| `api_key` | str | — | Bearer token if the endpoint requires auth |
| `timeout` | int | — | HTTP timeout in seconds (default 300) |

⚠ **API URL caveat**: the default `api_url` is a placeholder. ADMETlab3's public REST surface has changed across releases and may require API-key registration. Verify against current ADMETlab3 docs before relying on this in a real workflow.

---

## Outputs

| Field | Type | Description |
|---|---|---|
| `predictions` | list[dict] | One dict of endpoint → predicted value per input SMILES |
| `n_compounds` | int | Length of input list |
| `endpoints_returned` | list[str] | Endpoint names present in the response |
| `runtime_s` | float | Wall-clock time of the API call |
| `api_url` | str | Echo of the URL hit (helps debugging) |

Raises `RuntimeError("ADMETlab3 API unreachable: ...")` on network or HTTP failure.

---

## Implementation

`tools/admet/admetlab3.py` — `run_admetlab3()`

---

## Agent decision rules

- **Don't call on >1000 compounds** — service rate-limits + per-call latency makes this expensive at scale; pre-filter with `rdkit_props` and `filtering/ligand_filter` first
- **Cache responses** — for the same SMILES the prediction is deterministic; cache locally on (smiles, api_url) tuple to avoid re-hitting the API
- **Don't use as primary filter**: ADMETlab3 is one model — confirm hERG / hepatotox flags with at least one orthogonal predictor (e.g. `deeppurpose` Tox21, in-house assays)
- **Check `endpoints_returned`** — the API response shape varies between releases; if the agent expects a specific endpoint name, verify it's present before using
- **Upstream**: `rdkit_props` (cheaper triage) → `admetlab3` (deeper read on survivors)
- **Downstream**: `visualization/ligand_viz` to display property panels alongside structures

---

## Install

```bash
uv add requests  # already a core dependency
```

The ADMETlab3 service is a remote endpoint — no local model installation required. If the service requires a key, set `ADMETLAB3_API_KEY` in env and pass via `api_key=os.environ["ADMETLAB3_API_KEY"]`.

---

## References

- Fu L. et al. *Nucleic Acids Res.* 2024, gkae236 — ADMETlab 3.0
- Web service: https://admetlab3.scbdd.com
- License: free for academic use; commercial use requires permission
