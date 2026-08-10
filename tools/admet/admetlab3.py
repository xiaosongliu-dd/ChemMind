from __future__ import annotations

import time


def run_admetlab3(
    smiles: str | list[str],
    api_url: str = "https://admetlab3.scbdd.com/server/predict",
    api_key: str | None = None,
    timeout: int = 300,
) -> dict:
    """
    Predict ~70 ADMET endpoints via the ADMETlab3 web service. Sends SMILES
    via POST and parses a JSON response.

    NOTE on `api_url`: the default is a best-effort placeholder. ADMETlab3's
    public REST surface has changed across releases and may require an API
    key. Verify the current endpoint and request schema against the current
    ADMETlab3 documentation before relying on this in a production agent.

    Returns
    -------
    {
        "predictions":        list[dict],   # one dict of endpoints per SMILES
        "n_compounds":        int,
        "endpoints_returned": list[str],
        "runtime_s":          float,
        "api_url":            str,
    }
    """
    try:
        import requests
    except ImportError as e:
        raise RuntimeError(
            f"requests not installed (required for admetlab3): {e}"
        )

    smiles_list = [smiles] if isinstance(smiles, str) else list(smiles)

    headers = {"Content-Type": "application/json"}
    if api_key:
        headers["Authorization"] = f"Bearer {api_key}"

    t0 = time.perf_counter()
    try:
        resp = requests.post(
            api_url,
            json={"smiles": smiles_list},
            headers=headers,
            timeout=timeout,
        )
        resp.raise_for_status()
    except requests.RequestException as e:
        raise RuntimeError(f"ADMETlab3 API unreachable: {e}")

    data = resp.json()
    runtime = round(time.perf_counter() - t0, 1)

    predictions = data.get("predictions") or data.get("results") or []
    endpoints   = sorted(predictions[0].keys()) if predictions else []

    return {
        "predictions":        predictions,
        "n_compounds":        len(smiles_list),
        "endpoints_returned": endpoints,
        "runtime_s":          runtime,
        "api_url":            api_url,
    }
