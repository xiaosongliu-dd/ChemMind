from __future__ import annotations

import time
import xml.etree.ElementTree as ET
from typing import Any


_ENTREZ_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def run_pubmed(
    query: str,
    max_results: int = 20,
    *,
    sort: str = "relevance",
    api_key: str | None = None,
) -> dict[str, Any]:
    """
    Search PubMed via NCBI Entrez e-utilities and return article metadata + abstracts.

    Parameters
    ----------
    query : str
        PubMed search string, e.g. "CDK2 inhibitor crystal structure 2023".
    max_results : int
        Maximum articles to return (capped at 100).
    sort : str
        "relevance" (default) or "date" (most-recent first).
    api_key : str | None
        NCBI API key — raises rate limit from 3 → 10 req/s.

    Returns
    -------
    {
        articles    : list[{pmid, title, authors, journal, year, abstract, doi, pubmed_url}]
        n_results   : int
        total_found : int
        query       : str
        runtime_s   : float
    }
    """
    try:
        import requests
    except ImportError:
        raise RuntimeError("requests not installed: pip install requests")

    t0 = time.perf_counter()
    max_results = min(max_results, 100)

    common: dict[str, Any] = {"db": "pubmed", "retmode": "json"}
    if api_key:
        common["api_key"] = api_key

    # Step 1 — esearch: get PMIDs
    pmids, total_found = _esearch(requests, common, query, max_results, sort)

    if not pmids:
        return {
            "articles":    [],
            "n_results":   0,
            "total_found": total_found,
            "query":       query,
            "runtime_s":   round(time.perf_counter() - t0, 2),
        }

    # NCBI rate limit: 3 req/s without API key → wait before second call
    time.sleep(0.4)

    # Step 2 — efetch: retrieve full XML records
    articles = _efetch(requests, common, pmids)

    return {
        "articles":    articles,
        "n_results":   len(articles),
        "total_found": total_found,
        "query":       query,
        "runtime_s":   round(time.perf_counter() - t0, 2),
    }


def _esearch(
    requests,
    common: dict,
    query: str,
    max_results: int,
    sort: str,
) -> tuple[list[str], int]:
    params = {**common, "term": query, "retmax": max_results, "sort": sort}
    try:
        resp = requests.get(f"{_ENTREZ_BASE}/esearch.fcgi", params=params, timeout=30)
        resp.raise_for_status()
        data = resp.json()
    except Exception as e:
        raise RuntimeError(f"PubMed esearch failed: {e}")

    result = data.get("esearchresult", {})
    pmids       = result.get("idlist", [])
    total_found = int(result.get("count", 0))
    return pmids, total_found


def _efetch(requests, common: dict, pmids: list[str]) -> list[dict]:
    params = {
        **{k: v for k, v in common.items() if k != "retmode"},
        "id":      ",".join(pmids),
        "rettype": "abstract",
        "retmode": "xml",
    }
    try:
        resp = requests.get(f"{_ENTREZ_BASE}/efetch.fcgi", params=params, timeout=60)
        resp.raise_for_status()
        return _parse_xml(resp.text)
    except Exception as e:
        raise RuntimeError(f"PubMed efetch failed: {e}")


def _parse_xml(xml_text: str) -> list[dict]:
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return []

    articles = []
    for art in root.findall(".//PubmedArticle"):
        pmid = _text(art, ".//PMID") or ""
        title_el = art.find(".//ArticleTitle")
        title = "".join(title_el.itertext()) if title_el is not None else ""
        journal = _text(art, ".//Title") or ""
        year = _text(art, ".//PubDate/Year") or _text(art, ".//PubDate/MedlineDate", "")[:4]

        authors: list[str] = []
        for au in art.findall(".//Author"):
            last  = _text(au, "LastName")  or ""
            first = _text(au, "ForeName") or ""
            name  = f"{last} {first}".strip()
            if name:
                authors.append(name)

        abstract_parts: list[str] = []
        for ab in art.findall(".//AbstractText"):
            label = ab.get("Label", "")
            text  = "".join(ab.itertext())
            abstract_parts.append(f"{label}: {text}" if label else text)
        abstract = " ".join(abstract_parts)

        doi: str | None = None
        for aid in art.findall(".//ArticleId"):
            if aid.get("IdType") == "doi":
                doi = aid.text
                break

        articles.append({
            "pmid":        pmid,
            "title":       title.strip(),
            "authors":     authors[:6],
            "journal":     journal.strip(),
            "year":        year,
            "abstract":    abstract.strip(),
            "doi":         doi,
            "pubmed_url":  f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
        })

    return articles


def _text(el, xpath: str, default: str | None = None) -> str | None:
    found = el.find(xpath)
    return found.text if found is not None and found.text else default


def _find(el, xpath: str):
    return el.find(xpath)
