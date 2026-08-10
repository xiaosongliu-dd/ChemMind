# Skill: pubmed

**Tier:** L1 — atomic tool skill
**Category:** data
**Tool:** PubMed — NCBI Entrez e-utilities literature search
**GPU required:** No

---

## What this skill does

Searches PubMed for biomedical literature matching a query string. Returns article
titles, authors, journal names, publication years, and full abstracts. Use to retrieve
target biology context, known SAR from published studies, clinical precedent, or
mechanism-of-action details before designing experiments.

---

## When to invoke

- Before a virtual screen: gather known actives and binding modes from literature
- After a hit is found: look up its published activity profile
- Target characterisation: retrieve crystal structure reports, target validation papers
- Competitive intelligence: what scaffolds have been published for this target?

---

## Signature

```python
pubmed(
    query=<str>,            # PubMed search string
    max_results=20,         # 1–100
    sort="relevance",       # "relevance" or "date"
    api_key=None,           # NCBI API key (optional, raises rate limit)
)
→ {
    articles: [
        {
            pmid:       str,
            title:      str,
            authors:    list[str],   # first 6 authors
            journal:    str,
            year:       str,
            abstract:   str,
            doi:        str | None,
            pubmed_url: str,
        },
        ...
    ],
    n_results:   int,
    total_found: int,       # total PubMed hits (may exceed n_results)
    query:       str,
    runtime_s:   float,
}
```

---

## Query syntax tips

PubMed supports Boolean operators and field tags:
```
CDK2 inhibitor[Title] AND 2020:2024[DP]
EGFR kinase crystal structure[Title/Abstract]
"imatinib resistance"[MeSH Terms]
```

---

## Rate limits

Without API key: 3 requests/second (tool adds 0.4 s delay between esearch and efetch).
With API key: 10 requests/second.
Register at https://www.ncbi.nlm.nih.gov/account/ (free).

---

## References

- NCBI Entrez E-utilities: https://www.ncbi.nlm.nih.gov/books/NBK25501/
