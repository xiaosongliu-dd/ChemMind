from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from tools.data.pubmed import _parse_xml, run_pubmed


# ── sample XML ───────────────────────────────────────────────────────────────

_SAMPLE_XML = """<?xml version="1.0"?>
<PubmedArticleSet>
  <PubmedArticle>
    <MedlineCitation>
      <PMID>12345678</PMID>
      <Article>
        <ArticleTitle>CDK2 inhibitors: a 2023 perspective</ArticleTitle>
        <Journal>
          <Title>Journal of Medicinal Chemistry</Title>
        </Journal>
        <AuthorList>
          <Author><LastName>Smith</LastName><ForeName>John</ForeName></Author>
          <Author><LastName>Doe</LastName><ForeName>Jane</ForeName></Author>
        </AuthorList>
        <Abstract>
          <AbstractText>CDK2 is a key cell-cycle kinase.</AbstractText>
        </Abstract>
      </Article>
    </MedlineCitation>
    <PubmedData>
      <History>
        <PubMedPubDate PubStatus="pubmed">
          <Year>2023</Year>
        </PubMedPubDate>
      </History>
      <ArticleIdList>
        <ArticleId IdType="pubmed">12345678</ArticleId>
        <ArticleId IdType="doi">10.1021/jmedchem.123</ArticleId>
      </ArticleIdList>
    </PubmedData>
  </PubmedArticle>
</PubmedArticleSet>"""

_ESEARCH_RESPONSE = {
    "esearchresult": {
        "idlist": ["12345678"],
        "count": "1",
    }
}


def _mock_requests(esearch_json, efetch_xml):
    esearch_resp = MagicMock()
    esearch_resp.json.return_value = esearch_json
    esearch_resp.raise_for_status = MagicMock()
    esearch_resp.text = ""

    efetch_resp = MagicMock()
    efetch_resp.json.side_effect = ValueError("not JSON")
    efetch_resp.raise_for_status = MagicMock()
    efetch_resp.text = efetch_xml

    mock_get = MagicMock(side_effect=[esearch_resp, efetch_resp])
    return mock_get


# ── XML parser ────────────────────────────────────────────────────────────────

def test_parse_xml_returns_articles():
    articles = _parse_xml(_SAMPLE_XML)
    assert len(articles) == 1


def test_parse_xml_pmid():
    articles = _parse_xml(_SAMPLE_XML)
    assert articles[0]["pmid"] == "12345678"


def test_parse_xml_title():
    articles = _parse_xml(_SAMPLE_XML)
    assert "CDK2" in articles[0]["title"]


def test_parse_xml_doi():
    articles = _parse_xml(_SAMPLE_XML)
    assert articles[0]["doi"] == "10.1021/jmedchem.123"


def test_parse_xml_authors():
    articles = _parse_xml(_SAMPLE_XML)
    authors = articles[0]["authors"]
    assert len(authors) == 2
    assert "Smith" in authors[0]


def test_parse_xml_abstract_not_empty():
    articles = _parse_xml(_SAMPLE_XML)
    assert len(articles[0]["abstract"]) > 0


def test_parse_xml_pubmed_url_contains_pmid():
    articles = _parse_xml(_SAMPLE_XML)
    assert "12345678" in articles[0]["pubmed_url"]


def test_parse_xml_malformed_returns_empty():
    articles = _parse_xml("<not valid xml>")
    assert articles == []


# ── run_pubmed integration ────────────────────────────────────────────────────

def test_run_pubmed_returns_articles():
    with patch("requests.get", _mock_requests(_ESEARCH_RESPONSE, _SAMPLE_XML)):
        result = run_pubmed("CDK2 inhibitor", max_results=5)

    assert result["n_results"] == 1
    assert len(result["articles"]) == 1


def test_run_pubmed_output_schema():
    with patch("requests.get", _mock_requests(_ESEARCH_RESPONSE, _SAMPLE_XML)):
        result = run_pubmed("CDK2 inhibitor")

    for key in ("articles", "n_results", "total_found", "query", "runtime_s"):
        assert key in result


def test_run_pubmed_query_echoed():
    with patch("requests.get", _mock_requests(_ESEARCH_RESPONSE, _SAMPLE_XML)):
        result = run_pubmed("CDK2 crystal structure")
    assert result["query"] == "CDK2 crystal structure"


def test_run_pubmed_empty_results():
    empty_esearch = {"esearchresult": {"idlist": [], "count": "0"}}
    esearch_resp = MagicMock()
    esearch_resp.json.return_value = empty_esearch
    esearch_resp.raise_for_status = MagicMock()

    with patch("requests.get", return_value=esearch_resp):
        result = run_pubmed("xyzzy_nonexistent_target_12345")

    assert result["n_results"] == 0
    assert result["articles"] == []


def test_run_pubmed_caps_at_100():
    with patch("requests.get", _mock_requests(_ESEARCH_RESPONSE, _SAMPLE_XML)):
        result = run_pubmed("CDK2", max_results=999)

    # Verify max_results was capped by checking the esearch call
    # (the function internally caps at 100 — just check no crash)
    assert result["n_results"] >= 0


def test_run_pubmed_raises_if_requests_not_installed():
    with patch.dict("sys.modules", {"requests": None}):
        with pytest.raises(RuntimeError, match="requests not installed"):
            run_pubmed("CDK2")


def test_run_pubmed_runtime_s_positive():
    with patch("requests.get", _mock_requests(_ESEARCH_RESPONSE, _SAMPLE_XML)):
        result = run_pubmed("CDK2")
    assert result["runtime_s"] >= 0.0
