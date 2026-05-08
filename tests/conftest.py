from __future__ import annotations

import csv
from pathlib import Path
from unittest.mock import MagicMock

import pytest

FIXTURE_DIR = Path(__file__).parent / "test_target_folder"

# ABL1 kinase domain fragment — short enough to be fast, representative enough
# to exercise FASTA-accepting tools
ABL1_FASTA = (
    "MGQTYITLPSGKELSLFCNLREVLAQDLPGKHSWLQPKAAIPQLLPSIMDNTQFDDQHIK"
    "ELIPQTIPKPVFQKLEELATPPQQLQTCPNAQQKAMSLQKQNLQNIQMQTLHQLLNQQQTQ"
)


@pytest.fixture(scope="session")
def inhibitor_smiles() -> dict[str, str]:
    """All ABL1 inhibitors from the fixture CSV: {name: smiles}."""
    path = FIXTURE_DIR / "inhibitor_info.csv"
    with open(path) as fh:
        return {row["name"]: row["smiles"] for row in csv.DictReader(fh)}


@pytest.fixture(scope="session")
def imatinib(inhibitor_smiles) -> str:
    return inhibitor_smiles["Imatinib"]


@pytest.fixture(scope="session")
def all_smiles(inhibitor_smiles) -> list[str]:
    return list(inhibitor_smiles.values())


@pytest.fixture(scope="session")
def abl1_fasta() -> str:
    return ABL1_FASTA


@pytest.fixture
def mock_receptor_pdb(tmp_path) -> str:
    """Minimal valid PDB — enough for path/read tests."""
    pdb = tmp_path / "abl1_receptor.pdb"
    pdb.write_text(
        "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00  0.00           C\n"
        "ATOM      2  CA  GLY A   2       4.000   5.000   6.000  1.00  0.00           C\n"
        "ATOM      3  CA  LYS A   3       7.000   8.000   9.000  1.00  0.00           C\n"
        "END\n"
    )
    return str(pdb)


@pytest.fixture
def mock_nim_response():
    """Factory for fake NIM HTTP responses."""
    def _make(payload: dict):
        resp = MagicMock()
        resp.json.return_value = payload
        resp.raise_for_status = MagicMock()
        return resp
    return _make
