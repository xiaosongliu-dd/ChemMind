"""
ChemMind — agent/memory.py
===========================
Persistent molecule memory backed by ChromaDB.

Stores scored molecules with their SMILES, properties, and provenance
so the agent can retrieve similar hits across iterations and sessions.

Usage
-----
    mem = MoleculeMemory()
    mem.store({"smiles": "CCO", "affinity_kcal": -9.2, "qed": 0.81}, source_tool="boltz2")
    hits = mem.search("CCO", n=5)           # similarity search
    recent = mem.get_recent(n=20)           # last 20 stored
"""

from __future__ import annotations

import hashlib
import json
import logging
import time
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

_COLLECTION_NAME = "chemmind_molecules"


class MoleculeMemory:
    """
    Molecule store with semantic similarity search via ChromaDB.

    Falls back to an in-memory list if ChromaDB is not installed,
    so the agent runs without errors even before the full environment
    is set up.
    """

    def __init__(self, persist_dir: str = ".chemmind_memory"):
        self._persist_dir   = persist_dir
        self._client        = None
        self._collection    = None
        self._fallback: list[dict] = []     # used when ChromaDB unavailable
        self._init_chromadb()

    # ── Public API ────────────────────────────────────────────────────────────

    def store(self, molecule: dict[str, Any], source_tool: str = "unknown") -> str:
        """
        Persist a molecule record.

        Parameters
        ----------
        molecule : dict
            Must contain "smiles". Other keys (affinity, qed, sa_score, …) are stored
            as metadata and returned on retrieval.
        source_tool : str
            Which ChemMind tool produced this molecule.

        Returns
        -------
        str
            Unique molecule ID (SHA-256 of SMILES).
        """
        smiles = molecule.get("smiles", "")
        mol_id = _smiles_id(smiles)

        record = {
            **molecule,
            "mol_id":      mol_id,
            "source_tool": source_tool,
            "timestamp":   time.time(),
        }

        if self._collection is not None:
            try:
                self._collection.upsert(
                    ids=[mol_id],
                    documents=[smiles],
                    metadatas=[{k: str(v) for k, v in record.items()}],
                )
                logger.debug("Stored molecule %s via ChromaDB", mol_id[:8])
                return mol_id
            except Exception as exc:
                logger.warning("ChromaDB store failed, using fallback: %s", exc)

        # Fallback: in-memory list (deduped by mol_id)
        if not any(m.get("mol_id") == mol_id for m in self._fallback):
            self._fallback.append(record)
        return mol_id

    def search(self, query_smiles: str, n: int = 10) -> list[dict]:
        """
        Retrieve the n most similar molecules to a query SMILES.
        Uses ChromaDB cosine similarity on SMILES strings.
        Falls back to last-n if ChromaDB is unavailable.
        """
        if self._collection is not None:
            try:
                results = self._collection.query(
                    query_texts=[query_smiles],
                    n_results=min(n, self._collection.count()),
                )
                metas = results.get("metadatas", [[]])[0]
                return [_parse_meta(m) for m in metas]
            except Exception as exc:
                logger.warning("ChromaDB search failed: %s", exc)

        return self._fallback[-n:]

    def get_recent(self, n: int = 20) -> list[dict]:
        """Return the n most recently stored molecules."""
        if self._collection is not None:
            try:
                total = self._collection.count()
                if total == 0:
                    return []
                results = self._collection.get(
                    limit=min(n, total),
                    include=["metadatas"],
                )
                metas = results.get("metadatas", [])
                return [_parse_meta(m) for m in metas]
            except Exception as exc:
                logger.warning("ChromaDB get_recent failed: %s", exc)

        return self._fallback[-n:]

    def count(self) -> int:
        """Total molecules stored."""
        if self._collection is not None:
            try:
                return self._collection.count()
            except Exception:
                pass
        return len(self._fallback)

    def clear(self) -> None:
        """Wipe all stored molecules (useful between campaigns)."""
        if self._collection is not None:
            try:
                self._collection.delete(where={"source_tool": {"$ne": "__never__"}})
            except Exception as exc:
                logger.warning("ChromaDB clear failed: %s", exc)
        self._fallback.clear()

    # ── Internals ─────────────────────────────────────────────────────────────

    def _init_chromadb(self) -> None:
        try:
            import chromadb
            Path(self._persist_dir).mkdir(parents=True, exist_ok=True)
            self._client     = chromadb.PersistentClient(path=self._persist_dir)
            self._collection = self._client.get_or_create_collection(
                name=_COLLECTION_NAME,
                metadata={"hnsw:space": "cosine"},
            )
            logger.info(
                "MoleculeMemory: ChromaDB ready at %s (%d molecules)",
                self._persist_dir, self._collection.count(),
            )
        except ImportError:
            logger.warning(
                "chromadb not installed — using in-memory fallback. "
                "Run: uv add chromadb"
            )
        except Exception as exc:
            logger.warning("ChromaDB init failed (%s) — using in-memory fallback", exc)


# ── Helpers ───────────────────────────────────────────────────────────────────

def _smiles_id(smiles: str) -> str:
    """Deterministic SHA-256 ID for a SMILES string."""
    return hashlib.sha256(smiles.encode()).hexdigest()[:32]


def _parse_meta(meta: dict) -> dict:
    """Try to coerce numeric string values back to floats/ints."""
    out = {}
    for k, v in meta.items():
        try:
            out[k] = int(v) if "." not in v else float(v)
        except (ValueError, TypeError):
            out[k] = v
    return out
