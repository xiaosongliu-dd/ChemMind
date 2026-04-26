"""
ChemMind — agent/skill_loader.py
=================================
Progressive skill loading: reads L1/L2/L3 markdown files from disk only when
the agent first needs them, keeping the LLM context window lean.

Skill resolution order for a key like "L1/boltz2":
  1. skills/L1/structure/boltz2.md
  2. skills/L1/biologics/boltz2.md
  3. skills/L1/boltz2.md          (flat fallback)
  4. None                          (not found — caller handles gracefully)
"""

from __future__ import annotations

import logging
from functools import lru_cache
from pathlib import Path

logger = logging.getLogger(__name__)


class SkillLoader:
    """
    Loads ChemMind skill markdown files on demand.

    Parameters
    ----------
    root : Path
        Absolute path to the skills/ directory (contains L1/, L2/, L3/).
    """

    def __init__(self, root: Path):
        self.root = Path(root)
        self._cache: dict[str, str] = {}

    def load(self, key: str) -> str | None:
        """
        Load a skill by key and return its markdown content.

        Parameters
        ----------
        key : str
            Skill identifier, e.g. "L1/boltz2", "L2/antibody_design", "L3/planning".

        Returns
        -------
        str | None
            Markdown string, or None if the skill file is not found.
        """
        if key in self._cache:
            return self._cache[key]

        path = self._resolve(key)
        if path is None:
            logger.warning("Skill not found: %s (searched under %s)", key, self.root)
            return None

        content = path.read_text(encoding="utf-8")
        self._cache[key] = content
        logger.debug("Skill loaded: %s (%d chars)", key, len(content))
        return content

    def list_available(self, tier: str | None = None) -> list[str]:
        """
        List all available skill keys.

        Parameters
        ----------
        tier : str | None
            Filter by tier, e.g. "L1", "L2", "L3". None returns all.
        """
        pattern = f"{tier}/**/*.md" if tier else "**/*.md"
        paths   = sorted(self.root.glob(pattern))
        return [self._path_to_key(p) for p in paths]

    def reload(self, key: str) -> str | None:
        """Force re-read from disk (useful during development)."""
        self._cache.pop(key, None)
        return self.load(key)

    # ── Internals ─────────────────────────────────────────────────────────────

    def _resolve(self, key: str) -> Path | None:
        """
        Try to find the markdown file for a skill key.
        Supports both flat (L1/boltz2) and subdirectory (L1/structure/boltz2) keys.
        """
        # If the key already has a full relative path with subdirectory, try direct
        direct = self.root / f"{key}.md"
        if direct.exists():
            return direct

        # Otherwise search recursively under the tier directory
        parts = key.split("/", 1)
        if len(parts) == 2:
            tier, name = parts
            for candidate in (self.root / tier).rglob(f"{name}.md"):
                return candidate

        return None

    def _path_to_key(self, path: Path) -> str:
        """Convert an absolute path back to a skill key string."""
        rel = path.relative_to(self.root)
        return str(rel.with_suffix(""))
