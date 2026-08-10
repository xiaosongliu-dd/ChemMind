"""
tests/test_agent.py
===================
Unit tests for ChemMindAgent core logic:
  - SkillLoader (load, cache, missing, reload, list)
  - MoleculeMemory (fallback mode — no ChromaDB)
  - _flatten_response / _parse_response
  - _execute_tool retry logic
  - agent.run() end-to-end with mocked _llm_call
"""

from __future__ import annotations

import json
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

from agent.core import ChemMindAgent
from agent.memory import MoleculeMemory
from agent.skill_loader import SkillLoader


# ═══════════════════════════════════════════════════════════════════
# SkillLoader
# ═══════════════════════════════════════════════════════════════════

@pytest.fixture
def skill_root(tmp_path):
    """Minimal skills/ tree for SkillLoader tests."""
    (tmp_path / "L1" / "structure").mkdir(parents=True)
    (tmp_path / "L1" / "docking").mkdir(parents=True)
    (tmp_path / "L2").mkdir()
    (tmp_path / "L1" / "structure" / "boltz2.md").write_text("# Boltz-2")
    (tmp_path / "L1" / "docking" / "gnina.md").write_text("# GNINA")
    (tmp_path / "L2" / "lead_opt.md").write_text("# Lead Opt")
    return tmp_path


def test_skill_loader_load_returns_content(skill_root):
    loader = SkillLoader(skill_root)
    assert loader.load("L1/boltz2") == "# Boltz-2"


def test_skill_loader_caches_on_second_call(skill_root):
    loader = SkillLoader(skill_root)
    loader.load("L1/boltz2")
    (skill_root / "L1" / "structure" / "boltz2.md").write_text("# Modified")
    assert loader.load("L1/boltz2") == "# Boltz-2"


def test_skill_loader_returns_none_for_missing(skill_root):
    loader = SkillLoader(skill_root)
    assert loader.load("L1/nonexistent") is None


def test_skill_loader_reload_re_reads_disk(skill_root):
    loader = SkillLoader(skill_root)
    loader.load("L1/boltz2")
    (skill_root / "L1" / "structure" / "boltz2.md").write_text("# Updated")
    assert loader.reload("L1/boltz2") == "# Updated"


def test_skill_loader_list_all(skill_root):
    loader = SkillLoader(skill_root)
    keys = loader.list_available()
    assert "L1/structure/boltz2" in keys
    assert "L1/docking/gnina" in keys
    assert "L2/lead_opt" in keys


def test_skill_loader_list_filtered_by_tier(skill_root):
    loader = SkillLoader(skill_root)
    l1_keys = loader.list_available(tier="L1")
    assert all(k.startswith("L1/") for k in l1_keys)
    assert "L2/lead_opt" not in l1_keys


def test_skill_loader_direct_subdirectory_key(skill_root):
    loader = SkillLoader(skill_root)
    assert loader.load("L1/structure/boltz2") == "# Boltz-2"


# ═══════════════════════════════════════════════════════════════════
# MoleculeMemory — fallback mode (no ChromaDB)
# ═══════════════════════════════════════════════════════════════════

@pytest.fixture
def mem():
    with patch.dict("sys.modules", {"chromadb": None}):
        return MoleculeMemory(persist_dir="/tmp/_test_chemmind_mem")


def test_memory_store_returns_id(mem):
    mol_id = mem.store({"smiles": "CCO"})
    assert isinstance(mol_id, str) and len(mol_id) == 32


def test_memory_store_deduplicates(mem):
    mem.store({"smiles": "CCO"})
    mem.store({"smiles": "CCO"})
    assert mem.count() == 1


def test_memory_store_different_smiles(mem):
    mem.store({"smiles": "CCO"})
    mem.store({"smiles": "CCCO"})
    assert mem.count() == 2


def test_memory_get_recent(mem):
    mem.store({"smiles": "CCO"})
    mem.store({"smiles": "CCCO"})
    assert len(mem.get_recent(n=10)) == 2


def test_memory_search_returns_last_n_fallback(mem):
    for i in range(5):
        mem.store({"smiles": "C" * (i + 1) + "O"})
    assert len(mem.search("CCO", n=3)) == 3


def test_memory_count_starts_at_zero(mem):
    assert mem.count() == 0


def test_memory_count_after_store(mem):
    mem.store({"smiles": "C"})
    assert mem.count() == 1


def test_memory_clear(mem):
    mem.store({"smiles": "CCO"})
    mem.clear()
    assert mem.count() == 0


def test_memory_source_tool_stored(mem):
    mem.store({"smiles": "CCO"}, source_tool="boltz2")
    assert mem.get_recent()[0]["source_tool"] == "boltz2"


# ═══════════════════════════════════════════════════════════════════
# Agent fixture
# ═══════════════════════════════════════════════════════════════════

@pytest.fixture
def agent(tmp_path):
    """ChemMindAgent with real skill/prompt files but ChromaDB mocked out."""
    with patch.dict("sys.modules", {"chromadb": None}):
        a = ChemMindAgent(
            memory_persist_dir=str(tmp_path / ".mem"),
            verbose=False,
        )
    return a


# ═══════════════════════════════════════════════════════════════════
# _flatten_response
# ═══════════════════════════════════════════════════════════════════

def _block(type_, **kw):
    return SimpleNamespace(type=type_, **kw)


def test_flatten_text_only(agent):
    resp = SimpleNamespace(content=[_block("text", text="Hello")])
    assert agent._flatten_response(resp) == "Hello"


def test_flatten_tool_only(agent):
    resp = SimpleNamespace(content=[
        _block("tool_use", name="gnina", input={"ligand": "CCO"})
    ])
    data = json.loads(agent._flatten_response(resp))
    assert data["action"] == "gnina"
    assert data["action_input"] == {"ligand": "CCO"}


def test_flatten_mixed_tool_comes_first(agent):
    """When text and tool_use coexist, tool JSON must be the first line."""
    resp = SimpleNamespace(content=[
        _block("text", text="Let me run GNINA."),
        _block("tool_use", name="gnina", input={"ligand": "CCO"}),
    ])
    first_line = agent._flatten_response(resp).splitlines()[0]
    assert json.loads(first_line)["action"] == "gnina"


def test_flatten_multiple_text_blocks(agent):
    resp = SimpleNamespace(content=[
        _block("text", text="Part 1"),
        _block("text", text="Part 2"),
    ])
    assert agent._flatten_response(resp) == "Part 1\nPart 2"


def test_flatten_empty_content(agent):
    assert agent._flatten_response(SimpleNamespace(content=[])) == ""


# ═══════════════════════════════════════════════════════════════════
# _parse_response
# ═══════════════════════════════════════════════════════════════════

def test_parse_pure_json_tool_call(agent):
    raw = json.dumps({"action": "gnina", "action_input": {"ligand": "CCO"}, "thought": "Docking"})
    parsed = agent._parse_response(raw)
    assert parsed["type"] == "tool_call"
    assert parsed["action"] == "gnina"
    assert parsed["action_input"] == {"ligand": "CCO"}
    assert parsed["thought"] == "Docking"


def test_parse_json_first_line_with_trailing_text(agent):
    """Mixed flatten output: JSON on line 1, prose on line 2 — must detect as tool call."""
    tool_json = json.dumps({"action": "gnina", "action_input": {"ligand": "CCO"}})
    raw = f"{tool_json}\nI will now run GNINA docking."
    parsed = agent._parse_response(raw)
    assert parsed["type"] == "tool_call"
    assert parsed["action"] == "gnina"


def test_parse_final_answer_explicit(agent):
    raw = "Final Answer: The optimal molecule is CC1=CC=CC=C1."
    parsed = agent._parse_response(raw)
    assert parsed["type"] == "final_answer"
    assert "optimal molecule" in parsed["content"]


def test_parse_plain_text_fallback(agent):
    raw = "I need to think about this more carefully."
    parsed = agent._parse_response(raw)
    assert parsed["type"] == "final_answer"
    assert parsed["content"] == raw


def test_parse_empty_string(agent):
    assert agent._parse_response("")["type"] == "final_answer"


def test_parse_missing_thought_defaults_empty(agent):
    raw = json.dumps({"action": "gnina", "action_input": {}})
    assert agent._parse_response(raw)["thought"] == ""


def test_parse_missing_action_input_defaults_empty_dict(agent):
    raw = json.dumps({"action": "gnina"})
    assert agent._parse_response(raw)["action_input"] == {}


def test_parse_final_answer_strips_prefix(agent):
    raw = "Final Answer: Done."
    assert agent._parse_response(raw)["content"] == "Done."


# ═══════════════════════════════════════════════════════════════════
# _execute_tool
# ═══════════════════════════════════════════════════════════════════

def test_execute_tool_success(agent):
    agent.tool_registry.call = MagicMock(return_value={"score": 9.5})
    obs, success, elapsed = agent._execute_tool("gnina", {"ligand": "CCO"})
    assert success is True
    assert json.loads(obs)["score"] == 9.5


def test_execute_tool_retries_on_transient_failure(agent):
    agent.tool_registry.call = MagicMock(
        side_effect=[RuntimeError("transient"), {"score": 7.0}]
    )
    with patch("agent.core.time.sleep"):
        obs, success, _ = agent._execute_tool("gnina", {})
    assert success is True
    assert agent.tool_registry.call.call_count == 2


def test_execute_tool_returns_error_after_max_retries(agent):
    agent.tool_registry.call = MagicMock(side_effect=RuntimeError("always fails"))
    with patch("agent.core.time.sleep"):
        obs, success, _ = agent._execute_tool("gnina", {})
    assert success is False
    assert "ERROR" in obs
    assert agent.tool_registry.call.call_count == 3  # MAX_TOOL_RETRIES


def test_execute_tool_elapsed_is_non_negative(agent):
    agent.tool_registry.call = MagicMock(return_value={})
    _, _, elapsed = agent._execute_tool("gnina", {})
    assert elapsed >= 0.0


# ═══════════════════════════════════════════════════════════════════
# agent.run() end-to-end
# ═══════════════════════════════════════════════════════════════════

def test_run_immediate_final_answer(agent):
    with patch.object(agent, "_llm_call", return_value="Final Answer: Done."), \
         patch.object(agent, "_run_critic", return_value=""):
        result = agent.run("Design an EGFR inhibitor")
    assert "Done." in result.answer
    assert result.steps == []


def test_run_single_tool_then_answer(agent):
    tool_call = json.dumps({"action": "gnina", "action_input": {"ligand": "CCO"}})
    agent.tool_registry.call = MagicMock(return_value={"score": 9.0})
    with patch.object(agent, "_llm_call", side_effect=[tool_call, "Final Answer: Best is CCO."]), \
         patch.object(agent, "_run_critic", return_value=""):
        result = agent.run("Dock CCO to ABL1")
    assert len(result.steps) == 1
    assert result.steps[0].action == "gnina"
    assert result.steps[0].success is True


def test_run_l3_planning_in_skills_used(agent):
    with patch.object(agent, "_llm_call", return_value="Final Answer: done"), \
         patch.object(agent, "_run_critic", return_value=""):
        result = agent.run("test")
    assert "L3/planning" in result.skills_used


def test_run_lazy_l1_skill_loaded(agent):
    tool_call = json.dumps({"action": "gnina", "action_input": {}})
    agent.tool_registry.call = MagicMock(return_value={})
    agent.skill_loader.load = MagicMock(return_value="# skill")
    with patch.object(agent, "_llm_call", side_effect=[tool_call, "Final Answer: done"]), \
         patch.object(agent, "_run_critic", return_value=""):
        result = agent.run("test lazy loading")
    assert "L1/gnina" in result.skills_used


def test_run_max_steps_reached(agent):
    agent.tool_registry.call = MagicMock(return_value={})
    tool_call = json.dumps({"action": "gnina", "action_input": {}})
    with patch.object(agent, "_llm_call", return_value=tool_call), \
         patch.object(agent, "_run_critic", return_value=""), \
         patch("agent.core.MAX_STEPS", 3):
        result = agent.run("infinite loop test")
    assert "Max steps" in result.answer
    assert len(result.steps) == 3


def test_run_total_time_recorded(agent):
    with patch.object(agent, "_llm_call", return_value="Final Answer: done"), \
         patch.object(agent, "_run_critic", return_value=""):
        result = agent.run("quick test")
    assert result.total_time_s >= 0.0


def test_run_failed_tool_recorded_in_steps(agent):
    tool_call = json.dumps({"action": "gnina", "action_input": {}})
    agent.tool_registry.call = MagicMock(side_effect=RuntimeError("fail"))
    with patch.object(agent, "_llm_call", side_effect=[tool_call, "Final Answer: done"]), \
         patch.object(agent, "_run_critic", return_value=""), \
         patch("agent.core.time.sleep"):
        result = agent.run("test error handling")
    assert result.steps[0].success is False
    assert "ERROR" in result.steps[0].observation


def test_run_molecule_persistence(agent):
    tool_call = json.dumps({"action": "gnina", "action_input": {}})
    tool_result = {"molecules": [{"smiles": "CCO", "score": 9.0}]}
    agent.tool_registry.call = MagicMock(return_value=tool_result)
    agent.memory.store = MagicMock(return_value="abc123")
    with patch.object(agent, "_llm_call", side_effect=[tool_call, "Final Answer: done"]), \
         patch.object(agent, "_run_critic", return_value=""):
        agent.run("test persistence")
    agent.memory.store.assert_called_once()
