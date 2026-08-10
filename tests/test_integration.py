"""
tests/test_integration.py
=========================
End-to-end smoke tests for the ChemMindAgent ReAct loop.
All external calls (LLM + tools) are mocked so these tests run offline
and verify that the agent correctly chains: task → tool call → observation
→ final answer.
"""
from __future__ import annotations

import json
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

from agent.core import ChemMindAgent


# ── helpers ───────────────────────────────────────────────────────────────────

IMATINIB = "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cnnc2)n1"


def _make_response(tool_name: str, tool_input: dict):
    """Build a mock Anthropic response that calls one tool then answers."""
    return json.dumps({"action": tool_name, "action_input": tool_input})


def _final_answer(text: str):
    return f"Final Answer: {text}"


def _fake_llm_sequence(*responses):
    """Return a side_effect list for _llm_call that plays through responses."""
    return list(responses)


# ── fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture
def agent(tmp_path):
    return ChemMindAgent(
        model="claude-sonnet-4-6",
        verbose=False,
        memory_persist_dir=str(tmp_path / "memory"),
    )


# ── skill loader integration ──────────────────────────────────────────────────

def test_agent_loads_l3_planning_on_start(agent, tmp_path):
    """L3/planning skill is injected into the first system message."""
    responses = [_final_answer("Done.")]
    with patch.object(agent, "_llm_call", side_effect=responses):
        with patch.object(agent, "_run_critic", return_value=""):
            result = agent.run("test task")
    assert "L3/planning" in result.skills_used


def test_agent_loads_l1_skill_on_first_tool_use(agent):
    """Calling a tool triggers progressive L1 skill injection."""
    tool_result = {"smiles": IMATINIB, "mw": 493.6, "qed": 0.55,
                   "lipinski_pass": True, "sa_score": 3.1}
    responses = [
        _make_response("rdkit_props", {"smiles": IMATINIB}),
        _final_answer("QED is 0.55"),
    ]
    with patch.object(agent, "_llm_call", side_effect=responses):
        with patch.object(agent, "_run_critic", return_value=""):
            with patch("tools.admet.rdkit_props.run_rdkit_props", return_value=tool_result):
                result = agent.run(f"What is the QED of {IMATINIB}?")

    assert any("rdkit_props" in s for s in result.skills_used)


# ── tool dispatch ─────────────────────────────────────────────────────────────

def test_single_tool_call_produces_step(agent):
    tool_result = {"records": [], "n_records": 0}
    responses = [
        _make_response("chembl", {"target_id": "CHEMBL301", "standard_type": "IC50"}),
        _final_answer("No records found."),
    ]
    with patch.object(agent, "_llm_call", side_effect=responses):
        with patch.object(agent, "_run_critic", return_value=""):
            with patch("tools.data.chembl.run_chembl", return_value=tool_result):
                result = agent.run("Get ChEMBL IC50 data for CHEMBL301")

    assert len(result.steps) == 1
    assert result.steps[0].action == "chembl"


def test_multi_step_task_produces_multiple_steps(agent):
    rdkit_result  = {"smiles": IMATINIB, "qed": 0.55, "lipinski_pass": True}
    admet_result  = {"predictions": [{"smiles": IMATINIB, "herg_risk": "Low"}]}
    responses = [
        _make_response("rdkit_props", {"smiles": IMATINIB}),
        _make_response("admetlab3", {"smiles": IMATINIB}),
        _final_answer("Imatinib passes Lipinski and has low hERG risk."),
    ]
    with patch.object(agent, "_llm_call", side_effect=responses):
        with patch.object(agent, "_run_critic", return_value=""):
            with patch("tools.admet.rdkit_props.run_rdkit_props", return_value=rdkit_result):
                with patch("tools.admet.admetlab3.run_admetlab3", return_value=admet_result):
                    result = agent.run("Characterise imatinib properties and ADMET.")

    assert len(result.steps) == 2


def test_unknown_tool_produces_error_observation(agent):
    responses = [
        _make_response("nonexistent_tool_xyz", {}),
        _final_answer("Could not complete: tool not found."),
    ]
    with patch.object(agent, "_llm_call", side_effect=responses):
        with patch.object(agent, "_run_critic", return_value=""):
            result = agent.run("Use a nonexistent tool.")

    assert result.steps[0].success is False
    assert "ERROR" in result.steps[0].observation


# ── tool failure + retry ──────────────────────────────────────────────────────

def test_tool_failure_retries_and_records_error(agent):
    responses = [
        _make_response("rdkit_props", {"smiles": "INVALID"}),
        _final_answer("Tool failed."),
    ]
    with patch.object(agent, "_llm_call", side_effect=responses):
        with patch.object(agent, "_run_critic", return_value=""):
            with patch("tools.admet.rdkit_props.run_rdkit_props",
                       side_effect=ValueError("Invalid SMILES")):
                result = agent.run("Compute props for INVALID smiles.")

    assert result.steps[0].success is False


# ── result schema ─────────────────────────────────────────────────────────────

def test_result_has_required_fields(agent):
    responses = [_final_answer("Done.")]
    with patch.object(agent, "_llm_call", side_effect=responses):
        with patch.object(agent, "_run_critic", return_value=""):
            result = agent.run("Simple task.")

    assert hasattr(result, "answer")
    assert hasattr(result, "steps")
    assert hasattr(result, "skills_used")
    assert hasattr(result, "total_time_s")


def test_total_time_s_is_positive(agent):
    responses = [_final_answer("Done.")]
    with patch.object(agent, "_llm_call", side_effect=responses):
        with patch.object(agent, "_run_critic", return_value=""):
            result = agent.run("Simple task.")

    assert result.total_time_s >= 0.0


def test_final_answer_in_result(agent):
    responses = [_final_answer("42 is the answer.")]
    with patch.object(agent, "_llm_call", side_effect=responses):
        with patch.object(agent, "_run_critic", return_value=""):
            result = agent.run("What is the answer?")

    assert "42 is the answer." in result.answer


# ── max steps guard ───────────────────────────────────────────────────────────

def test_max_steps_reached_returns_partial(agent):
    """If the LLM never sends a final answer, the loop stops at MAX_STEPS."""
    tool_result = {}
    # Always return a tool call — never a final answer
    responses = [
        _make_response("rdkit_props", {"smiles": IMATINIB})
    ] * 50   # more than MAX_STEPS (40)

    with patch.object(agent, "_llm_call", side_effect=responses):
        with patch.object(agent, "_run_critic", return_value=""):
            with patch("tools.admet.rdkit_props.run_rdkit_props", return_value=tool_result):
                result = agent.run("Loop forever.")

    assert "Max steps reached" in result.answer
    assert len(result.steps) == 40
