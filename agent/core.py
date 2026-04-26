"""
ChemMind — agent/core.py
========================
ReAct-style agentic loop for LLM-driven drug design.

Architecture
------------
- Loads L1/L2/L3 skill markdown files on demand (progressive context loading)
- Dispatches tool calls via ToolRegistry
- Persists molecule + result memory via ChromaDB
- Recovers from tool failures with bounded retry logic

Usage
-----
    from agent.core import ChemMindAgent
    agent = ChemMindAgent()
    result = agent.run("Design a potent, selective EGFR inhibitor with good oral bioavailability")
"""

from __future__ import annotations

import json
import logging
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from agent.memory import MoleculeMemory
from agent.skill_loader import SkillLoader
from agent.tools import ToolRegistry

logger = logging.getLogger(__name__)

# ── Constants ───────────────────────────────────────────────────────────────────

MAX_STEPS        = 40        # max ReAct steps before forced stop
MAX_TOOL_RETRIES = 3         # retry budget per tool call
SKILLS_ROOT      = Path(__file__).parent.parent / "skills"
PROMPTS_ROOT     = Path(__file__).parent / "prompts"


# ── Data classes ────────────────────────────────────────────────────────────────

@dataclass
class Step:
    """Single ReAct step: thought → action → observation."""
    step_num:     int
    thought:      str
    action:       str
    action_input: dict[str, Any]
    observation:  str
    success:      bool
    elapsed_s:    float


@dataclass
class AgentResult:
    """Final output returned to caller."""
    answer:       str
    steps:        list[Step] = field(default_factory=list)
    molecules:    list[dict] = field(default_factory=list)
    skills_used:  list[str]  = field(default_factory=list)
    total_time_s: float      = 0.0


# ── System prompt ────────────────────────────────────────────────────────────────

def _build_system_prompt(skill_context: str) -> str:
    base = (PROMPTS_ROOT / "system.md").read_text()
    return f"{base}\n\n# Loaded skills\n\n{skill_context}"


# ── Agent ────────────────────────────────────────────────────────────────────────

class ChemMindAgent:
    """
    LLM drug design agent with three-tier hierarchical skill architecture.

    Tier 1 (L1) — atomic tool skills: one wrapper per tool (Boltz-2, GNINA, ...)
    Tier 2 (L2) — workflow skills: compose L1s into domain tasks
    Tier 3 (L3) — orchestration: planning + critic loaded here at init

    Skills are markdown files under skills/L{1,2,3}/**/*.md.
    The agent loads them progressively — only what the current task needs —
    to keep the context window tight.
    """

    def __init__(
        self,
        model: str = "claude-sonnet-4-6",
        temperature: float = 0.2,
        memory_persist_dir: str = ".chemmind_memory",
        verbose: bool = True,
    ):
        self.model         = model
        self.temperature   = temperature
        self.verbose       = verbose

        self.skill_loader  = SkillLoader(SKILLS_ROOT)
        self.tool_registry = ToolRegistry()
        self.memory        = MoleculeMemory(persist_dir=memory_persist_dir)

        self._steps: list[Step]    = []
        self._skills_used: set[str] = set()

        logger.info("ChemMindAgent initialised | model=%s", model)

    # ── Public API ───────────────────────────────────────────────────────────────

    def run(self, task: str) -> AgentResult:
        """
        Run the ReAct loop for a drug design task.

        Parameters
        ----------
        task : str
            Natural language goal, e.g.
            "Optimise this hit against CDK2: CC1=CN=C(N)N=C1 — improve selectivity"

        Returns
        -------
        AgentResult
            Final answer, all intermediate steps, and scored molecule list.
        """
        t0 = time.perf_counter()
        self._steps        = []
        self._skills_used  = set()

        # Always load L3 planning skill first so the agent can decompose the task
        plan_skill = self.skill_loader.load("L3/planning")
        self._skills_used.add("L3/planning")

        messages = self._init_messages(task, plan_skill)

        if self.verbose:
            print(f"\n{'═' * 60}\nChemMind ▸ {task}\n{'═' * 60}")

        # ── ReAct loop ───────────────────────────────────────────────────────────
        answer = "No answer produced."
        for step_num in range(1, MAX_STEPS + 1):
            response = self._llm_call(messages)
            parsed   = self._parse_response(response)

            if parsed["type"] == "final_answer":
                answer = parsed["content"]
                if self.verbose:
                    print(f"\n✓ Final answer after {step_num - 1} steps")
                break

            # ── Tool call ─────────────────────────────────────────────────────────
            tool_name  = parsed["action"]
            tool_input = parsed["action_input"]

            # Lazy-load the matching L1 skill into context on first use
            self._maybe_load_skill(tool_name, messages)

            # Execute with exponential-backoff retry
            observation, success, elapsed = self._execute_tool(tool_name, tool_input)

            step = Step(
                step_num=step_num,
                thought=parsed.get("thought", ""),
                action=tool_name,
                action_input=tool_input,
                observation=observation,
                success=success,
                elapsed_s=elapsed,
            )
            self._steps.append(step)

            if self.verbose:
                status = "✓" if success else "✗"
                print(f"  [{step_num:02d}] {status} {tool_name} ({elapsed:.1f}s)")

            messages.append({"role": "assistant", "content": response})
            messages.append({
                "role": "user",
                "content": f"Observation: {observation}",
            })

            # Persist any molecules that came back from this tool
            self._persist_molecules(tool_name, observation)

        else:
            answer = "Max steps reached. Returning partial results."
            logger.warning("Max steps (%d) reached for task: %s", MAX_STEPS, task)

        # Final critic pass
        critic_eval = self._run_critic(task, answer)

        return AgentResult(
            answer=f"{answer}\n\n{critic_eval}",
            steps=self._steps,
            molecules=self.memory.get_recent(n=20),
            skills_used=sorted(self._skills_used),
            total_time_s=round(time.perf_counter() - t0, 2),
        )

    # ── Internals ────────────────────────────────────────────────────────────────

    def _init_messages(self, task: str, plan_skill: str) -> list[dict]:
        system = _build_system_prompt(plan_skill)
        return [
            {"role": "system", "content": system},
            {"role": "user",   "content": task},
        ]

    def _llm_call(self, messages: list[dict]) -> str:
        """
        Call the LLM. Swap this block for any provider — Anthropic, OpenAI, Ollama.
        Returns a raw string that _parse_response will interpret.
        """
        try:
            import anthropic
            client = anthropic.Anthropic()

            system_msg = next(
                (m["content"] for m in messages if m["role"] == "system"), ""
            )
            user_msgs = [m for m in messages if m["role"] != "system"]

            resp = client.messages.create(
                model=self.model,
                max_tokens=4096,
                temperature=self.temperature,
                system=system_msg,
                messages=user_msgs,
                tools=self.tool_registry.anthropic_tool_specs(),
            )
            return self._flatten_response(resp)

        except ImportError:
            raise RuntimeError(
                "anthropic package not installed. Run: uv add anthropic"
            )

    def _flatten_response(self, resp) -> str:
        """Collapse Anthropic content blocks (text + tool_use) to one string."""
        parts = []
        for block in resp.content:
            if block.type == "text":
                parts.append(block.text)
            elif block.type == "tool_use":
                parts.append(json.dumps({
                    "action":       block.name,
                    "action_input": block.input,
                }))
        return "\n".join(parts)

    def _parse_response(self, response: str) -> dict:
        """
        Parse LLM output into one of two shapes:
          {"type": "tool_call",    "action": str, "action_input": dict, "thought": str}
          {"type": "final_answer", "content": str}
        """
        text = response.strip()

        if text.startswith("{"):
            try:
                data = json.loads(text)
                if "action" in data:
                    return {
                        "type":         "tool_call",
                        "thought":      data.get("thought", ""),
                        "action":       data["action"],
                        "action_input": data.get("action_input", {}),
                    }
            except json.JSONDecodeError:
                pass

        if "Final Answer:" in text:
            return {
                "type":    "final_answer",
                "content": text.split("Final Answer:", 1)[1].strip(),
            }

        return {"type": "final_answer", "content": text}

    def _maybe_load_skill(self, tool_name: str, messages: list[dict]) -> None:
        """
        Progressively inject the L1 skill markdown for a tool into the message
        context the first time that tool is called.
        Keeps context window lean — no unused skill docs loaded.
        """
        skill_key = f"L1/{tool_name}"
        if skill_key not in self._skills_used:
            skill_md = self.skill_loader.load(skill_key)
            if skill_md:
                self._skills_used.add(skill_key)
                messages.append({
                    "role":    "user",
                    "content": f"[Skill loaded: {skill_key}]\n\n{skill_md}",
                })

    def _execute_tool(
        self, tool_name: str, tool_input: dict
    ) -> tuple[str, bool, float]:
        """
        Execute a registered tool with exponential-backoff retry.
        Returns (observation_string, success_bool, elapsed_seconds).
        """
        for attempt in range(1, MAX_TOOL_RETRIES + 1):
            t0 = time.perf_counter()
            try:
                result  = self.tool_registry.call(tool_name, tool_input)
                elapsed = time.perf_counter() - t0
                return json.dumps(result, indent=2), True, round(elapsed, 2)
            except Exception as exc:
                elapsed = time.perf_counter() - t0
                logger.warning(
                    "Tool %s attempt %d/%d failed: %s",
                    tool_name, attempt, MAX_TOOL_RETRIES, exc,
                )
                if attempt == MAX_TOOL_RETRIES:
                    return (
                        f"ERROR after {MAX_TOOL_RETRIES} attempts: {exc}",
                        False,
                        round(elapsed, 2),
                    )
                time.sleep(2 ** attempt)

    def _persist_molecules(self, tool_name: str, observation: str) -> None:
        """Pull any SMILES candidates out of a tool response and store in memory."""
        try:
            data = json.loads(observation)
            mols = data.get("molecules") or data.get("candidates") or []
            for mol in mols:
                if isinstance(mol, dict) and "smiles" in mol:
                    self.memory.store(mol, source_tool=tool_name)
        except (json.JSONDecodeError, AttributeError):
            pass

    def _run_critic(self, task: str, answer: str) -> str:
        """Load the L3 critic skill and evaluate whether the task was completed."""
        critic_skill = self.skill_loader.load("L3/critic")
        self._skills_used.add("L3/critic")
        try:
            import anthropic
            client = anthropic.Anthropic()
            resp = client.messages.create(
                model=self.model,
                max_tokens=1024,
                temperature=0.1,
                system=critic_skill,
                messages=[{
                    "role":    "user",
                    "content": (
                        f"Original task: {task}\n\n"
                        f"Agent answer: {answer}\n\n"
                        "Evaluate: Did the agent complete the task? "
                        "What is missing or should be done next?"
                    ),
                }],
            )
            return "\n\n---\n**Critic evaluation:**\n" + resp.content[0].text
        except Exception as exc:
            logger.warning("Critic evaluation failed: %s", exc)
            return ""
