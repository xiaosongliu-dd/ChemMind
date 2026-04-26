# ChemMind

**An autonomous LLM drug design agent with a three-tier hierarchical skill architecture.**

ChemMind wraps the best open-source computational chemistry tools — Boltz-2, RFantibody, GNINA, OpenFE, ADMETlab3, and 20+ more — behind a unified agentic loop that a computational chemist can drive with plain English.

```
pip install chemmind
chemmind run "Design a selective CDK2 inhibitor with oral bioavailability — start from PDB 1HCL"
```

---

## Architecture

```
L3 — Agent orchestration     Planning skill · Critic skill
         │
L2 — Workflow skills         target_prep · hit_gen_sbdd · lead_opt
         │                   antibody_design · peptide_design · fep_campaign
L1 — Atomic tool skills      boltz2 · rfantibody · gnina · openmm · admetlab3 · ...
```

Skills are markdown files (`skills/L{1,2,3}/**/*.md`) loaded **on demand** — only the tools relevant to the current task enter the context window.

---

## Quickstart

```bash
git clone https://github.com/yourname/chemmind
cd chemmind
uv sync                          # install core deps
export ANTHROPIC_API_KEY="..."
export NVIDIA_API_KEY="..."      # for Boltz-2 NIM

python -c "
from agent.core import ChemMindAgent
agent = ChemMindAgent()
result = agent.run('Find a potent EGFR inhibitor — target is FASTA sequence: MRPSGTAGAALLALLAALCPASRAL...')
print(result.answer)
for mol in result.molecules[:5]:
    print(mol)
"
```

---

## Tool coverage

| Category | Tools |
|---|---|
| Structure & folding | Boltz-2, ESMFold, ColabFold, P2Rank, fpocket |
| Antibody design | RFantibody, AbDiffuser, IgFold, ABodyBuilder3, ProteinMPNN |
| Peptide design | RFdiffusion (peptide), BindCraft |
| Small molecule generation | DiffSBDD, REINVENT 4, MolGPT |
| Docking & scoring | GNINA, AutoDock-GPU, DiffDock |
| MD & FEP | OpenMM, OpenFE, MDAnalysis |
| ADMET | ADMETlab3, RDKit, DeepPurpose |
| Synthesis & data | ASKCOS, ChEMBL, RCSB PDB |

---

## Inspiration

ChemMind's skill architecture is inspired by MolClaw and OpenClaw — three-tier hierarchical skill loading with progressive context injection and `uvx`-isolated tool execution.

---

## License

MIT
