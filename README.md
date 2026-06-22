# Sigma Protocol — AI-Driven De Novo Drug Discovery Pipeline

**Cheon Woon-ki** (Independent Researcher) · Busan, Republic of Korea

> *"The role of the human researcher is not coding or wet-lab execution, but epistemic navigation: identifying where the system is failing, why, and how to restructure the search space."*

---

## Overview

Sigma Protocol is an iterative, AI-assisted pipeline for de novo small-molecule drug discovery. It operates on a **Lab-in-the-Loop** principle: each generation's failure data is systematically converted into structural rules that constrain the next generation's design space.

**No GPU. No institutional resources. No wet lab. Personal CPU only.**

Three completed projects:

| Project | Target | Generations | Final Lead | Key Result |
|---------|--------|-------------|------------|------------|
| **Epsilon-Gate** | AChE (Alzheimer's) | Gen0–Gen12 | direct_O_pyrrolidine | All 5 criteria passed simultaneously |
| **Cancer-Gate** | GRM4/mGluR4 (Breast cancer ADC) | Gen0–Gen15 | CG12-004 | All 6 criteria passed; FULL RMSD 0.14Å |
| **Battery-Gate** | Solid electrolyte optimization | Gen0–Gen8 | Li6PS4Cl2 | 16.29 mS/cm (matches experimental ~17 mS/cm) |

---

## Core Methodology

### Lab-in-the-Loop

Unlike high-throughput virtual screening, Sigma Protocol treats **failure as signal**. Each generation that fails at a specific filter produces a structural rule that eliminates that failure mode in the next generation.

```
Design (SMILES)
    ↓
Physicochemical Filter (MW, LogP, QED, Lipinski)
    ↓
ADMET Filter (hERG, DILI, HIA, BBB)
    ↓
Synthetic Accessibility (SA Score)
    ↓
Molecular Docking (AutoDock Vina)
    ↓
MD Stability (OpenMM, OpenFF Sage 2.2.0)
    ↓
Selectivity Screen (reverse docking)
    ↓
Failure Analysis → SAR Rules → Next Generation
```

### Reverse Inference

Rather than starting from scratch (cold-start problem), Cancer-Gate mines **known binders that failed for non-binding reasons** (pharmacokinetic failure, off-target toxicity) as structural seeds. This bypasses the cold-start problem while preserving validated binding cores.

### Tool Annotation System (tool_cards)

Inspired by **Amazon S3 Annotations** architecture: each tool carries a dedicated YAML annotation file that records version constraints, known conflicts, and execution history — without modifying the tool itself.

```yaml
# example: tool_cards/admet_ai.yaml
tool: admet_ai
verified_version: "latest"
known_issues:
  - pandas >= 2.2 conflict with RDKit PandasTools → pin pandas < 2.2
run_count: 47          # auto-incremented by update_tool_cards()
last_run_result: "OK"
used_in: [Epsilon-Gate, Cancer-Gate]
```

`run_count` and `last_run_result` are **automatically updated** on every `env_agent` run — no manual input required.

### Autonomous Agent Self-Exclusion (Log-Verified + Human-Confirmed)

During Gen0, WHY Agent and QoA Agent were included in the pipeline.
After 3 invocations (`why_log.json`: timestamps 2026-06-18 18:50, 18:51, 18:59),
Claude Code autonomously stopped calling them — no further calls appear across Gen1–Gen15.

When the researcher noticed this during a progress check, a deep-dive analysis
was conducted with Claude to understand why. Claude's reasoning: the rule-based
agents were pre-constraining its inference space before it could reason from raw data.
Based on this analysis, the researcher formally deleted both agents.

This is recorded as a two-stage verification event:
- **Log-Verified:** `why_log.json` confirms exactly 3 calls, then silence
- **Human-Confirmed:** researcher's judgment aligned with Claude's self-exclusion

> The system identified its own inefficiency before the human did.

---

### Decision Log

Every directional decision is recorded in a structured log:

```
[DECISION-CG-003]
Situation: CG14 — need to isolate 4-OH effect vs. F-reduction effect
Decision:  Include CG14-003 (CHF, no 4-OH) as control
Rationale: Without control, cannot confirm which variable drives Ki change
Result:    CG14-003 Ki=2.07μM → CF2 indispensability confirmed
LLM Insight: Experimental design must precede LLM operation
```

---

## Cancer-Gate Results

**Target:** mGluR4 TM allosteric pocket (PDB: 7E9H) as ADC linker trigger for breast cancer

### Established SAR Rules

```
[ESSENTIAL]
CF2-cyclopropane    → Walsh orbital lock; removal causes 6× Ki regression (CG14-003 control)
3,4-diF-Ph          → W773 π-stacking anchor; position-invariant
COOH                → ADC conjugation handle + hERG shielding
3,3-gem-diF pip     → N-adjacent gem-fluorine; optimal pocket fill

[OPTIMAL MODIFICATIONS]
pip 2-Me            → steric complementarity; Ki maintained
pip 4-OH            → hERG −41%; MD anchor effect (FULL RMSD 1.73→0.14Å)

[BLACKLIST]
tetrazole           → hERG ≥ 0.82 (Gen6, all candidates)
gem-diMe pip        → hERG ≥ 0.344 regardless of F count (Gen11)
pip 4-F / 4-Cl      → DILI ≥ 0.40 (Gen13, Gen15, confirmed)
```

### Final Lead: CG12-004

| Metric | Value | Threshold | Status |
|--------|-------|-----------|--------|
| Ki (mGluR4) | 0.580 μM | < 1.0 μM | ✅ |
| hERG | 0.154 | < 0.30 | ✅ ★ best |
| DILI | 0.312 | < 0.40 | ✅ |
| FULL RMSD (MD) | 0.14 Å | < 2.0 | ✅ ★ best |
| SI (mGluR5) | 19.3× | ≥ 10× | ✅ |
| Lipinski | PASS | PASS | ✅ |

> SMILES withheld for IP protection. Available upon reasonable research inquiry.

### Key Discovery: 4-OH Dual Function

The pip 4-hydroxyl group was introduced to reduce hERG. Unexpectedly, it simultaneously:
1. Reduced hERG by 41% (0.259 → 0.154) — electrostatic shielding of hERG channel
2. Improved MD stability 12× (FULL RMSD 1.73 → 0.14 Å) — H-bond anchor in GRM4 TM pocket

This dual function was not designed; it emerged from the data.

---

## Epsilon-Gate Results

**Target:** AChE catalytic anionic site (PDB: 1EVE) for Alzheimer's disease

**Final Lead:** direct_O_pyrrolidine (Gen11)

| Metric | Value | Status |
|--------|-------|--------|
| ΔG | −4.071 kcal/mol | ✅ |
| RMSD (MD) | 2.28 Å | ✅ |
| hERG | 0.14 | ✅ |
| DILI | 0.025 | ✅ |
| BBB | 0.90 | ✅ |

```
SMILES: CN1CCCC1OC12CC3CC(CC(O)(C3)C1)C2
```

---

## Pipeline Architecture

```
sigma_protocol/
├── cancer_gate/
│   ├── cancer_gate_pipeline.py      # Main pipeline (66KB)
│   ├── cancer_gate_gen*_report.json # Per-generation results (Gen1–Gen15)
│   ├── cancer_gate_md.py            # MD simulation runner
│   ├── validate_gen*.py             # Per-generation validation scripts
│   ├── md/
│   │   └── cancer_gate_md_results.json
│   ├── selectivity/
│   │   └── selectivity_results.json
│   └── docking/7E9H/
├── agents/
│   ├── env_agent.py                 # Environment check + tool_cards auto-update
│   └── orchestrator.py
├── tool_cards/
│   ├── rdkit.yaml
│   ├── admet_ai.yaml
│   ├── autodock_vina.yaml
│   ├── openmm.yaml
│   └── sa_score.yaml
├── logs/
│   └── why_log.json
├── decision_log_template.md         # Structured decision recording
├── FINAL_REPORT.txt
└── run_pipeline.py
```

---

## Tools

| Tool | Version | Role |
|------|---------|------|
| RDKit | 2026.03.3 | SMILES, physicochemical |
| ADMET-AI | latest | hERG, DILI, HIA, BBB |
| AutoDock Vina | 1.2.7 | Docking (exhaustiveness=48) |
| OpenMM | 8.1.1 | MD simulation |
| OpenFF Sage | 2.2.0 | Force field (SMIRNOFF) |

**Hardware:** Personal PC, CPU only.

---

## Research Manuscripts

Full manuscripts for Epsilon-Gate and Cancer-Gate are available in the [`/preprints`](./preprints) directory.

---

## Author

**Cheon Woon-ki (천운기)** · Independent Researcher · Busan, Republic of Korea

No institutional affiliation. No wet lab. No GPU.

Not a developer. All pipelines designed and operated exclusively through natural language prompting. No code written by hand.

---

## License

Code: MIT  
Data (generation reports, SAR rules): CC BY 4.0  
CG12-004 SMILES: withheld for IP protection
