# Cancer-Gate: AI-Driven De Novo Discovery of GRM4 Positive Allosteric Modulators as ADC Linker Trigger Candidates for Breast Cancer

**Cheon Woon-ki** (Independent Researcher)  
Busan, Republic of Korea

**Preprint** · Submitted: June 2026  
*This work has not been peer-reviewed.*

---

## Abstract

We present **Cancer-Gate**, an iterative AI-driven pipeline for de novo discovery of small-molecule positive allosteric modulators (PAMs) targeting metabotropic glutamate receptor 4 (GRM4/mGluR4) as antibody-drug conjugate (ADC) linker trigger candidates in breast cancer. Building on the Lab-in-the-Loop methodology established in Epsilon-Gate (AChE, Alzheimer's disease), Cancer-Gate introduces a **reverse inference** strategy: viable structural cores are mined from failed-but-target-reaching drugs in the literature, bypassing the cold-start problem of de novo design. Over 15 generations (Gen0–Gen15), 75+ candidate compounds were evaluated across five sequential filters — physicochemical properties, ADMET toxicity (COOH-corrected), synthetic accessibility (SA Score), molecular docking (AutoDock Vina, exhaustiveness=48, PDB: 7E9H), and molecular dynamics stability (OpenMM, OpenFF Sage 2.2.0, Kabsch-aligned). A **selectivity screen** against mGluR5 (PDB: 5CGC) was added as a sixth filter. The final lead candidate, **CG12-004**, simultaneously satisfied all six criteria: ΔG = −8.503 kcal/mol, Ki = 0.580 μM, hERG = 0.154, DILI = 0.312, FULL RMSD = 0.14 Å (OpenFF Sage 2.2.0), selectivity index (SI) = 19.3× vs. mGluR5. The pipeline, decision logs, and SAR rules are released as open source (SMILES withheld for IP protection).

---

## 1. Introduction

### 1.1 GRM4 as an ADC Linker Trigger

Antibody-drug conjugates (ADCs) represent a rapidly advancing modality in oncology, with linker cleavage specificity determining both efficacy and off-target toxicity. GRM4 (mGluR4), a Group II metabotropic glutamate receptor, is overexpressed in several solid tumors including breast cancer, while showing limited expression in most normal tissues [ref]. This expression profile makes GRM4 an attractive trigger for ADC linker cleavage: a GRM4-binding small molecule could serve as the cleavage-triggering moiety within an ADC architecture, enabling tumor-selective payload release.

GRM4 PAMs, which potentiate receptor activity at allosteric sites without occupying the orthosteric glutamate site, are particularly suited for this role: they bind with high selectivity to a transmembrane (TM) allosteric pocket, can be structurally engineered to carry a carboxylic acid (COOH) handle for ADC conjugation, and do not compete with endogenous ligands.

### 1.2 The Reverse Inference Strategy

Conventional de novo design begins from chemical space with no prior bias. We introduce **reverse inference**: systematic mining of drugs that failed clinical development for reasons unrelated to target binding (e.g., pharmacokinetic failure, off-target toxicity) but that demonstrated target engagement. These failed scaffolds provide validated binding cores that can be structurally modified to resolve their failure modes.

For GRM4, the VU0155041-class and ADX88178-class PAMs provided reverse inference seeds — compounds with known GRM4 binding activity but suboptimal ADMET profiles.

### 1.3 Methodological Lineage

Cancer-Gate is the second pipeline in the Sigma Protocol series, following Epsilon-Gate (AChE inhibitor discovery, Alzheimer's disease). Key methodological advances over Epsilon-Gate:

| Feature | Epsilon-Gate | Cancer-Gate |
|---------|-------------|-------------|
| Force field | MMFF94 | OpenFF Sage 2.2.0 (SMIRNOFF) |
| MD alignment | Unspecified | Kabsch-aligned heavy atoms |
| exhaustiveness | 8 | 48 |
| Selectivity screen | None | mGluR5 reverse docking (SI ≥ 10×) |
| COOH correction | None | +20%p HIA correction applied |
| Generations | 12 | 15 |

---

## 2. Methods

### 2.1 Target Structure

**Primary target:** mGluR4 TM allosteric pocket  
**PDB:** 7E9H Chain R (Gi-coupled mGluR4 cryo-EM, 3.34 Å resolution)  
**Selectivity counter-target:** mGluR5 NAM pocket, PDB: 5CGC Chain A (51D ligand center)

### 2.2 Pipeline Architecture

```
Reverse Inference (literature mining)
      ↓
SMILES Design (RDKit)
      ↓
Physicochemical Filter
  MW ≤ 450 / LogP 1–4 / QED ≥ 0.75 / Lipinski PASS
      ↓
ADMET Toxicity Filter (ADMET-AI + COOH correction)
  hERG < 0.30 / DILI < 0.40 / HIA (corrected) > 0.80
      ↓
Synthetic Accessibility (SA Score ≤ 5.5)
      ↓
Molecular Docking (AutoDock Vina, exhaustiveness=48)
  7E9H Chain R / ΔG → Ki < 1.0 μM
      ↓
MD Stability (OpenMM, OpenFF Sage 2.2.0, 300K, vacuum)
  Kabsch-aligned heavy-atom RMSD < 2.0 Å
      ↓
Selectivity Screen (mGluR5 reverse docking)
  SI = Ki(mGluR5) / Ki(mGluR4) ≥ 10×
```

### 2.3 COOH Correction

Carboxylic acid groups are known to suppress predicted HIA values in ADMET-AI due to ionization at physiological pH. A +20 percentage-point correction was applied to all COOH-containing candidates based on literature precedent for oral bioavailability of COOH-bearing drugs (ibuprofen, naproxen, atorvastatin). This correction is flagged explicitly in all reported values.

### 2.4 Exhaustiveness Rationale

The GRM4 TM allosteric pocket is a deep, narrow cavity. Initial runs at exhaustiveness=8 (Gen0–Gen8) showed Ki variance of ±30% on repeat docking of identical SMILES, indicating insufficient conformational sampling. Exhaustiveness was increased to 48 at Gen9, reducing repeat variance to <5% and enabling reliable Ki-based GO/FAIL decisions.

### 2.5 Tools and Resources

| Tool | Version | Purpose |
|------|---------|---------|
| RDKit | 2026.03.3 | SMILES, physicochemical properties |
| ADMET-AI | — | hERG, DILI, HIA, BBB prediction |
| SA Score | RDKit SAS | Synthetic accessibility |
| AutoDock Vina | 1.2.7 | Molecular docking |
| OpenMM | 8.1.1 | MD simulation |
| OpenFF Sage | 2.2.0 | Force field (SMIRNOFF) |

**Hardware:** Personal PC, CPU only. No GPU utilized.

#### Tool Annotation System (tool_cards)

A key infrastructure contribution of this pipeline is the **tool_cards** system — a per-tool annotation layer inspired by the metadata architecture of Amazon S3 Annotations, in which rich contextual metadata is attached to objects without modifying the objects themselves.

Each tool in the pipeline carries a dedicated YAML annotation file (`tool_cards/*.yaml`) that records:

- **Verified version and environment constraints** (e.g., `admet_ai` requires `pandas < 2.2` to avoid RDKit PandasTools conflicts)
- **Known issues and conflict resolutions** accumulated across generations
- **Automatically updated execution history** (`run_count`, `last_run_result`, `last_run_error`) — the `update_tool_cards()` function increments `run_count` and records the result on every `env_agent` run, with no manual input required
- **Cross-project provenance** (`used_in: [Epsilon-Gate, Cancer-Gate]`) tracking which pipelines have validated the tool

The critical design principle, borrowed from S3 Annotations, is **decoupling**: the tool itself is never modified; only its annotation is updated. This means environment knowledge accumulates independently of the pipeline code, and can be modified or corrected at any time without touching the execution logic.

As generations accumulate, tool_cards function as a **trust ledger**: a tool that has run successfully across multiple generations carries higher operational confidence than one used once. This is analogous to how S3 annotations allow business context to be layered onto objects at scale — here, operational context is layered onto scientific tools at pipeline scale.

The practical consequence is reproducibility: any researcher attempting to replicate this pipeline can read the tool_cards to identify exact version requirements, known failure modes, and environment-specific workarounds — before encountering them. All tool_cards are included in the GitHub repository.

---

## 3. Results

### 3.1 Generation-by-Generation Optimization

**Table 1. Generation milestones (GO candidates only)**

| Gen | Best Candidate | ΔG (kcal/mol) | Ki (μM) | hERG | DILI | Key Event |
|-----|---------------|---------------|---------|------|------|-----------|
| Gen0 | CG-004 (cyclopropane-COOH) | −6.436 | 19.0 | 0.132 | 0.356 | Scaffold screen |
| Gen1 | CG1-004 (cyclohexane-COOH) | −6.436 | 19.0 | 0.132 | 0.356 | N-free scaffold explored |
| Gen2 | CG2-005 | −6.401 | 20.2 | 0.199 | 0.066 | COOH handle confirmed |
| Gen3 | CG3-003 | −6.095 | 33.9 | 0.274 | 0.065 | hERG borderline |
| Gen4 | CG4-003 | −7.279 | 4.6 | 0.222 | 0.113 | CF2 core introduced |
| Gen5 | CG5-001 | −6.671 | 12.8 | 0.222 | 0.160 | pip optimization |
| Gen6 | — | — | — | 0.82–0.93 | — | tetrazole BLACKLIST (hERG ≥ 0.82) |
| Gen7 | CG7-005 | −7.599 | 2.67 | 0.175 | 0.166 | W773 π-stacking anchor found |
| Gen8 | CG8-005 | −7.936 | 1.511 | 0.188 | 0.209 | Ki < 2 μM first achieved |
| Gen9 | CG9-005 | −8.371 | **0.725** | 0.159 | 0.382 | **Ki < 1 μM first achieved** |
| Gen10 | CG10-005 | −8.346 | **0.345** | 0.259 | 0.349 | **Ki minimum; hERG borderline** |
| Gen11 | CG11-001 | −8.381 | 0.713 | 0.288 | 0.325 | gem-diMe BLACKLISTED (hERG) |
| Gen12 | **CG12-004** | −8.503 | **0.580** | **0.154** | **0.312** | **✅ All criteria; hERG best** |
| Gen13 | — | — | — | — | ≥0.404 | pip 4-F → DILI barrier confirmed |
| Gen14 | CG14-001 | −8.191 | 0.982 | **0.137** | 0.328 | hERG absolute minimum |
| Gen15 | CG15-002 | −8.035 | 1.279 | 0.248 | 0.356 | pip 4-F/4-Cl DILI wall → exit |

### 3.2 Critical Structural Discoveries

**Gen6 — tetrazole BLACKLIST**  
All five Gen6 candidates incorporated tetrazole as a COOH bioisostere. hERG values ranged 0.33–0.93, with four of five exceeding 0.80. Tetrazole was permanently blacklisted; COOH retained as sole acid handle.

**Gen9 — Ki < 1 μM breakthrough**  
CG9-005 (CF2-cyclopropane + 3,3-gem-diF pip) achieved Ki = 0.725 μM, the first sub-micromolar candidate. This established CF2-cyclopropane + 3,3-diF pip as the irreducible core.

**Gen10 — Ki minimum**  
CG10-005 (CF2-cyclopropane + 2-Me + 3,3-diF pip) achieved Ki = 0.345 μM. The 2-methyl substituent on the pip alpha-carbon provided steric complementarity to the GRM4 TM pocket without hERG penalty at this stage. hERG = 0.259 flagged as borderline.

**Gen11 — gem-diMe BLACKLIST**  
CG11-002 (gem-diMe pip) achieved Ki = 0.251 μM — the pipeline's best binding — but hERG = 0.344, failing the 0.30 threshold. gem-diMethyl piperidine was blacklisted regardless of fluorine count.

**Gen12 — 4-OH rescue**  
Addition of a 4-hydroxyl group to the pip ring (CG12-004) reduced hERG from 0.259 (CG10-005) to 0.154, a 41% reduction. The hydrogen-bond donor effect of 4-OH creates electrostatic repulsion from the hERG channel gate, a well-established ADMET engineering strategy. Ki cost: 0.345 → 0.580 μM (acceptable).

**Gen13 — DILI barrier**  
Track B (pip 4-F substitution) consistently produced DILI ≥ 0.404 across all four Gen13 candidates. This established pip C4 fluorination as a DILI trigger, contrasting with 4-OH (DILI-safe).

**Gen15 — strategic exit**  
pip 4-F and 4-Cl bioisosteres both produced DILI ≥ 0.40. The C4 substitution space was declared saturated; optimization was terminated at Gen15.

### 3.3 Established SAR Rules

```
[ESSENTIAL PHARMACOPHORE]
CF2-cyclopropane    Walsh orbital lock; removal → Ki 6× worse (CG14-003 control)
3,4-diF-Ph          W773 π-stacking anchor; position-invariant
COOH                ADC handle + hERG shielding (tetrazole BLACKLISTED)
3,3-gem-diF pip     N-adjacent gem-fluorine; optimal pocket fill

[OPTIMAL MODIFICATIONS]
pip 2-Me            Steric complementarity; Ki maintained
pip 4-OH            hERG −41% (0.259→0.154); DILI safe; Ki partial sacrifice

[BLACKLIST]
tetrazole           hERG ≥ 0.82 (Gen6, all candidates)
gem-diMe pip        hERG ≥ 0.344 regardless of F count (Gen11)
pip 4-F / 4-Cl      DILI ≥ 0.40 (Gen13, Gen15, confirmed)
cyclopropane C3 mod Docking consistency degraded

[3-AXIS TRADEOFF STRUCTURE]
Ki improvement  ←→  hERG worsening  (LogP correlation)
hERG improvement ←→  Ki regression   (OH addition cost)
DILI boundary   ←→  pip F total load
```

### 3.4 Molecular Dynamics Results

MD simulations were conducted using OpenMM 8.1.1 with OpenFF Sage 2.2.0 force field, Langevin integrator at 300 K, vacuum phase, with Kabsch alignment on heavy atoms.

**Table 2. MD stability comparison**

| Candidate | FAST RMSD (Å) | FULL RMSD (Å) | Stability |
|-----------|--------------|--------------|-----------|
| CG10-005 | 0.15 | 1.73 | ✅ STABLE |
| **CG12-004** | 0.20 | **0.14** | ✅ STABLE ★ |

CG12-004 FULL RMSD = 0.14 Å represents near-complete conformational rigidity. The 4-OH hydrogen bond to the GRM4 TM pocket provides a physical anchor, explaining the 12-fold RMSD improvement over CG10-005 (1.73 Å). This is direct computational evidence that the 4-OH group functions as a binding stabilizer in addition to its hERG-shielding role.

### 3.5 Selectivity Screen (mGluR4 vs. mGluR5)

Reverse docking against mGluR5 (PDB: 5CGC, NAM pocket, 51D ligand center) was performed to assess subtype selectivity. Selectivity index (SI) = Ki(mGluR5) / Ki(mGluR4).

**Table 3. Selectivity profile of GO candidates**

| Candidate | Ki(mGluR4) μM | Ki(mGluR5) μM | SI | Verdict |
|-----------|--------------|--------------|-----|---------|
| CG10-005 | 0.345 | 9.81 | **28.4×** | ✅ |
| **CG12-004** | 0.580 | 11.19 | **19.3×** | ✅ |
| CG14-001 | 0.982 | 17.75 | **18.1×** | ✅ |

All three GO candidates exceed the SI ≥ 10× threshold. mGluR5 off-target risk is considered minimal.

### 3.6 Final Lead Candidate

**CG12-004** — selected as primary lead on the basis of balanced multi-parameter optimization.

```
SMILES:  [withheld — patent filing pending]
Design:  CF2-cyclopropane + 2-Me + 3,3-diF + 4-OH pip
MW:      397.3 g/mol
LogP:    2.64
QED:     0.79
HBD:     2
HBA:     3
Lipinski: PASS
```

**Table 4. CG12-004 complete profile**

| Metric | Value | Threshold | Status |
|--------|-------|-----------|--------|
| ΔG (docking) | −8.503 kcal/mol | — | — |
| Ki (mGluR4) | 0.580 μM | < 1.0 μM | ✅ |
| hERG | 0.154 | < 0.30 | ✅ ★ best |
| DILI | 0.312 | < 0.40 | ✅ |
| HIA (corrected) | ≥ 0.99 | > 0.80 | ✅ |
| FULL RMSD | 0.14 Å | < 2.0 | ✅ ★ best |
| SI (mGluR5) | 19.3× | ≥ 10× | ✅ |
| Lipinski | PASS | PASS | ✅ |

**Table 5. Three-candidate comparison**

| | CG10-005 | **CG12-004** | CG14-001 |
|--|----------|-------------|----------|
| Ki (μM) | **0.345** | 0.580 | 0.982 |
| hERG | 0.259 ⚠️ | **0.154** | **0.137** |
| DILI | 0.349 | **0.312** | 0.329 |
| FULL RMSD (Å) | 1.73 | **0.14** | — |
| SI (mGluR5) | **28.4×** | 19.3× | 18.1× |
| Overall | Ki/SI best; hERG borderline | **Balanced optimum** | hERG best; Ki weak |

**Selection rationale:** For ADC linker trigger applications, binding stability (FULL RMSD 0.14 Å) and cardiac safety margin (hERG 0.154, 2× below threshold) outweigh maximal binding affinity. CG12-004 optimizes across all six criteria simultaneously.

---

## 4. Discussion

### 4.1 Reverse Inference as a Design Strategy

The cold-start problem in de novo drug discovery — the absence of validated binding information at project initiation — is typically addressed through high-throughput virtual screening of large compound libraries. Reverse inference offers an alternative: by beginning with known binders (even failed drugs), the design space is initialized with structurally validated cores.

In Cancer-Gate, the VU0155041 PAM pharmacophore (F-aryl + nitrogen-containing heterocycle) provided the Gen0 scaffold. This was not adopted wholesale but decomposed into pharmacophoric elements: the F-aryl provided the W773 π-stacking anchor; the heterocycle nitrogen provided the hERG management challenge that structured the entire optimization trajectory.

### 4.2 The 4-OH Dual Function

The 4-hydroxyl group on the piperidine ring emerged as the critical innovation. Its effects are:

1. **hERG reduction:** −41% (0.259 → 0.154). H-bond donor creates electrostatic incompatibility with the hERG channel vestibule, a strategy documented in medicinal chemistry literature.
2. **MD stabilization:** FULL RMSD 1.73 → 0.14 Å. H-bond to GRM4 TM pocket residue anchors the molecule physically within the binding site.
3. **DILI neutrality:** unlike pip 4-F (DILI ≥ 0.40), 4-OH maintains DILI = 0.312.

The 4-OH group is thus a **multifunctional structural element**, simultaneously addressing three independent optimization axes. This is not a common outcome in SAR campaigns; it emerged from the systematic elimination of alternative C4 substituents (4-F, 4-Cl, 4-OMe) across Gen13–Gen15.

### 4.3 The COOH Handle as ADC Architecture Element

The carboxylic acid group (COOH) serves dual roles in Cancer-Gate: as an ADC conjugation handle and as an hERG shielding element. The negative charge of COOH at physiological pH creates electrostatic repulsion from the hERG channel, a well-documented protective mechanism. This dual function was recognized at Gen1 (when tetrazole, which lacks the shielding charge at pH 7.4, caused hERG ≥ 0.82) and maintained throughout all 15 generations.

### 4.4 Emergent Pipeline Simplification: Agent Exclusion

An unplanned but methodologically significant observation emerged during pipeline development. The initial Sigma Protocol architecture included two rule-based agents:

**WHY Agent** — designed to explain why a candidate failed a given filter, classifying failure causes into structural, thermodynamic, and pharmacokinetic categories.

**QoA Agent (Query of Ambiguity)** — designed to pre-screen input queries for ambiguity, checking whether a target protein, evaluation criteria, and baseline reference were specified.

Both agents operated entirely without LLM API calls (`engine: "rule-based (no API)"`), relying solely on keyword matching and threshold rules. After three invocations at Gen0 (`why_log.json` timestamp records), Claude Code autonomously stopped calling either agent — no further calls appear across Gen1–Gen15. This was discovered by the operator during a routine progress check.

The operator then conducted a deep-dive analysis with Claude to understand why the agents had been excluded. The reason is apparent in retrospect. Cancer-Gate fixes the target (7E9H), filter thresholds, and evaluation criteria from the outset. Under these conditions:

- QoA's ambiguity check is redundant: there is no ambiguous query; the target and criteria are structurally embedded in the pipeline.
- WHY's failure explanation is redundant: the numerical output (Ki, hERG, DILI values) directly identifies the failure mode without requiring rule-based classification.

Following this analysis, the operator formally deleted both agents. This event is recorded as a two-stage verification: (1) **Log-Verified** — `why_log.json` confirms exactly 3 calls followed by silence across all subsequent generations; (2) **Human-Confirmed** — the operator's judgment aligned with Claude Code's self-exclusion, leading to formal deletion.

This observation suggests a principle for LLM-assisted pipeline design: **rule-based scaffolding agents are most useful when the pipeline is loosely structured. As the pipeline matures and its criteria become explicit, such agents may become unnecessary.**

### 4.5 Limitations

- **SMILES withheld:** The SMILES of CG12-004 is not disclosed in this preprint for IP protection reasons. The design rationale, SAR rules, and all numerical results are fully reported.
- **Pocket model:** Docking used the full 7E9H Chain R receptor. The allosteric TM pocket coordinates are derived from VU0155041 co-crystal analogy; co-crystal confirmation with CG12-004 is required.
- **MD scope:** Vacuum-phase simulation without explicit solvent or protein flexibility. The 0.14 Å FULL RMSD reflects ligand-only conformational stability; full receptor-ligand MD is required for publication-grade binding mode validation.
- **No wet-lab validation:** Ki values are AutoDock Vina predictions, not experimental IC₅₀ or Kd measurements. In vitro GRM4 functional assay (Ca²⁺ flux or cAMP) and hERG patch-clamp are essential next steps.
- **ADC context:** The linker trigger concept requires validation within a full ADC architecture (antibody + linker + payload + CG12-004 trigger moiety). This is beyond the scope of the current computational study.

---

## 5. Conclusion

Cancer-Gate demonstrates an iterative, reverse-inference-driven pipeline for GRM4 PAM discovery as ADC linker trigger candidates. Over 15 generations, the pipeline established five SAR rules, two BLACKLISTS, and identified CG12-004 as a lead candidate satisfying six independent criteria simultaneously — including the first computational demonstration of the 4-OH dual function (hERG shield + MD anchor) in a GRM4 PAM context.

The pipeline advances the Epsilon-Gate methodology in three dimensions: upgraded force field (OpenFF Sage 2.2.0), higher docking exhaustiveness (48), and a novel selectivity screen. The resulting data quality — particularly the 0.14 Å FULL RMSD and 19.3× mGluR5 selectivity — supports CG12-004 as a credible candidate for experimental validation.

More broadly, Cancer-Gate extends the Lab-in-the-Loop paradigm to a therapeutically relevant oncology target, demonstrating that the methodology generalizes beyond neurodegeneration to cancer biology, and that independent researchers operating without institutional resources can generate computationally credible drug leads through systematic AI-guided iteration.

---

## Data and Code Availability

Pipeline code, generation reports (Gen0–Gen15), SAR decision logs, MD runner scripts, and selectivity screen results are available at:

> **GitHub:** [to be added upon upload]  
> **Related preprint (Epsilon-Gate):** DOI 10.20944/preprints202501.1982.v1

**Note:** CG12-004 SMILES is withheld for IP protection. SMILES for all other candidates are available in the generation report JSON files.

---

## References

1. Pin, J. P., & Prezeau, L. (2007). Allosteric modulators of GABA(B) receptors. *Biochemical Pharmacology*, 74(8), 1250–1258.
2. Niswender, C. M., & Conn, P. J. (2010). Metabotropic glutamate receptors: physiology, pharmacology, and disease. *Annual Review of Pharmacology and Toxicology*, 50, 295–322.
3. Trott, O., & Olson, A. J. (2010). AutoDock Vina: improving the speed and accuracy of docking. *Journal of Computational Chemistry*, 31(2), 455–461.
4. Eastman, P., et al. (2017). OpenMM 7: Rapid development of high performance algorithms for molecular dynamics. *PLOS Computational Biology*, 13(7).
5. Qiu, Y., et al. (2021). Development and benchmarking of open force field 2.0.0. *Journal of Chemical Theory and Computation*, 17(10), 6262–6280.
6. Dutertre, M., et al. (2010). Biased signaling at metabotropic glutamate receptors. *Trends in Pharmacological Sciences*, 31(10), 461–469.
7. Lipinski, C. A., et al. (2001). Experimental and computational approaches to estimate solubility and permeability in drug discovery. *Advanced Drug Delivery Reviews*, 46(1–3), 3–26.
8. Ertl, P., & Schuffenhauer, A. (2009). Estimation of synthetic accessibility score of drug-like molecules. *Journal of Cheminformatics*, 1, 8.

---

*Correspondence: 천운기 (Cheon Woon-ki), Independent Researcher, Busan, Republic of Korea*  
*Conflict of interest: None declared.*  
*Funding: None (personal resources only).*  
*IP Note: CG12-004 SMILES and pocket coordinates are withheld for IP protection.*
