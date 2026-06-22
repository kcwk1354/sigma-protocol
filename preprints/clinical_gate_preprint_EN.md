# Structural Patterns of Clinical Trial Termination: A Token-Based Analysis of 10,000 Terminated Interventional Studies from ClinicalTrials.gov

**Woonki Cheon (천운기)**  
Independent Researcher, Busan, Republic of Korea  
Sigma Protocol Project

*Preprint — Not peer reviewed*  
*Date: June 21, 2026*

---

## Abstract

Clinical trial failure accounts for approximately 90% of drug development costs, yet the structural patterns underlying termination decisions remain poorly characterized in machine-readable form. We developed Clinical-Gate, a natural language tokenization pipeline that transforms unstructured termination descriptions from ClinicalTrials.gov into structured failure tokens. Applying this pipeline to 10,000 terminated interventional studies, we identified seven failure categories: ENROLLMENT, BUSINESS, EXTERNAL, SAFETY, EFFICACY, REGULATORY, and UNKNOWN. Among classifiable cases (n=6,471), ENROLLMENT failure (29.8%) was the dominant termination cause — nearly seven times more frequent than EFFICACY failure (4.3%). Phase-stratified analysis revealed three counter-intuitive findings: (1) PHASE4 studies exhibited the highest ENROLLMENT failure rate (57.3%), exceeding PHASE1 (32.2%) by 25 percentage points; (2) BUSINESS termination in PHASE1 (33.4%) matched the rate of SAFETY termination (15.4%), indicating that sponsor withdrawal is as consequential as toxicity in early-phase trials; and (3) mean trial duration at termination converged across all failure types (32–36 months), suggesting a universal "year-three crisis" independent of failure cause. These findings challenge conventional assumptions about clinical trial failure and provide a structured dataset for predictive modeling and reverse inference in drug design.

**Keywords:** clinical trial failure, tokenization, ClinicalTrials.gov, drug development, failure pattern analysis, Lab-in-the-Loop

---

## 1. Introduction

Drug development is characterized by high attrition. Of compounds entering Phase 1, fewer than 10% ultimately receive regulatory approval. The cost of this attrition — estimated at $1–2 billion per approved drug — is driven primarily by late-stage failures that could theoretically have been anticipated earlier in the development process.

Despite the scale of this problem, the structural patterns of clinical trial termination remain largely uncharacterized in a form amenable to computational analysis. ClinicalTrials.gov, the world's largest clinical trial registry with over 500,000 registered studies, contains termination descriptions in free-text form. These descriptions — "Lack of efficacy at interim analysis," "Insufficient enrollment," "Business decision" — encode rich information about why trials fail, but are not structured for pattern recognition or predictive modeling.

The challenge is not data availability but data form. The information exists; it has simply not been translated into the token space that language models and machine learning systems can process.

This paper introduces Clinical-Gate, a pipeline that performs this translation. By applying rule-based keyword classification to the WhyStopped field of terminated interventional studies, we convert unstructured termination descriptions into categorical failure tokens. We then apply this pipeline to 10,000 terminated studies and report the first large-scale structural analysis of clinical trial failure patterns stratified by trial phase, sponsor type, and trial duration.

Our analysis is motivated by two downstream applications: (1) predictive modeling of trial success probability given mechanism, indication, and phase; and (2) reverse inference — extracting viable mechanisms from failed trials to inform drug design, consistent with the methodology developed in prior Sigma Protocol projects (Epsilon-Gate, Cancer-Gate).

---

## 2. Methods

### 2.1 Data Source

Data were obtained from ClinicalTrials.gov via the version 2 REST API (https://clinicaltrials.gov/api/v2/studies). All records with `filter.overallStatus=TERMINATED` and study type INTERVENTIONAL were retrieved using paginated requests (pageSize=1000) with a 0.2-second inter-request delay to comply with rate limits. Data collection was performed on June 21, 2026.

Fields retrieved: NCTId, Phase, OverallStatus, WhyStopped, PrimaryOutcomeMeasure, Condition, InterventionName, InterventionType, LeadSponsorClass, EnrollmentCount, EnrollmentType, StartDate, CompletionDate.

### 2.2 Tokenization Pipeline

We developed a three-stage tokenization pipeline implemented as the Clinical-Gate tool suite.

**Stage 1: Data Fetch (clinicaltrials_api)**  
Raw study data retrieved from ClinicalTrials.gov API v2 and stored as JSON.

**Stage 2: Text Normalization (text_normalizer)**  
The WhyStopped field was classified into one of seven failure token categories using rule-based keyword matching (case-insensitive). Classification was iteratively refined across four generations (Gen0–Gen3) using a Lab-in-the-Loop approach: each generation's UNKNOWN residuals were analyzed for classifiable patterns, and keyword rules were expanded accordingly.

The seven failure token categories and their defining keywords are:

- **EFFICACY**: efficacy, endpoint, did not meet, insufficient effect, no benefit, MTD, insufficient clinical response, futility, stopped for futility, data did not support continuation
- **SAFETY**: safety, adverse, toxicity, death, serious, side effect, device deficiency, device failure
- **ENROLLMENT**: enrollment, accrual, recruitment, slow accrual, difficult to enroll, inability to recruit, lack of inclusion, unable to recruit, slow inclusion rates
- **REGULATORY**: FDA, IND, regulatory, hold, clinical hold, IRB, approval expired, compliance, protocol compliance, protocol deviation
- **BUSINESS**: funding, sponsor, business, financial, lack of funding, finance, company management, manufacturer, strategic decision, costs of study, decided not to pursue further development
- **EXTERNAL**: pandemic, covid, investigator, site closure, drug shortage, product expired, understaffing, competing studies, study discontinued, study halted, replaced with another study, production halt
- **UNKNOWN**: WhyStopped field absent, blank, or unclassifiable by above rules

A confidence score was assigned to each classification (single keyword match: 0.6; two or more matches: 0.8; no match: 0.0). Records with confidence below 0.5 were reassigned to UNKNOWN.

**Stage 3: Feature Tokenization**  
Additional structured tokens were derived: phase_token (PHASE1–PHASE4, NA), intervention_type (DRUG, BIOLOGICAL, DEVICE, PROCEDURE, OTHER), endpoint_category (SURVIVAL, RESPONSE_RATE, BIOMARKER, SAFETY, QOL, OTHER), sponsor_type (INDUSTRY, NIH, OTHER_GOV, ACADEMIC, OTHER), enrollment_ratio (actual/anticipated enrollment), and duration_months (CompletionDate minus StartDate).

### 2.3 Pipeline Validation

The tokenization pipeline was validated through a preflight analysis (n=50 TERMINATED studies) prior to full-scale processing. WhyStopped field presence was confirmed in 90% of TERMINATED records (null rate: 10%), establishing the feasibility of keyword-based classification. UNKNOWN rate converged to approximately 35–37% across all batches after Gen3 keyword refinement, consistent with an irreducible floor of null/blank records (~12%) and informationally neutral statements (~20–25%).

### 2.4 Analysis

Cross-tabulation of phase_token by failure_token was performed on the 6,471 classifiable records (UNKNOWN excluded). Sponsor type analysis compared INDUSTRY versus non-INDUSTRY sponsors. Duration analysis computed mean and median months to termination by failure category, excluding records with missing date fields.

---

## 3. Results

### 3.1 Dataset Overview

Of 10,000 terminated interventional studies retrieved, 6,471 (64.7%) were classified into one of six failure token categories. The remaining 3,529 (35.3%) were assigned UNKNOWN, attributable to absent WhyStopped entries (approximately 12%) and informationally neutral termination statements (approximately 23%).

Among classifiable records, ENROLLMENT was the dominant failure category (n=2,976; 29.8% of total), followed by BUSINESS (n=1,406; 14.1%), EXTERNAL (n=863; 8.6%), SAFETY (n=505; 5.1%), EFFICACY (n=427; 4.3%), and REGULATORY (n=294; 2.9%).

### 3.2 Phase-Stratified Failure Patterns

Phase-stratified analysis revealed distinct failure signatures by trial phase (Table 1).

**Table 1. failure_token distribution by phase_token (% within phase, classifiable records only)**

| | ENROLLMENT | BUSINESS | EXTERNAL | SAFETY | EFFICACY | REGULATORY | n |
|---------|-----------|---------|---------|--------|---------|-----------|---|
| PHASE1 | 32.2% | 33.4% | 8.6% | 15.4% | 6.1% | 4.4% | 1,121 |
| PHASE2 | 48.7% | 20.7% | 7.1% | 9.6% | 10.5% | 3.4% | 1,625 |
| PHASE3 | 40.0% | 19.4% | 6.5% | 11.7% | 14.9% | 7.5% | 643 |
| PHASE4 | 57.3% | 17.9% | 14.3% | 2.6% | 4.4% | 3.6% | 504 |
| NA | 49.6% | 18.7% | 20.8% | 3.4% | 2.7% | 4.8% | 2,578 |

**EFFICACY failure** increased monotonically from PHASE1 (6.1%) to PHASE3 (14.9%), consistent with the expectation that efficacy endpoints are most rigorously tested in late-phase trials. PHASE4 showed a marked decrease (4.4%), reflecting that post-marketing studies involve drugs with established efficacy.

**SAFETY failure** was highest in PHASE1 (15.4%) and lowest in PHASE4 (2.6%), confirming the role of early-phase trials as toxicity screening stages. The PHASE1 SAFETY rate was 5.9 times higher than PHASE4.

**ENROLLMENT failure** was highest in PHASE4 (57.3%), substantially exceeding PHASE1 (32.2%) and PHASE3 (40.0%). This pattern is discussed in Section 4.

**BUSINESS failure** was highest in PHASE1 (33.4%), exceeding PHASE2 (20.7%) by 13 percentage points.

### 3.3 Sponsor Type Analysis

INDUSTRY-sponsored trials (n=2,219) exhibited a BUSINESS failure rate of 33.2%, compared to 5.0–16.2% for non-INDUSTRY sponsors (Table 2). Conversely, non-INDUSTRY sponsors showed ENROLLMENT failure rates of 47.9–54.8%, compared to 29.3% for INDUSTRY.

**Table 2. failure_token distribution by sponsor_type (selected categories)**

| Sponsor | n | BUSINESS | ENROLLMENT | EFFICACY |
|---------|---|---------|-----------|---------|
| INDUSTRY | 2,219 | 33.2% | 29.3% | 11.4% |
| OTHER | 3,973 | 16.2% | 54.8% | 4.0% |
| NIH | 120 | 5.0% | 48.3% | 9.2% |
| OTHER_GOV | 94 | 11.7% | 47.9% | 4.3% |

### 3.4 Trial Duration at Termination

Mean duration to termination was 32–36 months across all failure categories (Table 3), a narrower range than expected given the mechanistic differences between failure types.

**Table 3. Duration at termination by failure_token**

| Failure Type | Mean (months) | Median (months) | n |
|-------------|--------------|----------------|---|
| EFFICACY | 32.8 | 27.8 | 424 |
| SAFETY | 32.2 | 25.2 | 501 |
| ENROLLMENT | 35.4 | 29.4 | 2,927 |
| REGULATORY | 36.3 | 28.0 | 291 |
| BUSINESS | 33.0 | 24.8 | 1,397 |
| EXTERNAL | 32.3 | 25.0 | 840 |

REGULATORY terminations had the longest mean duration (36.3 months), consistent with the protracted nature of regulatory disputes. SAFETY terminations had the shortest median duration (25.2 months), suggesting that safety signals tend to emerge within the first two years of a trial.

---

## 4. Discussion

### 4.1 Enrollment as the Primary Driver of Trial Termination

The most striking finding of this analysis is the dominance of ENROLLMENT failure (29.8% of all terminated studies, 46.0% of classifiable terminations). This substantially exceeds EFFICACY failure (4.3%), inverting the common narrative that clinical trials fail primarily because drugs do not work.

This inversion has important implications. The drug development community has invested heavily in improving target identification, biomarker selection, and efficacy prediction. Yet our data suggest that operational factors — the ability to recruit sufficient patients — represent a greater source of attrition than biological efficacy. Investment in enrollment strategy, site selection, and patient identification may yield higher returns than equivalent investment in preclinical efficacy optimization.

### 4.2 Counter-Intuitive Finding 1: PHASE4 Enrollment Crisis

PHASE4 studies exhibited the highest ENROLLMENT failure rate (57.3%), exceeding PHASE1 (32.2%) by 25 percentage points. This is counter-intuitive: PHASE4 studies involve drugs that have already received regulatory approval, suggesting that the patient population exists and the drug's safety and efficacy are established.

Two mechanisms may explain this finding. First, PHASE4 studies often target narrow subpopulations for label expansion, making recruitment inherently difficult. Second, the proliferation of competing trials for approved drugs creates enrollment competition that early-phase studies do not face. Whatever the mechanism, this finding suggests that enrollment risk does not diminish with drug maturity.

### 4.3 Counter-Intuitive Finding 2: Business Failure Rivals Safety in PHASE1

BUSINESS termination in PHASE1 (33.4%) was more than twice the rate of SAFETY termination (15.4%). This challenges the framing of early-phase trials as primarily safety screening exercises. For a substantial proportion of PHASE1 terminations, the limiting factor is not biological — it is financial. Sponsor withdrawal, budget constraints, and reprioritization decisions terminate more early-phase trials than toxicity findings.

This finding has implications for how early-phase trial portfolios are managed. Risk models that focus on toxicity probability without accounting for sponsor commitment and financial sustainability may systematically underestimate true attrition risk.

### 4.4 Counter-Intuitive Finding 3: The Year-Three Crisis

The convergence of mean termination duration across all failure types (32–36 months) suggests that clinical trials across all failure categories tend to reach a critical decision point at approximately three years. This "year-three crisis" appears to be structural rather than driven by any specific failure mechanism.

One interpretation is that three years represents the natural horizon for initial resource commitments in clinical research. Regardless of whether a trial is failing due to poor enrollment, safety signals, or efficacy concerns, the financial and institutional structures that support trial continuation tend to be stress-tested at approximately this interval.

### 4.5 Limitations

Several limitations should be noted. First, keyword-based classification introduces misclassification error, particularly for ambiguous termination descriptions. The 35.3% UNKNOWN rate represents the ceiling of this limitation. Second, the analysis is restricted to TERMINATED studies from ClinicalTrials.gov, which is predominantly US-centric and may not fully represent global trial termination patterns. Third, WhyStopped field entries are self-reported by study sponsors and may not accurately reflect the primary termination cause.

### 4.6 Future Directions

The tokenized dataset (tokenized_cumulative.csv) provides a foundation for two downstream applications. Predictive modeling — using phase, sponsor type, indication, and mechanism to predict failure probability — represents the near-term application. Reverse inference — identifying viable mechanisms from EFFICACY-failed trials — connects Clinical-Gate to the drug design pipeline established in prior Sigma Protocol projects.

Full-scale processing of all TERMINATED studies on ClinicalTrials.gov (estimated 50,000–70,000 records) is ongoing.

---

## 5. Conclusion

Clinical-Gate demonstrates that large-scale structural analysis of clinical trial failure patterns is achievable through natural language tokenization of publicly available registry data. Analysis of 10,000 terminated interventional studies reveals that enrollment failure, not efficacy failure, is the dominant cause of clinical trial termination — and that this pattern intensifies rather than diminishes in late-phase trials. The convergence of trial duration across failure types suggests a structural "year-three crisis" that transcends failure mechanism. These findings challenge prevailing assumptions about the drivers of drug development attrition and establish a structured dataset for predictive and generative applications in drug development.

---

## Data Availability

The tokenized dataset and pipeline code will be made available upon publication. Raw data are publicly accessible at https://clinicaltrials.gov/data-api/api.

## Acknowledgments

This work was conducted independently without institutional affiliation or external funding. All computational pipeline development was performed using Claude Code (Anthropic) via natural language prompting without direct programming.

## Competing Interests

The author declares no competing interests.

---

*Correspondence: Woonki Cheon, Busan, Republic of Korea*  
*Sigma Protocol Project — Clinical-Gate*
