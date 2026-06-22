# Clinical-Gate
**임상시험 실패 패턴 언어 토큰화 파이프라인**

*Sigma Protocol 세 번째 프로젝트 — 천운기 (독립 연구자, 부산)*

---

## 개요

Clinical-Gate는 ClinicalTrials.gov의 40만 건 이상 임상시험 데이터에서 **실패 패턴을 구조화된 언어 토큰으로 치환**하는 파이프라인이다.

기존에 임상시험 실패 데이터는 텍스트 요약 형태로만 존재했다. 구조화된 실패 패턴 토큰은 없었다. 이 프로젝트는 그 공백을 채운다.

```
비정형 텍스트 ("Study stopped due to lack of efficacy")
        ↓
구조화 토큰 (failure_token: EFFICACY, confidence: 0.87)
        ↓
패턴 분석 / 신약 설계 역추론
```

---

## 왜 이 프로젝트인가

### 수치로 보는 이유

신약 개발 비용의 **90%가 임상시험 실패 비용**이다. Phase 2 성공률은 약 30%, Phase 3는 약 60%. 실패 원인의 구조적 패턴이 있음에도 이를 언어 모델이 처리 가능한 형태로 치환한 시도는 거의 없었다.

### 세 가지 선정 기준 충족

Clinical-Gate는 Sigma Protocol의 프로젝트 선정 기준 3가지를 모두 충족한다.

**데이터 풍부** — ClinicalTrials.gov 40만 건 이상, 공개 API, 무료 접근  
**기준값 존재** — Phase 1/2/3 통과/실패 이진값 + 중단 사유 텍스트  
**미래 쓰임** — 실패 예측 모델 수요는 신약 개발 비용 구조상 폭발적

### 언어화 공백

| 도메인 | 언어화 상태 |
|--------|------------|
| 단백질 서열 | AlphaFold, ESM — 완료 |
| 소분자 구조 | SMILES, SELFIES — 완료 |
| 임상시험 실패 패턴 | **거의 안 됨** ← 여기 |

---

## 두 가지 활용 방향

파이프라인을 처음부터 두 방향 모두 지원하도록 설계했다.

### 예측 방향
실패 패턴 토큰 → 다음 임상 성공 확률 예측  
입력: 신약 후보의 메커니즘, 타깃, 적응증  
출력: Phase별 실패 확률 + 주요 위험 요인

### 역추론 방향
실패 원인 → viable mechanism 추출 → 신약 설계 시드  
Sigma Protocol의 핵심 방법론인 **역추론(Reverse Inference)** 적용  
실패한 약물에서 살릴 수 있는 코어를 추출해 다음 설계 공간 제약

---

## 토큰화 스키마

임상시험 1건 = 아래 토큰 벡터 1개

### 식별자 블록
| 토큰 | 타입 | 값 |
|------|------|----|
| `nct_id` | enum | NCT 고유 식별자 |
| `phase_token` | enum | PHASE1 \| PHASE2 \| PHASE3 \| PHASE4 \| NA |
| `status_token` | enum | COMPLETED \| TERMINATED \| WITHDRAWN \| SUSPENDED |

### 실패 유형 블록 ← 핵심
| 토큰 | 타입 | 값 |
|------|------|----|
| `failure_token` | derived | EFFICACY \| SAFETY \| ENROLLMENT \| REGULATORY \| BUSINESS \| UNKNOWN |
| `failure_confidence` | numeric | 0.0~1.0 (0.5 미만 → UNKNOWN 강제) |
| `why_stopped_raw` | text | WhyStopped 원문 — NLP 소스 |

### 개입 블록
| 토큰 | 타입 | 값 |
|------|------|----|
| `intervention_type` | enum | DRUG \| BIOLOGICAL \| DEVICE \| PROCEDURE \| OTHER |
| `drug_class_token` | derived | ATC 코드 매핑 |
| `mechanism_token` | derived | MOA 분류 |

### 엔드포인트 블록
| 토큰 | 타입 | 값 |
|------|------|----|
| `endpoint_category` | derived | SURVIVAL \| RESPONSE_RATE \| BIOMARKER \| SAFETY \| QOL \| OTHER |
| `endpoint_miss` | enum | True \| False \| NULL |
| `primary_endpoint_raw` | text | PrimaryOutcomeMeasure 원문 |

### 메타 블록
| 토큰 | 타입 | 값 |
|------|------|----|
| `condition_token` | derived | MeSH 코드 매핑 |
| `sponsor_type` | enum | INDUSTRY \| NIH \| OTHER_GOV \| ACADEMIC |
| `enrollment_ratio` | numeric | 실제등록 / 목표등록 |
| `duration_months` | numeric | StartDate → CompletionDate |

---

## 아키텍처

### 파이프라인 흐름

```
ClinicalTrials.gov API v2 (40만 건)
        ↓ clinicaltrials_api
raw_studies.json
        ↓ text_normalizer
tokenized.csv
        ↓ pattern_extractor
failure_pattern_df  →  예측 방향
viability_score     →  역추론 방향
```

### tool_cards 2-레이어 구조

Sigma Protocol의 tool_cards 시스템을 Clinical-Gate에 맞게 확장했다.

**Layer 1 — 항상 로드 (LLM 컨텍스트 상시 주입)**  
7개 이하 필드. LLM이 도구 선택 판단에 쓰는 최소 정보만.

**Layer 2 — 온디맨드 (`load_tool_detail()` 호출 시)**  
오류 발생 또는 설계 단계에서만 로드. 설치 절차, 스키마 이력, 상세 known_issues.

### Clinical-Gate 신규 필드 (Cancer-Gate 대비)

| 필드 | 위치 | 의미 |
|------|------|------|
| `status` | Layer 1 | OK \| RATE_LIMITED \| SCHEMA_CHANGED |
| `miss_rate` | Layer 1 | 해당 API 필드 누락률 |
| `last_fetch_count` | Layer 1 | 마지막 실행 레코드 수 |

---

## Sigma Protocol 방법론 적용

### Lab-in-the-Loop
Cancer-Gate에서 정립한 방법론을 그대로 계승.  
실패 데이터 → 구조 규칙 → 다음 세대 설계 공간 제약.  
Gen0 결과의 UNKNOWN 비율이 50% 초과이면 Gen1에서 키워드 룰 보완.

### 역추론 (Reverse Inference)
Cancer-Gate의 핵심 발견 방법을 임상 도메인에 적용.  
실패한 임상시험에서 살릴 수 있는 메커니즘을 추출 → 신약 설계 시드.

### Decision Log
판단 과정을 구조화 기록.  
`failure_token` 분류 근거, `confidence_threshold` 설정 이유 등 모든 판단 기록.

### 토큰 낭비 최소화 원칙
tool_cards 설계 시 "LLM이 이 정보를 읽고 다음 판단을 내릴 수 있는가"를 기준으로  
Layer 1/2를 분리. Layer 1 최대 7필드 유지.

---

## 완료된 프로젝트 (Sigma Protocol)

| 프로젝트 | 타깃 | 주요 결과 |
|---------|------|----------|
| Epsilon-Gate | AChE (알츠하이머) | direct_O_pyrrolidine — 5개 기준 동시 통과 |
| Cancer-Gate | GRM4 mGluR4 PAM | CG12-004 — RMSD 0.14Å, SI=19.5x |
| Battery-Gate | 고체 전해질 소재 | Li6PS4Cl2 — 16.29 mS/cm (실험값 ~17 독립 재현) |
| **Clinical-Gate** | **임상 실패 패턴** | **진행 중** |

---

## 데이터 소스

- **ClinicalTrials.gov API v2** — https://clinicaltrials.gov/api/v2/studies
- 등록 건수: 500,000건 이상 (2026년 기준)
- 갱신 주기: 월~금 매일 (ET 9AM)
- 라이선스: 공개 데이터, 무료 접근

---

## 현재 상태

- [x] 프로젝트 설계 완료
- [x] tool_cards 3종 설계 완료 (clinicaltrials_api, text_normalizer, pattern_extractor)
- [x] 토큰화 스키마 확정
- [x] Gen0 실행 — 100건 fetch + 토큰화 (UNKNOWN 92% → 키워드 룰 필요 확인)
- [x] Gen1 — 키워드 룰 보완 (UNKNOWN 92% → 41%)
- [x] Gen2 — 추가 보완 (UNKNOWN 41% → 21%), 전체 처리 준비 완료
- [x] 전체 데이터셋 1만 건 처리 + 분석 (6,471건 유효 / 3,529건 UNKNOWN 제외)
- [ ] 예측 모델 / 역추론 파이프라인

---

## 주요 발견 (Analysis001 — 1만 건 기준)

6,471건 유효 데이터에서 도출된 카운터인튜이티브 패턴 3개.

**1. PHASE4 모집 실패율이 PHASE1보다 높다**
```
PHASE1 ENROLLMENT 실패: 32.2%
PHASE4 ENROLLMENT 실패: 57.3%  ← 25%p 높음
```
일반적 예상(초기 임상이 어렵다)과 반대. 후기 임상일수록 적응증이 협소해지고 경쟁 시험이 증가하기 때문으로 추정.

**2. PHASE1 BUSINESS 실패율이 독성만큼 높다**
```
PHASE1 SAFETY  실패: 15.4%
PHASE1 BUSINESS 실패: 33.4%  ← 2배 이상
```
초기 임상에서 스폰서 철수·자금 부족이 독성 문제만큼 중요한 종료 원인.

**3. 실패 유형과 무관하게 "3년차 위기"가 존재한다**
```
EFFICACY    평균 임상 기간: 32.8개월
SAFETY      평균 임상 기간: 32.2개월
ENROLLMENT  평균 임상 기간: 35.4개월
REGULATORY  평균 임상 기간: 36.3개월
```
실패 원인과 무관하게 모든 유형이 32~36개월에 수렴. 원인이 아니라 타이밍이 구조적으로 고정되어 있음을 시사.

**Phase별 failure_token 분포 (유효 6,471건)**

|  | ENROLLMENT | BUSINESS | SAFETY | EFFICACY | REGULATORY |
|--|-----------|---------|--------|----------|-----------|
| PHASE1 | 32.2% | 33.4% | 15.4% | 6.1% | 4.4% |
| PHASE2 | 48.7% | 20.7% | 9.6% | 10.5% | 3.4% |
| PHASE3 | 40.0% | 19.4% | 11.7% | 14.9% | 7.5% |
| PHASE4 | 57.3% | 17.9% | 2.6% | 4.4% | 3.6% |

---

*독립 연구자 천운기 | Sigma Protocol | 부산, 2026*
