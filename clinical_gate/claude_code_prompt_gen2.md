# Clinical-Gate Gen2 실행 프롬프트
# Claude Code에 그대로 붙여넣기
# Gen1 교훈 반영: UNKNOWN 41건 원문 패턴 기반 키워드 룰 대폭 보완

---

## Gen1에서 발견된 문제
- UNKNOWN 41% (목표 30% 미만)
- UNKNOWN 원문 분석 결과 추가 가능한 키워드 패턴 9종 확인
- 키워드 룰 보완으로 30% 미만 달성 가능

## 목표
TERMINATED 건 100건 재처리.
보완된 키워드 룰 적용 후 UNKNOWN 30% 미만 달성.

## 작업 디렉토리
C:\sigma_protocol\clinical_gate\
디렉토리가 없으면 새로 생성한다.

## 단계별 지시

### Step 1 — 기존 데이터 재사용
raw_studies_gen1.json을 그대로 사용한다.
API를 다시 호출하지 않는다. 토큰화만 재실행한다.

### Step 2 — 토큰화 재실행 (보완된 키워드 룰 적용)
raw_studies_gen1.json을 읽어 각 레코드에 대해 토큰을 생성한다.

**failure_token** ← 핵심. 아래 키워드 룰 전체 적용
WhyStopped 텍스트를 읽고 다음 중 하나로 분류:
- EFFICACY
- SAFETY
- ENROLLMENT
- REGULATORY
- BUSINESS
- EXTERNAL
- UNKNOWN

**failure_confidence**
분류 신뢰도 0.0~1.0. 0.5 미만이면 UNKNOWN 강제.
규칙 기반 키워드 매칭 (대소문자 무시):

EFFICACY 키워드:
efficacy, endpoint, did not meet, insufficient effect, no benefit,
MTD, insufficient clinical response, no clinical response,
lack of efficacy, no efficacy, poor efficacy

SAFETY 키워드:
safety, adverse, toxicity, death, serious, side effect,
device deficiency, device failure, toxic, harm, risk

ENROLLMENT 키워드:
enrollment, accrual, recruitment, slow accrual, low accrual,
difficult to enroll, inability to recruit, poor enrollment,
under-enrollment, failed to enroll, insufficient enrollment

REGULATORY 키워드:
FDA, IND, regulatory, hold, clinical hold,
IRB, approval expired, protocol expired,
compliance, protocol compliance, protocol deviation,
GCP, violation

BUSINESS 키워드:
funding, sponsor, business, financial,
lack of funding, insufficient funds, finance, lack of finance,
company management, manufacturer, corporate, budget,
resource, prioritization

EXTERNAL 키워드:
pandemic, covid, investigator, site closure, natural disaster,
war, logistics, departure, drug shortage, product expired,
supplies unavailable, understaffing, competing studies,
PI closed, practice closed, PI and practice, supply issue,
manufacturing issue, material unavailable

**나머지 토큰은 Gen1과 동일하게 처리:**

phase_token: PHASE1 | PHASE2 | PHASE3 | PHASE4 | NA
status_token: TERMINATED
intervention_type: DRUG | BIOLOGICAL | DEVICE | PROCEDURE | OTHER
endpoint_category: SURVIVAL | RESPONSE_RATE | BIOMARKER | SAFETY | QOL | OTHER
endpoint_miss: PrimaryOutcomeMeasure 공백이면 True
sponsor_type: INDUSTRY | NIH | OTHER_GOV | ACADEMIC | OTHER
enrollment_ratio: 실제등록 / 목표등록. 계산 불가 시 NULL
duration_months: (CompletionDate - StartDate) 월 단위. 계산 불가 시 NULL
condition_token: Condition 첫 번째 값 소문자 정규화

### Step 3 — 저장
토큰화 결과를 tokenized_gen2.csv로 저장한다.
컬럼 순서:
nct_id, phase_token, status_token, failure_token, failure_confidence,
intervention_type, endpoint_category, endpoint_miss,
sponsor_type, enrollment_ratio, duration_months, condition_token,
why_stopped_raw, primary_endpoint_raw

### Step 4 — 품질 리포트
다음 수치를 콘솔에 출력하고 gen2_report.txt로도 저장한다:

- 전체 레코드 수
- failure_token 분포 (카테고리별 건수 및 %)
- UNKNOWN 비율 ← 핵심 지표
- Gen1 대비 UNKNOWN 감소량
- endpoint_miss 비율
- failure_confidence 평균
- UNKNOWN 잔존 건의 why_stopped_raw 원문 전체 출력

판단 기준:
```
UNKNOWN < 30% → "Gen2 성공. 전체 데이터셋 처리 준비 완료."
UNKNOWN 30~40% → "경고: 추가 키워드 보완 필요. UNKNOWN 원문 확인."
UNKNOWN > 40% → "키워드 룰 재설계 필요."
```

### Step 5 — tool_cards 자동 갱신
C:\sigma_protocol\clinical_gate\tool_cards\ 경로의 YAML 파일들을 업데이트한다:

text_normalizer.yaml:
- last_result: SUCCESS 또는 FAIL
- run_count: +1
- last_verified: 오늘 날짜
- last_run_error: 오류 있으면 기록

## 완료 조건
- tokenized_gen2.csv 생성 확인
- gen2_report.txt 생성 확인
- UNKNOWN 잔존 원문 출력 확인
- tool_cards 갱신 확인

## 주의사항
- API 재호출 없음. raw_studies_gen1.json 재사용
- 오류 발생 시 해당 레코드 UNKNOWN 처리 후 계속
- UNKNOWN 비율과 무관하게 끝까지 실행
