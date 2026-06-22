# Clinical-Gate Gen1 실행 프롬프트
# Claude Code에 그대로 붙여넣기
# Gen0 교훈 반영: TERMINATED 전용 fetch + filter.studyType 제거

---

## Gen0에서 발견된 문제 (반드시 확인)
1. filter.overallStatus에 COMPLETED를 포함하면 WhyStopped 없는 건이 대량 유입 → UNKNOWN 92%
2. filter.studyType=INTERVENTIONAL 파라미터가 API v2에서 400 에러 발생 → 제거

## 목표
TERMINATED 건 100건만 fetch해서 failure_token 실질 분류율을 확인한다.
목표: UNKNOWN 30% 미만

## 작업 디렉토리
C:\sigma_protocol\clinical_gate\
디렉토리가 없으면 새로 생성한다.

## 단계별 지시

### Step 1 — 환경 세팅
다음 패키지를 설치한다:
- requests >= 2.31
- pandas >= 2.1

### Step 2 — 데이터 fetch (clinicaltrials_api)
ClinicalTrials.gov API v2 엔드포인트:
https://clinicaltrials.gov/api/v2/studies

다음 조건으로 100건을 fetch한다:
- filter.overallStatus: TERMINATED        ← Gen0에서 수정. COMPLETED 제거
- pageSize: 100                            ← filter.studyType 제거 (v2 400 에러)
- format: json

가져올 필드:
- NCTId
- Phase
- OverallStatus
- WhyStopped
- PrimaryOutcomeMeasure
- Condition
- InterventionName
- InterventionType
- LeadSponsorClass
- EnrollmentCount
- EnrollmentType
- StartDate
- CompletionDate

결과를 raw_studies_gen1.json으로 저장한다.

### Step 3 — 토큰화 (text_normalizer)
raw_studies_gen1.json을 읽어 각 레코드에 대해 다음 토큰을 생성한다:

**phase_token**
Phase 값을 다음 중 하나로 변환:
PHASE1 | PHASE2 | PHASE3 | PHASE4 | NA

**status_token**
OverallStatus를 그대로 사용:
TERMINATED (이번 fetch는 전건 TERMINATED)

**failure_token** ← 핵심
WhyStopped 텍스트를 읽고 다음 중 하나로 분류:
- EFFICACY: 효능 부족, 1차 엔드포인트 미달
- SAFETY: 부작용, 독성, AE
- ENROLLMENT: 등록 부족, 모집 실패
- REGULATORY: FDA/규제기관 이슈
- BUSINESS: 자금, 스폰서 철수
- EXTERNAL: 외부 환경 요인 (pandemic, 연구자 이탈, 시설 폐쇄 등)
- UNKNOWN: WhyStopped 공백이거나 분류 불가

**failure_confidence**
분류 신뢰도 0.0~1.0. 0.5 미만이면 failure_token을 UNKNOWN으로 강제 변경.
규칙 기반 키워드 매칭으로 구현 (LLM 호출 없이):
- EFFICACY 키워드: efficacy, endpoint, did not meet, insufficient effect, no benefit
- SAFETY 키워드: safety, adverse, toxicity, death, serious, side effect
- ENROLLMENT 키워드: enrollment, accrual, recruitment, slow accrual, low accrual
- REGULATORY 키워드: FDA, IND, regulatory, hold, clinical hold
- BUSINESS 키워드: funding, sponsor, business, financial, lack of funding
- EXTERNAL 키워드: pandemic, covid, investigator, site closure, natural disaster, war, logistics, departure

**intervention_type**
InterventionType 열거값 그대로:
DRUG | BIOLOGICAL | DEVICE | PROCEDURE | OTHER

**endpoint_category**
PrimaryOutcomeMeasure 텍스트를 분류:
SURVIVAL | RESPONSE_RATE | BIOMARKER | SAFETY | QOL | OTHER

**endpoint_miss**
PrimaryOutcomeMeasure 필드가 비어있으면 True, 있으면 False

**sponsor_type**
LeadSponsorClass 값 매핑:
INDUSTRY | NIH | OTHER_GOV | ACADEMIC | OTHER

**enrollment_ratio**
EnrollmentCount(실제) / EnrollmentCount(목표).
EnrollmentType이 ACTUAL이면 분자, ANTICIPATED이면 분모.
계산 불가 시 NULL.

**duration_months**
(CompletionDate - StartDate) 월 단위. 계산 불가 시 NULL.

**condition_token**
Condition 필드 첫 번째 값을 소문자로 정규화.

### Step 4 — 저장
토큰화 결과를 tokenized_gen1.csv로 저장한다.
컬럼 순서:
nct_id, phase_token, status_token, failure_token, failure_confidence,
intervention_type, endpoint_category, endpoint_miss,
sponsor_type, enrollment_ratio, duration_months, condition_token,
why_stopped_raw, primary_endpoint_raw

### Step 5 — 품질 리포트
다음 수치를 콘솔에 출력하고 gen1_report.txt로도 저장한다:
- 전체 레코드 수
- failure_token 분포 (각 카테고리별 건수 및 %)
- UNKNOWN 비율 ← 핵심 지표. 30% 미만이면 성공
- endpoint_miss 비율
- failure_confidence 평균
- WhyStopped 공백 건수
- NULL 필드 현황 (컬럼별 누락률)

판단 기준을 리포트 마지막에 출력한다:
```
UNKNOWN < 30%  → "Gen1 성공. 전체 데이터셋 처리 준비 완료."
UNKNOWN 30~50% → "경고: 키워드 룰 추가 보완 필요. UNKNOWN 샘플 확인 요망."
UNKNOWN > 50%  → "키워드 룰 재설계 필요. UNKNOWN 원문 전체 출력."
```

UNKNOWN 건의 why_stopped_raw 원문을 전체 출력한다 (키워드 룰 보완용).

### Step 6 — tool_cards 자동 갱신
C:\sigma_protocol\clinical_gate\tool_cards\ 경로의 YAML 파일들을 업데이트한다:

clinicaltrials_api.yaml:
- last_result: SUCCESS 또는 FAIL
- last_fetch_count: fetch된 레코드 수
- miss_rate: endpoint_miss 비율
- status: OK
- run_count: +1
- last_verified: 오늘 날짜

text_normalizer.yaml:
- last_result: SUCCESS 또는 FAIL
- run_count: +1
- last_verified: 오늘 날짜
- last_run_error: 오류 있으면 기록

## 완료 조건
- tokenized_gen1.csv 생성 확인
- gen1_report.txt 생성 확인
- tool_cards YAML 갱신 확인
- UNKNOWN 원문 출력 확인

## 주의사항
- filter.studyType 파라미터 절대 사용하지 말 것 (v2 400 에러)
- 오류 발생 시 중단하지 말고 해당 레코드를 UNKNOWN으로 처리 후 계속
- UNKNOWN 비율이 높아도 중단하지 말고 원문 출력 후 완료
