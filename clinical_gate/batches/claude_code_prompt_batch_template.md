# Clinical-Gate 배치 처리 템플릿
# Claude Code에 그대로 붙여넣기
# 사용법: BATCH_NUMBER만 바꿔서 반복 실행

---

## 배치 설정
BATCH_NUMBER: 001              ← 실행 전 여기만 변경 (001, 002, 003 ...)
BATCH_SIZE: 1000
PAGE_TOKEN: null               ← 002부터는 이전 배치 last_page_token 값으로 교체

## 목표
TERMINATED 건을 1000건씩 배치로 fetch하고 토큰화한다.
배치마다 품질 리포트를 생성해 개선사항을 확인한다.

## 작업 디렉토리
C:\sigma_protocol\clinical_gate\batches\
디렉토리가 없으면 새로 생성한다.

## 단계별 지시

### Step 1 — 환경 세팅
다음 패키지를 설치한다:
- requests >= 2.31
- pandas >= 2.1

### Step 2 — 데이터 fetch
ClinicalTrials.gov API v2 엔드포인트:
https://clinicaltrials.gov/api/v2/studies

파라미터:
- filter.overallStatus: TERMINATED
- pageSize: 1000
- format: json
- pageToken: PAGE_TOKEN 값 (null이면 파라미터 생략 — 첫 배치)

주의:
- filter.studyType 파라미터 사용 금지 (v2 400 에러)
- 응답의 nextPageToken을 반드시 저장 — 다음 배치에 필요

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

결과를 raw_batch{BATCH_NUMBER}.json으로 저장한다.
nextPageToken을 next_page_token_{BATCH_NUMBER}.txt로 저장한다.

### Step 3 — 토큰화 (Gen2 확정 키워드 룰)

**failure_token**
WhyStopped 텍스트를 읽고 분류 (대소문자 무시):

EFFICACY 키워드:
efficacy, endpoint, did not meet, insufficient effect, no benefit,
MTD, insufficient clinical response, no clinical response,
lack of efficacy, no efficacy, poor efficacy,
futility, stopped for futility, interim futility, futility analysis,
data did not support continuation, data did not support

SAFETY 키워드:
safety, adverse, toxicity, death, serious, side effect,
device deficiency, device failure, toxic, harm, risk

ENROLLMENT 키워드:
enrollment, accrual, recruitment, slow accrual, low accrual,
difficult to enroll, inability to recruit, poor enrollment,
under-enrollment, failed to enroll, insufficient enrollment,
poor recruiting,
lack of inclusion, difficulty in recruiting, unable to recruit,
not interested in participating, small number of enrolled, low interest,
unable to meet adequate enrolment, slow inclusion rates, recruiting or enrolling halted,
enrolling halted, inclusion halted

REGULATORY 키워드:
FDA, IND, regulatory, hold, clinical hold,
IRB, approval expired, protocol expired,
compliance, protocol compliance, protocol deviation,
GCP, violation,
protocol design under review, protocol redesign, protocol revision

BUSINESS 키워드:
funding, sponsor, business, financial,
lack of funding, insufficient funds, finance, lack of finance,
company management, manufacturer, corporate, budget,
resource, prioritization, strategic decision,
costs of study, became prohibitive, cost prohibitive, too costly,
decided not to pursue further development, not pursue, no longer pursue

EXTERNAL 키워드:
pandemic, covid, investigator, site closure, natural disaster,
war, logistics, departure, drug shortage, product expired,
supplies unavailable, understaffing, competing studies,
PI closed, practice closed, PI and practice, supply issue,
manufacturing issue, material unavailable,
study discontinued, study halted, halted, discontinued,
replaced with another study, superseded, merged,
production halt, supply halt, manufacturing halt

UNKNOWN: WhyStopped 공백/null이거나 위 키워드 미해당

**failure_confidence**
키워드 매칭 1개: 0.6 / 2개 이상: 0.8 / 0개: 0.0
0.5 미만이면 UNKNOWN 강제

**나머지 토큰:**
phase_token: PHASE1 | PHASE2 | PHASE3 | PHASE4 | NA
status_token: TERMINATED
intervention_type: DRUG | BIOLOGICAL | DEVICE | PROCEDURE | OTHER
endpoint_category: SURVIVAL | RESPONSE_RATE | BIOMARKER | SAFETY | QOL | OTHER
endpoint_miss: PrimaryOutcomeMeasure 공백이면 True
sponsor_type: INDUSTRY | NIH | OTHER_GOV | ACADEMIC | OTHER
enrollment_ratio: 실제등록 / 목표등록. 계산 불가 시 NULL
duration_months: (CompletionDate - StartDate) 월 단위. 계산 불가 시 NULL
condition_token: Condition 첫 번째 값 소문자 정규화

### Step 4 — 저장
토큰화 결과를 tokenized_batch{BATCH_NUMBER}.csv로 저장한다.
컬럼 순서:
nct_id, phase_token, status_token, failure_token, failure_confidence,
intervention_type, endpoint_category, endpoint_miss,
sponsor_type, enrollment_ratio, duration_months, condition_token,
why_stopped_raw, primary_endpoint_raw

기존 tokenized_cumulative.csv가 있으면 이번 배치를 추가(append)한다.
없으면 새로 생성한다.

### Step 5 — 배치 품질 리포트
다음 수치를 콘솔에 출력하고 report_batch{BATCH_NUMBER}.txt로 저장한다:

[이번 배치]
- fetch 건수
- failure_token 분포 (카테고리별 건수 및 %)
- UNKNOWN 비율
- UNKNOWN 잔존 원문 샘플 10건 (키워드 룰 보완 확인용)
- endpoint_miss 비율
- failure_confidence 평균
- nextPageToken: (다음 배치에 붙여넣을 값)

[누적 현황]
- tokenized_cumulative.csv 총 레코드 수
- 누적 failure_token 분포
- 누적 UNKNOWN 비율

### Step 6 — tool_cards 자동 갱신
C:\sigma_protocol\clinical_gate\tool_cards\ 경로의 YAML 파일들을 업데이트한다:

clinicaltrials_api.yaml:
- last_result: SUCCESS 또는 FAIL
- last_fetch_count: 이번 배치 fetch 건수
- miss_rate: endpoint_miss 비율
- status: OK
- run_count: +1
- last_verified: 오늘 날짜

text_normalizer.yaml:
- last_result: SUCCESS 또는 FAIL
- run_count: +1
- last_verified: 오늘 날짜

## 완료 조건
- raw_batch{BATCH_NUMBER}.json 생성 확인
- tokenized_batch{BATCH_NUMBER}.csv 생성 확인
- tokenized_cumulative.csv 업데이트 확인
- next_page_token_{BATCH_NUMBER}.txt 생성 확인 ← 다음 배치에 필요
- report_batch{BATCH_NUMBER}.txt 생성 확인

## 다음 배치 실행 방법
1. report_batch{BATCH_NUMBER}.txt에서 nextPageToken 값 확인
2. 이 템플릿의 BATCH_NUMBER를 +1
3. PAGE_TOKEN을 nextPageToken 값으로 교체
4. Claude Code에 붙여넣고 실행

## 주의사항
- filter.studyType 파라미터 절대 사용 금지
- nextPageToken 없으면 전체 처리 완료
- 오류 시 해당 레코드 UNKNOWN 처리 후 계속
