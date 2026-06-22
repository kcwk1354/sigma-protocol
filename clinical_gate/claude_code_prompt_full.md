# Clinical-Gate 전체 데이터셋 처리 프롬프트
# Claude Code에 그대로 붙여넣기
# Gen2 키워드 룰 확정 적용 — TERMINATED 전건 처리

---

## 사전 확인
- Gen2 키워드 룰 확정 (UNKNOWN 21% — 현실적 하한선)
- filter.studyType 파라미터 사용 금지 (v2 400 에러)
- filter.overallStatus: TERMINATED 단독 사용

## 목표
ClinicalTrials.gov TERMINATED 전건을 페이지네이션으로 fetch하고
Gen2 키워드 룰로 토큰화해 전체 데이터셋을 생성한다.

## 작업 디렉토리
C:\sigma_protocol\clinical_gate\full\
디렉토리가 없으면 새로 생성한다.

## 단계별 지시

### Step 1 — 환경 세팅
다음 패키지를 설치한다:
- requests >= 2.31
- pandas >= 2.1

### Step 2 — 전체 데이터 fetch (페이지네이션)
ClinicalTrials.gov API v2 엔드포인트:
https://clinicaltrials.gov/api/v2/studies

**페이지네이션 방식:**
- pageSize: 1000 (최대값)
- 첫 요청 후 응답의 nextPageToken이 있으면 계속 fetch
- nextPageToken이 없으면 종료
- 요청 간 0.2초 대기 (rate limit 준수)
- 진행 상황을 100페이지마다 콘솔에 출력

파라미터:
- filter.overallStatus: TERMINATED
- pageSize: 1000
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

전체 원본을 raw_studies_full.json으로 저장한다.
fetch 완료 후 총 건수를 콘솔에 출력한다.

### Step 3 — 토큰화 (Gen2 확정 키워드 룰)

**failure_token**
WhyStopped 텍스트를 읽고 분류 (대소문자 무시):

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
under-enrollment, failed to enroll, insufficient enrollment,
poor recruiting

REGULATORY 키워드:
FDA, IND, regulatory, hold, clinical hold,
IRB, approval expired, protocol expired,
compliance, protocol compliance, protocol deviation,
GCP, violation

BUSINESS 키워드:
funding, sponsor, business, financial,
lack of funding, insufficient funds, finance, lack of finance,
company management, manufacturer, corporate, budget,
resource, prioritization, strategic decision

EXTERNAL 키워드:
pandemic, covid, investigator, site closure, natural disaster,
war, logistics, departure, drug shortage, product expired,
supplies unavailable, understaffing, competing studies,
PI closed, practice closed, PI and practice, supply issue,
manufacturing issue, material unavailable

UNKNOWN: WhyStopped 공백/null이거나 위 키워드 미해당

**failure_confidence**
0.0~1.0. 0.5 미만이면 UNKNOWN 강제.
키워드 매칭 1개: 0.6, 2개 이상: 0.8, 0개: 0.0

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
토큰화 결과를 tokenized_full.csv로 저장한다.
컬럼 순서:
nct_id, phase_token, status_token, failure_token, failure_confidence,
intervention_type, endpoint_category, endpoint_miss,
sponsor_type, enrollment_ratio, duration_months, condition_token,
why_stopped_raw, primary_endpoint_raw

1000건마다 중간 저장 (tokenized_full_checkpoint.csv 덮어쓰기).
중단되더라도 체크포인트에서 재개 가능하도록.

### Step 5 — 전체 품질 리포트
다음 수치를 콘솔에 출력하고 full_report.txt로 저장한다:

- 전체 fetch 건수
- 전체 토큰화 건수
- failure_token 분포 (카테고리별 건수 및 %)
- UNKNOWN 비율
- endpoint_miss 비율
- failure_confidence 평균
- sponsor_type 분포
- phase_token 분포
- intervention_type 분포
- NULL 필드 현황 (컬럼별 누락률)

### Step 6 — tool_cards 자동 갱신
C:\sigma_protocol\clinical_gate\tool_cards\ 경로의 YAML 파일들을 업데이트한다:

clinicaltrials_api.yaml:
- last_result: SUCCESS 또는 FAIL
- last_fetch_count: 전체 fetch 건수
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
- raw_studies_full.json 생성 확인
- tokenized_full.csv 생성 확인
- full_report.txt 생성 확인
- tool_cards 갱신 확인

## 주의사항
- filter.studyType 파라미터 절대 사용 금지 (v2 400 에러)
- 페이지네이션 중 네트워크 오류 시 해당 페이지 3회 재시도 후 스킵
- 토큰화 오류 시 해당 레코드 UNKNOWN 처리 후 계속
- 예상 소요 시간: fetch 20~40분, 토큰화 10~20분
