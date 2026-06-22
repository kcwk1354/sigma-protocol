# Clinical-Gate Analysis001 — Phase × failure_token 교차 분석
# Claude Code에 그대로 붙여넣기
# 목적: 10,000건 누적 데이터에서 Phase별 실패 패턴 추출

---

## 목표
tokenized_cumulative.csv를 읽어 Phase별 failure_token 분포를 분석한다.
UNKNOWN을 제외한 분류 가능 건만 대상으로 한다.

## 작업 디렉토리
C:\sigma_protocol\clinical_gate\analysis\
디렉토리가 없으면 새로 생성한다.

## 입력 파일
C:\sigma_protocol\clinical_gate\batches\tokenized_cumulative.csv

## 단계별 지시

### Step 1 — 환경 세팅
다음 패키지를 설치한다:
- pandas >= 2.1
- matplotlib
- seaborn

### Step 2 — 데이터 로드 및 전처리
tokenized_cumulative.csv를 읽는다.

전처리:
- UNKNOWN 제외 → 분류 가능 건만 분석 (failure_token != 'UNKNOWN')
- phase_token이 NA인 건 별도 집계
- PHASE1 / PHASE2 / PHASE3 / PHASE4 / NA 5개 그룹

### Step 3 — Phase × failure_token 교차 분석

**3-1. 전체 요약**
전체 10,000건 기준:
- 분류 가능 건수 / UNKNOWN 건수
- failure_token 전체 분포 (건수 + %)

**3-2. Phase별 failure_token 분포**
각 Phase별로 failure_token 분포를 계산한다.
행: Phase / 열: failure_token 카테고리
값: 건수와 % 모두 출력

**3-3. 핵심 비교 지표**
다음 수치를 별도로 추출한다:

EFFICACY 실패율 by Phase:
- PHASE1에서 EFFICACY 비율
- PHASE2에서 EFFICACY 비율
- PHASE3에서 EFFICACY 비율
→ Phase가 올라갈수록 EFFICACY 비율이 높아지는지 확인

ENROLLMENT 실패율 by Phase:
- PHASE1에서 ENROLLMENT 비율
- PHASE2에서 ENROLLMENT 비율
- PHASE3에서 ENROLLMENT 비율
→ 대규모 Phase일수록 모집 실패가 높은지 확인

SAFETY 실패율 by Phase:
- PHASE1에서 SAFETY 비율 (초기 독성 탐색 단계)
- PHASE2/3 대비 PHASE1 SAFETY 비율

**3-4. sponsor_type × failure_token 교차**
INDUSTRY vs ACADEMIC 스폰서별 실패 패턴 차이:
- INDUSTRY에서 BUSINESS 실패 비율
- ACADEMIC에서 ENROLLMENT 실패 비율
→ 스폰서 유형에 따라 실패 원인이 다른지 확인

**3-5. duration_months × failure_token**
실패 유형별 평균 임상 기간:
- EFFICACY 실패: 평균 몇 개월 후 중단
- ENROLLMENT 실패: 평균 몇 개월 후 중단
- SAFETY 실패: 평균 몇 개월 후 중단
→ 어떤 실패가 더 일찍/늦게 발생하는지

### Step 4 — 시각화
다음 차트를 생성하고 PNG로 저장한다:

chart1_phase_failure_heatmap.png
- Phase(행) × failure_token(열) 히트맵
- 값: 각 Phase 내 비율(%)
- 색상: 높을수록 진하게

chart2_failure_distribution.png
- failure_token별 전체 분포 막대 그래프
- UNKNOWN 포함/제외 두 가지 버전

chart3_efficacy_by_phase.png
- Phase별 EFFICACY 실패 비율 라인 차트
- ENROLLMENT 실패 비율도 함께 표시

chart4_sponsor_failure.png
- INDUSTRY vs ACADEMIC 실패 패턴 비교 막대 그래프

### Step 5 — 인사이트 리포트
analysis001_report.txt로 저장한다.

형식:
```
=== Clinical-Gate Analysis001 Report ===
날짜: YYYY-MM-DD
분석 대상: N건 (UNKNOWN 제외)

[핵심 발견]
1. ...
2. ...
3. ...

[Phase별 failure_token 분포 테이블]
...

[EFFICACY 실패율 by Phase]
PHASE1: X%
PHASE2: X%
PHASE3: X%
해석: ...

[ENROLLMENT 실패율 by Phase]
PHASE1: X%
PHASE2: X%
PHASE3: X%
해석: ...

[SAFETY 실패율 by Phase]
PHASE1: X%
PHASE2: X%
PHASE3: X%
해석: ...

[sponsor_type × failure_token]
...

[duration_months × failure_token]
EFFICACY 평균: X개월
ENROLLMENT 평균: X개월
SAFETY 평균: X개월
해석: ...

[카운터인튜이티브 발견]
(예상과 다른 결과가 있으면 명시)
```

## 완료 조건
- analysis001_report.txt 생성
- chart1~4 PNG 4종 생성
- 핵심 발견 3개 이상 도출

## 주의사항
- UNKNOWN은 분석에서 제외하되 전체 비율은 리포트에 명시
- NULL 값은 해당 분석에서 제외하고 건수 명시
- 수치는 소수점 1자리까지
