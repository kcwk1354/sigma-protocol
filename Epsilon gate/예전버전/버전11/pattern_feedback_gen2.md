# Epsilon-Gate Gen2 전체 파이프라인 패턴 보고서
# (Gen3 설계를 위한 컨텍스트 입력용)

## 파이프라인 3단계 결과 통합
물성 필터(Epsilon-Gate) → 도킹(AutoDock Vina) → MD 안정성(OpenMM 10ps)

---

## 단계별 발견 규칙

### 규칙 1: ether_ethyl 링커가 전 단계 압도적 우위
- QED 평균: ether_ethyl(0.916) > amide_ethyl(0.872) > amide_direct(0.867)
- 도킹 상위 5개 전부 ether_ethyl 조합
- MD Stable 유일 분자도 ether_ethyl
- → Gen3: ether_ethyl 링커를 기본으로 고정, 변형 탐색

### 규칙 2: 도킹 강도 ≠ MD 안정성 (핵심 발견)
- pyrimidyl: 도킹 1위(-3.844) → MD 불안정(RMSD 3.77Å) → 탈락
- fluorophenyl: 도킹 2위(-3.758) → MD 보통(RMSD 3.38Å) → 조건부
- pyridyl_4: 도킹 3위(-3.547) → MD 유일 안정(RMSD 2.97Å) → 최종 생존
- → 도킹 점수만으로 우열 판단 불가. 반드시 MD까지 봐야 함

### 규칙 3: pyridyl_4 꼬리가 안정성 최적
- 6원환 단일 질소(4번 위치) 구조가 가장 안정적 구조 제공
- pyridyl_2/3은 질소 위치 차이로 RMSD 미세하게 더 큼
- → Gen3: pyridyl_4 꼬리 유지, 여기에 치환기 추가 탐색

### 규칙 4: LogP 황금 구간 재확인
- MD 통과 분자 LogP: 3.115 (CNS 황금 구간 중심)
- 도킹 1위였던 pyrimidyl LogP: 2.510 (너무 극성 → 구조 불안정)
- → Gen3 목표 LogP: 3.0~3.5 유지 (Gen1 패턴 재확인)

### 규칙 5: 꼬리 극성이 양날의 검
- 질소 많을수록(pyrimidyl > pyridyl) 도킹 점수는 좋아짐
- 그러나 질소 많을수록 구조 유연성 증가 → MD 불안정
- → Gen3: 극성과 강직성 균형. 방향족 고리 추가로 구조 잠금

---

## Gen3 설계 지시

### 전략: "ether_ethyl + pyridyl_4 핵심 유지, 주변 변형"

**방향 1: 링커 길이 변형**
- ether_propyl (CH2 하나 추가) → 유연성↑, 더 깊은 포켓 접근 가능
- ether_direct (CH2 없이 직결) → 강직성↑, 안정성 테스트

**방향 2: 꼬리 치환기 추가 (pyridyl_4 기반)**
- methyl_pyridyl_4 (3-methyl-pyridyl) → LogP 미세 상승, 소수성 증가
- fluoro_pyridyl_4 (3-fluoro-pyridyl) → LogP 유지, 대사 안정성↑
- cyano_pyridyl_4 (3-cyano-pyridyl) → 수소결합 수용체 추가

**방향 3: 이중 꼬리 (새로운 시도)**
- bipyridyl → 두 개의 질소로 도킹 강화 + 고리 강직성으로 MD 안정 기대
- pyridyl_4 + OH → 수소결합 추가

### 제외 목록
- pyrimidyl 꼬리: 도킹 강하나 MD 불안정, 제외
- amide_direct 링커: QED 최저, 제외
- aminopyridyl 꼬리: QED 0.789로 최저, 제외
- chlorophenyl: Gen1부터 제외 유지

### 목표 수치
- LogP: 3.0~3.5 (pyridyl_4 기준 3.115 근방)
- QED: 0.90 이상
- MD RMSD 평균: 3.0Å 미만
- 도킹 ΔG: -3.5 이하 유지
