# Epsilon-Gate Gen3 전체 파이프라인 패턴 보고서
# (Gen4 설계를 위한 컨텍스트 입력용)

## 핵심 발견: hERG가 전원 탈락 원인

### Gen3 독성 데이터 요약
| 분자 | hERG | DILI | BBB |
|------|------|------|-----|
| ether_direct__cyano_pyridyl_4 | **0.317** | 0.411 | 0.828 |
| ether_direct__methyl_pyridyl_4 | 0.436 | 0.126 | 0.930 |
| ether_ethyl__pyridyl_4 | 0.708 | 0.199 | 0.933 |
| ether_direct__bipyridyl | 0.800 | 0.434 | 0.842 |

---

## 발견된 규칙

### 규칙 1: cyano 꼬리가 hERG 최저 (핵심 발견)
- cyano_pyridyl_4: hERG 0.317 → 기준(0.4) 근접, 최저
- 이유: CN기의 강한 전자흡인 효과 → 분자 극성↑ → hERG 채널 친화도↓
- → Gen4: cyano 계열 확장 탐색 (cyano 위치 변형, 이중 cyano)

### 규칙 2: ether_direct 링커가 hERG에 유리
- ether_direct__cyano: 0.317 vs ether_ethyl__cyano: 0.411
- 링커가 짧을수록 hERG 낮아짐 → 분자 유연성이 hERG에 영향
- → Gen4: ether_direct 링커 우선, 짧은 링커 집중

### 규칙 3: bipyridyl은 hERG 최악 (제외)
- hERG 0.800 → MD 안정성은 좋았지만 심장독성 위험
- 이중 질소 고리 구조 자체가 hERG 채널과 강하게 결합
- → Gen4: bipyridyl 제외

### 규칙 4: BBB는 전체적으로 양호
- 전 분자 BBB > 0.82 → CNS 통과 문제 없음
- 이 특성은 유지하면서 hERG만 낮추는 게 목표

### 규칙 5: DILI는 cyano 계열에서 경계
- cyano_pyridyl_4의 DILI 0.41 → 기준(0.4) 미세 초과
- → Gen4: DILI도 함께 모니터링 필요

---

## Gen4 설계 지시

### 전략: "hERG 낮추기 — 극성 강화 + 링커 최적화"

**방향 1: cyano 위치 변형 (ether_direct 고정)**
- 2-cyano-pyridyl (cyano를 2번 위치로)
- 3-cyano-pyridyl (현재 3번, 기준점 유지)
- dicyano-pyridyl (질소 + cyano 두 개) → 극성 극대화
- cyano-pyrimidyl (pyrimidyl + cyano 조합)

**방향 2: 새로운 극성 꼬리 탐색**
- sulfonamide 꼬리 → 극성 매우 강함, hERG 낮추는 효과 알려짐
- tetrazole 꼬리 → 생물학적 등가체, 극성↑
- carboxyl_pyridyl → COOH + 질소 조합

**방향 3: 링커 추가 변형**
- ether_direct 유지 (Gen3 최저 hERG)
- amino_direct → NH 링커, 수소결합 공여체 추가

### 제외 목록
- bipyridyl: hERG 0.80 → 완전 제외
- ether_ethyl + 소수성 꼬리 조합: hERG 높음
- fluorophenyl: hERG 0.59, 개선 여지 없음

### Gen4 목표 수치
- hERG < 0.4 (심장독성 안전)
- DILI < 0.4 (간손상 안전)
- BBB > 0.7 (CNS 통과 유지)
- SA Score < 5.5 (합성 가능)
- QED > 0.85
