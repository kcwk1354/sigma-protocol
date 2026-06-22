# Lab-in-the-Loop 1바퀴 실측 기록

## 루프 구조 (재학습 없이 컨텍스트만으로)
1세대 스크리닝 → 실패 데이터 생성
  → 패턴 문서로 구조화 (pattern_feedback_gen1.md)
  → 그 패턴을 다음 세대 설계에 반영
  → 2세대 스크리닝 (gen2_library.csv)

## 결과
| 지표 | Gen1 | Gen2 | 변화 |
|------|------|------|------|
| 통과율 | 91% | 95% | ↑ 개선됨 |
| 탈락 분자 | 3개 (전부 chlorophenyl) | 1개 | ↓ 실패 회피 성공 |
| 최고 QED | 0.925 | 0.925 | — (1차지표 천장) |
| 신규 후보 | - | pyrimidyl 계열 발견 | 패턴이 새 영역 제시 |

## 증명된 것
- 실패 데이터 → 패턴 → 개선 루프가 실제로 작동함
- 재학습(GPU) 없이 컨텍스트 엔지니어링만으로 다음 세대가 개선됨
- 개인 리소스로 무한히 돌릴 수 있는 루프

## 드러난 한계 (정직)
- 1차 물성 지표는 미세 구조 차이(질소 위치 등)를 구분 못 함
  → pyridyl_2/3/4 가 동일 QED/LogP 출력
- 미세 우열 판정은 도킹 ΔG / MD 필요 → 인프라 영역
- 즉 루프의 '거친 개선'은 개인이, '정밀 판정'은 인프라가

## 핵심 결론
설계 + 실패학습 루프 = 개인이 AI로 구동 가능 (오늘 실증)
정밀 검증(도킹 배치/MD/wet-lab) = 인프라 영역
→ 이 경계가 곧 CDMO 해자이자, 리소스가 필요한 지점

---

## Gen2 → Gen2+Docking 업데이트

### 추가된 파일
- `epsilon_gate_gen2_docking.py` — AutoDock Vina 도킹 단계 통합

### 파이프라인 변화
```
[기존]  RDKit 물성 필터 → gen2_library.csv
[추가]  → AutoDock Vina 도킹 (ΔG, kcal/mol)
        → gen2_docking_results.csv (도킹 순위)
```

### 도킹 설정
- 타겟: AChE (PDB: 1EVE, 알츠하이머)
- 결합 포켓: Center(2.0, -14.0, 20.0) / Box 20×20×20 Å
- 기준점: AP2601 ΔG = -10.81 kcal/mol
- 도킹 후보: QED 상위 5개

### 실행 환경
- 코랩에서 실제 PDB + prepare_receptor로 완전 실행
- 코드 내 코랩 설치 가이드 포함 (COLAB_SETUP 변수)

### 남은 한계
- 1ns MD 안정성 검증 (OpenMM) → Gen3에서 추가 예정
- 도킹 결과 기반 Gen3 패턴 피드백 루프 자동화 예정

---

## Gen3 루프 결과

### 설계 근거 (pattern_feedback_gen2.md)
- ether_ethyl 링커 고정 (Gen2 전 단계 우위)
- pyridyl_4 기반 유지 (MD 유일 안정)
- 신규 탐색: ether_propyl/direct 링커, cyano/methyl/fluoro/bipyridyl 꼬리

### 결과 비교

| 지표 | Gen2 | Gen3 | 변화 |
|------|------|------|------|
| 통과율 | 95% | 89% | ↓ (신규 구조 탐색 비용) |
| 최고 QED | 0.925 | 0.928 | ↑ 미세 개선 |
| 도킹 최고 ΔG | -3.844 | **-4.191** | ↑ 결합력 개선 |
| MD 안정 판정 | ⭐ 보통 | **⭐⭐ 안정** | ↑ 안정성 개선 |

### Gen3 최종 후보
**ether_direct__bipyridyl**
- 도킹 ΔG: -4.191 kcal/mol (Gen2 최고 대비 +0.347 개선)
- MD RMSD: 2.37Å (Gen2 최고 2.97Å 대비 개선)
- MD 판정: ⭐⭐ 안정 (Gen2는 ⭐ 보통이 최고)
- bipyridyl의 이중 고리 강직성이 도킹↑ + 안정성↑ 동시 달성

### 루프가 증명한 것
- Gen1 → Gen2: 물성 필터 통과율 개선 (91%→95%)
- Gen2 → Gen3: 도킹+MD 동시 개선 (-3.547/2.97Å → -4.191/2.37Å)
- 매 세대 실패 데이터가 다음 설계를 개선함

### Gen4 힌트 (다음 루프)
- bipyridyl 구조 추가 탐색 (methyl-bipyridyl, fluoro-bipyridyl)
- ether_direct 링커 + 다른 꼬리 조합 확대
- LogP 3.5~4.0 구간 탐색 (ether_direct는 LogP가 낮아 여유 있음)
