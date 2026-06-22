# Epsilon-Gate: AI-Driven Drug Discovery Engine

> **Gen0 수작업 독성 검사 → Gen12 병목 해결 엔진까지**  
> 비전공자가 오픈소스만으로 구축한 AChE 신약 후보 탐색 파이프라인

---

## 최종 후보 (Gen11)

```
direct_O_pyrrolidine
SMILES: CN1CCCC1OC12CC3CC(CC(O)(C3)C1)C2

물성:   MW=251.37  LogP=2.138  QED=0.816  SA=5.26
독성:   hERG=0.14 ✅  DILI=0.025 ✅  BBB=0.90 ✅
도킹:   ΔG=-4.071 kcal/mol ✅
MD:     RMSD=2.28Å ✅
```

**Gen2~Gen12 통틀어 처음으로 4지표를 동시에 통과한 분자.**

---

## 프로젝트 개요

- **타겟:** AChE (알츠하이머, PDB: 1EVE)
- **스캐폴드:** 아다만탄-OH
- **방법론:** Lab-in-the-Loop — 실패 데이터가 다음 설계를 만드는 자기개선 루프
- **리소스:** 개인 PC + 오픈소스만 사용 (GPU 없음)

---

## 파이프라인

```
화합물 설계 (RDKit SMILES)
    ↓
물성 필터 (Epsilon-Gate CNS)
  MW ≤ 450 / LogP 2~4 / TPSA ≤ 90 / QED ≥ 0.85
    ↓
독성 필터 (ADMET-AI)
  hERG < 0.4 / DILI < 0.4 / BBB > 0.5
    ↓
합성 용이성 (SA Score ≤ 5.5)
    ↓
도킹 (AutoDock Vina)
  AChE 결합 포켓 (TRP84/TYR121/PHE330/HIS440)
    ↓
MD 안정성 (OpenMM 10ps → 1ns)
  RMSD < 2.5Å
    ↓
병목 해결 엔진
  병목 감지 → 이전 최선값 유지 → 집중 최적화
```

---

## Gen0 → Gen12 학습 루프

| 세대 | 핵심 발견 | 결과 |
|------|-----------|------|
| Gen0 | ProTox-3 수작업 독성 검사 | "독성이 중요하다" |
| Gen1 | chlorophenyl LogP 초과 | F만 허용 |
| Gen2 | ether_ethyl 링커 최적 | 통과율 95% |
| Gen3 | bipyridyl 도킹↑ 안정↑ | hERG=0.80 탈락 |
| Gen4 | cyano+F → hERG 해결 | ΔG=-4.925 ✅ |
| Gen5 | F 두 개 → 도킹↓ 트레이드오프 | 패턴 확인 |
| Gen6 | 4-CN+2-F 위치 → RMSD 최저 | 2.26Å |
| Gen7 | 도네페질 파이페리딘 이식 | hERG=0.20 |
| Gen8 | N-메틸파이페리딘 하이브리드 | 균형점 |
| Gen9 | 고리 F → RMSD 악화 | 패턴 확인 |
| Gen10 | 병목 감지: MD가 주요 병목 | 엔진 구축 |
| Gen11 | 링커 최단화(CH2 제거) | **4지표 동시 달성** ✅ |
| Gen12 | N-방향족 = hERG↑ 발견 | 구조-활성 규칙 확립 |

---

## 전체 이정표

| 세대 | 분자 | ΔG | RMSD | hERG | DILI | 비고 |
|------|------|-----|------|------|------|------|
| Gen4 | cyano_F | -4.925 | 2.66Å | 0.39 | 0.30 | 도킹 최강 |
| Gen6 | 4cn_2f | -3.849 | 2.26Å | 0.40 | 0.29 | RMSD 최강 |
| Gen7 | adamantyl_pip | -4.467 | 3.48Å | 0.20 | 0.02 | 독성 최강 |
| Gen8 | Nmethyl_pip | -4.429 | 2.95Å | 0.29 | 0.015 | 균형점 |
| **Gen11** | **direct_O_pyrrolidine** | **-4.071** | **2.28Å** | **0.14** | **0.025** | **✅ 최종 후보** |

---

## 핵심 구조-활성 규칙 (Gen0~Gen12 학습)

```
1. chlorophenyl → 제외 (LogP 초과)
2. cyano+F 조합 → hERG 해결 핵심
3. F 위치: 4-CN+2-F = RMSD 최적 (activity cliff)
4. 파이페리딘 N-방향족 → hERG 급등 (Gen12 발견)
5. 고리 탄소 CN → hERG 유지
6. 링커 길이 = RMSD 결정 인자 (짧을수록 안정)
7. 아다만탄-O-고리(직결) = 강직성 + 독성 균형 최적
```

---

## 병목 해결 엔진 (epsilon_gate_bottleneck.py)

Gen10에서 구현한 핵심 방법론:

```python
# 병목 감지
bottleneck = detect_bottleneck({'물성': 0.85, '독성': 0.43, 'MD': 0.35})
# → MD 통과율 35% = 주요 병목

# 해결 전략
# 1. 병목 이전 최선값(독성) 유지
# 2. MD 병목 원인 분석 (유연성 → 링커 최단화)
# 3. 해당 부분만 집중 최적화
```

파레토 탐색보다 효율적: 수백 개 탐색 대신 10개 집중으로 병목 돌파.

---

## 파일 구조

```
├── epsilon_gate_gen2.py          # 물성 필터
├── epsilon_gate_gen2_docking.py  # 도킹 파이프라인
├── epsilon_gate_md.py            # MD 안정성
├── epsilon_gate_tox.py           # 독성 + SA Score
├── epsilon_gate_bottleneck.py    # 병목 해결 엔진 ← 핵심
├── epsilon_gate_gen3~12.py       # 세대별 설계
├── pattern_feedback_gen1~4.md    # 세대별 학습 문서
├── approved_drugs_analysis.csv   # 승인약 패턴 분석
├── LOOP_RECORD.md                # 전체 루프 기록
└── README.md
```

---

## 설치 및 실행

```bash
pip install rdkit vina meeko gemmi openmm openmmforcefields admet-ai pandas
```

```bash
# Gen11 후보 검증
python epsilon_gate_gen11.py

# 병목 해결 엔진
python epsilon_gate_bottleneck.py

# MD (N_STEPS 변경으로 스케일 조정)
python epsilon_gate_md.py --csv gen11_final_results.csv
```

---

## 다른 질병 적용

```python
# epsilon_gate_gen2_docking.py 내 3줄만 수정
PDB_ID        = "6LU7"        # 코로나 Mpro
POCKET_CENTER = (x, y, z)     # 해당 PDB 포켓 좌표
# → 나머지 파이프라인 동일
```

| 질병 | 타겟 | PDB |
|------|------|-----|
| 알츠하이머 | AChE | 1EVE ← 현재 |
| 파킨슨 | LRRK2 | 7LHW |
| 폐암 | EGFR | 1IEP |
| 코로나 | Mpro | 6LU7 |
| 당뇨 | DPP-4 | 2I78 |

---

## 한계 및 다음 단계

**현재 한계**
- Receptor: 포켓 부분구조(53원자) → 전체 1EVE 검증 필요
- MD: 진공 10ps → TIP3P 명시적 수분자 + 100ns
- wet-lab IC₅₀ 미측정

**다음 단계**
- [ ] 코랩에서 전체 1EVE 도킹 (실제 ΔG -8~-11 범위 확인)
- [ ] Gen11 후보 IC₅₀ CRO 측정
- [ ] Gen13: direct_O_pyrrolidine + CN 위치 최적화

---

## 참고 논문

1. **아다만탄 AChE 억제제 구조-활성 관계**  
   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6152672/

2. **F 도입으로 hERG 감소**  
   https://pubmed.ncbi.nlm.nih.gov/35397400/

3. **F 위치 변화와 activity cliff (최대 1300배)**  
   https://www.mdpi.com/1420-3049/28/2/490

4. **AChE 억제제 도킹 리뷰**  
   https://pmc.ncbi.nlm.nih.gov/articles/PMC9921523/
