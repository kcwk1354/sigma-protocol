# Epsilon-Gate: AI-Driven Drug Discovery Engine

> **비전공자가 AI로 설계한 알츠하이머 신약 후보 탐색 파이프라인**  
> Gen0 수작업 독성 검사 → Gen6 논문 근거 기반 F/CN 위치 최적화까지

---

## 프로젝트 개요

이 프로젝트는 **재학습(GPU) 없이 AI 컨텍스트 엔지니어링만으로** 신약 후보 분자를 설계·평가·개선하는 자동화 루프를 구현합니다.

- **타겟:** AChE (Acetylcholinesterase) — 알츠하이머 치료 표적
- **스캐폴드:** 아다만탄-OH (혈액뇌장벽 통과에 유리한 CNS 구조)
- **파이프라인:** 물성 필터 → 독성 필터 → 합성 용이성 → 도킹 → MD 안정성

---

## 핵심 아이디어

> **"파이프라인"이 아니라 "엔진"입니다.**

파이프라인은 설계된 경로를 흘러가는 것. 이 엔진은 **실패 데이터에서 패턴을 뽑아 다음 세대 설계에 반영**합니다. 돌릴수록 더 좋은 분자가 나옵니다.

```
실패 데이터 → 패턴 추출 → 설계 수정 → 재실행
      ↑__________________________________|
```

---

## 파이프라인 구조

```
화합물 설계 (RDKit SMILES 조합)
        ↓
물성 필터 (Epsilon-Gate CNS 기준)
  MW ≤ 450 / LogP 2.0~4.0 / TPSA ≤ 90
  HBD ≤ 3  / HBA ≤ 7     / QED ≥ 0.85
        ↓
독성 필터 (ADMET-AI, 104가지 엔드포인트)
  hERG < 0.4   심장독성
  DILI < 0.4   약물유발 간손상
  BBB  > 0.5   혈액뇌장벽 통과 (CNS 필수)
  AMES < 0.3   변이원성
        ↓
합성 용이성 (SA Score, RDKit)
  SA Score ≤ 5.5  (1=쉬움, 10=불가)
        ↓
도킹 (AutoDock Vina)
  AChE 결합 포켓 (1EVE, TRP84/TYR121/PHE330/HIS440/SER200)
  ΔG 낮을수록 강한 결합
        ↓
MD 안정성 (OpenMM)
  진공 10ps → 100ps/1ns (코랩 GPU)
  RMSD < 2.5Å 목표
        ↓
최종 후보 선정
```

---

## Gen0 → Gen6 학습 루프

### Gen0 — 수작업 시작
- ProTox-3 웹사이트에 SMILES 직접 입력, 스크린샷으로 독성 확인
- "독성이 왜 중요한지" 손으로 배움
- 이 경험이 이후 파이프라인 설계의 토대

### Gen1 → Gen2 — 물성 학습
- **발견:** chlorophenyl 꼬리 LogP 4.14 → CNS 상한 초과, 탈락
- **학습:** F(플루오로)는 OK, Cl(클로로)는 위험
- **결과:** 통과율 91% → 95%

### Gen2 → Gen3 — 도킹+MD 추가, 구조 학습
- **발견:** 도킹 1위 pyrimidyl이 MD에서 불안정 (RMSD 3.77Å)
- **학습:** 도킹 점수 ≠ MD 안정성, 반드시 둘 다 봐야 함
- **발견:** bipyridyl 이중 고리 강직성 → 도킹↑ + 안정성↑
- **결과:** 최고 후보 bipyridyl ΔG=-4.191, RMSD=2.37Å

### Gen3 → Gen4 — 독성 필터 추가, hERG 문제 발견
- **발견:** Gen3 전원 hERG 탈락 (bipyridyl hERG=0.80)
- **학습:** hERG = 심장 이온채널, 억제되면 부정맥 위험
- **전략:** cyano(CN) 극성기 추가 → hERG↓
- **결과:** cyano_fluorophenyl hERG=0.39 ✅, ΔG=-4.925 ✅

### Gen4 → Gen5 — F/CN 트레이드오프 발견
- **발견:** F 두 개(difluoro) → hERG↓ but 도킹도↓
- **학습:** F 극성화가 과하면 AChE 포켓 친화도도 떨어짐
- **결론:** 단일 F + CN이 균형점

### Gen5 → Gen6 — 논문 근거 기반 위치 최적화
- **논문 근거:**
  - [PMC6152672](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6152672/) — 아다만탄 AChE 억제제: 페닐 3번(메타) 전자흡인기 최적
  - [Grychowska 2022, EJMECH](https://pubmed.ncbi.nlm.nih.gov/35397400/) — F 위치에 따라 hERG + 타겟 결합 크게 달라짐
  - [Molecules 2023](https://www.mdpi.com/1420-3049/28/2/490) — F 위치 변화만으로 활성 최대 1300배 차이 (activity cliff)
- **가설:** 4-CN(파라) + 2-F(오르토) = 논문 기준 최적 조합
- **결과:** ether_direct__4cn_2f RMSD=2.26Å ⭐⭐ 안정 (Gen4 2.66Å 대비 개선)

---

## 세대별 핵심 지표

| 세대 | 통과율 | 주요 학습 | 최고 후보 | hERG | ΔG | RMSD |
|------|--------|-----------|-----------|------|-----|------|
| Gen1 | 91% | chlorophenyl 제거 | pyridyl_4 | 0.71 ❌ | — | — |
| Gen2 | 95% | ether_ethyl 링커 발견 | pyridyl_4 | 0.71 ❌ | -3.547 | 2.97Å |
| Gen3 | 89% | bipyridyl 강직성 | bipyridyl | 0.80 ❌ | -4.191 | 2.37Å |
| Gen4 | 71% | hERG 해결 (cyano+F) | cyano_F | **0.39 ✅** | **-4.925** | 2.66Å |
| Gen5 | 50% | F 트레이드오프 확인 | cyano_F | 0.39 ✅ | -3.757 | 2.95Å |
| Gen6 | 80% | 논문 근거 위치 최적화 | 4cn_2f | 0.40 ✅ | -3.849 | **2.26Å** |

---

## Gen6 최종 후보

```
ether_direct__4cn_2f
SMILES: OC12CC3(Oc4ccc(C#N)cc4F)CC(C1)CC(C2)C3

물성:   MW=301.36  LogP=3.55  QED=0.909
독성:   hERG=0.40  DILI=0.29  BBB=0.84  AMES=0.04
합성:   SA Score=4.77 (보통)
도킹:   ΔG=-3.849 kcal/mol (AChE 결합)
MD:     RMSD=2.26Å ⭐⭐ 안정
```

---

## 파일 구조

```
├── epsilon_gate_gen2.py          # Gen2 물성 필터
├── epsilon_gate_gen2_docking.py  # 도킹 파이프라인
├── epsilon_gate_md.py            # MD 안정성 검증
├── epsilon_gate_tox.py           # 독성 + SA Score 필터
├── epsilon_gate_gen3.py          # Gen3 설계
├── epsilon_gate_gen4.py          # Gen4 설계 (hERG 집중)
├── epsilon_gate_gen5.py          # Gen5 설계 (MD 안정성 집중)
├── epsilon_gate_gen6.py          # Gen6 설계 (논문 근거 위치 최적화)
├── pattern_feedback_gen1.md      # Gen1 실패 패턴
├── pattern_feedback_gen2.md      # Gen2 실패 패턴
├── pattern_feedback_gen3.md      # Gen3 실패 패턴
├── pattern_feedback_gen4.md      # Gen4 실패 패턴
├── LOOP_RECORD.md                # 전체 루프 기록
└── 1EVE_pocket.pdbqt             # AChE 결합 포켓 (하드코딩)
```

---

## 설치 및 실행

```bash
pip install rdkit vina meeko gemmi openmm openmmforcefields admet-ai pandas
```

```bash
# 물성 필터
python epsilon_gate_gen6.py

# 독성 + SA Score 필터
python epsilon_gate_tox.py --csv gen6_library.csv

# 도킹 (1EVE.pdb 필요 — 코랩 권장)
python epsilon_gate_gen2_docking.py

# MD 안정성 (10ps 기본, N_STEPS 변경으로 스케일 조정)
python epsilon_gate_md.py --csv gen6_final_results.csv
```

### MD 타임스케일 변경

```python
# epsilon_gate_md.py 내 한 줄만 수정
N_STEPS = 20_000     # 10ps  (기본, 이 환경)
N_STEPS = 200_000    # 100ps (코랩 CPU)
N_STEPS = 2_000_000  # 1ns   (코랩 GPU 권장)
```

---

## 다른 질병에 적용하기

타겟 단백질만 바꾸면 됩니다.

```python
# epsilon_gate_gen2_docking.py 내 3줄 수정
PDB_ID       = "6LU7"          # SARS-CoV-2 Mpro (코로나)
POCKET_CENTER = (x, y, z)      # 해당 PDB 결합 포켓 좌표
# → 나머지 파이프라인은 그대로
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
- Receptor: 포켓 부분구조(53원자) → 코랩 + 전체 1EVE로 절대값 검증 필요
- MD: 진공 10ps → 명시적 수분자(TIP3P) + 100ns 필요
- wet-lab IC₅₀ 실측 없음

**다음 단계**
- [ ] 코랩에서 전체 1EVE 도킹 (실제 ΔG -8~-11 kcal/mol 범위 확인)
- [ ] AiZynthFinder 합성 경로 탐색
- [ ] Gen7: 4cn_2f 기반 변형 탐색

---

## 참고 논문

1. **아다만탄 AChE 억제제 구조-활성 관계**  
   PMC6152672 — https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6152672/

2. **F 도입으로 hERG 감소: MAO-B 억제제 사례**  
   Grychowska et al. (2022), Eur J Med Chem  
   https://pubmed.ncbi.nlm.nih.gov/35397400/

3. **F 위치 변화와 activity cliff (최대 1300배)**  
   Molecules 2023, 28(2), 490  
   https://www.mdpi.com/1420-3049/28/2/490

4. **AChE 억제제 도킹 연구 리뷰 (2018-2022)**  
   Molecules 2023, 28(3), 1084  
   https://pmc.ncbi.nlm.nih.gov/articles/PMC9921523/

---

## 만든 사람

비전공자, 고졸, AI 경력 6개월.  
Gen0 수작업부터 Gen6 논문 근거 설계까지 — 학위 없이 AI로 구동한 신약 탐색 엔진.

> *"파이프라인이 아니라 엔진입니다. 실패가 다음 설계를 만듭니다."*
