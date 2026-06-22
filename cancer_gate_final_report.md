# Cancer-Gate 최종 리포트
**작성일:** 2026-06-19
**파이프라인:** Gen0 → Gen15 (GRM4 ADC 링커 트리거 타깃)

---

## 1. 프로젝트 개요

**목표:** GRM4를 ADC 링커 트리거로 활용하는 유방암 치료 후보물질 발굴
**방법론:** Epsilon-Gate (AChE) / Battery-Gate (Li6PS5Cl) 방법론 적용
**타깃:** mGluR4 (GRM4) PAM 포켓 — PDB 7E9H Chain R
**핵심 전략:** Reverse inference (실패 약물에서 viable core 채굴)

---

## 2. 파이프라인 구성

| 단계 | 도구 | 버전 |
|---|---|---|
| SMILES 설계 | RDKit | 2026.03.3 |
| 3D 최적화 | MMFF94s | — |
| ADMET 예측 | ADMET-AI | — |
| 분자 도킹 | AutoDock Vina | exhaustiveness=48 |
| SA Score | RDKit SAS | — |
| MD 시뮬레이션 | OpenMM + OpenFF Sage 2.2.0 | 8.1.1 |
| 선택성 스크리닝 | AutoDock Vina (역도킹) | 5CGC mGluR5 |

---

## 3. Gen0~Gen15 진행 요약

### GO 확정 화합물 전체

| 세대 | ID | Ki (μM) | DILI | hERG | LogP | 설계 |
|---|---|---|---|---|---|---|
| Gen9  | CG9-005  | 0.725 | 0.356 | 0.239 | — | CF2+3,3-diF pip |
| Gen10 | CG10-005 | 0.345 | 0.349 | 0.259 | 3.67 | CF2+2-Me+3,3-diF pip ★Ki최고 |
| Gen11 | CG11-001 | 0.713 | 0.325 | 0.287 | — | — |
| Gen12 | CG12-004 | 0.580 | 0.312 | 0.154 | 2.64 | CF2+2-Me+3,3-diF+4-OH pip ★최종선택 |
| Gen14 | CG14-001 | 0.982 | 0.328 | 0.137 | 2.35 | CF2+2-Me+3-monoF+5-OH pip |

### 경계선 FAIL (아까운 후보)

| 세대 | ID | Ki (μM) | DILI | hERG | 실패 원인 |
|---|---|---|---|---|---|
| Gen11 | CG11 | 0.251 | — | 0.350 | gem-diMe → hERG FAIL |
| Gen13 | CG13-B1 | 0.372 | 0.404 | 0.233 | 4-F pip → DILI 0.004 초과 |
| Gen13 | CG13-A1 | 0.478 | 0.407 | 0.126 | 5-OH → DILI 0.007 초과 |

---

## 4. 확정된 SAR (Structure-Activity Relationship)

### 핵심 구조 법칙

```
[CF2-cyclopropane 코어] — 절대 포기 불가
  └ CF2 → CHF 치환 시 Ki 6배 후퇴 (CG14-003 실증)
  └ Walsh 궤도 강화가 GRM4 포켓 정합성의 핵심

[pip 3,3-diF] — Ki의 핵심 기여자
  └ gem-diF > vicinal-diF (Ki 기준)
  └ 3,3-diF → 3-monoF: Ki 후퇴 (0.580 → 0.982)

[pip 4-OH] — hERG 보호 역할
  └ 4-OH 추가 시 hERG 0.259 → 0.154 (CG10-005 → CG12-004)
  └ H-bond donor 효과 — MD RMSD 0.14Å로 포켓 고정 물리적 확인

[pip 4-F] — DILI 장벽 원인
  └ 4-F 추가 시 DILI 0.349 → 0.404 (임계값 0.400 초과)
  └ Gen13-B1, Gen15-003 동일 패턴으로 확정

[pip C4 = DILI 최적 부착점]
  └ 4-OH: DILI OK / 5-OH, 6-OH: DILI FAIL
  └ C4가 DILI 임계 부착점으로 확정
```

### 3축 트레이드오프 구조

```
Ki 개선  ←→  hERG 악화  (LogP 상관)
hERG 개선 ←→  Ki 후퇴   (OH 추가 비용)
DILI 경계  ←→  pip F 총량 민감
```

---

## 5. MD 시뮬레이션 결과

### MD 결과 비교 (OpenFF Sage 2.2.0, 진공)

| 후보 | FAST RMSD | FULL RMSD | 안정성 |
|---|---|---|---|
| CG10-005 | 0.15Å | 1.73Å | ✅ STABLE |
| CG12-004 | 0.20Å | **0.14Å** | ✅ STABLE ★ |

CG12-004 FULL RMSD 0.14Å — CG10-005(1.73Å) 대비 12배 안정.
4-OH 수소결합이 GRM4 포켓 내에서 분자를 실제로 고정하는 물리적 증거.

---

## 6. 선택성 스크리닝 (mGluR4 vs mGluR5 역도킹)

**기준:** SI = Ki(mGluR5) / Ki(mGluR4) ≥ 10x → 선택적 PAM 판정

| 후보 | Ki(mGluR4) | Ki(mGluR5) | SI | 판정 |
|---|---|---|---|---|
| CG10-005 | 0.345μM | 9.843μM | **28.5x** | ✅ |
| CG12-004 | 0.580μM | 11.287μM | **19.5x** | ✅ |
| CG14-001 | 0.982μM | 17.749μM | **18.1x** | ✅ |

전 후보 SI ≥ 10x 통과 — mGluR5 교차반응성 리스크 없음 확인.

---

## 7. 종합 프로파일 및 최종 선택

### 최종 3종 비교

| 후보 | Ki(μM) | hERG | DILI | FULL RMSD | SI | 종합 |
|---|---|---|---|---|---|---|
| CG10-005 | **0.345** | 0.259⚠️ | 0.349 | 1.73Å | **28.5x** | Ki·SI 최고, hERG 경계 |
| **CG12-004** | 0.580 | **0.154** | **0.312** | **0.14Å** | 19.5x | **MD·hERG 최고, 균형 최적** |
| CG14-001 | 0.982 | **0.137** | 0.329 | — | 18.1x | hERG 최우수, Ki 약함 |

### ★ 최종 선택: CG12-004

```
SMILES: OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1C(C)C(F)(F)CC(O)C1
설계:   CF2-cyclopropane + 2-Me + 3,3-diF + 4-OH pip
```

| 지표 | 값 | 판정 |
|---|---|---|
| ΔG | -8.503 kcal/mol | — |
| Ki | 0.580 μM | ✅ |
| DILI | 0.312 | ✅ |
| hERG | 0.154 | ✅ 최고 |
| LogP | 2.64 | ✅ |
| MW | 397.3 | ✅ |
| FULL RMSD | 0.14Å | ✅ 최안정 |
| SI (mGluR5) | 19.5x | ✅ 선택적 |

**선택 근거:** ADC 링커 트리거 용도에서 결합 안정성이 Ki보다 우선. FULL RMSD 0.14Å + hERG 0.154 + SI 19.5x — 안전성·선택성·안정성 3축 동시 최적.

---

## 8. 파이프라인 완료 현황

### 완료
- ✅ Gen0~Gen15 SMILES 설계 및 도킹 (75개 이상)
- ✅ ADMET 필터링 (DILI / hERG / HIA)
- ✅ SAR 확정 (5개 구조 법칙)
- ✅ MD 시뮬레이션 (결합 안정성 RMSD)
- ✅ 역도킹 선택성 스크리닝 (mGluR5 SI)
- ✅ 최종 후보 확정: CG12-004

### 추가 검증 (GPU 환경 필요)
- ⬜ 전체 단백질+물 NPT MD — 더 엄밀한 결합 안정성
- ⬜ FEP (자유에너지 계산) — Ki 정밀 예측

### 다음 단계
1. **preprint 작성** — Epsilon-Gate 형식 준용, SMILES 비공개
2. **GPU MD** (선택) — Colab Pro or 로컬 NVIDIA, 논문급 완성도

---

## 9. 포트폴리오 가치

Cancer-Gate의 실질적 가치는 SMILES 구조 자체가 아니라 **파이프라인을 독립적으로 설계하고 실행했다는 사실**에 있다.

**구현 스택 (전부 오픈소스)**
RDKit / ADMET-AI / AutoDock Vina / OpenMM / OpenFF

**파이프라인 설계 역량**
- Gen0~Gen15 반복 최적화 — 75개 이상 화합물 스크리닝
- SAR 도출 → 병목 판단 → 방향 전환 전 과정 독립 설계
- MD 시뮬레이션, 역도킹 선택성 스크리닝까지 통합

**업계 맥락**
대부분의 AI 신약 스타트업은 Schrödinger / OpenEye 같은 상용 플랫폼을 구매하거나 CRO에 외주한다. 오픈소스만으로 동등한 파이프라인을 설계·실행할 수 있는 역량은 희소하다.

**다음 단계**
- Epsilon-Gate 프리프린트 (DOI 10.20944/preprints202501.1982.v1) + Cancer-Gate 리포트로 포트폴리오 구성
- 실험 검증 (in vitro IC50, 세포 독성) 확보 후 IP 전략 수립

---

*Cancer-Gate pipeline by 천운기 (천타프) — GRM4 ADC linker trigger approach*
*Methodology: Epsilon-Gate → Battery-Gate → Cancer-Gate iterative optimization*
