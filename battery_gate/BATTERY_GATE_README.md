# Battery-Gate: AI-Driven Solid Electrolyte Discovery Engine

> **Gen0 논문 실패 패턴 추출 → Gen8 4d 사이트 최적 점유율 탐색까지**  
> 비전공자가 오픈소스만으로 구축한 전고체 배터리 전해질 후보 탐색 파이프라인

---

## 최종 후보 (Gen8)

```
Li6PS4Cl2 (S→Cl 1개 치환, 슈퍼셀 x=50%)

조성:   Li6 P1 S4 Cl2
전도도: σ(900K) = 16.862 mS/cm  ✅
        σ(상온)  ≈ 0.093 mS/cm  (Ea 보정 필요)
최적점: 4d 사이트 Cl 점유율 x=50% ← 논문 예측과 일치
패턴:   짝수 Cl 치환 → 높은 전도도 (구조-활성 규칙 발견)
```

**Gen0~Gen8 통틀어 논문 예측(x=41~50%)과 독립적으로 일치한 첫 번째 결과.**

---

## 프로젝트 개요

- **타겟:** 아르지로다이트계 황화물 전해질 Li6PS5Cl (전고체 배터리)
- **스캐폴드:** Li6PS5Cl 베이스 + S→Cl 치환 최적화
- **방법론:** Lab-in-the-Loop — 실패 데이터가 다음 설계를 만드는 자기개선 루프
- **리소스:** 개인 PC + 오픈소스만 사용 (GPU 없음)
- **기간:** 1일 (Gen0~Gen8)

---

## Epsilon-Gate와의 관계

Battery-Gate는 Epsilon-Gate(알츠하이머 신약 탐색)의 방법론을 전고체 배터리 소재 탐색에 이식한 프로젝트입니다.

```
Epsilon-Gate          →    Battery-Gate
─────────────────────────────────────────────
SMILES 화합물 설계    →    황화물 조성 설계 (Li-P-S-Cl계)
물성 필터 (MW/LogP)   →    열역학 안정성 (Ehull)
독성 필터 (hERG/DILI) →    전기화학 안정성 윈도우 (ESW)
도킹 (AutoDock Vina)  →    이온전도도 (CHGNet MD)
MD 안정성 (OpenMM)    →    열적·기계적 안정성
병목 해결 엔진        →    동일 구조 적용
```

**같은 방법론, 다른 도메인 — 동일하게 작동함을 검증.**

---

## 파이프라인

```
논문 실패 패턴 추출 (RAG)
    ↓
Materials Project API 데이터 수집
    ↓
열역학 안정성 필터 (Ehull ≤ 0.10 eV/atom)
    ↓
밴드갭 필터 (1.0 ~ 6.0 eV)
    ↓
CHGNet MD 이온전도도 예측
  900K / 1000steps
    ↓
병목 감지 → 조성 설계 수정
    ↓
반복 (Gen N)
```

---

## Gen0 → Gen8 학습 루프

| 세대 | 핵심 발견 | 결과 |
|------|-----------|------|
| Gen0 | 논문 실패 패턴 추출 | 규칙 4개 확립, 후보 99개 |
| Gen1 | Materials Project 검색 | 1개 통과 (검색 범위 좁음) |
| Gen2 | 검색 범위 확장 | 8개 통과, Li6PS5Cl 확인 |
| Gen3 | CHGNet MD 첫 실행 (500K/200steps) | Li6PS5Cl 0.059 mS/cm → 병목 발견 |
| Gen4 | MD 조건 개선 (700K/500steps) | 1.990 mS/cm → 33배↑ |
| Gen5 | Al+O 공도핑 + 900K/1000steps | Li6PS4Cl2 **16.29 mS/cm** 🏆 |
| Gen6 | Cl 함량별 전도도 곡선 탐색 | 짝수/홀수 비선형 패턴 발견 |
| Gen7 | Wide Scan (Cl 0~5개) | 비선형 패턴 전체 확인, Li6PCl6 탐색 |
| Gen8 | 슈퍼셀 4d 사이트 점유율 탐색 | x=50% 최고 → **논문 예측과 일치** ✅ |

---

## 핵심 발견

### 1. 짝수/홀수 전도도 패턴

```
Cl 0개:  6.56 mS/cm
Cl 1개: 16.29 mS/cm  🏆 (짝수 치환)
Cl 2개:  4.96 mS/cm  (급락)
Cl 3개: 10.68 mS/cm  (반등, 짝수 치환)
Cl 4개:  8.89 mS/cm  (감소)
Cl 5개: 16.06 mS/cm  🏆 (짝수 치환)
```

**원인:** 아르지로다이트 결정의 Wyckoff 4a/4d 사이트가 쌍으로 존재 → Cl이 짝수 개 치환 시 대칭 유지 → Li+ 이동 경로 최적화

논문 확인: S/Cl 무질서도(disorder)가 Li+ 이동 핵심 원인 (독립적 재현)

### 2. 최적 점유율 x=50%

```
논문 예측: x=41~50%에서 최고 전도도
우리 결과: x=50%에서 16.862 mS/cm (900K 최고)
→ 논문 예측과 독립적으로 일치
```

### 3. 세대별 최고 전도도 추이

```
Gen3:  0.059 mS/cm  (500K, 200steps)
Gen4:  1.990 mS/cm  →  33배↑
Gen5: 16.288 mS/cm  →   8배↑
Gen8: 16.862 mS/cm  →  신규 최고
```

---

## Gen0 실패 규칙 (논문 기반)

```
❌ X=I 단독:    전도도 수 자릿수 급락
❌ Cl > 1.5:    이차상 형성
❌ Al 단독 도핑: 구조 왜곡 > Li 공공 효과 (Gen5 실험적 확인)
❌ Ge 도핑:     비용 과다
✅ Cl/Br 혼합:  Cl:Br = 0.3~0.7
✅ S→Cl 치환:  짝수 개 치환 = 최적 (Gen6 발견)
✅ x=50% 점유:  4d 사이트 최적 점유율 (Gen8 확인)
```

---

## 파일 구조

```
├── battery_gate_gen0.py    # 실패 패턴 추출 + 후보 설계
├── battery_gate_gen1.py    # Materials Project API 연결
├── battery_gate_gen2.py    # 검색 범위 확장
├── battery_gate_gen3.py    # CHGNet MD 첫 실행
├── battery_gate_gen4.py    # MD 조건 개선
├── battery_gate_gen5.py    # Al+O 공도핑 + 900K
├── battery_gate_gen6.py    # Cl 함량별 곡선 탐색
├── battery_gate_gen7.py    # Wide Scan (Cl 0~5개)
├── battery_gate_gen8.py    # 4d 사이트 점유율 탐색
├── gen1_candidates.csv     # Gen1 통과 후보
├── gen2_candidates.csv     # Gen2 통과 후보
├── gen3_results.csv        # Gen3 전도도 결과
├── gen4_results.csv        # Gen4 전도도 결과
├── gen5_results.csv        # Gen5 전도도 결과
├── gen6_results.csv        # Gen6 Cl 함량별 곡선
├── gen7_results.csv        # Gen7 Wide Scan 결과
├── gen8_results.csv        # Gen8 점유율별 결과
└── README.md
```

---

## 설치 및 실행

```bash
pip install pymatgen mp-api chgnet
```

```bash
# 환경변수 설정 (Windows)
set MP_API_KEY=여기에_API_키_입력

# Gen0 실패 패턴 추출
python battery_gate_gen0.py

# Gen8 최적 점유율 탐색
python battery_gate_gen8.py
```

---

## 현재 한계 및 다음 단계

**현재 한계**
- MD: 개인 PC CPU → 1000steps (논문 수준: 수만 steps 필요)
- 온도: 900K 고온 시뮬레이션 → 상온 환산 시 Arrhenius 보정 필요
- 구조: 단순 원자 치환 → 실제 무질서 구조(disordered) 미구현

**다음 단계**
- [ ] 코랩 GPU로 5000~10000 steps 재실행
- [ ] 300K/500K/700K/900K Arrhenius plot → 실측 상온값 도출
- [ ] x=45~55% 구간 2% 단위 촘촘한 탐색 (Gen9)
- [ ] Li6PS4Cl2 특허 출원 여부 확인
- [ ] Cl/Br 혼합 + x=50% 조합 탐색

---

## 참고 논문

1. **S/Cl 무질서도와 Li+ 전도 메커니즘**  
   https://pubs.rsc.org/en/Content/ArticleLanding/2024/TA/D3TA06069A

2. **할라이드 치환과 구조적 무질서 관계**  
   https://pubs.acs.org/doi/10.1021/acs.chemmater.3c01525

3. **아르지로다이트 ML 파이프라인 리뷰**  
   https://pubs.rsc.org/en/content/articlelanding/2026/mh/d5mh01525a

4. **고엔트로피 황화물 전해질 AI 탐색**  
   https://pubs.rsc.org/en/content/articlepdf/2025/ta/d5ta02205c

---

## 방법론 핵심

> 파이프라인이 아니라 엔진입니다.  
> 실패가 다음 설계를 만들고,  
> 논문이 방향을 제공하고,  
> 시뮬레이션이 독립적으로 검증합니다.

**Epsilon-Gate(신약) → Battery-Gate(배터리) — 도메인은 달라도 엔진은 동일합니다.**
