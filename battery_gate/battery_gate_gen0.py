"""
Battery-Gate Gen0
─────────────────
Epsilon-Gate 방법론의 배터리 버전
타겟: 아르지로다이트계 황화물 전해질 Li6PS5X (X=Cl/Br/I)
목표: 논문 실패 패턴 역산 → 최적 조성 공간 정의

Gen0 핵심 질문:
"왜 안 됐는가?"를 먼저 정리한다.
실패하지 않는 공간에서 후보를 찾는다.
"""

# ── 1. Gen0: 논문에서 추출한 실패 패턴 ──────────────────────────────────────

FAILURE_PATTERNS = {
    "X_substitution": {
        "I_fail": {
            "description": "X=I 단독 치환 → 이온전도도 수 자릿수 급락",
            "reason": "I- 이온 반경 과대 → Li+ 확산 경로 차단",
            "conductivity": "< 1e-6 S/cm (Li6PS5Cl 대비 수 자릿수 낮음)",
            "rule": "I 단독 사용 금지. Cl/Br 혼합만 허용"
        },
        "Br_optimal": {
            "description": "Cl/Br 혼합 → 최적 전도도 구간 존재",
            "best": "Li6PS5Cl0.5Br0.5 → 균형점",
            "rule": "Cl:Br = 0.3~0.7 범위 탐색"
        }
    },
    "doping": {
        "Cl_excess": {
            "description": "Cl 과잉 (>1.5) → 이차상(MgS 등) 형성",
            "threshold": "Cl content ≤ 1.5",
            "rule": "Cl 상한선 1.5 강제"
        },
        "Al_alone": {
            "description": "Al 단독 도핑 → 효과 제한적",
            "fix": "Al+Cl 공도핑 시 4.7배 향상 (7.29 mS/cm)",
            "rule": "단일 원소 도핑보다 공도핑 우선"
        },
        "O_doping": {
            "description": "O 도핑 최적점 x=0.25, 과잉 시 역효과",
            "best": "Li6PS4.75ClO0.25 → 4.7 mS/cm",
            "rule": "O 도핑량 0.1~0.3 범위"
        },
        "Mg_doping": {
            "description": "Mg 도핑: Cl 과잉과 결합 시 이차상",
            "rule": "Mg ≤ 0.125, Cl ≤ 1.5 동시 적용"
        }
    },
    "processing": {
        "pressure": {
            "description": "소결 압력 부족 → 계면 저항 급등",
            "minimum": "300 MPa 이상 권장",
            "rule": "합성 가능성 평가 시 압력 조건 포함"
        },
        "moisture": {
            "description": "수분 노출 → H2S 발생, 구조 붕괴",
            "rule": "드라이룸 합성 필수 → 합성 난이도 가산점"
        }
    }
}

# ── 2. 실패 패턴으로 정의한 탐색 공간 ─────────────────────────────────────────

SEARCH_SPACE = {
    "base": "Li6PS5X",
    "X_range": {
        "Cl": (0.3, 1.5),   # Cl 허용 범위
        "Br": (0.3, 0.7),   # Br 허용 범위
        "I":  (0.0, 0.0),   # I 금지
    },
    "dopants": {
        "Al": (0.05, 0.15),  # Al 도핑 범위
        "O":  (0.10, 0.30),  # O 도핑 범위
        "Mg": (0.05, 0.125), # Mg 도핑 범위
        "Si": (0.05, 0.15),  # Si 도핑 (신규 탐색)
        "Ge": (0.00, 0.00),  # Ge 금지 (비용)
    },
    "gate": {
        "ionic_conductivity_min": 3.0,   # mS/cm 최소 기준
        "Ehull_max": 0.05,               # eV/atom (열역학 안정성)
        "ESW_min": 2.0,                  # V (전기화학 안정성 윈도우)
        "synthesis_difficulty_max": 3,   # 1~5 척도
    }
}

# ── 3. Gen0 후보 조성 생성 ────────────────────────────────────────────────────

import itertools

def generate_gen1_candidates():
    """
    실패 패턴 역산으로 Gen1 후보 조성 생성
    Epsilon-Gate Gen1과 동일한 로직
    """
    candidates = []

    # 조합 1: Cl/Br 혼합 + Al/O 공도핑
    cl_range = [0.5, 0.75, 1.0, 1.25]
    br_range = [0.3, 0.5, 0.7]
    al_range = [0.0, 0.05, 0.10]
    o_range  = [0.0, 0.15, 0.25]

    for cl, br, al, o in itertools.product(cl_range, br_range, al_range, o_range):
        # 실패 규칙 적용
        if cl > 1.5: continue          # Cl 상한
        if cl + br > 1.8: continue     # 총 할로겐 과잉 방지
        if al > 0 and o > 0:           # 공도핑 우선
            label = f"Li6PS5Cl{cl}Br{br}_Al{al}O{o}"
        elif al > 0:
            label = f"Li6PS5Cl{cl}Br{br}_Al{al}"
        elif o > 0:
            label = f"Li6PS5Cl{cl}Br{br}_O{o}"
        else:
            label = f"Li6PS5Cl{cl}Br{br}_base"

        candidates.append({
            "id": label,
            "Cl": cl,
            "Br": br,
            "Al": al,
            "O":  o,
            "I":  0.0,
            "priority": "high" if (al > 0 and o > 0) else "medium"
        })

    return candidates

candidates = generate_gen1_candidates()

# ── 4. Gen0 결과 출력 ─────────────────────────────────────────────────────────

print("=" * 60)
print("Battery-Gate Gen0 — 실패 패턴 추출 완료")
print("=" * 60)

print("\n[실패 규칙 요약]")
print(f"  ❌ X=I 단독: 금지 (전도도 수 자릿수 급락)")
print(f"  ❌ Cl > 1.5: 금지 (이차상 형성)")
print(f"  ❌ Al 단독: 비권장 (공도핑 우선)")
print(f"  ❌ Ge 도핑: 금지 (비용 과다)")
print(f"  ✅ Cl/Br 혼합: Cl:Br = 0.3~0.7")
print(f"  ✅ Al+O 공도핑: 최우선 탐색")
print(f"  ✅ O 도핑: 0.10~0.30 범위")

print(f"\n[Gen1 후보 조성 수]: {len(candidates)}개")
print(f"  - High priority (공도핑): {sum(1 for c in candidates if c['priority']=='high')}개")
print(f"  - Medium priority (단일): {sum(1 for c in candidates if c['priority']=='medium')}개")

print(f"\n[Gen1 상위 후보 샘플 (High priority 10개)]")
high = [c for c in candidates if c['priority'] == 'high'][:10]
for i, c in enumerate(high, 1):
    print(f"  {i:2}. {c['id']}")
        
print(f"\n[Gate 기준]")
for k, v in SEARCH_SPACE['gate'].items():
    print(f"  {k}: {v}")

print("\n[다음 단계: Gen1]")
print("  Materials Project API → 각 후보 Ehull / ESW 계산")
print("  → 통과 후보만 CHGNet MD → 이온전도도 예측")
print("  → 병목 감지 → Gen2 설계")
print("=" * 60)
