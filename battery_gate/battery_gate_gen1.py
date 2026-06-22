"""
Battery-Gate Gen1
─────────────────
Materials Project API로 아르지로다이트 황화물 데이터 수집
Gen0 실패 패턴 기반 필터링 적용
목표: Ehull / 밴드갭 / 전기화학 안정성 필터 통과 후보 추출
"""

import os
import json
import csv
from mp_api.client import MPRester

# ── API 키 환경변수에서 로드 ───────────────────────────────────────────────────
API_KEY = os.environ.get("MP_API_KEY")
if not API_KEY:
    raise ValueError("MP_API_KEY 환경변수가 설정되지 않았습니다.")

# ── Gen0 Gate 기준 ────────────────────────────────────────────────────────────
GATE = {
    "ehull_max":   0.05,   # eV/atom 이하 → 열역학 안정
    "bandgap_min": 1.0,    # eV 이상 → 전자 절연체 (전해질 필수 조건)
    "bandgap_max": 6.0,    # eV 이하 → 너무 넓으면 합성 어려움
}

print("=" * 60)
print("Battery-Gate Gen1 — Materials Project 데이터 수집")
print("=" * 60)

with MPRester(API_KEY) as mpr:

    # ── 1. 아르지로다이트 구조 검색 ───────────────────────────────────────────
    print("\n[Step 1] Li-P-S-Cl 계열 검색 중...")

    # Li, P, S, Cl 포함 황화물 계열 검색
    results = mpr.materials.summary.search(
        elements=["Li", "P", "S", "Cl"],
        num_elements=(4, 4),          # 정확히 4원소
        energy_above_hull=(0, 0.1),   # 안정~준안정 범위
        fields=[
            "material_id",
            "formula_pretty",
            "energy_above_hull",
            "band_gap",
            "volume",
            "density",
            "nsites",
        ]
    )

    print(f"  검색 결과: {len(results)}개")

    # ── 2. 추가: Li-P-S-Br 계열 ──────────────────────────────────────────────
    print("\n[Step 2] Li-P-S-Br 계열 검색 중...")

    results_br = mpr.materials.summary.search(
        elements=["Li", "P", "S", "Br"],
        num_elements=(4, 4),
        energy_above_hull=(0, 0.1),
        fields=[
            "material_id",
            "formula_pretty",
            "energy_above_hull",
            "band_gap",
            "volume",
            "density",
            "nsites",
        ]
    )

    print(f"  검색 결과: {len(results_br)}개")

    # ── 3. 합치기 ─────────────────────────────────────────────────────────────
    all_results = list(results) + list(results_br)
    print(f"\n[Step 3] 전체 후보: {len(all_results)}개")

    # ── 4. Gen0 Gate 필터 적용 ────────────────────────────────────────────────
    print(f"\n[Step 4] Gate 필터 적용 중...")
    print(f"  기준: Ehull ≤ {GATE['ehull_max']} eV/atom")
    print(f"  기준: 밴드갭 {GATE['bandgap_min']} ~ {GATE['bandgap_max']} eV")

    passed = []
    failed_ehull = 0
    failed_bandgap = 0

    for r in all_results:
        ehull = r.energy_above_hull
        bgap  = r.band_gap if r.band_gap is not None else 0

        # Ehull 필터
        if ehull is None or ehull > GATE["ehull_max"]:
            failed_ehull += 1
            continue

        # 밴드갭 필터
        if bgap < GATE["bandgap_min"] or bgap > GATE["bandgap_max"]:
            failed_bandgap += 1
            continue

        passed.append({
            "material_id":  r.material_id,
            "formula":      r.formula_pretty,
            "ehull":        round(ehull, 4),
            "band_gap":     round(bgap, 3),
            "density":      round(r.density, 3) if r.density else None,
            "nsites":       r.nsites,
        })

    # ── 5. 결과 출력 ──────────────────────────────────────────────────────────
    print(f"\n[결과]")
    print(f"  전체:            {len(all_results)}개")
    print(f"  Ehull 탈락:      {failed_ehull}개")
    print(f"  밴드갭 탈락:     {failed_bandgap}개")
    print(f"  ✅ Gate 통과:    {len(passed)}개")

    if passed:
        # Ehull 기준 정렬
        passed.sort(key=lambda x: x["ehull"])

        print(f"\n[Gate 통과 후보 상위 20개]")
        print(f"{'#':>3} {'Formula':<30} {'Ehull':>8} {'BandGap':>9} {'Density':>9}")
        print("-" * 65)
        for i, p in enumerate(passed[:20], 1):
            print(f"{i:>3} {p['formula']:<30} {p['ehull']:>8.4f} {p['band_gap']:>9.3f} {p['density'] or 'N/A':>9}")

        # CSV 저장
        with open("gen1_candidates.csv", "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=passed[0].keys())
            writer.writeheader()
            writer.writerows(passed)
        print(f"\n✅ gen1_candidates.csv 저장 완료 ({len(passed)}개)")

    # ── 6. 병목 분석 (Gen0 → Gen1 학습) ─────────────────────────────────────
    print(f"\n[병목 분석]")
    total = len(all_results)
    if total > 0:
        ehull_rate  = (1 - failed_ehull/total) * 100
        bandgap_rate = (1 - failed_bandgap/(total - failed_ehull)) * 100 if (total - failed_ehull) > 0 else 0
        pass_rate   = len(passed)/total * 100
        print(f"  Ehull 통과율:   {ehull_rate:.1f}%")
        print(f"  밴드갭 통과율:  {bandgap_rate:.1f}%")
        print(f"  최종 통과율:    {pass_rate:.1f}%")

        if ehull_rate < bandgap_rate:
            print(f"\n  → 병목: 열역학 안정성 (Ehull)")
            print(f"     Gen2 방향: 더 안정한 조성 탐색 (Ehull 완화 또는 도핑 전략 수정)")
        else:
            print(f"\n  → 병목: 밴드갭 (전자 절연성)")
            print(f"     Gen2 방향: 전자 절연성 높은 원소 도핑 우선")

    print("\n[다음 단계: Gen2]")
    print("  통과 후보 → CHGNet MD → 이온전도도 예측")
    print("  병목 원인 → Gen2 조성 설계 수정")
    print("=" * 60)
