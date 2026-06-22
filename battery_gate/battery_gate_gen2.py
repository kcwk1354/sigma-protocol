"""
Battery-Gate Gen2
─────────────────
Gen1 병목 분석 결과: Ehull 탈락 67% → 검색 조건이 너무 좁았음
해결: 
  1. 원소 조합 확장 (4원소 → 3원소 포함, Ehull 범위 완화)
  2. 실제 아르지로다이트 구조 포함 (Li-P-S 베이스 + 할로겐)
  3. chemsys 방식으로 전환 (더 넓은 탐색)
"""

import os
import csv
from mp_api.client import MPRester

API_KEY = os.environ.get("MP_API_KEY")
if not API_KEY:
    raise ValueError("MP_API_KEY 환경변수가 설정되지 않았습니다.")

# ── Gen2 Gate (Ehull 완화) ────────────────────────────────────────────────────
GATE = {
    "ehull_max":   0.10,   # Gen1 0.05 → Gen2 0.10 완화 (준안정 포함)
    "bandgap_min": 1.0,
    "bandgap_max": 6.0,
    "density_max": 3.5,    # g/cm3 이하 (너무 무거우면 합성 어려움)
}

print("=" * 60)
print("Battery-Gate Gen2 — 확장 탐색")
print("=" * 60)
print(f"\n[Gen1 → Gen2 변경사항]")
print(f"  Ehull 기준 완화: 0.05 → 0.10 eV/atom")
print(f"  검색 방식: 4원소 고정 → chemsys 확장")
print(f"  탐색 계열 추가: Li-P-S-Cl-Br 혼합계")

all_results = []

with MPRester(API_KEY) as mpr:

    # ── 검색 1: Li-P-S-Cl (Ehull 완화) ───────────────────────────────────────
    print(f"\n[Step 1] Li-P-S-Cl 계열 (Ehull 완화)...")
    r1 = mpr.materials.summary.search(
        chemsys="Li-P-S-Cl",
        energy_above_hull=(0, 0.15),
        fields=["material_id","formula_pretty","energy_above_hull",
                "band_gap","density","nsites","elements"]
    )
    print(f"  결과: {len(r1)}개")
    all_results.extend(r1)

    # ── 검색 2: Li-P-S-Br ────────────────────────────────────────────────────
    print(f"\n[Step 2] Li-P-S-Br 계열...")
    r2 = mpr.materials.summary.search(
        chemsys="Li-P-S-Br",
        energy_above_hull=(0, 0.15),
        fields=["material_id","formula_pretty","energy_above_hull",
                "band_gap","density","nsites","elements"]
    )
    print(f"  결과: {len(r2)}개")
    all_results.extend(r2)

    # ── 검색 3: Li-S-Cl (단순계 — 아르지로다이트 유사 구조) ──────────────────
    print(f"\n[Step 3] Li-S-Cl 단순계...")
    r3 = mpr.materials.summary.search(
        chemsys="Li-S-Cl",
        energy_above_hull=(0, 0.10),
        fields=["material_id","formula_pretty","energy_above_hull",
                "band_gap","density","nsites","elements"]
    )
    print(f"  결과: {len(r3)}개")
    all_results.extend(r3)

    # ── 검색 4: Li-P-S (베이스 황화물) ───────────────────────────────────────
    print(f"\n[Step 4] Li-P-S 베이스계...")
    r4 = mpr.materials.summary.search(
        chemsys="Li-P-S",
        energy_above_hull=(0, 0.05),
        fields=["material_id","formula_pretty","energy_above_hull",
                "band_gap","density","nsites","elements"]
    )
    print(f"  결과: {len(r4)}개")
    all_results.extend(r4)

# ── 중복 제거 ─────────────────────────────────────────────────────────────────
seen = set()
unique = []
for r in all_results:
    if r.material_id not in seen:
        seen.add(r.material_id)
        unique.append(r)

print(f"\n[Step 5] 중복 제거 후: {len(unique)}개")

# ── Gate 필터 ─────────────────────────────────────────────────────────────────
print(f"\n[Gate 필터 적용]")
passed = []
failed = {"ehull": 0, "bandgap": 0, "density": 0}

for r in unique:
    ehull   = r.energy_above_hull or 999
    bgap    = r.band_gap or 0
    density = r.density or 999

    if ehull > GATE["ehull_max"]:
        failed["ehull"] += 1
        continue
    if bgap < GATE["bandgap_min"] or bgap > GATE["bandgap_max"]:
        failed["bandgap"] += 1
        continue
    if density > GATE["density_max"]:
        failed["density"] += 1
        continue

    passed.append({
        "material_id": r.material_id,
        "formula":     r.formula_pretty,
        "ehull":       round(ehull, 4),
        "band_gap":    round(bgap, 3),
        "density":     round(density, 3),
        "nsites":      r.nsites,
        "elements":    str([str(e) for e in r.elements]),
    })

# ── 결과 출력 ─────────────────────────────────────────────────────────────────
print(f"\n[Gen2 결과]")
print(f"  전체:        {len(unique)}개")
print(f"  Ehull 탈락:  {failed['ehull']}개")
print(f"  밴드갭 탈락: {failed['bandgap']}개")
print(f"  밀도 탈락:   {failed['density']}개")
print(f"  ✅ 통과:     {len(passed)}개")

if passed:
    passed.sort(key=lambda x: x["ehull"])

    print(f"\n[통과 후보 상위 20개]")
    print(f"{'#':>3} {'Formula':<28} {'Ehull':>7} {'BandGap':>8} {'Density':>8}")
    print("-" * 60)
    for i, p in enumerate(passed[:20], 1):
        print(f"{i:>3} {p['formula']:<28} {p['ehull']:>7.4f} {p['band_gap']:>8.3f} {p['density']:>8.3f}")

    with open("gen2_candidates.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=passed[0].keys())
        writer.writeheader()
        writer.writerows(passed)
    print(f"\n✅ gen2_candidates.csv 저장 완료 ({len(passed)}개)")

# ── 병목 분석 ─────────────────────────────────────────────────────────────────
print(f"\n[Gen1 → Gen2 병목 비교]")
print(f"  Gen1 통과: 1개 (검색 3개 중)")
print(f"  Gen2 통과: {len(passed)}개 (검색 {len(unique)}개 중)")

if len(passed) < 5:
    print(f"\n  → 병목: 여전히 데이터 부족")
    print(f"     Gen3 방향: Ehull 추가 완화 + 도핑 원소 확장")
elif len(passed) >= 5:
    print(f"\n  → 충분한 후보 확보")
    print(f"     Gen3 방향: CHGNet MD로 이온전도도 예측")

print(f"\n[다음 단계: Gen3]")
print(f"  통과 후보 상위 10개 → CHGNet 이온전도도 예측")
print("=" * 60)
