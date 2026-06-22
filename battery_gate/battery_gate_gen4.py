"""
Battery-Gate Gen4
─────────────────
Gen3 병목 분석:
  1. Li6PS5Cl/Br material_id 오류 → 재검색
  2. MD 200steps + 500K → Li+ 이동 불충분
  
Gen4 해결:
  1. chemsys로 Li6PS5Cl/Br 정확한 ID 재확인
  2. 온도 700K + 500steps로 상향
  3. Ge 포함 구조 명시적 제외
"""

import os
import csv
import numpy as np
from mp_api.client import MPRester
from chgnet.model import CHGNet
from chgnet.model.dynamics import MolecularDynamics

API_KEY = os.environ.get("MP_API_KEY")
if not API_KEY:
    raise ValueError("MP_API_KEY 환경변수가 설정되지 않았습니다.")

GATE = {
    "conductivity_min":    0.5,   # mS/cm
    "conductivity_target": 3.0,   # mS/cm
}

MD_TEMP  = 700   # K (Gen3 500K → Gen4 700K)
MD_STEPS = 500   # steps (Gen3 200 → Gen4 500)

print("=" * 60)
print("Battery-Gate Gen4 — ID 재확인 + MD 조건 개선")
print("=" * 60)
print(f"\n[Gen3 → Gen4 변경사항]")
print(f"  MD 온도:  500K → {MD_TEMP}K")
print(f"  MD steps: 200  → {MD_STEPS}")
print(f"  Ge 함유 구조 명시적 제외")

# ── Step 1: 정확한 material_id 재확인 ────────────────────────────────────────
print(f"\n[Step 1] Li6PS5Cl / Li6PS5Br 정확한 ID 재검색...")

targets = []

with MPRester(API_KEY) as mpr:

    # Li6PS5Cl 검색
    r_cl = mpr.materials.summary.search(
        chemsys="Li-P-S-Cl",
        energy_above_hull=(0, 0.15),
        fields=["material_id","formula_pretty","energy_above_hull",
                "band_gap","density","nsites","elements"]
    )
    for r in r_cl:
        elems = [str(e) for e in r.elements]
        if "Ge" in elems or "Ge" in r.formula_pretty:
            continue
        targets.append({
            "material_id": r.material_id,
            "formula": r.formula_pretty,
            "ehull": round(r.energy_above_hull, 4),
            "elements": elems,
        })
        print(f"  Li-P-S-Cl: {r.formula_pretty} ({r.material_id}) Ehull={r.energy_above_hull:.4f}")

    # Li6PS5Br 검색
    r_br = mpr.materials.summary.search(
        chemsys="Li-P-S-Br",
        energy_above_hull=(0, 0.15),
        fields=["material_id","formula_pretty","energy_above_hull",
                "band_gap","density","nsites","elements"]
    )
    for r in r_br:
        elems = [str(e) for e in r.elements]
        if "Ge" in elems:
            continue
        targets.append({
            "material_id": r.material_id,
            "formula": r.formula_pretty,
            "ehull": round(r.energy_above_hull, 4),
            "elements": elems,
        })
        print(f"  Li-P-S-Br: {r.formula_pretty} ({r.material_id}) Ehull={r.energy_above_hull:.4f}")

    # Gen3 검증된 Li3PS4 (Ge 없는 버전) 추가
    targets.append({
        "material_id": "mp-985583",
        "formula": "Li3PS4 (Ge-free)",
        "ehull": 0.0074,
        "elements": ["Li","P","S"],
    })

print(f"\n  총 타겟: {len(targets)}개 (Ge 제외)")

# ── Step 2: CHGNet MD 이온전도도 예측 ────────────────────────────────────────
print(f"\n[Step 2] CHGNet MD 실행 (T={MD_TEMP}K, steps={MD_STEPS})...")
model = CHGNet.load()

results = []

with MPRester(API_KEY) as mpr:
    for i, t in enumerate(targets, 1):
        mid     = t["material_id"]
        formula = t["formula"]
        print(f"\n  [{i}/{len(targets)}] {formula} ({mid})")

        try:
            structure = mpr.get_structure_by_material_id(mid)
            if structure is None:
                print(f"    ❌ 구조 없음")
                continue

            elems_in_struct = [str(s.specie) for s in structure]
            if "Ge" in elems_in_struct:
                print(f"    ⚠️ Ge 함유 → 제외 (Gen0 규칙)")
                continue

            print(f"    구조: {structure.formula} ({len(structure)} sites)")

            md = MolecularDynamics(
                atoms=structure,
                model=model,
                ensemble="nvt",
                temperature=MD_TEMP,
                timestep=2,
                trajectory="traj_gen4.traj",
                logfile="md_gen4.log",
                loginterval=10,
            )

            print(f"    MD 실행 중...")
            md.run(MD_STEPS)

            # MSD 계산
            from ase.io import read
            traj = read("traj_gen4.traj", index=":")

            li_indices = [j for j, site in enumerate(structure)
                         if str(site.specie) == "Li"]

            if len(li_indices) == 0:
                print(f"    ⚠️ Li 없음 — 스킵")
                continue

            positions = np.array([[atoms.positions[j] for j in li_indices]
                                  for atoms in traj])

            msd_values = []
            ref = positions[0]
            for pos in positions[1:]:
                disp = pos - ref
                msd  = np.mean(np.sum(disp**2, axis=-1))
                msd_values.append(msd)

            if len(msd_values) < 20:
                print(f"    ⚠️ 데이터 부족")
                continue

            timesteps = np.arange(1, len(msd_values)+1) * 2
            slope = np.polyfit(timesteps[-100:], msd_values[-100:], 1)[0]
            D     = slope / 6 * 1e-20 / 1e-15

            n_Li  = len(li_indices) / structure.volume * 1e24
            kB    = 1.38e-23
            e_    = 1.6e-19
            sigma = n_Li * e_**2 * D / (kB * MD_TEMP)
            sigma_mS = sigma * 1000

            status = "✅ PASS" if sigma_mS >= GATE["conductivity_target"] else \
                     "⚠️ LOW"  if sigma_mS >= GATE["conductivity_min"]    else \
                     "❌ FAIL"

            print(f"    D = {D:.2e} cm²/s")
            print(f"    σ = {sigma_mS:.4f} mS/cm  {status}")

            results.append({
                "material_id":        mid,
                "formula":            formula,
                "D_cm2s":             f"{D:.3e}",
                "conductivity_mS_cm": round(sigma_mS, 4),
                "status":             status,
                "temp_K":             MD_TEMP,
                "steps":              MD_STEPS,
            })

        except Exception as ex:
            print(f"    ❌ 오류: {ex}")

# ── 결과 정리 ─────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("[Gen4 최종 결과]")
print(f"{'#':>3} {'Formula':<25} {'전도도(mS/cm)':>14} {'판정'}")
print("-" * 58)

results.sort(key=lambda x: x["conductivity_mS_cm"], reverse=True)
for i, r in enumerate(results, 1):
    print(f"{i:>3} {r['formula']:<25} {r['conductivity_mS_cm']:>14} {r['status']}")

if results:
    with open("gen4_results.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)
    print(f"\n✅ gen4_results.csv 저장 완료")

# ── 병목 분석 ─────────────────────────────────────────────────────────────────
passed = [r for r in results if r["conductivity_mS_cm"] >= GATE["conductivity_target"]]
print(f"\n[병목 분석]")
print(f"  Gen3 PASS: 1개 (Ge 함유)")
print(f"  Gen4 PASS: {len(passed)}개 (Ge 제외)")

if len(passed) == 0:
    print(f"\n  → 병목: 순수 황화물계 전도도 부족")
    print(f"     Gen5 방향: Al+O 공도핑 조성 직접 구조 생성")
    print(f"                온도 900K + 1000steps 재검증")
else:
    print(f"\n  → 유망 후보 확보!")
    for r in passed:
        print(f"     ✅ {r['formula']}: {r['conductivity_mS_cm']} mS/cm")
    print(f"     Gen5 방향: 도핑 최적화로 전도도 추가 향상")

print(f"\n[Battery-Gate vs Epsilon-Gate 비교]")
print(f"  Epsilon-Gate Gen4: cyano+F 조합 → hERG 해결")
print(f"  Battery-Gate Gen4: Ge 제외 + MD 조건 개선 → 병목 재정의")
print("=" * 60)
