"""
Battery-Gate Gen3
─────────────────
Gen2 통과 후보 → CHGNet MD → Li+ 이온전도도 예측
타겟: Li6PS5Cl, Li6PS5Br, Li7PS6 등 8개 후보
방법: CHGNet MLFF 기반 MD 시뮬레이션 → MSD → 이온전도도 계산
"""

import os
import csv
import json
import numpy as np
from mp_api.client import MPRester
from chgnet.model import CHGNet
from chgnet.model.dynamics import MolecularDynamics
from pymatgen.core import Structure

API_KEY = os.environ.get("MP_API_KEY")
if not API_KEY:
    raise ValueError("MP_API_KEY 환경변수가 설정되지 않았습니다.")

# ── Gate 기준 ─────────────────────────────────────────────────────────────────
GATE = {
    "conductivity_min": 0.1,   # mS/cm 이상 (논문 기준 최소)
    "conductivity_target": 3.0, # mS/cm 이상 (목표)
}

# ── Gen2 통과 후보 material_id 목록 ──────────────────────────────────────────
# gen2_candidates.csv에서 상위 후보 직접 입력
GEN2_CANDIDATES = [
    {"material_id": "mp-985583",  "formula": "Li3PS4"},
    {"material_id": "mp-696128",  "formula": "Li3PS4"},
    {"material_id": "mp-1104821", "formula": "Li7P3S11"},
    {"material_id": "mp-1153395", "formula": "Li7PS6"},
    {"material_id": "mp-1040163", "formula": "Li6PS5Br"},
    {"material_id": "mp-985591",  "formula": "Li6PS5Cl"},
]

print("=" * 60)
print("Battery-Gate Gen3 — CHGNet MD 이온전도도 예측")
print("=" * 60)

# ── CHGNet 모델 로드 ──────────────────────────────────────────────────────────
print("\n[CHGNet 모델 로드 중...]")
model = CHGNet.load()
print("  ✅ CHGNet 로드 완료")

results = []

with MPRester(API_KEY) as mpr:
    for i, candidate in enumerate(GEN2_CANDIDATES, 1):
        mid = candidate["material_id"]
        formula = candidate["formula"]
        print(f"\n[{i}/{len(GEN2_CANDIDATES)}] {formula} ({mid})")

        try:
            # 구조 가져오기
            structure = mpr.get_structure_by_material_id(mid)
            print(f"  구조 로드: {structure.formula} ({len(structure)} sites)")

            # ── CHGNet MD 시뮬레이션 ──────────────────────────────────────────
            # 빠른 테스트: 500K, 200steps (실제는 1000K, 5000steps 권장)
            md = MolecularDynamics(
                atoms=structure,
                model=model,
                ensemble="nvt",
                temperature=500,        # K (높을수록 이온 이동 빠름)
                timestep=2,             # fs
                trajectory="traj.traj",
                logfile="md.log",
                loginterval=10,
            )

            print(f"  MD 시뮬레이션 실행 중 (200 steps)...")
            md.run(200)

            # ── MSD → 이온전도도 계산 ─────────────────────────────────────────
            # 간이 계산 (실제는 더 긴 시뮬레이션 필요)
            from ase.io import read
            traj = read("traj.traj", index=":")

            # Li 원자 MSD 계산
            li_indices = [j for j, site in enumerate(structure)
                         if str(site.specie) == "Li"]

            if len(li_indices) == 0:
                print(f"  ⚠️ Li 원자 없음 — 스킵")
                continue

            positions = np.array([[atoms.positions[j] for j in li_indices]
                                  for atoms in traj])

            # MSD 계산
            msd_values = []
            ref = positions[0]
            for pos in positions[1:]:
                disp = pos - ref
                msd = np.mean(np.sum(disp**2, axis=-1))
                msd_values.append(msd)

            if len(msd_values) < 10:
                print(f"  ⚠️ 데이터 부족 — 스킵")
                continue

            # 확산계수 D (단위: Å²/fs → cm²/s 변환)
            timesteps = np.arange(1, len(msd_values)+1) * 2  # fs
            slope = np.polyfit(timesteps[-50:], msd_values[-50:], 1)[0]
            D = slope / 6 * 1e-20 / 1e-15  # cm²/s

            # Nernst-Einstein 이온전도도 (간이)
            n_Li = len(li_indices) / structure.volume * 1e24  # /cm³
            T = 500  # K
            kB = 1.38e-23
            e = 1.6e-19
            sigma = n_Li * e**2 * D / (kB * T)  # S/cm
            sigma_mS = sigma * 1000  # mS/cm

            # Gate 판정
            passed = sigma_mS >= GATE["conductivity_min"]
            status = "✅ PASS" if sigma_mS >= GATE["conductivity_target"] else \
                     "⚠️ LOW" if passed else "❌ FAIL"

            print(f"  확산계수 D: {D:.2e} cm²/s")
            print(f"  이온전도도: {sigma_mS:.4f} mS/cm  {status}")

            results.append({
                "material_id": mid,
                "formula": formula,
                "D_cm2s": f"{D:.3e}",
                "conductivity_mS_cm": round(sigma_mS, 4),
                "status": status,
                "temperature_K": 500,
                "md_steps": 200,
            })

        except Exception as e:
            print(f"  ❌ 오류: {e}")
            results.append({
                "material_id": mid,
                "formula": formula,
                "D_cm2s": "ERROR",
                "conductivity_mS_cm": 0,
                "status": f"ERROR: {str(e)[:50]}",
                "temperature_K": 500,
                "md_steps": 200,
            })

# ── 결과 정리 ─────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("[Gen3 최종 결과]")
print(f"{'#':>3} {'Formula':<20} {'전도도(mS/cm)':>14} {'판정'}")
print("-" * 55)

results.sort(key=lambda x: float(x["conductivity_mS_cm"]) 
             if x["conductivity_mS_cm"] != "ERROR" else -1, reverse=True)

for i, r in enumerate(results, 1):
    print(f"{i:>3} {r['formula']:<20} {str(r['conductivity_mS_cm']):>14} {r['status']}")

# CSV 저장
if results:
    with open("gen3_results.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)
    print(f"\n✅ gen3_results.csv 저장 완료")

# ── 병목 분석 ─────────────────────────────────────────────────────────────────
passed_target = [r for r in results 
                 if isinstance(r["conductivity_mS_cm"], float) 
                 and r["conductivity_mS_cm"] >= GATE["conductivity_target"]]

print(f"\n[병목 분석]")
print(f"  목표 통과 (≥3.0 mS/cm): {len(passed_target)}개")

if len(passed_target) == 0:
    print(f"  → 병목: 이온전도도 부족")
    print(f"     Gen4 방향: 도핑 원소 추가 / 온도 최적화 / 구조 변형")
else:
    print(f"  → 유망 후보 발견!")
    for r in passed_target:
        print(f"     ✅ {r['formula']}: {r['conductivity_mS_cm']} mS/cm")
    print(f"     Gen4 방향: 해당 조성 도핑 최적화")

print(f"\n[주의사항]")
print(f"  MD 200steps는 빠른 스크리닝용입니다.")
print(f"  최종 검증은 5000steps 이상 권장합니다.")
print(f"\n[다음 단계: Gen4]")
print(f"  병목 해결 → 도핑 최적화 → 최종 후보 선별")
print("=" * 60)
