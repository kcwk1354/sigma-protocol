"""
Battery-Gate Gen6
─────────────────
Gen5 핵심 발견:
  - Li6PS4Cl2 (Cl 증량) → 16.29 mS/cm 압도적 1위
  - Al 도핑 → 오히려 전도도 하락 (구조 왜곡)

Gen6 목표:
  Cl 증량 최적점 탐색
  S→Cl 치환 개수: 1개(Gen5) → 0.5/1/1.5/2개 단계적 탐색
  임계점: 과잉 Cl → 이차상 형성 (Gen0 규칙: Cl ≤ 1.5)
  
구조-활성 규칙 확립:
  Epsilon-Gate Gen6: F 위치 하나가 activity cliff
  Battery-Gate Gen6: Cl 함량이 전도도 cliff
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
    "conductivity_target": 3.0,
    "conductivity_top":    10.0,  # Gen5 최고 16.29 → 10 이상을 최우수로
}

MD_TEMP  = 900
MD_STEPS = 1000

print("=" * 60)
print("Battery-Gate Gen6 — Cl 증량 최적점 탐색")
print("=" * 60)
print(f"\n[Gen5 핵심 학습]")
print(f"  ❌ Al 도핑:  1.41 mS/cm (구조 왜곡 > Li 공공 효과)")
print(f"  ✅ Cl 증량: 16.29 mS/cm (이동 경로 확장 압도적)")
print(f"\n[Gen6 전략]")
print(f"  Cl 함량 단계적 변화 → 최적점 + 임계점 동시 탐색")
print(f"  Gen0 규칙 재검증: Cl ≤ 1.5 임계점 실제 확인")

model = CHGNet.load()
print(f"\n✅ CHGNet 로드 완료")

def run_md(structure, label, temp=900, steps=1000):
    try:
        li_idx = [j for j, s in enumerate(structure) if str(s.specie) == "Li"]
        if not li_idx:
            return None, "Li 없음"

        md = MolecularDynamics(
            atoms=structure, model=model, ensemble="nvt",
            temperature=temp, timestep=2,
            trajectory=f"traj_{label}.traj",
            logfile=f"md_{label}.log", loginterval=20,
        )
        md.run(steps)

        from ase.io import read
        traj = read(f"traj_{label}.traj", index=":")
        pos  = np.array([[a.positions[j] for j in li_idx] for a in traj])

        msd = [np.mean(np.sum((p - pos[0])**2, axis=-1)) for p in pos[1:]]
        if len(msd) < 50:
            return None, "데이터 부족"

        ts    = np.arange(1, len(msd)+1) * 2
        slope = np.polyfit(ts[-200:], msd[-200:], 1)[0]
        D     = slope / 6 * 1e-20 / 1e-15

        n_Li  = len(li_idx) / structure.volume * 1e24
        sigma = n_Li * (1.6e-19)**2 * D / (1.38e-23 * temp) * 1000
        return round(sigma, 4), None

    except Exception as e:
        return None, str(e)[:60]

results = []

with MPRester(API_KEY) as mpr:
    base = mpr.get_structure_by_material_id("mp-985592")  # Li6PS5Cl
    print(f"\n베이스 구조: {base.formula} ({len(base)} sites)")

    s_sites = [i for i, s in enumerate(base) if str(s.specie) == "S"]
    print(f"S 사이트 수: {len(s_sites)}개 (치환 가능)")

    # ── 실험 설계 ─────────────────────────────────────────────────────────────
    experiments = [
        {"label": "Cl0",   "n_sub": 0, "desc": "Li6PS5Cl 베이스 (Cl 0개 추가)"},
        {"label": "Cl0.5", "n_sub": 0, "desc": "Li6PS5Cl + Br 혼합 (Cl/Br 비교)"},
        {"label": "Cl1",   "n_sub": 1, "desc": "Li6PS4Cl2 (S→Cl 1개, Gen5 최고)"},
        {"label": "Cl2",   "n_sub": 2, "desc": "Li6PS3Cl3 (S→Cl 2개, 과잉 탐색)"},
        {"label": "Cl3",   "n_sub": 3, "desc": "Li6PS2Cl4 (S→Cl 3개, 임계점)"},
    ]

    for i, exp in enumerate(experiments, 1):
        print(f"\n[{i}/{len(experiments)}] {exp['desc']}")

        struct = base.copy()
        n = exp["n_sub"]

        if n == 0 and exp["label"] == "Cl0.5":
            # Cl/Br 혼합 비교 — S 1개를 Br로
            struct.replace(s_sites[0], "Br")
            formula_label = "Li6PS4ClBr"
        elif n > 0:
            if n > len(s_sites):
                print(f"  ⚠️ S 사이트 부족 ({len(s_sites)}개) — 스킵")
                continue
            for k in range(n):
                struct.replace(s_sites[k], "Cl")
            formula_label = f"Li6PS{5-n}Cl{1+n}"
        else:
            formula_label = "Li6PS5Cl (베이스)"

        print(f"  구조: {struct.formula}")

        sigma, err = run_md(struct, exp["label"])
        if sigma:
            status = "🏆 TOP"  if sigma >= GATE["conductivity_top"]    else \
                     "✅ PASS" if sigma >= GATE["conductivity_target"]  else \
                     "⚠️ LOW"
            print(f"  σ = {sigma} mS/cm  {status}")
            results.append({
                "formula":            formula_label,
                "S_replaced":         n,
                "conductivity_mS_cm": sigma,
                "status":             status,
                "desc":               exp["desc"],
            })
        else:
            print(f"  ❌ {err}")

# ── 결과 정리 ─────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("[Gen6 최종 결과 — Cl 함량별 전도도]")
print(f"{'#':>3} {'Formula':<22} {'S치환':>6} {'전도도':>12} {'판정'}")
print("-" * 58)

results.sort(key=lambda x: x["S_replaced"])
for i, r in enumerate(results, 1):
    print(f"{i:>3} {r['formula']:<22} {r['S_replaced']:>6} {r['conductivity_mS_cm']:>12} {r['status']}")

if results:
    with open("gen6_results.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)
    print(f"\n✅ gen6_results.csv 저장 완료")

# ── 구조-활성 규칙 도출 ───────────────────────────────────────────────────────
print(f"\n[Gen6 구조-활성 규칙 도출]")
if len(results) >= 3:
    vals = [(r["S_replaced"], r["conductivity_mS_cm"]) for r in results]
    vals.sort()
    peak = max(vals, key=lambda x: x[1])
    print(f"  최적 Cl 치환 개수: S→Cl {peak[0]}개")
    print(f"  최고 전도도: {peak[1]} mS/cm")

    # 임계점 확인
    if len(vals) >= 3:
        for j in range(1, len(vals)):
            if vals[j][1] < vals[j-1][1] * 0.7:
                print(f"  ⚠️ 임계점 발견: S→Cl {vals[j][0]}개에서 전도도 급락")
                print(f"     → Gen0 규칙 'Cl 과잉 금지' 실험적 확인")
                break
        else:
            print(f"  → 임계점 미발견: Gen7에서 추가 탐색 필요")

print(f"\n[세대별 최고 전도도 추이]")
print(f"  Gen3:  0.059 mS/cm")
print(f"  Gen4:  1.990 mS/cm  →  33배↑")
print(f"  Gen5: 16.288 mS/cm  →   8배↑")
print(f"  Gen6: {max([r['conductivity_mS_cm'] for r in results], default=0):.3f} mS/cm")

print(f"\n[Epsilon-Gate vs Battery-Gate 패턴 비교]")
print(f"  Epsilon Gen6: F 위치 → activity cliff 발견")
print(f"  Battery Gen6: Cl 함량 → 전도도 cliff 탐색")

print(f"\n[다음 단계: Gen7]")
print(f"  최적 Cl 함량 고정 + Br 혼합 비율 최적화")
print(f"  → 논문 미탐색 조성 공간 진입")
print("=" * 60)
