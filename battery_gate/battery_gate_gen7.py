"""
Battery-Gate Gen7
─────────────────
Gen6 발견: 비선형 패턴
  Cl 1개: 16.29 🏆
  Cl 2개:  4.96 (급락)
  Cl 3개: 10.68 🏆 (반등)

Gen7 전략: Wide Scan
  Cl 0~5개 전체 범위 탐색
  - 구조 붕괴 지점 확인
  - 비선형 패턴 전체 그림 완성
  - 최적 구간 결정 → Gen8 촘촘한 탐색

핵심 질문:
  Cl 4개, 5개에서 어떻게 되는가?
  구조가 유지되는가, 무너지는가?
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
    "conductivity_top":    10.0,
}

MD_TEMP  = 900
MD_STEPS = 1000

print("=" * 60)
print("Battery-Gate Gen7 — Wide Scan (Cl 0~5개)")
print("=" * 60)
print(f"\n[Gen6 확인된 패턴]")
print(f"  Cl 0개:  6.56 mS/cm")
print(f"  Cl 1개: 16.29 mS/cm 🏆")
print(f"  Cl 2개:  4.96 mS/cm (급락)")
print(f"  Cl 3개: 10.68 mS/cm 🏆 (반등)")
print(f"  Cl 4개: ??? (미탐색)")
print(f"  Cl 5개: ??? (S 전부 치환 — 구조 붕괴 가능)")
print(f"\n[Gen7 목표]")
print(f"  전체 범위 완성 → 비선형 패턴 규명")
print(f"  구조 붕괴 지점 확인")

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
        return None, str(e)[:80]

# ── Gen6 확인값 + Gen7 신규 탐색 ─────────────────────────────────────────────
experiments = [
    {"n": 0, "label": "cl0", "desc": "Li6PS5Cl (베이스)"},
    {"n": 1, "label": "cl1", "desc": "Li6PS4Cl2 (Gen6 🏆)"},
    {"n": 2, "label": "cl2", "desc": "Li6PS3Cl3 (Gen6 급락)"},
    {"n": 3, "label": "cl3", "desc": "Li6PS2Cl4 (Gen6 반등)"},
    {"n": 4, "label": "cl4", "desc": "Li6PSCl5 (신규 탐색)"},
    {"n": 5, "label": "cl5", "desc": "Li6PCl6 (S 전부 치환 — 구조 붕괴?)"},
]

results = []

with MPRester(API_KEY) as mpr:
    base = mpr.get_structure_by_material_id("mp-985592")
    s_sites = [i for i, s in enumerate(base) if str(s.specie) == "S"]
    print(f"\n베이스: {base.formula} | S 사이트: {len(s_sites)}개")

    for exp in experiments:
        n = exp["n"]
        print(f"\n[{n+1}/6] {exp['desc']}")

        struct = base.copy()

        if n > len(s_sites):
            print(f"  ⚠️ S 사이트 부족 — 스킵")
            continue

        for k in range(n):
            struct.replace(s_sites[k], "Cl")

        formula = struct.formula
        print(f"  구조: {formula}")

        # 구조 붕괴 사전 체크 — S가 없으면 경고
        s_count = sum(1 for s in struct if str(s.specie) == "S")
        if s_count == 0:
            print(f"  ⚠️ S 원자 없음 — 아르지로다이트 구조 붕괴 가능")

        sigma, err = run_md(struct, exp["label"])

        if sigma is not None:
            status = "🏆 TOP"  if sigma >= GATE["conductivity_top"]    else \
                     "✅ PASS" if sigma >= GATE["conductivity_target"]  else \
                     "⚠️ LOW"
            print(f"  σ = {sigma} mS/cm  {status}")
            results.append({
                "n_Cl_added":         n,
                "formula":            formula,
                "S_remaining":        s_count,
                "conductivity_mS_cm": sigma,
                "status":             status,
                "desc":               exp["desc"],
            })
        else:
            print(f"  ❌ {err}")
            results.append({
                "n_Cl_added":         n,
                "formula":            formula,
                "S_remaining":        s_count,
                "conductivity_mS_cm": 0,
                "status":             f"ERROR",
                "desc":               exp["desc"],
            })

# ── 결과 정리 ─────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("[Gen7 Wide Scan 결과 — Cl 함량별 전도도 전체]")
print(f"{'Cl추가':>6} {'Formula':<20} {'S잔여':>6} {'전도도':>12} {'판정'}")
print("-" * 60)

for r in results:
    print(f"{r['n_Cl_added']:>6} {r['formula']:<20} {r['S_remaining']:>6} {r['conductivity_mS_cm']:>12} {r['status']}")

if results:
    with open("gen7_results.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)
    print(f"\n✅ gen7_results.csv 저장 완료")

# ── 패턴 분석 ─────────────────────────────────────────────────────────────────
print(f"\n[Gen7 패턴 분석]")
valid = [(r["n_Cl_added"], r["conductivity_mS_cm"])
         for r in results if r["conductivity_mS_cm"] > 0]

if len(valid) >= 4:
    peak = max(valid, key=lambda x: x[1])
    print(f"  전체 최고: Cl {peak[0]}개 → {peak[1]} mS/cm")

    # 비선형 패턴 확인
    vals = [v[1] for v in sorted(valid)]
    print(f"\n  전도도 추이:")
    for n, sigma in sorted(valid):
        bar = "█" * int(sigma / 2)
        print(f"    Cl {n}개: {sigma:>8.2f} mS/cm  {bar}")

    # 최적 구간 결정
    top_candidates = [(n, s) for n, s in valid if s >= GATE["conductivity_top"]]
    if top_candidates:
        ns = [n for n, s in top_candidates]
        print(f"\n  🏆 TOP 구간: Cl {min(ns)}~{max(ns)}개")
        print(f"     → Gen8: 이 구간 0.25개 단위로 촘촘하게 탐색")

print(f"\n[구조 안정성 분석]")
for r in results:
    if r["S_remaining"] == 0:
        if r["conductivity_mS_cm"] > 0:
            print(f"  S=0 (Li6PCl6): 구조 유지, σ={r['conductivity_mS_cm']} mS/cm")
        else:
            print(f"  S=0 (Li6PCl6): 구조 붕괴 확인 ← Gen0 규칙 재확인")

print(f"\n[세대별 최고 전도도 추이]")
print(f"  Gen3:  0.059 mS/cm")
print(f"  Gen4:  1.990 mS/cm  →  33배↑")
print(f"  Gen5: 16.288 mS/cm  →   8배↑")
print(f"  Gen6: 16.288 mS/cm")
print(f"  Gen7: {max([r['conductivity_mS_cm'] for r in results], default=0):.3f} mS/cm")

print(f"\n[다음 단계: Gen8]")
print(f"  TOP 구간 촘촘한 탐색 → 진짜 최적점 확정")
print(f"  → 논문 미탐색 조성 최종 확인")
print("=" * 60)
