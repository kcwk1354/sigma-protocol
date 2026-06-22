"""
Battery-Gate Gen8
─────────────────
Gen7 발견: 짝수/홀수 비선형 패턴
논문 확인: 4d 사이트 S/Cl 무질서도가 핵심
  - x=41% → 가장 안정적 구조
  - x=50% → 최고 전도도 (논문: 4.6 mS/cm)

Gen8 목표:
  슈퍼셀(2x1x1) 구조로 0.5개 단위 치환 가능하게
  Cl 1~2개 사이 구간 촘촘하게 탐색
  4d 사이트 점유율 25/37.5/50/62.5/75% 탐색
  
핵심 질문:
  우리 파이프라인에서 x=41~50% 최적점이 재현되는가?
  상온 환산 5 mS/cm 이상 달성 가능한가?
"""

import os
import csv
import numpy as np
from mp_api.client import MPRester
from chgnet.model import CHGNet
from chgnet.model.dynamics import MolecularDynamics
from pymatgen.core import Structure

API_KEY = os.environ.get("MP_API_KEY")
if not API_KEY:
    raise ValueError("MP_API_KEY 환경변수가 설정되지 않았습니다.")

MD_TEMP  = 900
MD_STEPS = 1000

print("=" * 60)
print("Battery-Gate Gen8 — 4d 사이트 최적 점유율 탐색")
print("=" * 60)
print(f"\n[논문 기반 타겟]")
print(f"  4d 사이트 S/Cl 무질서도 x=41~50% → 최고 전도도")
print(f"  논문 최고: x=50% → 4.6 mS/cm (상온)")
print(f"  우리 목표: 상온 환산 5 mS/cm 이상")
print(f"\n[Gen8 전략]")
print(f"  2x1x1 슈퍼셀 → 26 sites (유닛셀 13 x 2)")
print(f"  S 사이트 10개 → 0.5개 단위 치환 가능")
print(f"  Cl 치환 비율: 1/2/3/4/5개 (10%, 20%, 30%, 40%, 50%)")

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
        msd  = [np.mean(np.sum((p - pos[0])**2, axis=-1)) for p in pos[1:]]

        if len(msd) < 50:
            return None, "데이터 부족"

        ts    = np.arange(1, len(msd)+1) * 2
        slope = np.polyfit(ts[-200:], msd[-200:], 1)[0]
        D     = slope / 6 * 1e-20 / 1e-15
        n_Li  = len(li_idx) / structure.volume * 1e24
        sigma = n_Li * (1.6e-19)**2 * D / (1.38e-23 * temp) * 1000

        # 상온 환산 (Arrhenius, Ea~0.2eV 가정)
        Ea    = 0.2 * 1.6e-19  # J
        kB    = 1.38e-23
        T_rt  = 298  # K
        ratio = np.exp(-Ea/kB * (1/T_rt - 1/temp))
        sigma_rt = sigma * ratio

        return round(sigma, 4), round(sigma_rt, 4), None

    except Exception as e:
        return None, None, str(e)[:60]

results = []

with MPRester(API_KEY) as mpr:
    # 유닛셀 로드
    base_unit = mpr.get_structure_by_material_id("mp-985592")
    print(f"\n베이스 유닛셀: {base_unit.formula} ({len(base_unit)} sites)")

    # 2x1x1 슈퍼셀 생성
    base = base_unit.make_supercell([2, 1, 1])
    print(f"슈퍼셀 생성: {base.formula} ({len(base)} sites)")

    s_sites = [i for i, s in enumerate(base) if str(s.specie) == "S"]
    print(f"S 사이트: {len(s_sites)}개")

    # 탐색 실험 설계
    # 슈퍼셀에 S가 10개 → n개 치환 시 점유율 n/10
    experiments = [
        {"n": 1,  "label": "x10",  "x_pct": 10,   "desc": "Cl 10% (x=0.1)"},
        {"n": 2,  "label": "x20",  "x_pct": 20,   "desc": "Cl 20% (x=0.2)"},
        {"n": 3,  "label": "x30",  "x_pct": 30,   "desc": "Cl 30% (x=0.3)"},
        {"n": 4,  "label": "x40",  "x_pct": 40,   "desc": "Cl 40% ← 논문 최적점 근처"},
        {"n": 5,  "label": "x50",  "x_pct": 50,   "desc": "Cl 50% ← 논문 최고 전도도"},
        {"n": 6,  "label": "x60",  "x_pct": 60,   "desc": "Cl 60% (과잉 탐색)"},
        {"n": 7,  "label": "x70",  "x_pct": 70,   "desc": "Cl 70% (과잉 탐색)"},
    ]

    for exp in experiments:
        n     = exp["n"]
        label = exp["label"]
        x_pct = exp["x_pct"]
        desc  = exp["desc"]

        print(f"\n[{n}/7] {desc}")

        if n > len(s_sites):
            print(f"  ⚠️ S 사이트 부족 — 스킵")
            continue

        struct = base.copy()
        for k in range(n):
            struct.replace(s_sites[k], "Cl")

        cl_count = sum(1 for s in struct if str(s.specie) == "Cl")
        s_count  = sum(1 for s in struct if str(s.specie) == "S")
        print(f"  구조: {struct.formula} | Cl={cl_count} S={s_count}")

        result = run_md(struct, label)

        if result[0] is not None:
            sigma, sigma_rt, _ = result
            status = "🏆 TOP" if sigma_rt >= 5.0 else \
                     "✅ PASS" if sigma_rt >= 1.0 else "⚠️ LOW"
            print(f"  σ(900K)  = {sigma} mS/cm")
            print(f"  σ(상온)  ≈ {sigma_rt} mS/cm  {status}")

            results.append({
                "x_pct":       x_pct,
                "n_Cl":        n,
                "formula":     struct.formula,
                "sigma_900K":  sigma,
                "sigma_RT_est": sigma_rt,
                "status":      status,
                "desc":        desc,
            })
        else:
            print(f"  ❌ {result[2]}")

# ── 결과 정리 ─────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("[Gen8 결과 — 4d 사이트 점유율별 전도도]")
print(f"{'x%':>5} {'σ(900K)':>10} {'σ(상온)':>10} {'판정'}")
print("-" * 45)

for r in sorted(results, key=lambda x: x["x_pct"]):
    print(f"{r['x_pct']:>5}% {r['sigma_900K']:>10} {r['sigma_RT_est']:>10} {r['status']}")

if results:
    with open("gen8_results.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)
    print(f"\n✅ gen8_results.csv 저장 완료")

# ── 최적점 분석 ───────────────────────────────────────────────────────────────
print(f"\n[최적점 분석]")
if results:
    best = max(results, key=lambda x: x["sigma_RT_est"])
    print(f"  최고 상온 전도도: x={best['x_pct']}% → {best['sigma_RT_est']} mS/cm")
    print(f"  논문 예측 최적점: x=41~50%")

    if 35 <= best["x_pct"] <= 55:
        print(f"  ✅ 논문 예측과 일치 — 파이프라인 신뢰도 확인")
    else:
        print(f"  ⚠️ 논문 예측과 불일치 — 추가 분석 필요")

# ── 전도도 곡선 시각화 ────────────────────────────────────────────────────────
print(f"\n[전도도 곡선]")
print(f"  x%   σ(상온) mS/cm")
print(f"  {'─'*35}")
for r in sorted(results, key=lambda x: x["x_pct"]):
    bar   = "█" * int(r["sigma_RT_est"] * 5)
    arrow = " ← 논문 최적" if 35 <= r["x_pct"] <= 55 else ""
    print(f"  {r['x_pct']:>3}%  {r['sigma_RT_est']:>6.3f}  {bar}{arrow}")

# ── 세대별 추이 ───────────────────────────────────────────────────────────────
print(f"\n[세대별 최고 전도도 추이]")
print(f"  Gen3:  0.059 mS/cm (900K)")
print(f"  Gen4:  1.990 mS/cm (900K)")
print(f"  Gen5: 16.288 mS/cm (900K)")
print(f"  Gen7: 16.288 mS/cm (900K)")
if results:
    best_900 = max(results, key=lambda x: x["sigma_900K"])
    print(f"  Gen8: {best_900['sigma_900K']} mS/cm (900K)")
    print(f"        {best['sigma_RT_est']} mS/cm (상온 환산)")

print(f"\n[다음 단계]")
print(f"  최적점 확정 → Battery-Gate README 작성")
print(f"  → Epsilon-Gate + Battery-Gate 포트폴리오 완성")
print("=" * 60)
