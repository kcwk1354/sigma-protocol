"""
Epsilon-Gate 병목 해결 엔진
==============================
핵심 철학:
  병목 발생 시 → 이전 단계 최선값 유지
              → 병목 단계만 집중 최적화
              → 통과하면 다음 단계 진행

병목 감지 기준:
  물성 통과율 < 70%  → 물성 병목
  독성 통과율 < 40%  → 독성 병목  ← Gen4에서 경험 (20%)
  도킹 통과율 < 50%  → 도킹 병목
  MD 통과율   < 50%  → MD 병목    ← 현재 주요 병목 (25~67%)

현재 분석 결과:
  MD가 주요 병목 — 파이페리딘 유연성이 원인
  독성 최선값(Gen8 Nmethyl_pip) 유지 + MD 강직성 강화

Gen10 전략:
  파이페리딘 → benzopiperidine/isoquinoline 교체
  아다만탄-스파이로-파이페리딘 탐색
"""

import os, sys, warnings, tempfile
warnings.filterwarnings('ignore')
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import QED, Descriptors, Crippen, Lipinski, rdMolDescriptors

# ──────────────────────────────────────────
# 1. 병목 감지기
# ──────────────────────────────────────────
def detect_bottleneck(stage_rates):
    """
    각 단계 통과율 → 병목 단계 감지
    stage_rates: {'물성': 0.89, '독성': 0.43, '도킹': 1.0, 'MD': 0.33}
    """
    thresholds = {'물성': 0.70, '독성': 0.40, '도킹': 0.50, 'MD': 0.50}
    bottlenecks = []
    for stage, rate in stage_rates.items():
        if rate < thresholds.get(stage, 0.5):
            bottlenecks.append((stage, rate))
    # 가장 심각한 병목 반환
    if bottlenecks:
        return min(bottlenecks, key=lambda x: x[1])
    return None

# ──────────────────────────────────────────
# 2. 구조적 병목 원인 분석
# ──────────────────────────────────────────
def analyze_md_bottleneck(stable_smiles_list, unstable_smiles_list):
    """
    안정/불안정 분자 구조 비교 → 원인 패턴 추출
    """
    def features(smi):
        mol = Chem.MolFromSmiles(smi)
        if not mol: return None
        return {
            'RotBonds': Descriptors.NumRotatableBonds(mol),
            'ArRings':  rdMolDescriptors.CalcNumAromaticRings(mol),
            'MW':       Descriptors.MolWt(mol),
            'HasPip':   int(mol.HasSubstructMatch(Chem.MolFromSmarts('[NX3;R;!a]'))),
            'HasCN':    int(mol.HasSubstructMatch(Chem.MolFromSmarts('C#N'))),
            'HasF':     int(mol.HasSubstructMatch(Chem.MolFromSmarts('[F]'))),
        }

    sf = [f for s in stable_smiles_list   if (f := features(s))]
    uf = [f for s in unstable_smiles_list if (f := features(s))]

    if not sf or not uf: return {}

    causes = {}
    for k in sf[0]:
        sv = sum(f[k] for f in sf)/len(sf)
        uv = sum(f[k] for f in uf)/len(uf)
        if abs(uv - sv) > 0.2:
            causes[k] = {'stable': round(sv,2), 'unstable': round(uv,2),
                         'diff': round(uv-sv,2)}
    return causes

# ──────────────────────────────────────────
# 3. Gen10 분자 설계 — MD 병목 집중 해결
# ──────────────────────────────────────────

CNS_GATE = {
    "MW":   (None, 480), "LogP": (1.5, 4.5),
    "TPSA": (None, 95),  "HBD":  (None, 3),
    "HBA":  (None, 8),   "QED":  (0.75, None),
    "RotB": (None, 7),   "Lipinski_violations": (None, 0),
}

def evaluate(smiles, gate=CNS_GATE):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"status": "INVALID", "reasons": ["파싱실패"], "metrics": {}}
    mw   = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    hbd  = Lipinski.NumHDonors(mol)
    hba  = Lipinski.NumHAcceptors(mol)
    qed  = QED.qed(mol)
    rotb = Descriptors.NumRotatableBonds(mol)
    lip  = sum([mw>500, logp>5, hbd>5, hba>10])
    m = {"MW":round(mw,2),"LogP":round(logp,3),"TPSA":round(tpsa,2),
         "HBD":hbd,"HBA":hba,"QED":round(qed,3),
         "RotB":rotb,"Lipinski_violations":lip}
    reasons = []
    for k, v in gate.items():
        lo, hi = v; val = m[k]
        if lo is not None and val < lo: reasons.append(f"{k}={round(val,2)}<{lo}")
        if hi is not None and val > hi: reasons.append(f"{k}={round(val,2)}>{hi}")
    return {"status": "PASS" if not reasons else "FAIL",
            "reasons": reasons, "metrics": m}

# Gen10 후보: 파이페리딘 유연성 → 강직성 구조로 교체
# 독성 최선값(Gen8) 핵심 유지: 아다만탄-O-CH2-N(CH3)-[고리]
CANDIDATES_GEN10 = {
    # ── Gen8 기준점 (병목 이전 최선값)
    "ref_gen8_Nmethyl":
        "OC12CC3(OCC4CCN(C)CC4)CC(C1)CC(C2)C3",

    # ── 전략 1: benzopiperidine (파이페리딘 + 방향족 잠금)
    # 테트라히드로이소퀴놀린 (5+6 융합, 강직성↑)
    "tetrahydroisoquinoline":
        "OC12CC3(OCC4CNc5ccccc54)CC(C1)CC(C2)C3",

    # N-메틸 테트라히드로이소퀴놀린
    "Nmethyl_THIQ":
        "OC12CC3(OCC4CN(C)c5ccccc54)CC(C1)CC(C2)C3",

    # ── 전략 2: 스파이로 연결 (아다만탄-스파이로-파이페리딘)
    # 스파이로[아다만탄-파이페리딘] — 강직성 극대화
    "spiro_adamantyl_pip":
        "OC12CC3(CC4(CCNCC4)CC(C1)CC2)C3",

    # ── 전략 3: 피롤리딘 (5원환, 파이페리딘보다 강직)
    "pyrrolidine":
        "OC12CC3(OCC4CCCN4C)CC(C1)CC(C2)C3",

    # ── 전략 4: 모르폴린 (O+N 고리, 극성↑ hERG↓)
    "morpholine":
        "OC12CC3(OCC4CCOCC4)CC(C1)CC(C2)C3",

    # N-메틸 모르폴린
    "Nmethyl_morpholine":
        "OC12CC3(OCC4CCN(C)OC4)CC(C1)CC(C2)C3",

    # ── 전략 5: 이소퀴놀린 (완전 방향족, 강직성 최대)
    "isoquinoline":
        "OC12CC3(OCc4nccc5ccccc45)CC(C1)CC(C2)C3",

    # ── 전략 6: Gen8 + 방향족 고리 추가 잠금
    # N-메틸파이페리딘 + 페닐 (4번 위치 강직)
    "Nmethyl_pip_4phenyl":
        "OC12CC3(OCC4CCN(C)C(c5ccccc5)C4)CC(C1)CC(C2)C3",
}

def run_gen10():
    rows = []
    for mol_id, smi in CANDIDATES_GEN10.items():
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print(f"  파싱실패: {mol_id}")
            continue
        canon = Chem.MolToSmiles(mol)
        r = evaluate(canon)

        # RotBonds 체크 (MD 병목 원인)
        rot = Descriptors.NumRotatableBonds(mol)
        ar  = rdMolDescriptors.CalcNumAromaticRings(mol)

        row = {"ID": f"gen10__{mol_id}", "SMILES": canon,
               "STATUS": r["status"],
               "REASONS": " ; ".join(r["reasons"]) if r["reasons"] else "—",
               "RotBonds_raw": rot, "ArRings": ar}
        row.update(r["metrics"])
        rows.append(row)
    return pd.DataFrame(rows)

if __name__ == "__main__":
    print("\n" + "="*70)
    print("  Epsilon-Gate 병목 해결 엔진")
    print("="*70)

    # 병목 감지
    stage_rates = {'물성': 0.85, '독성': 0.43, '도킹': 1.0, 'MD': 0.35}
    bottleneck = detect_bottleneck(stage_rates)
    print(f"\n  🔍 병목 감지: {bottleneck[0]} (통과율 {bottleneck[1]*100:.0f}%)")
    print(f"  → 병목 이전 최선값 유지, 해당 단계 집중 최적화")

    print(f"\n  전략: 파이페리딘 유연성 → 강직성 구조 교체")
    print(f"  목표: RMSD < 2.5Å (현재 최고 2.70Å)")
    print(f"  독성 기준 유지: hERG < 0.4, DILI < 0.4")

    df = run_gen10()
    n     = len(df)
    npass = (df.STATUS == "PASS").sum()

    print(f"\n  생성 {n}개 / 통과 {npass}개 / 통과율 {npass/n*100:.0f}%")
    print("="*70)

    passed = df[df.STATUS=="PASS"].copy()

    sys.path.append(os.path.join(__import__('rdkit').RDConfig.RDContribDir, 'SA_Score'))
    import sascorer
    passed['SA'] = passed['SMILES'].apply(
        lambda s: round(sascorer.calculateScore(Chem.MolFromSmiles(s)),2)
        if Chem.MolFromSmiles(s) else None)

    print(f"\n  {'ID':<35} {'QED':>5} {'RotB':>5} {'ArR':>4} {'SA':>5}")
    print("  " + "-"*58)
    for _, r in passed.sort_values('RotBonds_raw').iterrows():
        print(f"  {r['ID']:<35} {r['QED']:>5.3f} {r['RotBonds_raw']:>5} "
              f"{r['ArRings']:>4} {r['SA']:>5.2f}")

    if (df.STATUS!="PASS").any():
        print(f"\n  [탈락]")
        print(df[df.STATUS!="PASS"][["ID","REASONS"]].to_string(index=False))

    df.to_csv("gen10_library.csv", index=False)
    print(f"\n  저장: gen10_library.csv")
    print(f"""
  📌 병목 해결 로직:
     Gen8 최선값(독성) → 유지
     파이페리딘 → 강직성 구조 교체 → MD 안정성 목표
     세 지표 동시 달성 가능성 가장 높은 Gen
""")
