"""
Epsilon-Gate Gen7 — 승인약 패턴 기반 설계 (Gen0 시드)
======================================================
출처: PubChem 공개 데이터 (public domain)
  - 도네페질 (Aricept, 1996 FDA 승인)
  - 리바스티그민 (Exelon, 2000 FDA 승인)
  - 갈란타민 (Razadyne, 2001 FDA 승인)

승인약 분석 결과:
  · 황금 LogP: 1.5~4.0 (평균 2.40)
  · 황금 SA:   2.5~4.7 (평균 3.34) → 현재 4.77보다 훨씬 낮음
  · 공통 모티프: 방향족 고리 + 질소 + 메톡시

도네페질 핵심 결합 모티프:
  · 파이페리딘 → TRP84 π-cation 결합
  · 디메톡시벤질 → PHE330/TYR121 π-π stacking
  · 인다논 코어 → 포켓 깊숙이 삽입

Gen7 전략:
  1. 아다만탄 + 파이페리딘 → 도네페질 결합 패턴 이식
  2. SA Score 낮추기 (목표 ≤ 4.0)
  3. Gen6 검증된 cyano+F 유지
"""

from rdkit import Chem
from rdkit.Chem import QED, Descriptors, Crippen, Lipinski
import pandas as pd, warnings
warnings.filterwarnings('ignore')

CNS_GATE = {
    "MW":   (None, 460), "LogP": (1.5, 4.2),   # 승인약 황금구간 반영
    "TPSA": (None, 90),  "HBD":  (None, 3),
    "HBA":  (None, 7),   "QED":  (0.80, None),  # 승인약 수준으로 완화
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

SCAFFOLD = "OC1(C2)CC3CC(C1)CC(C2)([*])C3"

# Gen7: 승인약 모티프 이식
# 전략별 SMILES 직접 설계
CANDIDATES = {
    # ── Gen6 기준점 유지
    "ref_gen6_4cn2f":
        "OC12CC3(Oc4ccc(C#N)cc4F)CC(C1)CC(C2)C3",

    # ── 전략 1: 파이페리딘 이식 (도네페질 핵심 모티프)
    # 아다만탄 + O-CH2-piperidine (TRP84 π-cation 결합 유도)
    "adamantyl_piperidine":
        "OC12CC3(OCC4CCNCC4)CC(C1)CC(C2)C3",

    # 아다만탄 + N-piperidine (직결)
    "adamantyl_N_piperidine":
        "OC12CC3(N4CCCCC4)CC(C1)CC(C2)C3",

    # ── 전략 2: 파이페리딘 + cyano (Gen6 성과 + 도네페질 모티프)
    "piperidine_cyano":
        "OC12CC3(OCC4CCN(Cc5ccc(C#N)cc5)CC4)CC(C1)CC(C2)C3",

    # ── 전략 3: 디메톡시벤질 (도네페질 π-π stacking 모티프)
    "dimethoxybenzyl":
        "OC12CC3(OCc4cc(OC)c(OC)cc4)CC(C1)CC(C2)C3",

    # 디메톡시벤질 + F (Gen4 hERG 해결책 결합)
    "dimethoxybenzyl_F":
        "OC12CC3(OCc4cc(OC)c(OC)c(F)c4)CC(C1)CC(C2)C3",

    # ── 전략 4: SA Score 낮추기 (단순 구조)
    # 메톡시페닐 (도네페질 메톡시 모티프, 단순 버전)
    "methoxyphenyl_simple":
        "OC12CC3(Oc4ccc(OC)cc4)CC(C1)CC(C2)C3",

    # 메톡시 + F
    "methoxyphenyl_F":
        "OC12CC3(Oc4ccc(OC)c(F)c4)CC(C1)CC(C2)C3",

    # ── 전략 5: 갈란타민 모티프 (OH + 메톡시 조합)
    "galantamine_motif":
        "OC12CC3(Oc4ccc(O)c(OC)c4)CC(C1)CC(C2)C3",
}

def run():
    rows = []
    for mol_id, smi in CANDIDATES.items():
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print(f"  파싱실패: {mol_id}")
            continue
        canon = Chem.MolToSmiles(mol)
        r = evaluate(canon)
        row = {"ID": f"gen7__{mol_id}", "SMILES": canon,
               "STATUS": r["status"],
               "REASONS": " ; ".join(r["reasons"]) if r["reasons"] else "—",
               "전략": mol_id.split('_')[0] if '_' in mol_id else mol_id}
        row.update(r["metrics"])
        rows.append(row)
    return pd.DataFrame(rows)

if __name__ == "__main__":
    import sys, os
    from rdkit import RDConfig
    sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
    import sascorer

    df = run()
    n     = len(df)
    npass = (df.STATUS == "PASS").sum()

    print(f"\nGen7: 승인약 패턴 기반 설계")
    print(f"출처: PubChem 공개 데이터 (도네페질/갈란타민/리바스티그민)")
    print(f"생성 {n}개 / 통과 {npass}개 / 통과율 {npass/n*100:.0f}%")
    print("="*80)

    # SA Score 추가
    passed = df[df.STATUS == "PASS"].copy()
    sa_scores = []
    for _, row in passed.iterrows():
        mol = Chem.MolFromSmiles(row['SMILES'])
        sa = round(sascorer.calculateScore(mol), 2) if mol else None
        sa_scores.append(sa)
    passed['SA_Score'] = sa_scores

    print("\n[통과 분자 — SA Score 순 (낮을수록 합성 쉬움)]")
    print(passed[["ID","MW","LogP","QED","SA_Score"]].sort_values('SA_Score').to_string(index=False))

    print(f"\n[승인약 기준점 비교]")
    print(f"  도네페질:  SA=2.52  LogP=4.01  QED=0.805")
    print(f"  Gen6 최고: SA=4.77  LogP=3.55  QED=0.909")

    if (df.STATUS != "PASS").any():
        print("\n[탈락]")
        print(df[df.STATUS!="PASS"][["ID","REASONS"]].to_string(index=False))

    df.to_csv("gen7_library.csv", index=False)
    print(f"\n저장: gen7_library.csv")
