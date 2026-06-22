"""
Epsilon-Gate Gen9 — 하이브리드 코어 + F/CN 위치 최적화
=========================================================
Gen8 학습:
  Nmethyl_pip_dimethoxy: hERG=0.29, DILI=0.015, ΔG=-4.429, RMSD=2.95Å
  → 코어(N-메틸파이페리딘) 유지, 꼬리만 교체

Gen6 교훈 재적용:
  4-CN + 2-F 위치 조합 → RMSD 2.26Å (최고 안정)
  논문 [PMC6152672]: 3번 메타 전자흡인기 = AChE 결합 최적
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6152672/

Gen9 핵심 가설:
  아다만탄-O-CH2-N(CH3)-piperidine 코어
  + F/CN 최적 위치 꼬리
  → 세 지표(hERG / 도킹 / RMSD) 동시 최적화
"""

from rdkit import Chem
from rdkit.Chem import QED, Descriptors, Crippen, Lipinski
import pandas as pd, warnings
warnings.filterwarnings('ignore')

CNS_GATE = {
    "MW":   (None, 480), "LogP": (1.5, 4.5),
    "TPSA": (None, 95),  "HBD":  (None, 3),
    "HBA":  (None, 8),   "QED":  (0.78, None),
    "RotB": (None, 9),   "Lipinski_violations": (None, 0),
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

# Gen9: N-메틸파이페리딘 코어 고정 + 꼬리 최적화
# 코어: 아다만탄-OH + O-CH2-piperidine(N-CH3)
# 꼬리: 파이페리딘 4번 위치에 다양한 치환기

CANDIDATES = {
    # ── Gen8 기준점
    "ref_Nmethyl_pip":
        "OC12CC3(OCC4CCN(C)CC4)CC(C1)CC(C2)C3",

    # ── Gen6 검증 조합: 4-CN + 2-F → 파이페리딘 4번 탄소에 적용
    "pip_4c_4cn2f":
        "OC12CC3(OCC4CCN(C)C(Cc5ccc(C#N)cc5F)C4)CC(C1)CC(C2)C3",

    # ── 파이페리딘 4C + 디메톡시벤질 + F (Gen7 최강 도킹 통합)
    "pip_4c_dimethoxy_F":
        "OC12CC3(OCC4CCN(C)C(Cc5cc(OC)c(OC)c(F)c5)C4)CC(C1)CC(C2)C3",

    # ── 파이페리딘 4C + 메톡시 (단순화, SA Score↓)
    "pip_4c_methoxy":
        "OC12CC3(OCC4CCN(C)C(Cc5ccc(OC)cc5)C4)CC(C1)CC(C2)C3",

    # ── 파이페리딘 4C + CN 직결 (Gen4 핵심 전자흡인기)
    "pip_4c_cn":
        "OC12CC3(OCC4CCN(C)C(c5ccc(C#N)cc5)C4)CC(C1)CC(C2)C3",

    # ── 파이페리딘 N-시아노페닐 (N-치환은 가능)
    "pip_N_cyanophenyl":
        "OC12CC3(OCC4CCN(c5ccc(C#N)cc5)CC4)CC(C1)CC(C2)C3",

    # ── 4,4-difluoro piperidine (고리 자체 F 도입)
    "difluoro_pip_core":
        "OC12CC3(OCC4CC(F)(F)CN(C)C4)CC(C1)CC(C2)C3",

    # ── N-메틸파이페리딘 단독 (기준점)
    "pip_core_only":
        "OC12CC3(OCC4CCN(C)CC4)CC(C1)CC(C2)C3",
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
        row = {"ID": f"gen9__{mol_id}", "SMILES": canon,
               "STATUS": r["status"],
               "REASONS": " ; ".join(r["reasons"]) if r["reasons"] else "—"}
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

    print(f"\nGen9: N-메틸파이페리딘 코어 + F/CN 위치 최적화")
    print(f"코어: Gen8 Nmethyl_pip (hERG=0.29, DILI=0.015) 유지")
    print(f"꼬리: Gen6 F/CN 최적 위치 적용")
    print(f"생성 {n}개 / 통과 {npass}개 / 통과율 {npass/n*100:.0f}%")
    print("="*80)

    passed = df[df.STATUS == "PASS"].copy()
    sa_list = []
    for _, row in passed.iterrows():
        mol = Chem.MolFromSmiles(row['SMILES'])
        sa_list.append(round(sascorer.calculateScore(mol), 2) if mol else None)
    passed['SA_Score'] = sa_list

    print("\n[통과 분자 — QED 순]")
    print(passed[["ID","MW","LogP","QED","SA_Score","HBA"]].sort_values('QED', ascending=False).to_string(index=False))

    if (df.STATUS != "PASS").any():
        print("\n[탈락]")
        print(df[df.STATUS!="PASS"][["ID","REASONS"]].to_string(index=False))

    print(f"\n[Gen8 기준점]")
    print(f"  Nmethyl_pip_dimethoxy: ΔG=-4.429  RMSD=2.95Å  hERG=0.29  DILI=0.015")

    df.to_csv("gen9_library.csv", index=False)
    print(f"\n저장: gen9_library.csv")
