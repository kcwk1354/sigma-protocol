"""
Epsilon-Gate Gen8 — 하이브리드 설계
======================================
Gen7 두 최고 후보의 강점 결합:

  dimethoxybenzyl_F    → 도킹 최강 (ΔG=-4.697) / π-π stacking
  adamantyl_piperidine → 독성 최강 (hERG=0.20, DILI=0.02) / π-cation

하이브리드 전략:
  파이페리딘(독성 해결) + 디메톡시벤질/F(도킹 강화) 조합
  → 두 특성을 한 분자에서 달성

참고 논문:
  [PMC6152672] 아다만탄 AChE 억제제 — peripheral anionic site 결합
  [Grychowska 2022] F 위치 → hERG 조절
  도네페질 구조: 파이페리딘 + 벤질 조합이 AChE 결합 핵심
"""

from rdkit import Chem
from rdkit.Chem import QED, Descriptors, Crippen, Lipinski
import pandas as pd, warnings
warnings.filterwarnings('ignore')

CNS_GATE = {
    "MW":   (None, 480), "LogP": (1.5, 4.5),
    "TPSA": (None, 95),  "HBD":  (None, 3),
    "HBA":  (None, 8),   "QED":  (0.78, None),
    "RotB": (None, 8),   "Lipinski_violations": (None, 0),
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

# Gen8 하이브리드 후보
CANDIDATES = {
    # ── 기준점 (Gen7 두 후보)
    "ref_dimethoxybenzyl_F":
        "OC12CC3(OCc4cc(OC)c(OC)c(F)c4)CC(C1)CC(C2)C3",
    "ref_adamantyl_pip":
        "OC12CC3(OCC4CCNCC4)CC(C1)CC(C2)C3",

    # ── 전략 1: 파이페리딘 N에 디메톡시벤질 직결 (도네페질 패턴)
    # 아다만탄-O-CH2-piperidine(N-CH2-dimethoxybenzyl)
    "pip_N_dimethoxybenzyl":
        "OC12CC3(OCC4CCN(Cc5cc(OC)c(OC)cc5)CC4)CC(C1)CC(C2)C3",

    # ── 전략 2: 파이페리딘 N에 F-벤질 (hERG 유지 + 도킹 강화)
    "pip_N_fluorobenzyl":
        "OC12CC3(OCC4CCN(Cc5ccc(F)cc5)CC4)CC(C1)CC(C2)C3",

    # ── 전략 3: 파이페리딘 N에 CN-F 페닐 (Gen6 검증 조합)
    "pip_N_cn_fluorobenzyl":
        "OC12CC3(OCC4CCN(Cc5ccc(C#N)cc5F)CC4)CC(C1)CC(C2)C3",

    # ── 전략 4: N-메틸파이페리딘 + 디메톡시벤질
    # 도네페질 정확한 패턴 재현
    "Nmethyl_pip_dimethoxy":
        "OC12CC3(OCC4CCN(C)CC4)CC(C1)CC(C2)C3",

    # ── 전략 5: 파이페리딘 + 메톡시페닐 (단순화, SA Score 낮추기)
    "pip_N_methoxybenzyl":
        "OC12CC3(OCC4CCN(Cc5ccc(OC)cc5)CC4)CC(C1)CC(C2)C3",

    # ── 전략 6: 직접 연결 (링커 없이)
    # 아다만탄-N-piperidine-CH2-dimethoxybenzyl
    "direct_pip_dimethoxy":
        "OC12CC3(N4CCC(Cc5cc(OC)c(OC)cc5)CC4)CC(C1)CC(C2)C3",
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
        row = {"ID": f"gen8__{mol_id}", "SMILES": canon,
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

    print(f"\nGen8: 하이브리드 설계")
    print(f"전략: dimethoxybenzyl_F(도킹) × adamantyl_piperidine(독성) 결합")
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

    print(f"\n[Gen7 기준점]")
    print(f"  dimethoxybenzyl_F:    ΔG=-4.697  hERG=0.36  DILI=0.31")
    print(f"  adamantyl_piperidine: ΔG=-4.467  hERG=0.20  DILI=0.02")

    df.to_csv("gen8_library.csv", index=False)
    print(f"\n저장: gen8_library.csv")
