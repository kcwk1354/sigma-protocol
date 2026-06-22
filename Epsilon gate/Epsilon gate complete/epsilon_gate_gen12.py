"""
Epsilon-Gate Gen12 — 도킹 강화 (Gen11 안정성 유지)
====================================================
Gen11 최고 후보: direct_O_pyrrolidine
  ΔG=-4.071 ✅  RMSD=2.28Å ✅  hERG=0.14 ✅  DILI=0.025 ✅
  → 네 지표 동시 달성, 단 ΔG 여유 있음 (-4.5 목표)

Gen12 전략:
  피롤리딘 코어 유지 (RMSD/독성 보장)
  + 도킹 강화 요소 추가:
    - 방향족 꼬리 (π-π stacking, TYR121/PHE330)
    - CN기 (수소결합, Gen4 검증)
    - F (할로겐 결합, Gen4 검증)
    - N-방향족 치환 (π-cation, TRP84)

핵심 제약:
  RotBonds ≤ 3 (Gen11 2에서 약간 허용)
  hERG < 0.4 유지 (Gen11 0.14 수준)
  RMSD < 2.5Å 목표 유지
"""

from rdkit import Chem
from rdkit.Chem import QED, Descriptors, Crippen, Lipinski
import pandas as pd, warnings
warnings.filterwarnings('ignore')

CNS_GATE = {
    "MW":   (None, 460), "LogP": (1.5, 4.2),
    "TPSA": (None, 95),  "HBD":  (None, 3),
    "HBA":  (None, 7),   "QED":  (0.78, None),
    "RotB": (None, 7),   "Lipinski_violations": (None, 0),
}

def evaluate(smiles, gate=CNS_GATE):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"status":"INVALID","reasons":["파싱실패"],"metrics":{}}
    mw   = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    hbd  = Lipinski.NumHDonors(mol)
    hba  = Lipinski.NumHAcceptors(mol)
    qed  = QED.qed(mol)
    rotb = Descriptors.NumRotatableBonds(mol)
    lip  = sum([mw>500,logp>5,hbd>5,hba>10])
    m = {"MW":round(mw,2),"LogP":round(logp,3),"TPSA":round(tpsa,2),
         "HBD":hbd,"HBA":hba,"QED":round(qed,3),
         "RotB":rotb,"Lipinski_violations":lip}
    reasons=[]
    for k,v in gate.items():
        lo,hi=v; val=m[k]
        if lo is not None and val<lo: reasons.append(f"{k}={round(val,2)}<{lo}")
        if hi is not None and val>hi: reasons.append(f"{k}={round(val,2)}>{hi}")
    return {"status":"PASS" if not reasons else "FAIL",
            "reasons":reasons,"metrics":m}

# Gen11 코어: 아다만탄-O-pyrrolidine(N-CH3), RotBonds=2
# 도킹 강화 방향: 피롤리딘 N 또는 아다만탄 주변에 방향족/CN/F 추가

CANDIDATES = {
    # ── 기준점
    "ref_gen11":
        "OC12CC3(OC4CCCN4C)CC(C1)CC(C2)C3",

    # ── 전략 1: 피롤리딘 N → 방향족 치환 (π-cation, TRP84 타겟)
    # N-페닐 피롤리딘
    "pyr_N_phenyl":
        "OC12CC3(OC4CCN(c5ccccc5)C4)CC(C1)CC(C2)C3",

    # N-시아노페닐 (Gen4 CN 검증 + π-cation)
    "pyr_N_cyanophenyl":
        "OC12CC3(OC4CCN(c5ccc(C#N)cc5)C4)CC(C1)CC(C2)C3",

    # N-플루오로페닐 (hERG 유지 + π-cation)
    "pyr_N_fluorophenyl":
        "OC12CC3(OC4CCN(c5ccc(F)cc5)C4)CC(C1)CC(C2)C3",

    # ── 전략 2: 피롤리딘 3번 탄소에 치환 (RotBonds 유지)
    # 3-CN 피롤리딘 (Gen4 CN 패턴)
    "pyr_3CN":
        "OC12CC3(OC4CCN(C)C4C#N)CC(C1)CC(C2)C3",

    # 3-F 피롤리딘 (할로겐 결합)
    "pyr_3F":
        "OC12CC3(OC4CCN(C)C4F)CC(C1)CC(C2)C3",

    # ── 전략 3: 아다만탄 OH → CN으로 교체 (도킹 포켓 상호작용 변화)
    "adamantyl_CN_pyr":
        "N#CC12CC3(OC4CCCN4C)CC(C1)CC(C2)C3",

    # ── 전략 4: 피롤리딘 + 메틸 추가 (입체장애로 포켓 밀착)
    # 2-메틸 피롤리딘 (키랄 중심)
    "pyr_2methyl":
        "OC12CC3(OC4CCCN4)CC(C1)CC(C2)C3",

    # N-메틸 + 3-OH 피롤리딘 (수소결합 추가)
    "pyr_3OH":
        "OC12CC3(OC4CCN(C)C4O)CC(C1)CC(C2)C3",
}

def run():
    rows=[]
    from rdkit.Chem import rdMolDescriptors
    for mol_id, smi in CANDIDATES.items():
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print(f"  파싱실패: {mol_id}"); continue
        canon = Chem.MolToSmiles(mol)
        r = evaluate(canon)
        rot = Descriptors.NumRotatableBonds(mol)
        ar  = rdMolDescriptors.CalcNumAromaticRings(mol)
        row = {"ID":f"gen12__{mol_id}","SMILES":canon,
               "STATUS":r["status"],
               "REASONS":" ; ".join(r["reasons"]) if r["reasons"] else "—",
               "RotBonds_raw":rot,"ArRings":ar}
        row.update(r["metrics"])
        rows.append(row)
    return pd.DataFrame(rows)

if __name__ == "__main__":
    import sys, os
    from rdkit import RDConfig
    from rdkit.Chem import rdMolDescriptors
    sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
    import sascorer

    df = run()
    n     = len(df)
    npass = (df.STATUS=="PASS").sum()

    print(f"\nGen12: 도킹 강화 (Gen11 안정성 유지)")
    print(f"코어: direct_O_pyrrolidine (RMSD=2.28, hERG=0.14)")
    print(f"목표: ΔG < -4.5 + 세 지표 유지")
    print(f"생성 {n}개 / 통과 {npass}개 / 통과율 {npass/n*100:.0f}%")
    print("="*75)

    passed = df[df.STATUS=="PASS"].copy()
    passed['SA'] = passed['SMILES'].apply(
        lambda s: round(sascorer.calculateScore(Chem.MolFromSmiles(s)),2)
        if Chem.MolFromSmiles(s) else None)

    print(f"\n  {'ID':<32} {'QED':>5} {'RotB':>5} {'SA':>5}  전략")
    print("  "+"-"*65)
    desc = {
        'ref_gen11':         'Gen11 기준점',
        'pyr_N_phenyl':      'N-페닐 π-cation',
        'pyr_N_cyanophenyl': 'N-CN페닐 π-cation+HBA',
        'pyr_N_fluorophenyl':'N-F페닐 π-cation+할로겐',
        'pyr_3CN':           '3-CN Gen4패턴',
        'pyr_3F':            '3-F 할로겐결합',
        'adamantyl_CN_pyr':  '아다만탄CN 교체',
        'pyr_2methyl':       '2-메틸 입체효과',
        'pyr_3OH':           '3-OH 수소결합추가',
    }
    for _,r in passed.sort_values('RotBonds_raw').iterrows():
        sid = r['ID'].replace('gen12__','')
        print(f"  {r['ID']:<32} {r['QED']:>5.3f} {r['RotBonds_raw']:>5} "
              f"{r['SA']:>5.2f}  {desc.get(sid,'')}")

    if (df.STATUS!="PASS").any():
        print(f"\n  [탈락]")
        print(df[df.STATUS!="PASS"][["ID","REASONS"]].to_string(index=False))

    df.to_csv("gen12_library.csv", index=False)
    print(f"\n  저장: gen12_library.csv")
