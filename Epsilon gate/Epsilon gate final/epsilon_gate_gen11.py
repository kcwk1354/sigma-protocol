"""
Epsilon-Gate Gen11 — 링커 최단화 설계
======================================
병목 해결 엔진이 알려준 패턴:
  RMSD를 낮추는 건 고리 종류가 아니라 링커 강직성
  ether_direct(짧은 링커) = RMSD 최저
  파이페리딘/모르폴린 = 독성 최저

Gen11 핵심 가설:
  링커를 CH2 하나 줄이거나 완전 직결
  → 유연성↓ → RMSD↓ → 세 지표 동시 달성

비교:
  Gen8:  아다만탄-O-CH2-N(CH3)-piperidine  RMSD=2.95
  Gen10: 아다만탄-O-CH2-morpholine          RMSD=4.12
  Gen11: 아다만탄-O-morpholine (CH2 제거)   RMSD=?
         아다만탄-N-morpholine (O도 제거)    RMSD=?
"""

from rdkit import Chem
from rdkit.Chem import QED, Descriptors, Crippen, Lipinski, rdMolDescriptors
import pandas as pd, warnings
warnings.filterwarnings('ignore')

CNS_GATE = {
    "MW":   (None, 450), "LogP": (1.5, 4.5),
    "TPSA": (None, 95),  "HBD":  (None, 3),
    "HBA":  (None, 8),   "QED":  (0.75, None),
    "RotB": (None, 6),   "Lipinski_violations": (None, 0),  # RotB 기준 강화
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

CANDIDATES = {
    # ── 기준점들
    "ref_gen8_Nmethyl":
        "OC12CC3(OCC4CCN(C)CC4)CC(C1)CC(C2)C3",       # CH2 있음, RMSD=2.95
    "ref_gen10_morpholine":
        "OC12CC3(OCC4CCOCC4)CC(C1)CC(C2)C3",           # CH2 있음, RMSD=4.12

    # ── 전략 1: CH2 제거 (O 직결)
    "direct_O_morpholine":
        "OC12CC3(OC4CCOCC4)CC(C1)CC(C2)C3",            # O-morpholine 직결
    "direct_O_Nmethylpip":
        "OC12CC3(OC4CCN(C)CC4)CC(C1)CC(C2)C3",         # O-N메틸파이페리딘 직결
    "direct_O_pyrrolidine":
        "OC12CC3(OC4CCCN4C)CC(C1)CC(C2)C3",            # O-N메틸피롤리딘 직결

    # ── 전략 2: O도 제거 (N 직결 — 최단)
    "direct_N_morpholine":
        "OC12CC3(N4CCOCC4)CC(C1)CC(C2)C3",             # N-morpholine 직결
    "direct_N_Nmethylpip":
        "OC12CC3(N4CCN(C)CC4)CC(C1)CC(C2)C3",          # N,N'-디메틸파이페라진
    "direct_N_pyrrolidine":
        "OC12CC3(N4CCCC4)CC(C1)CC(C2)C3",              # N-피롤리딘 직결

    # ── 전략 3: 스파이로 (링커 자체를 없앰 — 완전 강직)
    # 아다만탄과 파이페리딘을 탄소 하나로 공유
    "spiro_pip":
        "OC12CC3(C4(CCNCC4)CC(C1)CC2)C3",              # 스파이로[아다만탄-파이페리딘]
    "spiro_morpholine":
        "OC12CC3(C4(CCOCC4)CC(C1)CC2)C3",              # 스파이로[아다만탄-모르폴린]
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
        rot = Descriptors.NumRotatableBonds(mol)
        ar  = rdMolDescriptors.CalcNumAromaticRings(mol)
        row = {"ID": f"gen11__{mol_id}", "SMILES": canon,
               "STATUS": r["status"],
               "REASONS": " ; ".join(r["reasons"]) if r["reasons"] else "—",
               "RotBonds_raw": rot, "ArRings": ar}
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

    print(f"\nGen11: 링커 최단화 설계")
    print(f"가설: CH2 제거 → 유연성↓ → RMSD↓")
    print(f"목표: RMSD < 2.5Å + hERG < 0.4 + ΔG < -4.0 동시 달성")
    print(f"생성 {n}개 / 통과 {npass}개 / 통과율 {npass/n*100:.0f}%")
    print("="*75)

    passed = df[df.STATUS=="PASS"].copy()
    passed['SA'] = passed['SMILES'].apply(
        lambda s: round(sascorer.calculateScore(Chem.MolFromSmiles(s)),2)
        if Chem.MolFromSmiles(s) else None)

    print(f"\n  {'ID':<32} {'QED':>5} {'RotB':>5} {'SA':>5}  링커구조")
    print("  "+"-"*65)
    linker_desc = {
        'ref_gen8_Nmethyl':      'O-CH2-pip (기준)',
        'ref_gen10_morpholine':  'O-CH2-mor (기준)',
        'direct_O_morpholine':   'O-mor (CH2제거)',
        'direct_O_Nmethylpip':   'O-pip (CH2제거)',
        'direct_O_pyrrolidine':  'O-pyr (CH2제거)',
        'direct_N_morpholine':   'N-mor (O+CH2제거)',
        'direct_N_Nmethylpip':   'N-pip (O+CH2제거)',
        'direct_N_pyrrolidine':  'N-pyr (O+CH2제거)',
        'spiro_pip':             '스파이로-pip',
        'spiro_morpholine':      '스파이로-mor',
    }
    for _, r in passed.sort_values('RotBonds_raw').iterrows():
        short_id = r['ID'].replace('gen11__','')
        desc = linker_desc.get(short_id, '')
        print(f"  {r['ID']:<32} {r['QED']:>5.3f} {r['RotBonds_raw']:>5} "
              f"{r['SA']:>5.2f}  {desc}")

    if (df.STATUS!="PASS").any():
        print(f"\n  [탈락]")
        print(df[df.STATUS!="PASS"][["ID","REASONS"]].to_string(index=False))

    df.to_csv("gen11_library.csv", index=False)
    print(f"\n  저장: gen11_library.csv")
    print(f"""
  핵심 비교 포인트:
    RotBonds 낮을수록 → MD 안정성↑ 예상
    기준점(ref) RotBonds=3 → Gen11 목표 RotBonds≤2
""")
