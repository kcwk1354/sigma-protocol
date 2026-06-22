"""
Epsilon-Gate Gen2 — 1세대 실패 패턴을 반영한 다음 세대 설계
=============================================================
1세대 학습 규칙 적용:
  - chlorophenyl 제거 (LogP 초과 실패)
  - pyridyl 중심 확장 (최적 영역, 질소 위치 변형)
  - amide_ethyl 링커 우선 (가장 안정)
  - 목표 LogP 3.0~3.5
"""
from rdkit import Chem
from rdkit.Chem import QED, Descriptors, Crippen, Lipinski
import pandas as pd, itertools

CNS_GATE = {"name":"CNS","MW":(None,450),"LogP":(2.0,4.0),"TPSA":(None,90),
            "HBD":(None,3),"HBA":(None,7),"QED":(0.5,None),"RotB":(None,8),
            "Lipinski_violations":(None,0)}

def evaluate(smiles, gate):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"status":"INVALID","reasons":["파싱실패"],"metrics":{}}
    mw=Descriptors.MolWt(mol); logp=Crippen.MolLogP(mol); tpsa=Descriptors.TPSA(mol)
    hbd=Lipinski.NumHDonors(mol); hba=Lipinski.NumHAcceptors(mol); qed=QED.qed(mol)
    rotb=Descriptors.NumRotatableBonds(mol)
    lip=sum([mw>500,logp>5,hbd>5,hba>10])
    m={"MW":round(mw,2),"LogP":round(logp,3),"TPSA":round(tpsa,2),"HBD":hbd,
       "HBA":hba,"QED":round(qed,3),"RotB":rotb,"Lipinski_violations":lip}
    chk={"MW":mw,"LogP":logp,"TPSA":tpsa,"HBD":hbd,"HBA":hba,"QED":qed,
         "RotB":rotb,"Lipinski_violations":lip}
    rs=[]
    for k,b in gate.items():
        if k=="name": continue
        lo,hi=b; v=chk[k]
        if lo is not None and v<lo: rs.append(f"{k}={round(v,2)}<{lo}")
        if hi is not None and v>hi: rs.append(f"{k}={round(v,2)}>{hi}")
    return {"status":"PASS" if not rs else "FAIL","reasons":rs,"metrics":m}

SCAFFOLD = "OC1(C2)CC3CC(C1)CC(C2)([*])C3"

# Gen2: 학습된 규칙 반영
# - chlorophenyl 삭제됨
# - pyridyl 변형 3종 추가 (질소 위치)
# - amide_ethyl 우선 + 안정 링커만
LINKERS = {
    "amide_ethyl":  "C(=O)NCC{TAIL}",   # 최안정 (우선)
    "amide_direct": "C(=O)N{TAIL}",
    "ether_ethyl":  "OCC{TAIL}",
}
TAILS = {
    "pyridyl_4":     "C4=CC=NC=C4",        # 1세대 최적
    "pyridyl_3":     "C4=CC=CN=C4",        # 질소 위치 변형 (신규)
    "pyridyl_2":     "C4=NC=CC=C4",        # 질소 위치 변형 (신규)
    "fluorophenyl":  "C4=CC=C(F)C=C4",     # 할로겐은 F만 유지
    "methoxyphenyl": "C4=CC=C(OC)C=C4",    # AP2601 계열
    "pyrimidyl":     "C4=NC=CC=N4",        # 질소 2개 (신규, 극성 강화)
    "aminopyridyl":  "C4=CC=NC(N)=C4",     # 신규 변형
}

def run():
    rows=[]
    for (ln,lt),(tn,ts) in itertools.product(LINKERS.items(),TAILS.items()):
        smi = SCAFFOLD.replace("[*]", lt.replace("{TAIL}",ts))
        mol = Chem.MolFromSmiles(smi)
        canon = Chem.MolToSmiles(mol) if mol else smi
        r = evaluate(canon, CNS_GATE)
        row={"ID":f"{ln}__{tn}","SMILES":canon,"STATUS":r["status"],
             "REASONS":" ; ".join(r["reasons"]) if r["reasons"] else "—"}
        row.update(r["metrics"]); rows.append(row)
    return pd.DataFrame(rows)

if __name__=="__main__":
    df=run()
    n=len(df); npass=(df.STATUS=="PASS").sum()
    print(f"\nGen2: 학습 규칙 반영 설계")
    print(f"생성 {n}개 / 통과 {npass}개 / 통과율 {npass/n*100:.0f}%")
    print("="*100)
    passed=df[df.STATUS=="PASS"].sort_values("QED",ascending=False)
    print("\n[통과 분자 — QED 순]")
    print(passed[["ID","MW","LogP","TPSA","QED"]].to_string(index=False))
    if (df.STATUS!="PASS").any():
        print("\n[탈락]")
        print(df[df.STATUS!="PASS"][["ID","LogP","REASONS"]].to_string(index=False))

    # Gen1 vs Gen2 비교
    g1=pd.read_csv("v2_generated_library.csv")
    g1p=(g1.STATUS=="PASS").mean()*100
    print("\n"+"="*100)
    print(f"Gen1 통과율: {g1p:.0f}%  →  Gen2 통과율: {npass/n*100:.0f}%")
    print(f"Gen1 통과분자 최고 QED: {g1[g1.STATUS=='PASS']['QED'].max():.3f}")
    print(f"Gen2 통과분자 최고 QED: {passed['QED'].max():.3f}")
    df.to_csv("gen2_library.csv",index=False)
