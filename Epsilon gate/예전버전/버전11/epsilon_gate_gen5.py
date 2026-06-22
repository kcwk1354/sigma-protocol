"""
Epsilon-Gate Gen5 — MD 안정성 확보 집중 설계
==============================================
Gen4 학습 규칙:
  - cyano+F 조합 유지 (hERG+도킹 동시 해결 검증)
  - MD 안정성 개선: 고리 강직성 강화
  - ether_direct 링커 유지 (hERG 최적)
  - 목표: RMSD < 2.5Å (Gen4 2.66Å → 개선)
"""

from rdkit import Chem
from rdkit.Chem import QED, Descriptors, Crippen, Lipinski
import pandas as pd, itertools, warnings
warnings.filterwarnings('ignore')

CNS_GATE = {
    "MW":   (None, 460), "LogP": (2.0, 4.2),
    "TPSA": (None, 95),  "HBD":  (None, 3),
    "HBA":  (None, 7),   "QED":  (0.82, None),
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
    m = {"MW":round(mw,2), "LogP":round(logp,3), "TPSA":round(tpsa,2),
         "HBD":hbd, "HBA":hba, "QED":round(qed,3),
         "RotB":rotb, "Lipinski_violations":lip}
    reasons = []
    for k, v in gate.items():
        lo, hi = v
        val = m[k]
        if lo is not None and val < lo: reasons.append(f"{k}={round(val,2)}<{lo}")
        if hi is not None and val > hi: reasons.append(f"{k}={round(val,2)}>{hi}")
    return {"status": "PASS" if not reasons else "FAIL",
            "reasons": reasons, "metrics": m}

SCAFFOLD = "OC1(C2)CC3CC(C1)CC(C2)([*])C3"

# Gen5: ether_direct 고정 + 강직성 강화 꼬리
LINKERS = {
    "ether_direct":  "O{TAIL}",           # Gen4 hERG 최적, 유지
    "ether_methyl":  "OC{TAIL}",          # 메틸 추가 → 완충
}

TAILS = {
    # ── Gen4 기준점 (cyano+F 검증 조합)
    "cyano_fluorophenyl_4cn2f": "c1cc(F)cc(C#N)c1",   # Gen4 통과, 기준

    # ── 방향 1: 고리 강직성 강화 (MD 안정성 목표)
    "fluoronaphthalenyl":       "c1ccc2cc(F)ccc2c1",   # 나프탈렌 (고리 2개)
    "cyano_naphthalenyl":       "c1ccc2cc(C#N)ccc2c1", # 나프탈렌 + CN
    "benzofuranyl":             "c1ccc2occc2c1",        # 벤조푸란 (융합 고리)
    "indanyl":                  "C1CCc2ccccc21",        # 인단 (5+6 융합, 포화)

    # ── 방향 2: cyano+F 위치 변형
    "cf_ortho":    "c1ccccc1C#N",                       # 오르토 CN 단독
    "fcn_para":    "c1cc(F)ccc1C#N",                    # 파라 F + CN
    "difluoro_cn": "c1c(F)cc(C#N)cc1F",                # 이중 F + CN

    # ── 방향 3: sulfonyl 계열 확장 (Gen4 hERG 최저 0.304)
    "sulfonyl_cn_phenyl": "c1cc(C#N)ccc1S(N)(=O)=O",  # sulfonyl + CN
}

def run():
    rows = []
    for (ln, lt), (tn, ts) in itertools.product(LINKERS.items(), TAILS.items()):
        raw = SCAFFOLD.replace("[*]", lt.replace("{TAIL}", ts))
        mol = Chem.MolFromSmiles(raw)
        if mol is None:
            continue
        canon = Chem.MolToSmiles(mol)
        r = evaluate(canon)
        row = {"ID": f"{ln}__{tn}", "SMILES": canon,
               "STATUS": r["status"],
               "REASONS": " ; ".join(r["reasons"]) if r["reasons"] else "—"}
        row.update(r["metrics"])
        rows.append(row)
    return pd.DataFrame(rows)

if __name__ == "__main__":
    df = run()
    n     = len(df)
    npass = (df.STATUS == "PASS").sum()

    print(f"\nGen5: MD 안정성 집중 설계")
    print(f"생성 {n}개 / 통과 {npass}개 / 통과율 {npass/n*100:.0f}%")
    print("="*80)

    passed = df[df.STATUS == "PASS"].sort_values("QED", ascending=False)
    if len(passed):
        print("\n[통과 분자 — QED 순]")
        print(passed[["ID","MW","LogP","TPSA","QED","RotB"]].to_string(index=False))

    if (df.STATUS != "PASS").any():
        print("\n[탈락]")
        print(df[df.STATUS!="PASS"][["ID","REASONS"]].to_string(index=False))

    try:
        g4 = pd.read_csv("gen4_library.csv")
        print(f"\nGen4 통과율: {(g4.STATUS=='PASS').mean()*100:.0f}%  →  Gen5: {npass/n*100:.0f}%")
    except Exception:
        pass

    df.to_csv("gen5_library.csv", index=False)
    print(f"\n저장: gen5_library.csv")
