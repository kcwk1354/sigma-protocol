"""
Epsilon-Gate Gen4 — hERG 낮추기 집중 설계
==========================================
Gen3 학습 규칙:
  - cyano 꼬리가 hERG 최저 (0.317) → 확장 탐색
  - ether_direct 링커가 hERG에 유리 → 고정
  - bipyridyl 완전 제외 (hERG 0.80)
  - 극성 강화 방향: sulfonamide, tetrazole, dicyano
  목표: hERG < 0.4, DILI < 0.4, BBB > 0.7
"""

from rdkit import Chem
from rdkit.Chem import QED, Descriptors, Crippen, Lipinski
import pandas as pd, itertools, sys, os, warnings
warnings.filterwarnings('ignore')

CNS_GATE = {
    "MW":   (None, 450), "LogP": (2.0, 4.0),
    "TPSA": (None, 90),  "HBD":  (None, 3),
    "HBA":  (None, 7),   "QED":  (0.85, None),  # Gen4: QED 기준 상향
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

# ──────────────────────────────────────────
# Gen4: hERG 낮추기 집중 설계
# 핵심: ether_direct 링커 + 극성 꼬리
# ──────────────────────────────────────────
SCAFFOLD = "OC1(C2)CC3CC(C1)CC(C2)([*])C3"

LINKERS = {
    # Gen3 hERG 최저 링커 → 고정
    "ether_direct":   "O{TAIL}",
    # 새로운 시도: 아민 링커 (수소결합 공여체 추가)
    "amino_direct":   "N{TAIL}",
    # 비교군 유지
    "amide_ethyl":    "C(=O)NCC{TAIL}",
}

TAILS = {
    # ── Gen3 최저 hERG 기준점 유지
    "cyano_pyridyl_4":     "c1ccncc1C#N",       # Gen3 hERG 최저, 기준점

    # ── 방향 1: cyano 위치 변형
    "cyano_pyridyl_3":     "c1cnccc1C#N",       # cyano 3번 위치
    "cyano_pyridyl_2":     "c1ccccn1C#N",       # cyano 2번 위치 (오르토)

    # ── 방향 2: 극성 극대화 (이중 치환)
    "dicyano_phenyl":      "c1cc(C#N)cc(C#N)c1", # 이중 cyano → 극성 극대화
    "cyano_fluorophenyl":  "c1cc(F)cc(C#N)c1",   # F + CN 조합
    "cyano_pyrimidyl":     "c1ncncc1C#N",         # pyrimidyl + cyano

    # ── 방향 3: 새로운 극성 꼬리
    "trifluoromethyl_py":  "c1ccncc1C(F)(F)F",   # CF3 → 강한 전자흡인
    "sulfonyl_phenyl":     "c1ccc(S(=O)(=O)N)cc1", # sulfonamide
    "hydroxymethyl_py":    "c1ccncc1CO",           # CH2OH → 극성↑, 수소결합↑
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
        row = {
            "ID":      f"{ln}__{tn}",
            "SMILES":  canon,
            "STATUS":  r["status"],
            "REASONS": " ; ".join(r["reasons"]) if r["reasons"] else "—",
        }
        row.update(r["metrics"])
        rows.append(row)
    return pd.DataFrame(rows)

if __name__ == "__main__":
    df = run()
    n     = len(df)
    npass = (df.STATUS == "PASS").sum()

    print(f"\nGen4: hERG 낮추기 집중 설계")
    print(f"생성 {n}개 / 통과 {npass}개 / 통과율 {npass/n*100:.0f}%")
    print("="*80)

    passed = df[df.STATUS == "PASS"].sort_values("QED", ascending=False)
    if len(passed):
        print("\n[통과 분자 — QED 순]")
        print(passed[["ID","MW","LogP","TPSA","QED","HBA"]].to_string(index=False))
    else:
        print("\n통과 분자 없음 — 기준 완화 필요")

    if (df.STATUS != "PASS").any():
        print("\n[탈락]")
        print(df[df.STATUS!="PASS"][["ID","REASONS"]].to_string(index=False))

    # 이전 세대 비교
    try:
        g3 = pd.read_csv("gen3_library.csv")
        g3p = (g3.STATUS == "PASS").mean() * 100
        print(f"\nGen3 통과율: {g3p:.0f}%  →  Gen4 통과율: {npass/n*100:.0f}%")
        if len(passed):
            print(f"Gen3 최고 QED: {g3[g3.STATUS=='PASS']['QED'].max():.3f}")
            print(f"Gen4 최고 QED: {passed['QED'].max():.3f}")
    except Exception:
        pass

    df.to_csv("gen4_library.csv", index=False)
    print(f"\n저장: gen4_library.csv")
