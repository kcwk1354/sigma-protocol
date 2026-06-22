"""
Epsilon-Gate Gen6 — 논문 근거 기반 F/CN 위치 최적화
=====================================================
논문 근거:
  [1] PMC6152672 — 아다만탄 AChE 억제제:
      페닐 3번 위치(메타) 전자흡인기가 AChE 결합 최적
      결합 포켓: peripheral anionic site (Tyr124, Phe295)
  [2] Grychowska 2022 (EJMECH) — F 도입 hERG 감소 전략:
      F 위치에 따라 타겟 결합 활성 크게 달라짐
  [3] Molecules 2023 — F 위치 activity cliff:
      위치 변화만으로 활성 최대 1300배 차이

Gen5 트레이드오프 해석:
  F 두 개(difluoro) → hERG↓ but 도킹↓ (Gen5 실측 확인)
  논문 예측과 정확히 일치

Gen6 핵심 가설:
  3-CN + 4-F (메타 CN, 파라 F) = 논문 기준 최적 위치
  → AChE 결합 유지 (3번 전자흡인기) + hERG 억제 (F 극성)
"""

from rdkit import Chem
from rdkit.Chem import QED, Descriptors, Crippen, Lipinski
import pandas as pd, warnings
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

# Gen6: ether_direct 고정 + F/CN 위치 체계적 탐색
LINKER = "O{TAIL}"  # ether_direct 고정 (Gen4 hERG 최적 확인)

# 논문 [1] 근거: 3번 위치(메타) 전자흡인기 최적
# 논문 [2][3] 근거: F 위치가 결합 활성 결정적
TAILS = {
    # ── 기준점 (Gen4 통과)
    "ref_3cn_2f":    "c1cc(F)cc(C#N)c1",      # Gen4 기준 (2-F, 3-CN)

    # ── 논문 [1] 제안: 3-CN(메타) 고정 + F 위치 변형
    "3cn_4f":        "c1ccc(F)cc1C#N",         # 3-CN + 4-F (파라 F) ← 논문 최적 예측
    "3cn_5f":        "c1cc(F)c(C#N)cc1",       # 3-CN + 5-F (메타 F, CN 인접)
    "3cn_only":      "c1ccccc1C#N",            # 3-CN 단독 (F 없음, 순수 CN 효과)
    "4cn_3f":        "c1cc(F)cc(C#N)c1",       # 4-CN + 3-F (위치 교환)

    # ── 논문 [3] activity cliff 탐색: 동일 치환기, 다른 위치
    "2cn_4f":        "c1ccc(F)c(C#N)c1",       # 2-CN + 4-F (오르토 CN)
    "4cn_2f":        "c1ccc(C#N)cc1F",         # 4-CN + 2-F

    # ── 논문 [2] 근거: F를 고리 외부 위치로 (aliphatic F)
    "3cn_cf2":       "c1ccc(CC(F)F)cc1C#N",    # CF2 지방족 F + CN (새로운 시도)

    # ── Gen4 sulfonyl 계열 (hERG 0.304 최저) + 위치 최적화
    "sulfonyl_3cn":  "c1cc(C#N)ccc1S(N)(=O)=O", # sulfonyl + 3-CN
    "3f_sulfonyl":   "c1cc(F)ccc1S(N)(=O)=O",   # 3-F + sulfonyl
}

def run():
    rows = []
    for tn, ts in TAILS.items():
        raw = SCAFFOLD.replace("[*]", LINKER.replace("{TAIL}", ts))
        mol = Chem.MolFromSmiles(raw)
        if mol is None:
            print(f"  파싱실패: {tn}")
            continue
        canon = Chem.MolToSmiles(mol)
        r = evaluate(canon)
        row = {"ID": f"ether_direct__{tn}", "SMILES": canon,
               "STATUS": r["status"],
               "REASONS": " ; ".join(r["reasons"]) if r["reasons"] else "—",
               "논문근거": "PMC6152672+EJMECH2022+Molecules2023"}
        row.update(r["metrics"])
        rows.append(row)
    return pd.DataFrame(rows)

if __name__ == "__main__":
    df = run()
    n     = len(df)
    npass = (df.STATUS == "PASS").sum()

    print(f"\nGen6: 논문 근거 기반 F/CN 위치 최적화")
    print(f"근거: PMC6152672 | EJMECH 2022 | Molecules 2023")
    print(f"생성 {n}개 / 통과 {npass}개 / 통과율 {npass/n*100:.0f}%")
    print("="*80)

    passed = df[df.STATUS == "PASS"].sort_values("QED", ascending=False)
    if len(passed):
        print("\n[통과 분자 — QED 순]")
        print(passed[["ID","MW","LogP","TPSA","QED","HBA"]].to_string(index=False))

    if (df.STATUS != "PASS").any():
        print("\n[탈락]")
        print(df[df.STATUS!="PASS"][["ID","REASONS"]].to_string(index=False))

    df.to_csv("gen6_library.csv", index=False)
    print(f"\n저장: gen6_library.csv")
    print(f"\n핵심 가설: 3-CN(메타) + 4-F(파라) 조합")
    print(f"  → 논문 [1]: 3번 전자흡인기 = AChE 결합 최적")
    print(f"  → 논문 [2]: F 위치 = hERG 조절 결정적")
