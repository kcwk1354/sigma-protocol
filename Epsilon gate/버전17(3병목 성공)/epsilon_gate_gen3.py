"""
Epsilon-Gate Gen3 — Gen2 전체 파이프라인 패턴 반영
====================================================
Gen2 학습 규칙:
  - ether_ethyl 링커 고정 (전 단계 압도적 우위)
  - pyridyl_4 기반 유지 (MD 유일 안정)
  - 링커 길이 변형 탐색 (propyl, direct)
  - 꼬리 치환기 추가 (methyl, fluoro, cyano + pyridyl_4)
  - bipyridyl 신규 시도 (도킹 강화 + 강직성)
  - 목표 LogP 3.0~3.5
"""

from rdkit import Chem
from rdkit.Chem import QED, Descriptors, Crippen, Lipinski
import pandas as pd, itertools

CNS_GATE = {
    "name": "CNS",
    "MW":   (None, 450), "LogP": (2.0,  4.0),
    "TPSA": (None, 90),  "HBD":  (None, 3),
    "HBA":  (None, 7),   "QED":  (0.5,  None),
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
        if k == "name": continue
        lo, hi = v
        v = m[k]
        if lo is not None and v < lo: reasons.append(f"{k}={round(v,2)}<{lo}")
        if hi is not None and v > hi: reasons.append(f"{k}={round(v,2)}>{hi}")
    return {"status": "PASS" if not reasons else "FAIL",
            "reasons": reasons, "metrics": m}

# ──────────────────────────────────────────
# Gen3 분자 설계
# 규칙: ether_ethyl + pyridyl_4 핵심 유지, 주변 변형
# ──────────────────────────────────────────

# 스캐폴드: 아다만탄-OH (동일 유지)
SCAFFOLD = "OC1(C2)CC3CC(C1)CC(C2)([*])C3"

# Gen3 링커: ether_ethyl 중심 + 길이 변형
LINKERS = {
    # ── 기존 최적 유지
    "ether_ethyl":   "OCC{TAIL}",       # Gen2 승자, 기준점
    # ── 링커 길이 변형
    "ether_propyl":  "OCCC{TAIL}",      # CH2 하나 추가 → 더 깊은 포켓 접근
    "ether_direct":  "O{TAIL}",         # 직결 → 강직성↑ 안정성 테스트
    # ── amide 계열 1개 유지 (비교군)
    "amide_ethyl":   "C(=O)NCC{TAIL}",  # Gen2 안정 링커, 비교용
}

# Gen3 꼬리: pyridyl_4 기반 변형 집중
TAILS = {
    # ── Gen2 유일 MD 통과, 기준점
    "pyridyl_4":          "c1ccncc1",           # 기준점 유지
    # ── 방향 1: pyridyl_4에 치환기 추가
    "methyl_pyridyl_4":   "c1cc(C)ncc1",        # 3-methyl → LogP↑ 소수성↑
    "fluoro_pyridyl_4":   "c1cc(F)ncc1",        # 3-fluoro → 대사 안정성↑
    "cyano_pyridyl_4":    "c1cc(C#N)ncc1",      # 3-cyano → HBA 추가
    # ── 방향 2: 새로운 구조 (강직성 강화)
    "bipyridyl":          "c1ccnc(c1)-c1ccncc1", # 두 pyridyl 연결 → 강직↑
    "pyridyl_4_OH":       "c1cc(O)ncc1",         # 3-OH → 수소결합↑
    # ── 방향 3: 기존 좋은 꼬리 일부 유지 (비교)
    "fluorophenyl":       "c1ccc(F)cc1",         # Gen2 MD 2위, 비교군
}

def build_smiles(linker_template, tail_smiles):
    """SMILES 조합: 링커 + 꼬리"""
    # {TAIL} 자리에 꼬리 삽입
    if "{TAIL}" in linker_template:
        # 꼬리가 방향족이면 Cc1... 형태로 연결
        filled = linker_template.replace("{TAIL}", f"c1{tail_smiles[2:]}" 
                                         if tail_smiles.startswith("c1") 
                                         else tail_smiles)
        # 더 안전한 방법: 직접 SMILES 구성
        tail_part = tail_smiles
        filled = linker_template.replace("{TAIL}", tail_part)
    else:
        # ether_direct: O{TAIL}
        filled = linker_template.replace("{TAIL}", tail_smiles)
    return SCAFFOLD.replace("[*]", filled)

def run():
    rows = []
    for (ln, lt), (tn, ts) in itertools.product(LINKERS.items(), TAILS.items()):
        raw_smi = SCAFFOLD.replace("[*]", lt.replace("{TAIL}", ts))
        mol = Chem.MolFromSmiles(raw_smi)
        if mol is None:
            # SMILES 파싱 실패 시 스킵
            continue
        canon = Chem.MolToSmiles(mol)
        r = evaluate(canon)
        row = {
            "ID":     f"{ln}__{tn}",
            "SMILES": canon,
            "STATUS": r["status"],
            "REASONS": " ; ".join(r["reasons"]) if r["reasons"] else "—",
        }
        row.update(r["metrics"])
        rows.append(row)
    return pd.DataFrame(rows)

if __name__ == "__main__":
    df = run()
    n     = len(df)
    npass = (df.STATUS == "PASS").sum()

    print(f"\nGen3: 전 파이프라인 패턴 반영 설계")
    print(f"생성 {n}개 / 통과 {npass}개 / 통과율 {npass/n*100:.0f}%")
    print("="*90)

    passed = df[df.STATUS == "PASS"].sort_values("QED", ascending=False)
    print("\n[통과 분자 — QED 순]")
    print(passed[["ID","MW","LogP","TPSA","QED","HBA"]].to_string(index=False))

    if (df.STATUS != "PASS").any():
        print("\n[탈락]")
        print(df[df.STATUS != "PASS"][["ID","REASONS"]].to_string(index=False))

    # Gen2 vs Gen3 비교
    gen2 = pd.read_csv("gen2_library.csv")
    g2p  = (gen2.STATUS == "PASS").mean() * 100
    print(f"\n{'='*90}")
    print(f"Gen2 통과율: {g2p:.0f}%  →  Gen3 통과율: {npass/n*100:.0f}%")
    print(f"Gen2 최고 QED: {gen2[gen2.STATUS=='PASS']['QED'].max():.3f}")
    print(f"Gen3 최고 QED: {passed['QED'].max():.3f}")
    print(f"\nGen3 신규 탐색 영역:")
    new_linkers = [l for l in LINKERS if l not in ['ether_ethyl','amide_ethyl']]
    new_tails   = [t for t in TAILS   if t not in ['pyridyl_4','fluorophenyl']]
    print(f"  링커: {new_linkers}")
    print(f"  꼬리: {new_tails}")

    df.to_csv("gen3_library.csv", index=False)
    print(f"\n저장: gen3_library.csv")
