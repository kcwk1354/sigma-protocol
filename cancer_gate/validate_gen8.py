import sys
sys.stdout.reconfigure(encoding='utf-8')
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

# Gen8: Gen7 최우수 기능 조합 synergy test
# CG7-002 (C2-Me, Ki=3.70) × CG7-004 (4-F pip, Ki=3.61) × CG7-005 (4,4-diF pip, Ki=2.67)
# 4,4-gem-diF: pip N1CCC(F)(F)CC1 (4위치 gem-diF)

candidates = [
    (
        "CG8-001",
        "OC(=O)C1(c2ccc(F)c(F)c2)C(C)C1CN1CCC(F)(F)CC1",
        "C2-Me + 4,4-diF-pip",
        "CG7-002 + CG7-005 조합: C2-Me 소수성 × 4,4-diF pip 링 고정 synergy 검증",
    ),
    (
        "CG8-002",
        "OC(=O)C1(c2ccc(F)c(F)c2)C(F)C1CN1CCC(F)(F)CC1",
        "C2-F + 4,4-diF-pip",
        "CG7-001 + CG7-005 조합: C2-F 기하학 × 4,4-diF pip synergy",
    ),
    (
        "CG8-003",
        "OC(=O)C1(c2ccc(F)c(F)c2)C(C)C1CN1CCC(F)CC1",
        "C2-Me + 4-F-pip",
        "CG7-002 + CG7-004 조합: C2-Me × 4-F pip synergy",
    ),
    (
        "CG8-004",
        "OC(=O)C1(c2ccc(F)c(F)c2)C(F)C1CN1CCC(F)CC1",
        "C2-F + 4-F-pip",
        "CG7-001 + CG7-004 조합: C2-F × 4-F pip synergy",
    ),
    (
        "CG8-005",
        "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1CCCCC1",
        "C2-gem-diF (CF2 cyclopropane) + pip",
        "C2 위치 gem-diF: cyclopropane 링 전자 효과 최대화 + 기하학 변형 탐색",
    ),
]

print("=" * 72)
print("  Gen8 후보 RDKit 검증")
print("  전략: Gen7 최우수 (C2-Me/C2-F + 4,4-diF/4-F pip) 조합 synergy")
print("=" * 72)

all_valid = True
for cid, smi, name, rationale in candidates:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print(f"\n{cid}: INVALID — {smi}")
        all_valid = False
        continue

    mw   = round(Descriptors.MolWt(mol), 1)
    logp = round(Descriptors.MolLogP(mol), 2)
    sp3  = round(rdMolDescriptors.CalcFractionCSP3(mol), 3)
    hbd  = rdMolDescriptors.CalcNumHBD(mol)
    hba  = rdMolDescriptors.CalcNumHBA(mol)
    tpsa = round(Descriptors.TPSA(mol), 1)
    aroR = rdMolDescriptors.CalcNumAromaticRings(mol)
    ro5  = "PASS" if (mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10) else "FAIL"

    warn = []
    if logp > 4.5: warn.append("logP>4.5")
    if aroR > 1:   warn.append(f"방향족{aroR}개")
    flag = "  ⚠ " + " / ".join(warn) if warn else "  ✓"

    print(f"\n{cid}: VALID{flag}")
    print(f"  이름: {name}")
    print(f"  SMILES: {smi}")
    print(f"  MW={mw}  LogP={logp}  Fsp3={sp3}  HBD={hbd}  HBA={hba}  TPSA={tpsa}  Ro5={ro5}")
    print(f"  전략: {rationale}")

print("\n" + "=" * 72)
print(f"  결과: {'전체 VALID → gen8 실행 가능' if all_valid else 'INVALID 존재 → 수정 필요'}")
print("=" * 72)
