import sys
sys.stdout.reconfigure(encoding='utf-8')
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

# Gen6 전략: CG4-003 라세미 베이스 (OC(=O)C1(c2ccc(F)c(F)c2)CC1CN1CCCCC1)
# Option A: cyclopropane C3에 소수성 그룹 추가 → TM5 클램프 확장
# Option C: COOH → tetrazole bioisostere → 이온 접촉 개선
#
# cyclopropane 구조:
#   C1(COOH/tetrazole)(3,4-diF-Ph)  ← C_A (quaternary, C1)
#   C_B (CH2 bridge)
#   C3(CH2pip)(new group)            ← C_C (C3 position)
#
# 고유 ring-closure 번호 사용 (재사용 없음):
#   ring1=main cyclopropane, ring2=diF-Ph, ring3=tetrazole, ring4=pip/cyclopropyl, ring5=pip

candidates = [
    (
        "CG6-001",
        "C1(c2ccc(F)c(F)c2)(c3nnn[nH]3)CC1CN4CCCCC4",
        "Tetrazole bioisostere (COOH→tetrazole, C3 미변경)",
        "Option C: 순수 bioisostere 기준점 / 이온 접촉 개선",
    ),
    (
        "CG6-002",
        "C1(c2ccc(F)c(F)c2)(c3nnn[nH]3)CC1(C)CN4CCCCC4",
        "Tetrazole + C3-methyl",
        "Option A+C: 소수성 +CH3 → TM5 Val/Leu 클램프 접촉 / Fsp3 유지",
    ),
    (
        "CG6-003",
        "C1(c2ccc(F)c(F)c2)(c3nnn[nH]3)CC1(OC)CN4CCCCC4",
        "Tetrazole + C3-OMe",
        "Option A+C: -OMe H-bond acceptor → TM Ser/Thr 추가 접촉 / 극성 클램프",
    ),
    (
        "CG6-004",
        "C1(c2ccc(F)c(F)c2)(c3nnn[nH]3)CC1(C2CC2)CN4CCCCC4",
        "Tetrazole + C3-cyclopropyl",
        "Option A+C: cyclopropyl Walsh orbital → 추가 π-stacking / Fsp3 최대화",
    ),
    (
        "CG6-005",
        "OC(=O)C1(c2ccc(F)c(F)c2)CC1(OC)CN3CCCCC3",
        "COOH 유지 + C3-OMe (대조군)",
        "Option A only: COOH 기준 C3-OMe 효과 단독 검증 / CG6-003과 직접 비교",
    ),
]

print("=" * 72)
print("  Gen6 후보 RDKit 검증")
print("  베이스: CG4-003 (Ki=4.6μM) | 목표: Ki<1μM, DILI<0.4, hERG<0.3")
print("=" * 72)

all_valid = True
for cid, smi, name, rationale in candidates:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print(f"\n{cid}: INVALID — {smi}")
        all_valid = False
        continue

    mw    = round(Descriptors.MolWt(mol), 1)
    logp  = round(Descriptors.MolLogP(mol), 2)
    sp3   = round(rdMolDescriptors.CalcFractionCSP3(mol), 3)
    hbd   = rdMolDescriptors.CalcNumHBD(mol)
    hba   = rdMolDescriptors.CalcNumHBA(mol)
    tpsa  = round(Descriptors.TPSA(mol), 1)
    rings = rdMolDescriptors.CalcNumRings(mol)
    aroR  = rdMolDescriptors.CalcNumAromaticRings(mol)

    # Lipinski Ro5 check
    ro5 = "PASS" if (mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10) else "FAIL"
    # DILI pre-screen: flag high logP or naphthalene-like
    warn = []
    if logp > 4.5:
        warn.append("logP>4.5→DILI위험")
    if aroR > 1:
        warn.append(f"방향족{aroR}개→DILI확인필요")
    if mw > 450:
        warn.append("MW>450")

    flag = "  ⚠ " + " / ".join(warn) if warn else "  ✓"

    print(f"\n{cid}: VALID{flag}")
    print(f"  이름: {name}")
    print(f"  SMILES: {smi}")
    print(f"  MW={mw}  LogP={logp}  Fsp3={sp3}  HBD={hbd}  HBA={hba}  TPSA={tpsa}  Ro5={ro5}")
    print(f"  방향족링={aroR}  전체링={rings}")
    print(f"  전략: {rationale}")

print("\n" + "=" * 72)
if all_valid:
    print("  전체 5개 VALID → cancer_gate_pipeline.py gen6 실행 가능")
else:
    print("  일부 INVALID → SMILES 수정 필요")
print("=" * 72)
