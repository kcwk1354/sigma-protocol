import sys
sys.stdout.reconfigure(encoding='utf-8')
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

# Gen7 전략: CG4-003 라세미 베이스 (Ki=4.6μM), COOH 유지 필수
# tetrazole 블랙리스트 (Gen6: hERG=0.82~0.93 확인)
# C3 수식 블랙리스트 (Gen6: docking 악화 확인)
#
# 탐색 축:
#  1. cyclopropane C2 브릿지 수식 (C2-F / C2-Me → 아직 미탐색)
#  2. piperidine 링 내부 치환 (4-Me, 4-F, 3,3-diF)
#  3. piperidine N-adjacent 메틸 (2-Me: N 기하학 변경)
#
# 목표: DILI<0.4 / hERG<0.3 / Ki<1μM

candidates = [
    (
        "CG7-001",
        "OC(=O)C1(c2ccc(F)c(F)c2)C(F)C1CN1CCCCC1",
        "C2-불소화 cyclopropane + COOH",
        "C2-H→C2-F: 대사 안정화 + cyclopropane 기하학 미조정 → TM pocket 밀착도 변화",
    ),
    (
        "CG7-002",
        "OC(=O)C1(c2ccc(F)c(F)c2)C(C)C1CN1CCCCC1",
        "C2-methyl cyclopropane + COOH",
        "C2-H→C2-Me: 소수성 클램프 TM5 방향 탐색 / C2 수식 비교기준",
    ),
    (
        "CG7-003",
        "OC(=O)C1(c2ccc(F)c(F)c2)CC1CN1CCC(C)CC1",
        "4-methyl-piperidine + COOH",
        "pip-4-Me: TM7 Val775/Phe817 소수성 포켓 심층 접촉 / 링 편평화 효과",
    ),
    (
        "CG7-004",
        "OC(=O)C1(c2ccc(F)c(F)c2)CC1CN1CCC(F)CC1",
        "4-fluoro-piperidine + COOH",
        "pip-4-F: C-F 지방족 결합 → 전자 효과 + gauche 효과 → pip 링 구조 고정",
    ),
    (
        "CG7-005",
        "OC(=O)C1(c2ccc(F)c(F)c2)CC1CN1CCC(F)(F)CC1",
        "3,3-gem-difluoro-piperidine + COOH",
        "pip-3,3-diF: gem-diF 대사 차단 + 링 평면성 감소 → N lone pair 방향 최적화",
    ),
]

print("=" * 72)
print("  Gen7 후보 RDKit 검증")
print("  베이스: CG4-003 (Ki=4.6μM) | 블랙리스트: tetrazole / C3-수식")
print("  탐색: C2-브릿지 수식 / piperidine 링 내부 치환")
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
print(f"  결과: {'전체 VALID → gen7 실행 가능' if all_valid else 'INVALID 존재 → 수정 필요'}")
print("=" * 72)
