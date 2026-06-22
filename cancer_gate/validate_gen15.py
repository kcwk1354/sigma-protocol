from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

cands = [
    ("CG15-001", "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1CC(F)(F)C(F)CC1",
     "N-deMe+3,3-diF+4-F pip"),
    ("CG15-002", "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1C(C)C(F)C(F)CC1",
     "3,4-diF pip (gem→vicinal 재배치)"),
    ("CG15-003", "OC(=O)C1(c2cc(F)c(F)cc2)C(F)(F)C1CN1C(C)C(F)(F)C(F)CC1",
     "Ph-3,4-diF+pip-4-F"),
    ("CG15-004", "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1C(C)C(F)(F)C(Cl)CC1",
     "4-F→4-Cl (바이오이소스터)"),
]

print(f"{'ID':<12} {'Valid':<5} {'MW':>7} {'LogP':>6} {'HBD':>4} {'HBA':>4}  Note")
print("-" * 75)
all_ok = True
for cid, smi, note in cands:
    mol = Chem.MolFromSmiles(smi)
    if mol:
        mw = Descriptors.MolWt(mol)
        lp = Descriptors.MolLogP(mol)
        hbd = rdMolDescriptors.CalcNumHBD(mol)
        hba = rdMolDescriptors.CalcNumHBA(mol)
        print(f"{cid:<12} {'OK':<5} {mw:>7.1f} {lp:>6.2f} {hbd:>4} {hba:>4}  {note}")
    else:
        print(f"{cid:<12} FAIL                              {note}")
        all_ok = False

print()
print("--- 참조 (DILI 경계선) ---")
print("CG10-005  LogP=3.67  Ki=0.345uM  DILI=0.349 OK  hERG=0.259")
print("CG13-B1   LogP=3.62  Ki=0.372uM  DILI=0.404 FAIL hERG=0.233  ← 목표: 이 구조에서 DILI만 내리기")
print("CG12-004  LogP=2.64  Ki=0.580uM  DILI=0.312 OK  hERG=0.154")
print()
print("Gen15 목표: Ki < 0.400uM + DILI < 0.400 + hERG < 0.300")
print()
print("Result:", "ALL VALID" if all_ok else "SOME INVALID")
