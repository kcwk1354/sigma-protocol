from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

cands = [
    ("CG14-001", "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1C(C)C(F)CC(O)C1",
     "3-monoF+4-OH pip"),
    ("CG14-002", "OC(=O)C1(c2ccc(F)c(F)c2)[C@@H](F)C1CN1C(C)C(F)(F)CC(O)C1",
     "CHF-cyclopropane+4-OH pip"),
    ("CG14-003", "OC(=O)C1(c2ccc(F)c(F)c2)[C@@H](F)C1CN1C(C)C(F)(F)CCC1",
     "CHF-cyclopropane only (대조군)"),
]

print(f"{'ID':<12} {'Valid':<5} {'MW':>7} {'LogP':>6} {'HBD':>4} {'HBA':>4}  Note")
print("-" * 72)
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
print("--- 참조 ---")
print("CG10-005  MW=381.3  LogP=3.67  HBD=1  Ki=0.345uM  hERG=0.259  (현재 최고 Ki)")
print("CG12-004  MW=397.3  LogP=2.64  HBD=2  Ki=0.580uM  hERG=0.154  (4-OH, hERG최우수)")
print("CG13-B1   MW=399.3  LogP=3.62  HBD=1  Ki=0.372uM  hERG=0.233  (4-F, DILI FAIL)")
print()
print("Result:", "ALL VALID" if all_ok else "SOME INVALID")
