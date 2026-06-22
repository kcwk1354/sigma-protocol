"""Gen10 SMILES validation — RDKit"""
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

candidates = [
    ("CG10-001", "OC(=O)C1(c2cccc(C(F)(F)F)c2)C(F)(F)C1CN1CC(F)(F)CCC1",
     "A1: 3-CF3-Ph"),
    ("CG10-002", "OC(=O)C1(c2ccc(Cl)c(F)c2)C(F)(F)C1CN1CC(F)(F)CCC1",
     "A2: 3-F,4-Cl-Ph"),
    ("CG10-003", "OC(=O)C1(c2cc(F)c(F)c(F)c2)C(F)(F)C1CN1CC(F)(F)CCC1",
     "A3: 3,4,5-triF-Ph"),
    ("CG10-004", "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1(C)CN1CC(F)(F)CCC1",
     "B1: C3-Me cyclopropane"),
    ("CG10-005", "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1C(C)C(F)(F)CCC1",
     "C1: 2-Me-3,3-diF-pip"),
]

print(f"{'ID':<12} {'Valid':<6} {'MW':>7} {'HBD':>5} {'HBA':>5} {'Note'}")
print("-" * 70)
all_valid = True
for cid, smi, note in candidates:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print(f"{cid:<12} {'FAIL':<6} {'N/A':>7} {'N/A':>5} {'N/A':>5} {note}")
        all_valid = False
    else:
        mw = Descriptors.MolWt(mol)
        hbd = rdMolDescriptors.CalcNumHBD(mol)
        hba = rdMolDescriptors.CalcNumHBA(mol)
        print(f"{cid:<12} {'OK':<6} {mw:>7.1f} {hbd:>5} {hba:>5} {note}")

print()
print(f"Result: {'ALL VALID' if all_valid else 'SOME INVALID'}")
