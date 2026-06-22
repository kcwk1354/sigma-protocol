from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

cands = [
    ("CG12-001", "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1C(C)(C)C(F)CCC1",
     "gem-diMe + 3-monoF pip"),
    ("CG12-002", "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1C(C)(C)CCCC1",
     "gem-diMe + plain pip"),
    ("CG12-003", "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1C(C)(C)CC(F)(F)CC1",
     "gem-diMe + 4,4-diF pip"),
    ("CG12-004", "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1C(C)C(F)(F)C(O)CC1",
     "2-Me + 3,3-diF + 4-OH pip"),
    ("CG12-005", "OC(=O)C1(c2ccc(F)c(F)c2)C(F)C1CN1C(C)(C)C(F)(F)CCC1",
     "gem-diMe + 3,3-diF + CF-cyclopropane"),
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
print("CG10-005 ref: MW=381.3  LogP=3.67  (hERG=0.259)")
print("CG11-002 ref: MW=395.3  LogP=~4.1  (hERG=0.344, FAIL)")
print()
print("Result:", "ALL VALID" if all_ok else "SOME INVALID")
