from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

cands = [
    # Track A: CG12-004 scaffold + OH position shift
    ("CG13-A1", "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1C(C)C(F)(F)CC(O)C1",
     "A: 5-OH  (C4→C5)"),
    ("CG13-A2", "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1C(C)C(F)(F)CCC1O",
     "A: 6-OH  (C4→C6, alpha-OH)"),
    # Track B: CG10-005 base + 4-position substituent (vs 4-OH=CG12-004)
    ("CG13-B1", "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1C(C)C(F)(F)C(F)CC1",
     "B: 4-F   (F bioisostere, no H-bond donor)"),
    ("CG13-B2", "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1C(C)C(F)(F)C(OC)CC1",
     "B: 4-OMe (methoxy, polar but no donor)"),
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
print("CG10-005  MW=381.3  LogP=3.67  HBD=1  Ki=0.345μM  hERG=0.259")
print("CG12-004  MW=397.3  LogP=2.64  HBD=2  Ki=0.580μM  hERG=0.154  (4-OH)")
print()
print("Result:", "ALL VALID" if all_ok else "SOME INVALID")
