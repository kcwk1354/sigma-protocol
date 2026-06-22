import sys
sys.stdout.reconfigure(encoding='utf-8')
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

candidates = [
    ("CG9-001", "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1CCC(F)(F)CC1",
     "CF2-cyclopropane + 4,4-diF pip (CG8-005 × CG8-001 핵심 조합)"),
    ("CG9-002", "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1CCC(F)CC1",
     "CF2-cyclopropane + 4-F pip"),
    ("CG9-003", "OC(=O)C1(c2ccc(F)c(F)c2)C(CC)C1CN1CCC(F)(F)CC1",
     "C2-Et + 4,4-diF pip (CG8-001 C2-Me→Et 확장)"),
    ("CG9-004", "OC(=O)C1(c2cc(F)c(F)c(F)c2)C(F)(F)C1CN1CCCCC1",
     "2,3,4-triF-Ph + CF2-cyclopropane + pip (phenyl ring 추가 F)"),
    ("CG9-005", "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1CC(F)(F)CCC1",
     "CF2-cyclopropane + 3,3-diF pip (pip N-adjacent gem-diF)"),
]

print("=" * 70)
all_valid = True
for cid, smi, note in candidates:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print(f"{cid}: INVALID — {smi}"); all_valid = False; continue
    mw   = round(Descriptors.MolWt(mol), 1)
    logp = round(Descriptors.MolLogP(mol), 2)
    sp3  = round(rdMolDescriptors.CalcFractionCSP3(mol), 3)
    hbd  = rdMolDescriptors.CalcNumHBD(mol)
    hba  = rdMolDescriptors.CalcNumHBA(mol)
    tpsa = round(Descriptors.TPSA(mol), 1)
    aroR = rdMolDescriptors.CalcNumAromaticRings(mol)
    ro5  = "PASS" if mw<=500 and logp<=5 and hbd<=5 and hba<=10 else "FAIL"
    warn = ("  ⚠ logP>4.5" if logp>4.5 else "") + (f"  ⚠ aroR={aroR}" if aroR>1 else "")
    print(f"{cid}: VALID{warn}  MW={mw}  LogP={logp}  Fsp3={sp3}  HBD={hbd}  HBA={hba}  Ro5={ro5}")
    print(f"       {note}")
print("=" * 70)
print("ALL VALID" if all_valid else "INVALID EXISTS")
