import sys
sys.stdout.reconfigure(encoding='utf-8')
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

candidates = [
    ("CG4-001", "OC(=O)C1(c2ccc(F)c(F)c2)CC1CCN1CCOCC1",
     "cyclopropane-COOH/diF-Ph + morpholine tail"),
    ("CG4-002", "OC(=O)C1(c2ccc(F)c(F)c2)CC1OC",
     "cyclopropane-COOH/diF-Ph + 2-OMe"),
    ("CG4-003", "OC(=O)C1(c2ccc(F)c(F)c2)CC1CN1CCCCC1",
     "cyclopropane-COOH/diF-Ph + piperidinyl-CH2"),
    ("CG4-004", "OC(=O)C1(c2ccc(F)c(F)c2)CC1c1ccncc1",
     "cyclopropane-COOH/diF-Ph + pyridyl"),
    ("CG4-005", "OC(=O)C1(c2ccc(F)c(F)c2)CC1(C)C",
     "cyclopropane-COOH/diF-Ph + gem-dimethyl"),
]

for cid, smi, note in candidates:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print(f"{cid}: INVALID")
        continue
    mw   = round(Descriptors.MolWt(mol), 1)
    logp = round(Descriptors.MolLogP(mol), 2)
    sp3  = round(rdMolDescriptors.CalcFractionCSP3(mol), 3)
    rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    hbd  = rdMolDescriptors.CalcNumHBD(mol)
    hba  = rdMolDescriptors.CalcNumHBA(mol)
    print(f"{cid}: VALID  MW={mw}  LogP={logp}  Fsp3={sp3}  ArRings={rings}  HBD={hbd}  HBA={hba} | {note}")
