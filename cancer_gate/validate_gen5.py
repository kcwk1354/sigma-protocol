import sys
sys.stdout.reconfigure(encoding='utf-8')
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

# CG4-003 base: OC(=O)C1(c2ccc(F)c(F)c2)CC1CN1CCCCC1
# cyclopropane C1(COOH)(aryl) - C2H2 - C3H(CH2N)
# Strategy:
#  CG5-001/002: cis/trans explicit stereochemistry
#  CG5-003:     piperidine N-oxide (N+ O- strong H-bond)
#  CG5-004:     piperazine (NH H-bond donor replaces piperidine)
#  CG5-005:     acetylamide (NHCO H-bond pair, removes piperidine ring)

candidates = [
    ("CG5-001",
     "OC(=O)[C@@]1(c2ccc(F)c(F)c2)C[C@H]1CN1CCCCC1",
     "CG4-003 cis (1S,3R): aryl/CH2N same face → tighter TM fit"),
    ("CG5-002",
     "OC(=O)[C@@]1(c2ccc(F)c(F)c2)C[C@@H]1CN1CCCCC1",
     "CG4-003 trans (1S,3S): aryl/CH2N opposite face → extended reach"),
    ("CG5-003",
     "OC(=O)C1(c2ccc(F)c(F)c2)CC1C[N+]1(CCCCC1)[O-]",
     "N-oxide of piperidine: strong H-bond acceptor/donor + polarity"),
    ("CG5-004",
     "OC(=O)C1(c2ccc(F)c(F)c2)CC1CN1CCNCC1",
     "Piperidine → Piperazine: NH H-bond donor + ring polarity"),
    ("CG5-005",
     "OC(=O)C1(c2ccc(F)c(F)c2)CC1CNC(C)=O",
     "Piperidine → acetylamide: NHCO H-bond pair, MW reduction"),
]

for cid, smi, note in candidates:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print(f"{cid}: INVALID — {smi}")
        continue
    mw   = round(Descriptors.MolWt(mol), 1)
    logp = round(Descriptors.MolLogP(mol), 2)
    sp3  = round(rdMolDescriptors.CalcFractionCSP3(mol), 3)
    rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    hbd  = rdMolDescriptors.CalcNumHBD(mol)
    hba  = rdMolDescriptors.CalcNumHBA(mol)
    tpsa = round(Descriptors.TPSA(mol), 1)
    print(f"{cid}: VALID  MW={mw}  LogP={logp}  Fsp3={sp3}  HBD={hbd}  HBA={hba}  TPSA={tpsa}")
    print(f"       {note}")
    print()
