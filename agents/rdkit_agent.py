# Sigma Protocol v3.0 - RDKit Agent
# 구조화학 계산 (로컬 RDKit) - API 호출 없음

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from rdkit import Chem
from rdkit.Chem import Descriptors, QED, rdMolDescriptors
from calibration.bias_record import get_rdkit_bias


def run_rdkit(smiles: str) -> dict:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"error": "유효하지 않은 SMILES", "smiles": smiles}

    bias = get_rdkit_bias()

    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    qed_val = QED.qed(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    tpsa = Descriptors.TPSA(mol)
    rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    rings = rdMolDescriptors.CalcNumRings(mol)
    aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    formula = rdMolDescriptors.CalcMolFormula(mol)
    heavy_atoms = mol.GetNumHeavyAtoms()

    # Lipinski Rule of 5 위반 횟수
    ro5_violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])

    # Veber 규칙 (경구 생체이용률)
    veber_pass = rot_bonds <= 10 and tpsa <= 140

    # QED 경계 구간 경고 (캘리브레이션 적용)
    warn_lo, warn_hi = bias["qed_warn_range"]
    qed_warn = warn_lo <= qed_val <= warn_hi

    return {
        "smiles": smiles,
        "formula": formula,
        "heavy_atoms": heavy_atoms,
        "molecular_weight": round(mw, 3),
        "logP": round(logp, 3),
        "QED": round(qed_val, 4),
        "qed_warn": qed_warn,
        "qed_warn_message": bias["warn_message"] if qed_warn else None,
        "HBD": hbd,
        "HBA": hba,
        "TPSA": round(tpsa, 3),
        "rotatable_bonds": rot_bonds,
        "rings": rings,
        "aromatic_rings": aromatic_rings,
        "lipinski_ro5_violations": ro5_violations,
        "lipinski_pass": ro5_violations <= 1,
        "veber_pass": veber_pass,
        "calibration": "RDKit 보정값 없음 (QED 0.48~0.52 경계 경고만 적용)",
    }
