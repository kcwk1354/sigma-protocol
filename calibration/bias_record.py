# Sigma Protocol v3.0 - Calibration Bias Records
# 소스별 캘리브레이션 보정값 정의 및 적용 함수

BIAS_RECORDS = {
    "rdkit": {
        # 보정값 없음 - 단순 경고 구간만 정의
        "correction": None,
        "qed_warn_range": (0.48, 0.52),
        "warn_message": "QED 경계 구간 (0.48~0.52): 약효 가능성 판단 보류 권고",
    },
    "docking": {
        # ΔG < -8 kcal/mol 구간에서 30% 과소추정
        "underestimation_threshold": -8.0,
        "underestimation_factor": 0.30,
        "note": "ΔG < -8 구간 30% 과소추정 → corrected = raw / (1 - 0.30)",
    },
    "admet": {
        # COOH 구조 존재 시 경구 생체이용률(F%) 15~25%p 과소추정
        # hERG Safe 예측은 신뢰도 높음 (False Negative 낮음)
        "cooh_smarts": "[CX3](=O)[OX2H1]",
        "cooh_bioavailability_correction_low": 15,   # %p
        "cooh_bioavailability_correction_high": 25,  # %p
        "herg_safe_confidence": "high",
        "note": "COOH 구조: HIA F% +15~25%p 과소추정. hERG Safe 예측은 신뢰.",
    },
}


def get_rdkit_bias() -> dict:
    return BIAS_RECORDS["rdkit"]


def get_docking_bias() -> dict:
    return BIAS_RECORDS["docking"]


def get_admet_bias() -> dict:
    return BIAS_RECORDS["admet"]


def apply_docking_correction(raw_dg: float) -> dict:
    bias = get_docking_bias()
    if raw_dg < bias["underestimation_threshold"]:
        corrected = raw_dg / (1.0 - bias["underestimation_factor"])
        return {
            "raw_dG": raw_dg,
            "corrected_dG": round(corrected, 3),
            "correction_applied": True,
            "note": bias["note"],
        }
    return {
        "raw_dG": raw_dg,
        "corrected_dG": raw_dg,
        "correction_applied": False,
        "note": "임계값(ΔG < -8) 미달 - 보정 없음",
    }


def apply_admet_correction(predictions: dict, smiles: str) -> dict:
    from rdkit import Chem

    bias = get_admet_bias()
    mol = Chem.MolFromSmiles(smiles)
    cooh_pattern = Chem.MolFromSmarts(bias["cooh_smarts"])
    has_cooh = bool(mol and mol.HasSubstructMatch(cooh_pattern))

    result = dict(predictions)
    result["cooh_detected"] = has_cooh

    if has_cooh:
        low = bias["cooh_bioavailability_correction_low"] / 100.0
        high = bias["cooh_bioavailability_correction_high"] / 100.0
        mid_correction = (low + high) / 2.0

        # HIA_Hou (0~1 확률값)에 보정 적용
        if "HIA_Hou" in predictions:
            raw_val = float(predictions["HIA_Hou"])
            result["HIA_Hou_corrected"] = round(min(1.0, raw_val + mid_correction), 4)

        result["bioavailability_correction_applied"] = True
        result["correction_note"] = (
            f"COOH 구조 감지: F% +{bias['cooh_bioavailability_correction_low']}"
            f"~{bias['cooh_bioavailability_correction_high']}%p 보정 적용"
        )
    else:
        result["bioavailability_correction_applied"] = False
        result["correction_note"] = "COOH 미감지 - 보정 없음"

    return result
