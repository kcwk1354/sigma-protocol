# Sigma Protocol v3.0 - ADMET Agent
# 독성학/약동학 계산 (로컬 ADMET-AI) - API 호출 없음

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from calibration.bias_record import apply_admet_correction

# 주요 지표 목록 (표시 및 추출용)
KEY_METRICS = [
    "HIA_Hou",               # 경구 흡수율 (Human Intestinal Absorption)
    "Caco2_Wang",            # Caco-2 투과성
    "BBB_Martins",           # 혈뇌장벽 통과
    "hERG",                  # 심장 독성 (hERG 채널 차단)
    "hERG_at_1uM",
    "hERG_at_10uM",
    "AMES",                  # 변이원성
    "DILI",                  # 약물유발 간독성
    "ClinTox",               # 임상 독성
    "Skin_Reaction",         # 피부 반응
    "Lipophilicity_AstraZeneca",
    "Solubility_AqSolDB",    # 수용해도
    "CYP3A4_Substrate_CarbonMangels",
    "CYP2D6_Substrate_CarbonMangels",
    "CYP3A4_Inhibitor_Veith",
    "CYP2C9_Inhibitor_Veith",
]

_model_cache = None


def _get_model():
    global _model_cache
    if _model_cache is None:
        from admet_ai import ADMETModel
        _model_cache = ADMETModel()
    return _model_cache


def run_admet(smiles: str) -> dict:
    try:
        model = _get_model()
        preds = model.predict(smiles=smiles)

        # DataFrame 또는 dict 처리
        if hasattr(preds, "iloc"):
            raw = preds.iloc[0].to_dict()
        elif isinstance(preds, dict):
            raw = preds
        else:
            raw = {}

        # numpy 타입 → Python float 변환
        raw_serialized = _serialize(raw)

        # 캘리브레이션 보정 적용 (COOH 구조 감지 시 HIA 보정)
        corrected = apply_admet_correction(raw_serialized, smiles)

        return {
            "smiles": smiles,
            "raw_predictions": raw_serialized,
            "corrected_predictions": corrected,
            "key_metrics": _extract_key_metrics(corrected),
            "calibration_applied": corrected.get("bioavailability_correction_applied", False),
            "calibration_note": corrected.get("correction_note", ""),
        }

    except Exception as e:
        return {
            "error": str(e),
            "smiles": smiles,
            "hint": "ADMET-AI 모델 로딩 실패. 첫 실행 시 가중치 다운로드가 필요할 수 있습니다.",
        }


def _serialize(d: dict) -> dict:
    result = {}
    for k, v in d.items():
        if hasattr(v, "item"):          # numpy scalar
            result[k] = round(v.item(), 4)
        elif isinstance(v, float):
            result[k] = round(v, 4)
        else:
            result[k] = v
    return result


def _extract_key_metrics(preds: dict) -> dict:
    return {k: preds[k] for k in KEY_METRICS if k in preds}
