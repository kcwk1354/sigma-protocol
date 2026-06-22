# Sigma Protocol v3.0 - Orchestrator
# 변경: qoa_result, why_result 파라미터 제거
# 점수 체계: RDKit 25점 + Docking 30점 + ADMET 35점 + 캘리브레이션 10점 = 100점

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))


def run_orchestration(
    smiles: str,
    meta_result: dict,
    rdkit_result: dict,
    docking_result: dict,
    admet_result: dict,
) -> dict:
    km   = admet_result.get("key_metrics", {})
    corr = admet_result.get("corrected_predictions", {})

    # ── RDKit 점수 (25점) ────────────────────────────────
    rdkit_score = 25
    ro5_v = rdkit_result.get("lipinski_ro5_violations", 0)
    rdkit_score -= ro5_v * 4
    if not rdkit_result.get("veber_pass", True):
        rdkit_score -= 4
    qed = rdkit_result.get("QED", 0.5)
    if qed < 0.3:   rdkit_score -= 6
    elif qed < 0.5: rdkit_score -= 2
    rdkit_score = max(0, rdkit_score)

    # ── Docking 점수 (30점) ──────────────────────────────
    dg = docking_result.get("corrected_dG_kcal_mol", -5.0)
    if   dg <= -10: docking_score = 30
    elif dg <= -8:  docking_score = 25
    elif dg <= -6:  docking_score = 18
    elif dg <= -4:  docking_score = 10
    else:           docking_score = 4

    # ── ADMET 점수 (35점) ────────────────────────────────
    admet_score = 35
    herg = _f(km.get("hERG") or km.get("hERG_at_10uM"))
    ames = _f(km.get("AMES"))
    dili = _f(km.get("DILI"))
    hia  = _f(corr.get("HIA_Hou_corrected") or km.get("HIA_Hou"))
    bbb  = _f(km.get("BBB_Martins"))
    clintox = _f(km.get("ClinTox"))

    critical_toxicity = []
    if herg    is not None and herg > 0.5:    admet_score -= 12; critical_toxicity.append(f"hERG({herg:.3f})")
    if ames    is not None and ames > 0.5:    admet_score -= 12; critical_toxicity.append(f"AMES({ames:.3f})")
    if dili    is not None and dili > 0.5:    admet_score -= 10; critical_toxicity.append(f"DILI({dili:.3f})")
    if hia     is not None and hia  < 0.4:    admet_score -= 5
    if bbb     is not None and bbb  < 0.2:    admet_score -= 3
    if clintox is not None and clintox > 0.6: admet_score -= 5;  critical_toxicity.append(f"ClinTox({clintox:.3f})")
    admet_score = max(0, admet_score)

    # ── 캘리브레이션 보너스 (10점) ───────────────────────
    calib_bonus  = 0
    calib_applied = []
    if admet_result.get("calibration_applied"):
        calib_bonus += 5; calib_applied.append(admet_result.get("calibration_note", "ADMET COOH 보정"))
    if docking_result.get("correction_applied"):
        calib_bonus += 5; calib_applied.append("Docking ΔG 보정")

    # ── 최종 점수 / 판정 ─────────────────────────────────
    final_score = min(100, rdkit_score + docking_score + admet_score + calib_bonus)

    if final_score >= 70:
        verdict    = "Go"
        risk_level = "low" if final_score >= 80 else "medium"
    elif final_score >= 40:
        verdict    = "Conditional"
        risk_level = "medium" if final_score >= 55 else "high"
    else:
        verdict    = "No-Go"
        risk_level = "critical" if critical_toxicity else "high"

    if len(critical_toxicity) >= 2:
        verdict     = "No-Go"
        risk_level  = "critical"
        final_score = min(final_score, 39)

    confidence = _confidence(final_score, km)
    strengths, weaknesses = _sw(rdkit_result, dg, herg, ames, dili, hia, bbb, qed, calib_applied)
    next_steps = _next_steps(verdict, critical_toxicity, rdkit_result, hia, dg)

    return {
        "final_score": final_score,
        "score_breakdown": {
            "rdkit_score":       rdkit_score,
            "docking_score":     docking_score,
            "admet_score":       admet_score,
            "calibration_bonus": calib_bonus,
        },
        "verdict":            verdict,
        "verdict_confidence": confidence,
        "strengths":          strengths,
        "weaknesses":         weaknesses,
        "calibration_adjustments_applied": calib_applied or ["없음"],
        "next_steps":         next_steps,
        "risk_level":         risk_level,
        "overall_summary": (
            f"{smiles[:30]}{'...' if len(smiles)>30 else ''} → "
            f"{verdict} (점수: {final_score}/100, 신뢰도: {confidence:.0%}, "
            f"리스크: {risk_level.upper()}). "
            + (f"치명적 독성: {', '.join(critical_toxicity)}." if critical_toxicity
               else "치명적 독성 미검출.")
        ),
        "engine": "rule-based (no API)",
    }


def _f(v):
    if v is None: return None
    try:    return float(v)
    except: return None


def _confidence(score, km):
    key_count  = sum(1 for k in ["hERG","AMES","DILI","HIA_Hou","BBB_Martins"] if km.get(k) is not None)
    completeness = key_count / 5.0
    extremity    = abs(score - 50) / 50.0
    return round(min(0.95, 0.5 + completeness * 0.3 + extremity * 0.2), 2)


def _sw(rdkit_result, dg, herg, ames, dili, hia, bbb, qed, calib_applied):
    strengths, weaknesses = [], []
    if rdkit_result.get("lipinski_pass"):
        strengths.append("Lipinski Ro5 준수")
    else:
        weaknesses.append(f"Lipinski {rdkit_result.get('lipinski_ro5_violations')}개 위반")
    if qed >= 0.6:  strengths.append(f"QED={qed:.3f} — 높은 약물 유사성")
    elif qed < 0.4: weaknesses.append(f"QED={qed:.3f} — 낮은 약물 유사성")
    if dg <= -8:    strengths.append(f"결합에너지 강함 (ΔG={dg})")
    elif dg > -5:   weaknesses.append(f"결합에너지 약함 (ΔG={dg})")
    if herg is not None:
        if herg <= 0.3: strengths.append(f"hERG 안전 ({herg:.3f})")
        elif herg > 0.5: weaknesses.append(f"hERG 차단 위험 ({herg:.3f})")
    if ames is not None and ames > 0.5: weaknesses.append(f"AMES 변이원성 ({ames:.3f})")
    if dili is not None and dili > 0.5: weaknesses.append(f"DILI 간독성 ({dili:.3f})")
    if hia  is not None and hia >= 0.7: strengths.append(f"경구 흡수율 우수 (HIA={hia:.3f})")
    elif hia is not None and hia < 0.4: weaknesses.append(f"경구 흡수율 낮음 (HIA={hia:.3f})")
    if calib_applied: strengths.append("캘리브레이션 보정 적용")
    return strengths[:4], weaknesses[:4]


def _next_steps(verdict, critical_toxicity, rdkit_result, hia, dg):
    steps = []
    if verdict == "No-Go":
        if critical_toxicity:
            steps.append(f"독성 우선 해결: {', '.join(critical_toxicity[:2])} scaffold 변경")
        steps.append("유사 구조 라이브러리에서 독성 회피 analog 탐색")
    elif verdict == "Conditional":
        if not rdkit_result.get("lipinski_pass"):
            steps.append("MW·LogP 최적화 (scaffold hopping)")
        if hia is not None and hia < 0.5:
            steps.append("경구 흡수 개선: 프로드러그 또는 제형 최적화")
        if dg > -6:
            steps.append("H-bond 최적화로 결합에너지 개선")
        steps.append("최적화 후 ADMET 재평가")
    else:
        steps.extend([
            "in vitro 효능 실험 (IC50/Ki 측정)",
            "대사 안정성 실험 (CYP assay)",
            "선택성 패널 스크리닝",
            "동물 PK 실험 계획 수립",
        ])
    return steps[:4]
