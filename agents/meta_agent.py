# Sigma Protocol v3.0 - Meta Agent
# 4문항 모호성 체크: 순수 Python 룰베이스 (API 호출 없음)

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))


# 목표 타겟 관련 키워드
_TARGET_KEYWORDS = [
    "cdk2", "egfr", "vegfr", "her2", "akt", "pi3k", "kras", "mtor",
    "jak", "raf", "braf", "parp", "bcl", "protease", "kinase", "receptor",
    "단백질", "타겟", "target", "inhibitor", "억제", "길항", "agonist",
]
# 평가 기준 키워드
_CRITERIA_KEYWORDS = [
    "독성", "흡수", "admet", "약효", "효능", "ic50", "ki", "kd", "ec50",
    "bioavailability", "생체이용률", "대사", "안전성", "평가", "분석",
    "선택성", "specificity", "selectivity", "toxicity",
]
# 비교 대상 키워드
_BASELINE_KEYWORDS = [
    "대비", "비교", "baseline", "reference", "기준", "대조", "vs",
    "compared", "known drug", "기존 약물", "lead compound",
]


def _answer(flag: bool, partial_flag: bool = False) -> tuple[str, float]:
    if flag:
        return "Yes", 0.0
    if partial_flag:
        return "Partial", 0.5
    return "No", 1.0


def check_ambiguity(query: str, smiles: str) -> dict:
    query_lower = query.lower()

    # Q1. 목표 단백질/타겟이 명시되었는가?
    full_target = any(kw in query_lower for kw in _TARGET_KEYWORDS[:10])
    partial_target = any(kw in query_lower for kw in _TARGET_KEYWORDS[10:])
    q1_ans, q1_score = _answer(full_target, partial_target)
    q1_note = (
        "타겟 단백질 명시됨" if q1_ans == "Yes"
        else "일반적 평가 목적만 언급 (타겟 미명시)" if q1_ans == "No"
        else "평가 목적은 있으나 구체적 타겟 미명시"
    )

    # Q2. 평가 기준이 충분히 구체적인가?
    full_criteria = sum(1 for kw in _CRITERIA_KEYWORDS if kw in query_lower) >= 2
    partial_criteria = any(kw in query_lower for kw in _CRITERIA_KEYWORDS)
    q2_ans, q2_score = _answer(full_criteria, partial_criteria)
    q2_note = (
        "복수 평가 기준 명시" if q2_ans == "Yes"
        else "평가 기준 미명시 — 종합 평가로 해석" if q2_ans == "No"
        else "단일 평가 기준만 언급"
    )

    # Q3. SMILES가 유효한가? (RDKit 검증)
    q3_ans, q3_score, q3_note = _validate_smiles(smiles)

    # Q4. 비교 대상(baseline)이 존재하는가?
    has_baseline = any(kw in query_lower for kw in _BASELINE_KEYWORDS)
    q4_ans, q4_score = _answer(has_baseline)
    q4_note = (
        "비교 대상 언급됨" if has_baseline
        else "baseline 미명시 — 절대 평가 수행"
    )

    ambiguity_score = round((q1_score + q2_score + q3_score + q4_score) / 4.0, 3)

    # SMILES 무효만 halt — 나머지는 clarify/proceed
    if q3_ans == "No":
        proceed = "halt"
    elif ambiguity_score <= 0.40:
        proceed = "proceed"
    else:
        proceed = "clarify"  # 타겟·baseline 미명시는 clarify로 계속 진행

    return {
        "q1_target":   {"answer": q1_ans, "note": q1_note},
        "q2_criteria": {"answer": q2_ans, "note": q2_note},
        "q3_smiles":   {"answer": q3_ans, "note": q3_note},
        "q4_baseline": {"answer": q4_ans, "note": q4_note},
        "ambiguity_score": ambiguity_score,
        "proceed_recommendation": proceed,
        "summary": (
            f"모호성 점수 {ambiguity_score:.2f} — "
            f"타겟:{q1_ans} / 기준:{q2_ans} / SMILES:{q3_ans} / baseline:{q4_ans}. "
            f"권고: {proceed}."
        ),
        "engine": "rule-based (no API)",
    }


def _validate_smiles(smiles: str) -> tuple[str, float, str]:
    if not smiles or not smiles.strip():
        return "No", 1.0, "SMILES 비어 있음"
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return "No", 1.0, "RDKit 파싱 실패 — 유효하지 않은 SMILES"
        return "Yes", 0.0, f"RDKit 검증 통과 (원자 {mol.GetNumHeavyAtoms()}개)"
    except ImportError:
        # RDKit 없으면 문법 간단 체크
        depth = 0
        for c in smiles:
            if c == '(':
                depth += 1
            elif c == ')':
                depth -= 1
            if depth < 0:
                return "No", 1.0, "괄호 불균형"
        if depth != 0:
            return "No", 1.0, f"미닫힌 괄호 {depth}개"
        return "Partial", 0.5, "RDKit 미설치 — 문법 간단 검증만 수행"
