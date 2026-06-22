"""
Epsilon-Gate 독성 필터 — ADMET-AI 기반
========================================
Gen0에서 수작업(ProTox-3 웹)으로 하던 독성 검사를 자동화.
파이프라인에서 도킹 전에 실행하여 독성 분자를 미리 제거.

설치:
    pip install admet-ai

파이프라인 위치:
    물성 필터 → [독성 필터] → 도킹 → MD → 최종 후보

핵심 독성 기준 (CNS 신약 기준):
    AMES          < 0.3   변이원성 없음 (발암 위험)
    hERG          < 0.4   심장독성 낮음
    DILI          < 0.4   약물유발 간손상 낮음
    BBB_Martins   > 0.5   뇌혈관장벽 통과 (CNS 필수)
    Bioavailability_Ma > 0.5  경구 생체이용률
    Carcinogens   < 0.2   발암성 없음
"""

import os, warnings
import pandas as pd
warnings.filterwarnings('ignore')

# ──────────────────────────────────────────
# 독성 게이트 기준값
# ──────────────────────────────────────────
TOX_GATE = {
    # 항목명              (방향,  기준값,  설명)
    "AMES":              ("<",   0.3,   "변이원성 — 낮을수록 안전"),
    "hERG":              ("<",   0.4,   "심장독성 — 낮을수록 안전"),
    "DILI":              ("<",   0.4,   "약물유발 간손상 — 낮을수록 안전"),
    "BBB_Martins":       (">",   0.5,   "뇌혈관장벽 통과 — 높을수록 CNS 유리"),
    "Bioavailability_Ma":(">",   0.5,   "경구 생체이용률 — 높을수록 좋음"),
    "Carcinogens_Lagunin":("<",  0.2,   "발암성 — 낮을수록 안전"),
}

# ──────────────────────────────────────────
# 단일 분자 독성 예측
# ──────────────────────────────────────────
def predict_tox(smiles, model=None):
    """
    SMILES → 독성 프로파일 딕셔너리
    model: ADMETModel 인스턴스 (없으면 자동 로드)
    """
    if model is None:
        from admet_ai import ADMETModel
        model = ADMETModel()

    result = model.predict(smiles=smiles)

    profile = {}
    for key in TOX_GATE:
        val = result.get(key, None)
        if val is not None:
            profile[key] = round(float(val[0]) if hasattr(val,'__getitem__') else float(val), 3)
        else:
            profile[key] = None
    return profile

# ──────────────────────────────────────────
# 독성 게이트 통과 판정
# ──────────────────────────────────────────
def tox_gate(profile):
    """
    독성 프로파일 → PASS/FAIL + 탈락 이유
    """
    reasons = []
    for key, (direction, threshold, desc) in TOX_GATE.items():
        val = profile.get(key)
        if val is None:
            continue
        if direction == "<" and val >= threshold:
            reasons.append(f"{key}={val:.2f}≥{threshold} ({desc})")
        elif direction == ">" and val <= threshold:
            reasons.append(f"{key}={val:.2f}≤{threshold} ({desc})")
    return "PASS" if not reasons else "FAIL", reasons

# ──────────────────────────────────────────
# 독성 등급
# ──────────────────────────────────────────
def tox_grade(profile):
    """위험 점수 합산 → 등급"""
    score = 0
    score += profile.get("AMES", 0) * 2          # 변이원성 가중치 높음
    score += profile.get("hERG", 0) * 2           # 심장독성 가중치 높음
    score += profile.get("DILI", 0)
    score += profile.get("Carcinogens_Lagunin", 0) * 2
    score -= profile.get("BBB_Martins", 0) * 0.5  # BBB는 높을수록 좋음
    score -= profile.get("Bioavailability_Ma", 0) * 0.5

    if score < 0.3:   return "🟢 안전"
    elif score < 0.6: return "🟡 주의"
    else:             return "🔴 위험"

# ──────────────────────────────────────────
# 배치 실행
# ──────────────────────────────────────────
def run_tox_filter(csv_path="gen3_library.csv", top_n=None):
    """
    라이브러리 CSV → 독성 필터 → 결과 저장
    """
    from admet_ai import ADMETModel

    print("\n" + "="*68)
    print("  Epsilon-Gate 독성 필터 (ADMET-AI)")
    print("  Gen0 ProTox-3 수작업 → 자동화 버전")
    print("="*68)

    df = pd.read_csv(csv_path)
    passed = df[df['STATUS'] == 'PASS']
    if top_n:
        passed = passed.nlargest(top_n, 'QED')

    print(f"\n  대상: {len(passed)}개 분자\n")
    print(f"  {'ID':<35} {'AMES':>5} {'hERG':>5} {'BBB':>5} {'DILI':>5} {'BA':>5}  {'등급':<8} {'판정'}")
    print("  " + "-"*80)

    print("  모델 로딩 중...", end="\r")
    model = ADMETModel()
    print("  모델 로딩 완료.  ")

    results = []
    for _, row in passed.iterrows():
        profile = predict_tox(row['SMILES'], model)
        status, reasons = tox_gate(profile)
        grade  = tox_grade(profile)

        ames = profile.get('AMES', 0)
        herg = profile.get('hERG', 0)
        bbb  = profile.get('BBB_Martins', 0)
        dili = profile.get('DILI', 0)
        ba   = profile.get('Bioavailability_Ma', 0)

        print(f"  {row['ID']:<35} {ames:>5.2f} {herg:>5.2f} {bbb:>5.2f} "
              f"{dili:>5.2f} {ba:>5.2f}  {grade:<8} {status}")

        results.append({
            'ID':           row['ID'],
            'SMILES':       row['SMILES'],
            'QED':          row['QED'],
            'LogP':         row['LogP'],
            'AMES':         ames,
            'hERG':         herg,
            'BBB_Martins':  bbb,
            'DILI':         dili,
            'Bioavailability': ba,
            'Carcinogens':  profile.get('Carcinogens_Lagunin', 0),
            'Tox_Grade':    grade,
            'Tox_Status':   status,
            'Tox_Reasons':  " | ".join(reasons) if reasons else "—",
        })

    res_df = pd.DataFrame(results)
    safe   = res_df[res_df['Tox_Status'] == 'PASS']
    failed = res_df[res_df['Tox_Status'] == 'FAIL']

    print(f"\n{'='*68}")
    print(f"  결과: {len(safe)}개 통과 / {len(failed)}개 탈락")

    if len(safe) > 0:
        print(f"\n  ✅ 독성 통과 분자:")
        for _, r in safe.sort_values('QED', ascending=False).iterrows():
            print(f"     {r['ID']:<35} {r['Tox_Grade']}")

    if len(failed) > 0:
        print(f"\n  ❌ 독성 탈락 분자:")
        for _, r in failed.iterrows():
            print(f"     {r['ID']:<35} → {r['Tox_Reasons']}")

    out = os.path.join(os.path.dirname(csv_path), 'gen3_tox_results.csv')
    res_df.to_csv(out, index=False)
    print(f"\n  저장: {out}")

    # 파이프라인 연계: 독성 통과 분자만 도킹으로 넘김
    safe_csv = os.path.join(os.path.dirname(csv_path), 'gen3_tox_passed.csv')
    safe.to_csv(safe_csv, index=False)
    print(f"  도킹 입력용: {safe_csv}")
    print(f"\n  ⚠️  기준 해석")
    print(f"     AMES < 0.3  → 변이원성 없음")
    print(f"     hERG < 0.4  → 심장독성 낮음")
    print(f"     BBB  > 0.5  → 뇌혈관장벽 통과 (CNS 필수)")
    print(f"     DILI < 0.4  → 간손상 위험 낮음")
    print(f"     BA   > 0.5  → 경구 투여 가능\n")

    return res_df

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--csv',   default='gen3_library.csv')
    parser.add_argument('--top_n', type=int, default=None)
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    run_tox_filter(
        csv_path=os.path.join(script_dir, args.csv),
        top_n=args.top_n
    )
