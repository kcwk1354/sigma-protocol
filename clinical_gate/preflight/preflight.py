import requests
import json
import datetime

OUTPUT_DIR = r"C:\sigma_protocol\clinical_gate\preflight"

def fetch_studies():
    url = "https://clinicaltrials.gov/api/v2/studies"
    params = {
        "filter.overallStatus": "TERMINATED",
        "filter.studyType": "INTERVENTIONAL",
        "pageSize": 50,
        "format": "json",
        "fields": "NCTId,Phase,OverallStatus,WhyStopped,PrimaryOutcomeMeasure,Condition,InterventionType"
    }
    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    return resp.json()

def extract_field(study, *keys):
    try:
        node = study
        for k in keys:
            node = node[k]
        return node
    except (KeyError, TypeError):
        return None

def main():
    print("Fetching TERMINATED studies from ClinicalTrials.gov...")
    raw = fetch_studies()

    with open(f"{OUTPUT_DIR}/preflight_raw.json", "w", encoding="utf-8") as f:
        json.dump(raw, f, ensure_ascii=False, indent=2)
    print(f"Raw data saved. Total studies in response: {len(raw.get('studies', []))}")

    studies = raw.get("studies", [])
    total = len(studies)

    records = []
    for s in studies:
        try:
            nct_id = extract_field(s, "protocolSection", "identificationModule", "nctId")
            phase_list = extract_field(s, "protocolSection", "designModule", "phases")
            phase = phase_list[0] if phase_list else None
            overall_status = extract_field(s, "protocolSection", "statusModule", "overallStatus")
            why_stopped = extract_field(s, "protocolSection", "statusModule", "whyStopped")
            primary_outcomes = extract_field(s, "protocolSection", "outcomesModule", "primaryOutcomes")
            primary_outcome = primary_outcomes[0].get("measure") if primary_outcomes else None
            conditions = extract_field(s, "protocolSection", "conditionsModule", "conditions")
            condition = conditions[0] if conditions else None
            interventions = extract_field(s, "protocolSection", "armsInterventionsModule", "interventions")
            intervention_type = interventions[0].get("type") if interventions else None

            records.append({
                "NCTId": nct_id,
                "Phase": phase,
                "OverallStatus": overall_status,
                "WhyStopped": why_stopped,
                "PrimaryOutcomeMeasure": primary_outcome,
                "Condition": condition,
                "InterventionType": intervention_type
            })
        except Exception as e:
            records.append({
                "NCTId": None, "Phase": None, "OverallStatus": None,
                "WhyStopped": None, "PrimaryOutcomeMeasure": None,
                "Condition": None, "InterventionType": None
            })

    # WhyStopped analysis
    ws_exists = [r for r in records if r["WhyStopped"] is not None]
    ws_blank = [r for r in records if r["WhyStopped"] == ""]
    ws_null = [r for r in records if r["WhyStopped"] is None]
    ws_effective_blank = len(ws_blank) + len(ws_null)
    ws_blank_rate = ws_effective_blank / total * 100 if total > 0 else 0

    # PrimaryOutcome analysis
    po_exists = [r for r in records if r["PrimaryOutcomeMeasure"] is not None and r["PrimaryOutcomeMeasure"] != ""]
    po_blank = [r for r in records if not r["PrimaryOutcomeMeasure"]]

    # Phase distribution
    phase_counts = {"PHASE1": 0, "PHASE2": 0, "PHASE3": 0, "other": 0}
    for r in records:
        p = r["Phase"]
        if p in phase_counts:
            phase_counts[p] += 1
        else:
            phase_counts["other"] += 1

    # Judgment
    if ws_blank_rate < 30:
        judgment = "Gen0 진행 가능. 설계 유지."
    elif ws_blank_rate <= 60:
        judgment = "경고: failure_token 보완 필요. Gen0 전에 대안 검토."
    else:
        judgment = "설계 변경 필요. WhyStopped 단독 의존 불가."

    # WhyStopped samples
    ws_samples = [r for r in records if r["WhyStopped"] and r["WhyStopped"] != ""][:10]

    # Console output
    print("\n" + "="*50)
    print("=== Clinical-Gate Preflight Report ===")
    print(f"날짜: {datetime.date.today()}")
    print(f"fetch 건수: {total}")
    print()
    print("[WhyStopped]")
    print(f"  존재: {len(ws_exists)}건 ({len(ws_exists)/total*100:.1f}%)")
    print(f"  공백(''): {len(ws_blank)}건 ({len(ws_blank)/total*100:.1f}%)")
    print(f"  완전 누락(null): {len(ws_null)}건 ({len(ws_null)/total*100:.1f}%)")
    print(f"  실질 공백률: {ws_blank_rate:.1f}%")
    print(f"  판단: {judgment}")
    print()
    print("[PrimaryOutcomeMeasure]")
    print(f"  존재: {len(po_exists)}건 ({len(po_exists)/total*100:.1f}%)")
    print(f"  공백/누락: {len(po_blank)}건 ({len(po_blank)/total*100:.1f}%)")
    print()
    print("[Phase 분포]")
    print(f"  PHASE1: {phase_counts['PHASE1']}건")
    print(f"  PHASE2: {phase_counts['PHASE2']}건")
    print(f"  PHASE3: {phase_counts['PHASE3']}건")
    print(f"  기타: {phase_counts['other']}건")
    print()
    print("[WhyStopped 샘플 10건]")
    for r in ws_samples:
        print(f"  {r['NCTId']} | \"{r['WhyStopped']}\"")

    # Save report
    report = f"""=== Clinical-Gate Preflight Report ===
날짜: {datetime.date.today()}
fetch 건수: {total}

[WhyStopped]
존재: {len(ws_exists)}건 ({len(ws_exists)/total*100:.1f}%)
공백: {len(ws_blank)}건 ({len(ws_blank)/total*100:.1f}%)
실질 공백률: {ws_blank_rate:.1f}%
판단: {judgment}

[PrimaryOutcomeMeasure]
존재: {len(po_exists)}건 ({len(po_exists)/total*100:.1f}%)
공백: {len(po_blank)}건 ({len(po_blank)/total*100:.1f}%)

[Phase 분포]
PHASE1: {phase_counts['PHASE1']}건
PHASE2: {phase_counts['PHASE2']}건
PHASE3: {phase_counts['PHASE3']}건
기타: {phase_counts['other']}건

[WhyStopped 샘플 10건]
""" + "\n".join(f"{r['NCTId']} | \"{r['WhyStopped']}\"" for r in ws_samples)

    with open(f"{OUTPUT_DIR}/preflight_report.txt", "w", encoding="utf-8") as f:
        f.write(report)
    print("\npreflight_report.txt 저장 완료.")

if __name__ == "__main__":
    main()
