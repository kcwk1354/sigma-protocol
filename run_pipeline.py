# Sigma Protocol v3.0 - 메인 실행 파일
# 사용법: python run_pipeline.py [SMILES] [QUERY] [docking_mode]
# 변경: QoA(Step2), WHY(Step6) 제거 — Claude가 원본 데이터 직접 추론

import sys
import os

if sys.stdout.encoding and sys.stdout.encoding.lower() not in ("utf-8", "utf8"):
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")

import json
from pathlib import Path
from datetime import datetime

sys.path.insert(0, str(Path(__file__).parent))

from agents.meta_agent import check_ambiguity
from agents.rdkit_agent import run_rdkit
from agents.docking_agent import run_docking
from agents.admet_agent import run_admet
from agents.orchestrator import run_orchestration

DEFAULT_SMILES = "C1CCN(C1)c1nc2ccccc2n1CC(=O)O"
DEFAULT_QUERY  = "이 화합물의 신약 개발 가능성을 종합 평가하라."
LOG_PATH = Path(__file__).parent / "logs" / "pipeline_log.json"


class _JsonEncoder(json.JSONEncoder):
    def default(self, obj):
        if hasattr(obj, "item"):   return obj.item()
        if hasattr(obj, "tolist"): return obj.tolist()
        return super().default(obj)


def _save_log(data: dict) -> None:
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    logs = []
    if LOG_PATH.exists():
        try:
            logs = json.loads(LOG_PATH.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            logs = []
    logs.append(data)
    LOG_PATH.write_text(
        json.dumps(logs, ensure_ascii=False, indent=2, cls=_JsonEncoder),
        encoding="utf-8",
    )
    print(f"[LOG] 저장: {LOG_PATH}")


def _section(title: str) -> None:
    print(f"\n{'─'*60}")
    print(f"  {title}")
    print(f"{'─'*60}")


def run_pipeline(
    smiles: str = DEFAULT_SMILES,
    query:  str = DEFAULT_QUERY,
    docking_mode: str = "simulate",
) -> dict:
    print(f"\n{'='*60}")
    print(f"  Sigma Protocol v3.0")
    print(f"{'='*60}")
    print(f"  SMILES : {smiles}")
    print(f"  Query  : {query}")
    print(f"  Docking: {docking_mode} 모드")
    print(f"{'='*60}")

    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    results: dict = {
        "run_id":       run_id,
        "smiles":       smiles,
        "query":        query,
        "docking_mode": docking_mode,
        "timestamp":    datetime.now().isoformat(),
    }

    # ── Step 1: Meta Agent ─────────────────────────────
    _section("[1/4] Meta Agent — 모호성 체크")
    meta = check_ambiguity(query, smiles)
    results["meta"] = meta
    print(f"  모호성 점수: {meta.get('ambiguity_score', 'N/A')}")
    print(f"  진행 권고 : {meta.get('proceed_recommendation', 'N/A')}")
    print(f"  요약      : {meta.get('summary', 'N/A')}")

    if meta.get("proceed_recommendation") == "halt":
        print("\n[HALT] SMILES 무효 또는 치명적 모호성 — 파이프라인 종료.")
        results["halt_reason"] = meta.get("summary", "모호성 과다")
        _save_log(results)
        return results

    # ── Step 2: RDKit ──────────────────────────────────
    _section("[2/4] RDKit Agent — 구조화학 계산")
    rdkit = run_rdkit(smiles)
    results["rdkit"] = rdkit
    if "error" in rdkit:
        print(f"  [ERROR] {rdkit['error']}")
    else:
        print(f"  MW       : {rdkit.get('molecular_weight')} g/mol")
        print(f"  LogP     : {rdkit.get('logP')}")
        print(f"  QED      : {rdkit.get('QED')}")
        print(f"  TPSA     : {rdkit.get('TPSA')} Å²")
        print(f"  HBD/HBA  : {rdkit.get('HBD')} / {rdkit.get('HBA')}")
        print(f"  Lipinski : {'PASS' if rdkit.get('lipinski_pass') else 'FAIL'} "
              f"(위반 {rdkit.get('lipinski_ro5_violations')}개)")
        print(f"  Veber    : {'PASS' if rdkit.get('veber_pass') else 'FAIL'}")

    # ── Step 3: Docking ────────────────────────────────
    _section(f"[3/4] Docking Agent ({docking_mode} 모드)")
    docking = run_docking(smiles, mode=docking_mode)
    results["docking"] = docking
    if "error" in docking and "vina_fallback" not in docking:
        print(f"  [ERROR] {docking['error']}")
    else:
        print(f"  Raw ΔG      : {docking.get('raw_dG_kcal_mol')} kcal/mol")
        print(f"  Corrected ΔG: {docking.get('corrected_dG_kcal_mol')} kcal/mol")
        print(f"  추정 Ki     : {docking.get('estimated_Ki_uM')} μM")

    # ── Step 4: ADMET ──────────────────────────────────
    _section("[4/4] ADMET Agent — 독성/약동학")
    admet = run_admet(smiles)
    results["admet"] = admet
    if "error" in admet:
        print(f"  [ERROR] {admet['error']}")
    else:
        km = admet.get("key_metrics", {})
        print(f"  HIA   : {km.get('HIA_Hou', 'N/A')}")
        print(f"  BBB   : {km.get('BBB_Martins', 'N/A')}")
        print(f"  hERG  : {km.get('hERG', 'N/A')}")
        print(f"  AMES  : {km.get('AMES', 'N/A')}")
        print(f"  DILI  : {km.get('DILI', 'N/A')}")

    # ── Orchestrator ───────────────────────────────────
    _section("최종 취합 판정")
    orch = run_orchestration(smiles, meta, rdkit, docking, admet)
    results["orchestrator"] = orch

    print(f"\n{'='*60}")
    print(f"  ■ 최종 판정 : {orch.get('verdict', 'N/A')}")
    print(f"  ■ 종합 점수 : {orch.get('final_score', 'N/A')} / 100")
    sb = orch.get("score_breakdown", {})
    if sb:
        print(f"     RDKit {sb.get('rdkit_score',0):>3}pt | "
              f"Docking {sb.get('docking_score',0):>3}pt | "
              f"ADMET {sb.get('admet_score',0):>3}pt | "
              f"Calib +{sb.get('calibration_bonus',0):>2}pt")
    print(f"  ■ 신뢰도    : {orch.get('verdict_confidence', 'N/A')}")
    print(f"  ■ 리스크    : {orch.get('risk_level', 'N/A').upper()}")
    print(f"{'='*60}")

    for s in orch.get("strengths",  []): print(f"    + {s}")
    for w in orch.get("weaknesses", []): print(f"    - {w}")
    for i, step in enumerate(orch.get("next_steps", []), 1):
        print(f"    {i}. {step}")

    print(f"\n  [요약] {orch.get('overall_summary', 'N/A')}\n")

    _save_log(results)
    return results


if __name__ == "__main__":
    _smiles = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_SMILES
    _query  = sys.argv[2] if len(sys.argv) > 2 else DEFAULT_QUERY
    _mode   = sys.argv[3] if len(sys.argv) > 3 else "simulate"
    run_pipeline(smiles=_smiles, query=_query, docking_mode=_mode)
