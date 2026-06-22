# Sigma Protocol v3.0 - Orchestrator env_agent 패치
# 변경: qoa_result, why_result 파라미터 제거

from env_agent import run_env_check


def run_orchestration_with_env(
    smiles: str,
    meta_result: dict,
    rdkit_result: dict,
    docking_result: dict,
    admet_result: dict,
) -> dict:
    # ── Step 0: 환경 점검 ────────────────────────────────
    env_result = run_env_check()

    if env_result["proceed"] == "halt":
        return {
            "final_score":    0,
            "verdict":        "Env-Halt",
            "risk_level":     "critical",
            "env_check":      env_result,
            "overall_summary": "환경 점검 실패 — 파이프라인 중단.\n" + env_result["summary"],
            "engine":         "env_agent (halt)",
        }

    # ── Step 1: 기존 orchestration ───────────────────────
    from orchestrator import run_orchestration
    orch_result = run_orchestration(
        smiles         = smiles,
        meta_result    = meta_result,
        rdkit_result   = rdkit_result,
        docking_result = docking_result,
        admet_result   = admet_result,
    )

    # ── Step 2: env 경고 병기 ────────────────────────────
    orch_result["env_check"] = {
        "status":    env_result["status"],
        "warnings":  env_result["warnings"],
        "conflicts": env_result["conflicts"],
    }
    if env_result["warnings"]:
        orch_result["overall_summary"] += (
            f"\n[ENV] {env_result['status']}: "
            + " / ".join(env_result["warnings"])
        )

    return orch_result


# ── 사용법 ───────────────────────────────────────────────
# 기존: result = run_orchestration(smiles, meta, rdkit, docking, admet)
# 변경: result = run_orchestration_with_env(smiles, meta, rdkit, docking, admet)
