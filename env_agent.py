# Sigma Protocol v3.0 - Environment Agent
# 파이프라인 실행 전 환경 자동 점검
# - tool_cards/*.yaml 스캔
# - 경로/설치 상태 확인
# - 충돌 감지 및 리포트
# API 호출 없음 (순수 Python 룰베이스)

import sys
import importlib
import subprocess
from pathlib import Path

# ── 경로 고정 (백업 원본 유지) ───────────────────────────
SIGMA_ROOT   = Path("C:/sigma_protocol")
VINA_EXE     = SIGMA_ROOT / "vina" / "vina.exe"
TOOL_CARD_DIR = Path(__file__).parent / "tool_cards"

# ── 툴별 점검 스펙 ───────────────────────────────────────
TOOL_CHECKS = {
    "rdkit": {
        "import": "rdkit",
        "version_attr": "rdkit.__version__",
        "card": "rdkit.yaml",
        "required": True,
    },
    "admet_ai": {
        "import": "admet_ai",
        "version_attr": None,
        "card": "admet_ai.yaml",
        "required": True,
    },
    "sa_score": {
        "import": "rdkit",          # RDKit 내장
        "version_attr": None,
        "card": "sa_score.yaml",
        "required": True,
    },
    "autodock_vina": {
        "import": None,             # 바이너리 직접 확인
        "exe_path": VINA_EXE,
        "card": "autodock_vina.yaml",
        "required": False,          # simulate 모드 폴백 가능
    },
    "openmm": {
        "import": "openmm",
        "version_attr": "openmm.__version__",
        "card": "openmm.yaml",
        "required": False,          # MD 단계에서만 필요
    },
}


def run_env_check() -> dict:
    """
    파이프라인 실행 전 전체 환경 점검.
    반환값은 orchestrator가 그대로 수신.
    """
    results = {}
    critical_failures = []
    warnings = []

    # ── 1. 각 툴 점검 ────────────────────────────────────
    for tool_name, spec in TOOL_CHECKS.items():
        result = _check_tool(tool_name, spec)
        results[tool_name] = result

        if not result["ok"]:
            if spec.get("required"):
                critical_failures.append(
                    f"[CRITICAL] {tool_name}: {result['issue']}"
                )
            else:
                warnings.append(
                    f"[WARN] {tool_name}: {result['issue']} "
                    f"→ {result.get('fallback', '폴백 없음')}"
                )

    # ── 2. Sigma 루트 경로 점검 ──────────────────────────
    path_check = _check_sigma_path()
    if not path_check["ok"]:
        warnings.append(f"[WARN] sigma_root: {path_check['issue']}")

    # ── 3. 충돌 감지 ─────────────────────────────────────
    conflicts = _detect_conflicts()

    # ── 4. 판정 ──────────────────────────────────────────
    if critical_failures:
        proceed = "halt"
        status  = "CRITICAL"
    elif warnings:
        proceed = "proceed_with_warnings"
        status  = "WARNING"
    else:
        proceed = "proceed"
        status  = "OK"

    final = {
        "status": status,
        "proceed": proceed,
        "tool_results": results,
        "critical_failures": critical_failures,
        "warnings": warnings,
        "conflicts": conflicts,
        "path_check": path_check,
        "summary": _make_summary(status, results, critical_failures, warnings),
        "engine": "rule-based (no API)",
    }

    # ── tool_card 자동 갱신 (실행마다 누적) ─────────────
    update_tool_cards(final)

    return final


# ── 툴 개별 점검 ─────────────────────────────────────────

def _check_tool(tool_name: str, spec: dict) -> dict:
    # 바이너리 직접 확인 (AutoDock Vina)
    if spec.get("exe_path"):
        exe = spec["exe_path"]
        if exe.exists():
            return {
                "ok": True,
                "tool": tool_name,
                "path": str(exe),
                "note": "바이너리 확인됨",
            }
        else:
            return {
                "ok": False,
                "tool": tool_name,
                "issue": f"바이너리 없음: {exe}",
                "fallback": "simulate 모드 자동 전환 (docking_agent.py)",
            }

    # Python import 확인
    if spec.get("import"):
        try:
            mod = importlib.import_module(spec["import"])
            version = "unknown"
            if spec.get("version_attr"):
                try:
                    version = eval(spec["version_attr"])
                except Exception:
                    pass
            return {
                "ok": True,
                "tool": tool_name,
                "version": version,
                "note": "import 성공",
            }
        except ImportError as e:
            return {
                "ok": False,
                "tool": tool_name,
                "issue": f"import 실패: {e}",
                "fallback": spec.get("fallback", "설치 필요"),
            }

    return {"ok": True, "tool": tool_name, "note": "점검 불필요"}


def _check_sigma_path() -> dict:
    if SIGMA_ROOT.exists():
        return {
            "ok": True,
            "path": str(SIGMA_ROOT),
            "note": "sigma_root 존재 확인",
        }
    return {
        "ok": False,
        "path": str(SIGMA_ROOT),
        "issue": f"sigma_root 없음: {SIGMA_ROOT} — 백업 경로 확인 필요",
    }


def _detect_conflicts() -> list:
    """
    알려진 충돌 패턴 감지.
    현재 등록: pandas 버전 충돌 (RDKit + ADMET-AI)
    """
    conflicts = []
    try:
        import pandas as pd
        version = tuple(int(x) for x in pd.__version__.split(".")[:2])
        if version >= (2, 2):
            conflicts.append(
                f"pandas {pd.__version__} ≥ 2.2 감지 — "
                "RDKit PandasTools 충돌 가능. "
                "권고: pip install 'pandas>=2.0,<2.2'"
            )
    except ImportError:
        pass
    return conflicts


# ── 요약 메시지 ──────────────────────────────────────────

def _make_summary(status, results, critical_failures, warnings) -> str:
    ok_tools  = [t for t, r in results.items() if r["ok"]]
    bad_tools = [t for t, r in results.items() if not r["ok"]]

    lines = [
        f"[ENV_AGENT] 상태: {status}",
        f"  정상: {', '.join(ok_tools) if ok_tools else '없음'}",
    ]
    if bad_tools:
        lines.append(f"  문제: {', '.join(bad_tools)}")
    if critical_failures:
        lines.append(f"  치명 오류: {len(critical_failures)}건 → 파이프라인 중단")
    if warnings:
        lines.append(f"  경고: {len(warnings)}건 → 폴백 모드로 진행")
    return "\n".join(lines)


# ── tool_card 자동 갱신 ──────────────────────────────────

def update_tool_cards(env_result: dict) -> None:
    """
    env_agent 실행 결과를 각 tool_card yaml에 자동 기록.
    last_verified / last_run_result / last_run_error / run_count 갱신.
    """
    from datetime import date
    today = str(date.today())

    for tool_name, result in env_result["tool_results"].items():
        card_name = TOOL_CHECKS.get(tool_name, {}).get("card")
        if not card_name:
            continue

        card_path = TOOL_CARD_DIR / card_name
        if not card_path.exists():
            continue

        try:
            text = card_path.read_text(encoding="utf-8")

            # run_count 추출 및 증가
            run_count = 0
            for line in text.splitlines():
                if line.startswith("run_count:"):
                    try:
                        run_count = int(line.split(":")[1].strip())
                    except ValueError:
                        pass
                    break

            run_count += 1
            ok = result.get("ok", False)
            error_msg = result.get("issue", "null") if not ok else "null"

            # 갱신할 필드 처리
            new_lines = []
            fields_updated = {"last_verified": False, "last_run_result": False,
                              "last_run_error": False, "run_count": False}

            for line in text.splitlines():
                if line.startswith("last_verified:"):
                    new_lines.append(f'last_verified: "{today}"')
                    fields_updated["last_verified"] = True
                elif line.startswith("last_run_result:"):
                    new_lines.append(f'last_run_result: "{"OK" if ok else "FAIL"}"')
                    fields_updated["last_run_result"] = True
                elif line.startswith("last_run_error:"):
                    new_lines.append(f"last_run_error: {error_msg}")
                    fields_updated["last_run_error"] = True
                elif line.startswith("run_count:"):
                    new_lines.append(f"run_count: {run_count}")
                    fields_updated["run_count"] = True
                else:
                    new_lines.append(line)

            # 없는 필드는 파일 끝에 추가
            if not fields_updated["last_run_result"]:
                new_lines.append(f'last_run_result: "{"OK" if ok else "FAIL"}"')
            if not fields_updated["last_run_error"]:
                new_lines.append(f"last_run_error: {error_msg}")
            if not fields_updated["run_count"]:
                new_lines.append(f"run_count: {run_count}")

            card_path.write_text("\n".join(new_lines) + "\n", encoding="utf-8")

        except Exception:
            pass  # 카드 갱신 실패해도 파이프라인 중단 안 함


# ── 단독 실행 시 점검 결과 출력 ──────────────────────────

if __name__ == "__main__":
    import json
    result = run_env_check()
    print(result["summary"])
    print()
    if result["conflicts"]:
        print("충돌 감지:")
        for c in result["conflicts"]:
            print(f"  ⚠ {c}")
    if result["critical_failures"]:
        print("\n치명 오류:")
        for f in result["critical_failures"]:
            print(f"  ✗ {f}")
        sys.exit(1)
    if result["warnings"]:
        print("\n경고:")
        for w in result["warnings"]:
            print(f"  △ {w}")
    print(f"\n→ 진행 여부: {result['proceed']}")
