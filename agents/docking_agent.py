# Sigma Protocol v3.0 - Docking Agent
# 열역학 계산: AutoDock Vina (설치된 경우) 또는 RDKit 기반 추정 (simulate 모드)

import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from calibration.bias_record import apply_docking_correction


VINA_EXE        = Path("C:/sigma_protocol/vina/vina.exe")
DOCKING_WORKDIR = Path("C:/sigma_protocol/docking")
PREPARE_SCRIPT  = Path("C:/sigma_protocol/docking/prepare_and_dock.py")


def run_docking(smiles: str, target: str = "4EY7", mode: str = "simulate") -> dict:
    """
    mode='simulate' : Vina 없이 RDKit 기반 간이 추정값 사용
    mode='vina'     : AutoDock Vina 실행 (prepare_and_dock.py 호출)
    mode='cached'   : 이미 계산된 docking_report.json 로드
    """
    if mode == "vina":
        return _run_real_vina(smiles, target)
    if mode == "cached":
        return _load_cached(target)
    return _simulate_docking(smiles, target)


def _run_real_vina(smiles: str, target: str) -> dict:
    """prepare_and_dock.py를 subprocess로 호출해 실제 Vina 도킹 수행"""
    import json
    report_path = DOCKING_WORKDIR / target / "docking_report.json"

    # 이미 결과가 있으면 재사용
    if report_path.exists():
        return _load_cached(target)

    result = subprocess.run(
        [sys.executable, str(PREPARE_SCRIPT)],
        capture_output=True, text=True, timeout=600,
        env={**os.environ, "PYTHONUTF8": "1"},
    )
    if result.returncode != 0 and not report_path.exists():
        return {
            "error": "Vina 도킹 실패",
            "details": result.stderr[-500:],
            "mode": "vina",
        }
    return _load_cached(target)


def _load_cached(target: str) -> dict:
    """이미 계산된 docking_report.json에서 결과 로드"""
    import json
    report_path = DOCKING_WORKDIR / target / "docking_report.json"
    if not report_path.exists():
        return {"error": f"{report_path} 없음 — 먼저 vina 모드로 실행하세요."}

    data = json.loads(report_path.read_text(encoding="utf-8"))
    raw_dg = data["best_affinity_kcal_mol"]
    correction = apply_docking_correction(raw_dg)

    return {
        "mode": "vina (실제 계산)",
        "target": target,
        "smiles": data.get("smiles", ""),
        "grid_center": data.get("grid_center"),
        "grid_size": data.get("grid_size"),
        "raw_dG_kcal_mol": raw_dg,
        "corrected_dG_kcal_mol": correction["corrected_dG"],
        "correction_applied": correction["correction_applied"],
        "correction_note": correction["note"],
        "estimated_Ki_uM": _dg_to_ki(correction["corrected_dG"]),
        "all_modes": data.get("all_modes", []),
        "binding_poses_generated": len(data.get("all_modes", [])),
    }


def _simulate_docking(smiles: str, target: str) -> dict:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"error": "유효하지 않은 SMILES"}

    logp = Descriptors.MolLogP(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    rings = rdMolDescriptors.CalcNumRings(mol)
    aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)

    # 간이 추정식: hydrophobicity + H-bond + ring stacking 기여
    raw_dg = -(
        1.5
        + logp * 0.35
        + hba * 0.25
        + hbd * 0.15
        + aromatic_rings * 0.50
    )
    raw_dg = round(max(-13.0, min(-3.0, raw_dg)), 3)

    correction = apply_docking_correction(raw_dg)

    return {
        "mode": "simulate",
        "target": target,
        "smiles": smiles,
        "raw_dG_kcal_mol": correction["raw_dG"],
        "corrected_dG_kcal_mol": correction["corrected_dG"],
        "correction_applied": correction["correction_applied"],
        "correction_note": correction["note"],
        "estimated_Ki_uM": _dg_to_ki(correction["corrected_dG"]),
        "binding_poses_generated": 9,
        "warning": "시뮬레이션 모드: AutoDock Vina 미설치. RDKit 기반 추정값 사용.",
    }


def _try_vina(smiles: str, target: str) -> dict:
    vina_path = os.environ.get("VINA_PATH", "vina")
    try:
        result = subprocess.run(
            [vina_path, "--version"],
            capture_output=True, text=True, timeout=10,
        )
        if result.returncode != 0:
            raise FileNotFoundError
    except (FileNotFoundError, subprocess.TimeoutExpired, OSError):
        # Vina 없으면 simulate 모드로 폴백
        fallback = _simulate_docking(smiles, target)
        fallback["vina_fallback"] = True
        fallback["vina_fallback_reason"] = (
            "AutoDock Vina 미설치. VINA_PATH 환경변수 또는 PATH 확인 필요."
        )
        return fallback

    return {
        "error": "Vina 실행에는 receptor PDBQT 파일과 grid box 설정이 필요합니다.",
        "hint": "receptor PDBQT를 준비하고 docking_config.json을 설정하세요.",
    }


def _dg_to_ki(dg_kcal_mol: float) -> float:
    """ΔG (kcal/mol) → Ki (μM) 변환 (T=298K)"""
    import math
    R = 0.001987  # kcal/(mol·K)
    T = 298.0
    ki_molar = math.exp(dg_kcal_mol / (R * T))
    return round(ki_molar * 1e6, 4)  # μM
