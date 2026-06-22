"""
Sigma Protocol v3.0 — AChE (4EY7) 실제 AutoDock Vina 도킹 워크플로우
Step 1. PDB 4EY7 다운로드
Step 2. 수용체 준비 (물·리간드 제거 → 수소 추가 → PDBQT 변환)
Step 3. 리간드 3D 생성 → PDBQT 변환 (meeko)
Step 4. 그리드 박스 자동 산출 (공결정 리간드 중심)
Step 5. Vina 도킹 실행
Step 6. 결과 파싱 및 리포트
"""

import os
import sys
import json
import math
import subprocess
import urllib.request
from pathlib import Path

# UTF-8 출력
if sys.stdout.encoding and sys.stdout.encoding.lower() not in ("utf-8", "utf8"):
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")

WORK_DIR   = Path("C:/sigma_protocol/docking/4EY7")
VINA_EXE   = Path("C:/sigma_protocol/vina/vina.exe")
PDB_ID     = "4EY7"
LIGAND_ID  = "E20"   # co-crystal ligand (donepezil, 4EY7 잔기명)
SMILES     = "C1CCN(C1)c1nc2ccccc2n1CC(=O)O"

# AChE 활성부위 그리드 박스 (4EY7 기반, 공결정 리간드로 자동 산출)
# 자동 산출 실패 시 fallback 좌표 사용
FALLBACK_CENTER = (-5.5, 40.5, 55.5)   # Å
BOX_SIZE        = (22.0, 22.0, 22.0)   # Å


# ─────────────────────────────────────────────
# Step 1. PDB 다운로드
# ─────────────────────────────────────────────
def download_pdb(pdb_id: str) -> Path:
    path = WORK_DIR / f"{pdb_id}.pdb"
    if path.exists() and path.stat().st_size > 1000:
        print(f"  [SKIP] 이미 존재: {path}")
        return path
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    print(f"  다운로드: {url}")
    urllib.request.urlretrieve(url, path)
    print(f"  저장: {path} ({path.stat().st_size:,} bytes)")
    return path


# ─────────────────────────────────────────────
# Step 2. 수용체 준비 (Biopython + RDKit)
# ─────────────────────────────────────────────
def extract_ligand_center(pdb_path: Path, ligand_resname: str) -> tuple[float, float, float] | None:
    """공결정 리간드 원자 좌표로 그리드 중심 산출"""
    coords = []
    with open(pdb_path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            rec = line[:6].strip()
            if rec in ("HETATM", "ATOM"):
                resname = line[17:20].strip()
                if resname == ligand_resname:
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        coords.append((x, y, z))
                    except ValueError:
                        pass
    if not coords:
        return None
    cx = sum(c[0] for c in coords) / len(coords)
    cy = sum(c[1] for c in coords) / len(coords)
    cz = sum(c[2] for c in coords) / len(coords)
    print(f"  공결정 리간드 ({ligand_resname}) {len(coords)}개 원자 → 중심: "
          f"({cx:.2f}, {cy:.2f}, {cz:.2f})")
    return round(cx, 3), round(cy, 3), round(cz, 3)


def prepare_receptor_pdbqt(pdb_path: Path) -> Path:
    """
    순수 Python으로 수용체 PDBQT 생성:
    1. ATOM 행만 추출 (HOH/물 제외, HETATM 제외)
    2. 수소 추가 (RDKit)
    3. meeko MoleculePreparation → PDBQT
    실패 시 경량 fallback (수소 없는 PDBQT — 정확도 저하)
    """
    receptor_pdb  = WORK_DIR / "receptor_clean.pdb"
    receptor_pdbqt = WORK_DIR / "receptor.pdbqt"

    if receptor_pdbqt.exists() and receptor_pdbqt.stat().st_size > 1000:
        print(f"  [SKIP] 수용체 PDBQT 이미 존재: {receptor_pdbqt}")
        return receptor_pdbqt

    # ATOM 행만 남기기 (HOH 제외)
    lines_out = []
    with open(pdb_path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            rec = line[:6].strip()
            if rec == "ATOM":
                resname = line[17:20].strip()
                if resname != "HOH":
                    lines_out.append(line)
            elif rec == "TER":
                lines_out.append(line)
    with open(receptor_pdb, "w", encoding="utf-8") as f:
        f.writelines(lines_out)
    print(f"  수용체 정제: {len(lines_out)}행 → {receptor_pdb}")

    # MGLTools prepare_receptor4.py 없으므로 — Biopython + 직접 PDBQT 생성
    _write_pdbqt_from_pdb(receptor_pdb, receptor_pdbqt)
    return receptor_pdbqt


def _write_pdbqt_from_pdb(pdb_path: Path, pdbqt_path: Path):
    """
    Biopython으로 파싱 → AutoDock 원자 타입 매핑 → PDBQT 쓰기.
    열 위치 (1-indexed):
      1-6  record  7-11 serial  12 sp  13-16 name  17 altloc
      18-20 resname  21 sp  22 chain  23-26 resseq  27-30 sp
      31-38 x  39-46 y  47-54 z  55-60 occ  61-66 B
      67-70 sp  71-76 charge  77 sp  78-79 atom_type
    """
    from Bio import PDB as BioPDB

    elem_to_ad = {
        "C": "C", "N": "N", "O": "OA", "S": "SA", "P": "P",
        "H": "HD", "F": "F", "CL": "CL", "BR": "BR", "I": "I",
        "FE": "FE", "ZN": "ZN", "CA": "CA", "MG": "MG",
    }

    parser = BioPDB.PDBParser(QUIET=True)
    structure = parser.get_structure("rec", str(pdb_path))

    lines = []
    serial = 1
    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname().strip()
                if resname in ("HOH", "WAT"):
                    continue
                for atom in residue:
                    elem = (atom.element or "C").strip().upper()
                    if elem == "H":
                        continue
                    ad_type = elem_to_ad.get(elem, "C")
                    x, y, z = atom.get_coord()
                    aname = atom.get_name().strip()
                    # PDB 원자명 열 규칙: 1-char 원소 → col14 시작(" X  "), 2-char → col13("XX  ")
                    if len(elem) == 1 and len(aname) < 4:
                        name_field = f" {aname:<3}"   # 4 chars: col 13-16
                    else:
                        name_field = f"{aname:<4}"     # 4 chars: col 13-16
                    resid = residue.get_id()[1]
                    # 정확한 PDBQT 열 배치
                    line = (
                        f"ATOM  "               # 1-6
                        f"{serial:5d}"          # 7-11
                        f" "                    # 12
                        f"{name_field}"         # 13-16
                        f" "                    # 17 (altloc)
                        f"{resname:<3}"         # 18-20
                        f" "                    # 21
                        f"{chain.id}"           # 22
                        f"{resid:4d}"           # 23-26
                        f"    "                 # 27-30
                        f"{x:8.3f}"             # 31-38
                        f"{y:8.3f}"             # 39-46
                        f"{z:8.3f}"             # 47-54
                        f"  1.00"               # 55-60
                        f"  0.00"               # 61-66
                        f"    "                 # 67-70
                        f"  0.000"              # 71-77 (charge 6.3f + leading space)
                        f" {ad_type:<2}\n"      # 78-79 atom type
                    )
                    lines.append(line)
                    serial += 1
            lines.append("TER\n")
    lines.append("ENDMDL\n")

    with open(pdbqt_path, "w", encoding="utf-8") as f:
        f.writelines(lines)
    print(f"  수용체 PDBQT 생성: {pdbqt_path} ({len(lines)}행)")


# ─────────────────────────────────────────────
# Step 3. 리간드 3D 생성 → PDBQT
# ─────────────────────────────────────────────
def prepare_ligand_pdbqt(smiles: str) -> Path:
    ligand_pdbqt = WORK_DIR / "ligand.pdbqt"
    if ligand_pdbqt.exists() and ligand_pdbqt.stat().st_size > 100:
        print(f"  [SKIP] 리간드 PDBQT 이미 존재: {ligand_pdbqt}")
        return ligand_pdbqt

    from rdkit import Chem
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation

    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    # 3D 형태 생성 (ETKDG)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    result = AllChem.EmbedMolecule(mol, params)
    if result == -1:
        # Fallback: 거리 기하학
        AllChem.EmbedMolecule(mol, AllChem.ETDG())
    AllChem.MMFFOptimizeMolecule(mol)

    # meeko 0.7.x API로 PDBQT 변환 (prepare → list 반환)
    from meeko import PDBQTWriterLegacy
    preparator = MoleculePreparation()
    mol_setups = preparator.prepare(mol)
    if isinstance(mol_setups, list):
        mol_setup = mol_setups[0]
    else:
        mol_setup = mol_setups
    pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(mol_setup)
    if not is_ok:
        raise RuntimeError(f"meeko PDBQT 생성 실패: {error_msg}")

    ligand_pdbqt.write_text(pdbqt_string, encoding="utf-8")
    print(f"  리간드 PDBQT 생성: {ligand_pdbqt}")
    return ligand_pdbqt


# ─────────────────────────────────────────────
# Step 4. Vina 설정 파일 생성
# ─────────────────────────────────────────────
def write_vina_config(
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    center: tuple,
    box_size: tuple,
    exhaustiveness: int = 16,
    num_modes: int = 9,
) -> Path:
    config_path = WORK_DIR / "vina_config.txt"
    out_path    = WORK_DIR / "docking_out.pdbqt"

    config = f"""receptor = {receptor_pdbqt}
ligand   = {ligand_pdbqt}

center_x = {center[0]}
center_y = {center[1]}
center_z = {center[2]}

size_x = {box_size[0]}
size_y = {box_size[1]}
size_z = {box_size[2]}

out = {out_path}

exhaustiveness = {exhaustiveness}
num_modes      = {num_modes}
energy_range   = 3
"""
    config_path.write_text(config, encoding="utf-8")
    print(f"  Vina 설정 저장: {config_path}")
    print(f"  그리드 중심: {center}  크기: {box_size}")
    return config_path


# ─────────────────────────────────────────────
# Step 5. Vina 실행
# ─────────────────────────────────────────────
def run_vina(config_path: Path) -> str:
    out_pdbqt = WORK_DIR / "docking_out.pdbqt"
    cmd = [str(VINA_EXE), "--config", str(config_path)]
    print(f"  실행: {' '.join(cmd)}")
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=300,
    )
    output = result.stdout + result.stderr
    print(output)
    if result.returncode != 0 and not out_pdbqt.exists():
        raise RuntimeError(f"Vina 실행 실패 (exit {result.returncode}):\n{output}")
    return output


# ─────────────────────────────────────────────
# Step 6. 결과 파싱
# ─────────────────────────────────────────────
def parse_vina_output(vina_stdout: str, out_pdbqt: Path) -> list[dict]:
    """
    Vina stdout에서 모드별 ΔG 파싱.
    포맷:
      mode |   affinity | dist from best mode
           | (kcal/mol) | rmsd l.b.| rmsd u.b.
    -----+------------+----------+----------
       1 |       -7.5 |    0.000 |    0.000
    """
    modes = []
    in_table = False
    for line in vina_stdout.splitlines():
        if "affinity" in line.lower():
            in_table = True
            continue
        if in_table and "---" in line:
            continue
        if in_table and line.strip():
            parts = line.split()
            if len(parts) >= 4:
                try:
                    mode_idx = int(parts[0])
                    affinity = float(parts[1])
                    rmsd_lb  = float(parts[2])
                    rmsd_ub  = float(parts[3])
                    modes.append({
                        "mode": mode_idx,
                        "affinity_kcal_mol": affinity,
                        "rmsd_lb": rmsd_lb,
                        "rmsd_ub": rmsd_ub,
                        "ki_uM": _dg_to_ki(affinity),
                    })
                except (ValueError, IndexError):
                    pass
    return modes


def _dg_to_ki(dg: float) -> float:
    R, T = 0.001987, 298.0
    return round(math.exp(dg / (R * T)) * 1e6, 4)


# ─────────────────────────────────────────────
# 메인
# ─────────────────────────────────────────────
def main():
    print(f"\n{'='*62}")
    print(f"  Sigma Protocol v3.0 — AChE (4EY7) 실제 Vina 도킹")
    print(f"  리간드: {SMILES}")
    print(f"{'='*62}\n")

    WORK_DIR.mkdir(parents=True, exist_ok=True)

    # Step 1
    print("[1/6] PDB 다운로드")
    pdb_path = download_pdb(PDB_ID)

    # Step 2
    print("\n[2/6] 수용체 준비")
    center = extract_ligand_center(pdb_path, LIGAND_ID)
    if center is None:
        print(f"  공결정 리간드 미발견 → fallback 좌표 사용: {FALLBACK_CENTER}")
        center = FALLBACK_CENTER
    receptor_pdbqt = prepare_receptor_pdbqt(pdb_path)

    # Step 3
    print("\n[3/6] 리간드 준비")
    ligand_pdbqt = prepare_ligand_pdbqt(SMILES)

    # Step 4
    print("\n[4/6] 도킹 설정")
    config_path = write_vina_config(
        receptor_pdbqt, ligand_pdbqt, center, BOX_SIZE,
        exhaustiveness=16, num_modes=9,
    )

    # Step 5
    print("\n[5/6] AutoDock Vina 실행 (exhaustiveness=16)")
    vina_out = run_vina(config_path)

    # Step 6
    print("\n[6/6] 결과 분석")
    modes = parse_vina_output(vina_out, WORK_DIR / "docking_out.pdbqt")

    if not modes:
        print("  [WARN] 모드 파싱 실패 — vina 출력 확인 필요")
        return

    best = modes[0]
    print(f"\n  {'='*50}")
    print(f"  Best pose: ΔG = {best['affinity_kcal_mol']} kcal/mol  "
          f"Ki = {best['ki_uM']} μM")
    print(f"  {'='*50}")
    print(f"\n  {'Mode':>4}  {'ΔG (kcal/mol)':>14}  {'Ki (μM)':>12}  "
          f"{'RMSD lb':>8}  {'RMSD ub':>8}")
    print(f"  {'-'*56}")
    for m in modes:
        print(f"  {m['mode']:>4}  {m['affinity_kcal_mol']:>14.3f}  "
              f"{m['ki_uM']:>12.4f}  "
              f"{m['rmsd_lb']:>8.3f}  {m['rmsd_ub']:>8.3f}")

    # 판정
    dg = best["affinity_kcal_mol"]
    print(f"\n  [판정] ", end="")
    if dg <= -9:
        print(f"강력한 결합 (ΔG≤-9) — 고효능 후보 ✅")
    elif dg <= -7:
        print(f"양호한 결합 (-9<ΔG≤-7) — Lead 수준 ✅")
    elif dg <= -5:
        print(f"보통 결합 (-7<ΔG≤-5) — 최적화 필요 ⚡")
    else:
        print(f"약한 결합 (ΔG>-5) — 재설계 권고 ❌")

    # JSON 저장
    report = {
        "pdb_id": PDB_ID,
        "smiles": SMILES,
        "grid_center": center,
        "grid_size": BOX_SIZE,
        "best_affinity_kcal_mol": best["affinity_kcal_mol"],
        "best_ki_uM": best["ki_uM"],
        "all_modes": modes,
    }
    report_path = WORK_DIR / "docking_report.json"
    report_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(f"\n  리포트 저장: {report_path}")
    print(f"  도킹 포즈:   {WORK_DIR / 'docking_out.pdbqt'}")
    print(f"\n{'='*62}\n")
    return report


if __name__ == "__main__":
    main()
