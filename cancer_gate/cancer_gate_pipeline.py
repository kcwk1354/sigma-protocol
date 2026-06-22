"""
Cancer-Gate Pipeline v1.0
타깃: GRM4 (mGluR4) | ADC 링커 트리거 후보물질 발굴
PDB:  7E9H (Gi-bound mGlu4 cryo-EM, TM 알로스테릭 포켓)

Epsilon-Gate 교훈 반영:
  - DILI < 0.4  (간독성 사전 차단)
  - hERG < 0.3  (IC50 > 30μM 근사)
  - Ki 목표 < 1μM  (ΔG ≤ -8.2 kcal/mol)
"""

import sys, os, json, math, subprocess, urllib.request
from pathlib import Path
from dataclasses import dataclass, field, asdict
from typing import Optional

if sys.stdout.encoding and sys.stdout.encoding.lower() not in ("utf-8", "utf8"):
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")

# ─────────────────────────────────────────────
# 0. 설정
# ─────────────────────────────────────────────
BASE_DIR   = Path("C:/sigma_protocol/cancer_gate")
DOCK_DIR   = BASE_DIR / "docking" / "7E9H"
VINA_EXE   = Path("C:/sigma_protocol/vina/vina.exe")
PDB_ID     = "7E9H"

# ADC 링커 트리거 필터 (Epsilon-Gate 교훈)
DILI_THRESHOLD  = 0.40   # < 0.4  PASS
HERG_THRESHOLD  = 0.30   # < 0.3  PASS (IC50 > ~30μM 근사)
KI_TARGET_UM    = 1.0    # < 1μM  목표
DG_TARGET       = -8.2   # kcal/mol (Ki 1μM 근사)

# Gen0 후보 5종 (GRM4 PAM 스캐폴드 + ADC 컨쥬게이션 핸들)
GEN0_CANDIDATES = [
    {
        "id": "CG-001",
        "smiles": "OC(=O)C1(CC1)c1ccc(OC(F)(F)F)cc1",
        "name": "1-(4-trifluoromethoxyphenyl)cyclopropane-1-COOH",
        "scaffold": "Cyclopropane (VU0155041형 PAM)",
        "adc_handle": "COOH",
        "rationale": "GRM4 TM 알로스테릭 포켓의 소수성 영역 + cyclopropane torsional rigidity",
    },
    {
        "id": "CG-002",
        "smiles": "OC(=O)c1ccc(C(=O)N2CCC(F)(F)CC2)cc1",
        "name": "4-(4,4-difluoropiperidine-1-carbonyl)benzoic acid",
        "scaffold": "Difluoropiperidine-benzoyl",
        "adc_handle": "COOH",
        "rationale": "difluoropiperidine: metabolic stability + sp3 풍부화 / benzoic acid linker 핸들",
    },
    {
        "id": "CG-003",
        "smiles": "OC(=O)c1ccc2cc(NC(=O)c3ccc(F)cc3)ccc2c1",
        "name": "6-(4-fluorobenzamido)naphthalene-2-carboxylic acid",
        "scaffold": "Naphthalene biaryl amide (ADX88178형)",
        "adc_handle": "COOH",
        "rationale": "ADX88178 core scaffold: GRM4 PAM 검증 scaffold + fluorophenyl 선택성",
    },
    {
        "id": "CG-004",
        "smiles": "Nc1ccc(CC2(C(=O)Nc3ccc(F)cc3)CC2)cc1",
        "name": "1-(4-aminobenzyl)-N-(4-fluorophenyl)cyclopropane-1-carboxamide",
        "scaffold": "Cyclopropane + 아닐린 (NH2 ADC 핸들)",
        "adc_handle": "NH2",
        "rationale": "Primary amine → mAb lysine 컨쥬게이션 / cyclopropane PAM scaffold 유지",
    },
    {
        "id": "CG-005",
        "smiles": "OC(=O)c1ccc(-c2nc3cc(F)ccc3[nH]2)cc1",
        "name": "4-(5-fluorobenzimidazol-2-yl)benzoic acid",
        "scaffold": "Benzimidazole-benzoate",
        "adc_handle": "COOH",
        "rationale": "벤즈이미다졸 NH: H-bond donor to GRM4 Glu pocket / COOH linker 핸들",
    },
]

# Gen1 후보 5종 (Gen0 블랙리스트 적용: 나프탈렌/아닐린/CF3O 제거, 방향족 max 1개, Fsp3 > 0.4)
GEN1_CANDIDATES = [
    {
        "id": "CG1-001",
        "smiles": "OC(=O)C1CCN(Cc2ccc(F)cc2)CC1",
        "name": "1-(4-fluorobenzyl)piperidine-4-carboxylic acid",
        "scaffold": "Piperidine-4-COOH + 4-F-benzyl",
        "adc_handle": "COOH",
        "rationale": "VU0155041형 PAM pharmacophore: F-benzyl TM3/5 소수성 포켓 / piperidine sp3 풍부화",
    },
    {
        "id": "CG1-002",
        "smiles": "OC(=O)CN1CCC(Cc2ccc(F)cc2)CC1",
        "name": "2-(4-(4-fluorobenzyl)piperidin-1-yl)acetic acid",
        "scaffold": "N-glycinate-piperidine + 4-F-benzyl",
        "adc_handle": "COOH",
        "rationale": "CG1-001 대비 N-acetic acid 링커: 유연한 COOH ADC 핸들 / 같은 F-benzyl 약물단",
    },
    {
        "id": "CG1-003",
        "smiles": "OC(=O)C1CN(Cc2ccc(F)cc2)C1",
        "name": "1-(4-fluorobenzyl)azetidine-3-carboxylic acid",
        "scaffold": "Azetidine-3-COOH + 4-F-benzyl",
        "adc_handle": "COOH",
        "rationale": "4원환 축소 링: TM 포켓 tight fit 탐색 / Fsp3 최대화",
    },
    {
        "id": "CG1-004",
        "smiles": "OC(=O)C1CCC(Cc2ccc(F)cc2)CC1",
        "name": "4-(4-fluorobenzyl)cyclohexane-1-carboxylic acid",
        "scaffold": "Cyclohexane-COOH + 4-F-benzyl (N-free)",
        "adc_handle": "COOH",
        "rationale": "질소 제거 → hERG 최소화 / 순수 탄화수소 scaffold = 낮은 DILI 기대",
    },
    {
        "id": "CG1-005",
        "smiles": "OC(=O)C1CCN(Cc2ccncc2)CC1",
        "name": "1-(pyridin-4-ylmethyl)piperidine-4-carboxylic acid",
        "scaffold": "Piperidine-4-COOH + 4-pyridylmethyl",
        "adc_handle": "COOH",
        "rationale": "F-phenyl -> pyridine: N H-bond to TM Ser / 더 낮은 DILI 예상",
    },
]


# ─────────────────────────────────────────────
# 1. 수용체 준비
# ─────────────────────────────────────────────
def download_and_prepare_receptor() -> tuple[Path, tuple]:
    pdb_path = DOCK_DIR / f"{PDB_ID}.pdb"

    # PDB 다운로드
    if not pdb_path.exists() or pdb_path.stat().st_size < 1000:
        url = f"https://files.rcsb.org/download/{PDB_ID}.pdb"
        print(f"  다운로드: {url}")
        urllib.request.urlretrieve(url, pdb_path)
        print(f"  저장: {pdb_path} ({pdb_path.stat().st_size:,} bytes)")

    # 체인 A (mGlu4)만 추출 + 물 제거
    receptor_pdb   = DOCK_DIR / "receptor_R.pdb"
    receptor_pdbqt = DOCK_DIR / "receptor_R.pdbqt"

    lines_R = []
    hetatm_resnames = set()
    with open(pdb_path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            rec = line[:6].strip()
            if rec == "ATOM" and line[21] == "R":   # Chain R = mGlu4 receptor
                if line[17:20].strip() != "HOH":
                    lines_R.append(line)
            elif rec == "HETATM" and line[21] == "R":
                rn = line[17:20].strip()
                hetatm_resnames.add(rn)
                if rn not in ("HOH", "WAT"):
                    lines_R.append(line)
    with open(receptor_pdb, "w", encoding="utf-8") as f:
        f.writelines(lines_R)
    print(f"  체인 R 추출: {len(lines_R)}행 | HETATM 잔기: {hetatm_resnames}")

    # 리간드 중심으로 그리드 박스 산출
    center = _find_grid_center(pdb_path, hetatm_resnames - {"HOH", "WAT"})

    # PDBQT 변환 (캐시)
    if not receptor_pdbqt.exists() or receptor_pdbqt.stat().st_size < 1000:
        _write_pdbqt(receptor_pdb, receptor_pdbqt)
    else:
        print(f"  [SKIP] 수용체 PDBQT 이미 존재")

    return receptor_pdbqt, center


def _find_grid_center(pdb_path: Path, lig_resnames: set) -> tuple:
    """공결정 소분자 리간드로 그리드 중심 산출. 없으면 7E9H TM3-7 Calpha centroid fallback.
    SEP/TPO 같은 수정 아미노산은 리간드로 취급하지 않음."""
    _MODIFIED_AA = {"SEP", "TPO", "MSE", "HYP", "MLZ", "PTR", "NEP", "CGU"}
    small_mol_resnames = lig_resnames - _MODIFIED_AA
    coords = []
    if small_mol_resnames:
        target_resname = sorted(small_mol_resnames)[0]
        with open(pdb_path, encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line[:6].strip() == "HETATM" and line[17:20].strip() == target_resname:
                    try:
                        coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                    except ValueError:
                        pass
    if coords:
        cx = round(sum(c[0] for c in coords) / len(coords), 3)
        cy = round(sum(c[1] for c in coords) / len(coords), 3)
        cz = round(sum(c[2] for c in coords) / len(coords), 3)
        print(f"  그리드 중심 (소분자 리간드 기반): ({cx}, {cy}, {cz})")
        return cx, cy, cz
    # 7E9H Chain R TM3-7 Calpha centroid (GRM4 PAM 알로스테릭 포켓)
    fallback = (211.697, 188.139, 202.871)
    print(f"  그리드 중심 (7E9H ChainR TM3-7 centroid): {fallback}")
    return fallback


def _write_pdbqt(pdb_path: Path, pdbqt_path: Path):
    from Bio import PDB as BioPDB
    elem_to_ad = {
        "C": "C", "N": "N", "O": "OA", "S": "SA", "P": "P",
        "F": "F", "CL": "CL", "BR": "BR", "I": "I",
        "FE": "FE", "ZN": "ZN", "CA": "CA", "MG": "MG",
    }
    parser = BioPDB.PDBParser(QUIET=True)
    structure = parser.get_structure("rec", str(pdb_path))
    lines, serial = [], 1
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname().strip() in ("HOH", "WAT"):
                    continue
                for atom in residue:
                    elem = (atom.element or "C").strip().upper()
                    if elem == "H":
                        continue
                    ad_type = elem_to_ad.get(elem, "C")
                    x, y, z = atom.get_coord()
                    aname = atom.get_name().strip()
                    name_field = f" {aname:<3}" if len(elem) == 1 and len(aname) < 4 else f"{aname:<4}"
                    resid = residue.get_id()[1]
                    line = (
                        f"ATOM  {serial:5d} {name_field} "
                        f"{residue.get_resname():<3} {chain.id}{resid:4d}    "
                        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00    "
                        f"  0.000 {ad_type:<2}\n"
                    )
                    lines.append(line)
                    serial += 1
            lines.append("TER\n")
    lines.append("ENDMDL\n")
    pdbqt_path.write_text("".join(lines), encoding="utf-8")
    print(f"  수용체 PDBQT: {pdbqt_path} ({len(lines)}행)")


# ─────────────────────────────────────────────
# 2. 리간드 준비
# ─────────────────────────────────────────────
def prepare_ligand(cand: dict) -> Optional[Path]:
    cid = cand["id"]
    smiles = cand["smiles"]
    path = DOCK_DIR / f"ligand_{cid}.pdbqt"
    if path.exists():
        return path
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from meeko import MoleculePreparation, PDBQTWriterLegacy

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"  [{cid}] SMILES 파싱 실패")
            return None
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        if AllChem.EmbedMolecule(mol, params) == -1:
            AllChem.EmbedMolecule(mol, AllChem.ETDG())
        AllChem.MMFFOptimizeMolecule(mol)
        preparator = MoleculePreparation()
        setups = preparator.prepare(mol)
        setup = setups[0] if isinstance(setups, list) else setups
        pdbqt_str, ok, err = PDBQTWriterLegacy.write_string(setup)
        if not ok:
            print(f"  [{cid}] meeko 실패: {err}")
            return None
        path.write_text(pdbqt_str, encoding="utf-8")
        return path
    except Exception as e:
        print(f"  [{cid}] 리간드 준비 오류: {e}")
        return None


# ─────────────────────────────────────────────
# 3. Vina 도킹
# ─────────────────────────────────────────────
def run_vina_for(cand: dict, receptor_pdbqt: Path, center: tuple) -> dict:
    cid = cand["id"]
    ligand_pdbqt = prepare_ligand(cand)
    if ligand_pdbqt is None:
        return {"error": "리간드 준비 실패", "id": cid}

    out_pdbqt  = DOCK_DIR / f"out_{cid}.pdbqt"
    config_path = DOCK_DIR / f"config_{cid}.txt"
    gen = cand["id"][:3]
    gen_num = int(cand["id"][2:].split('-')[0]) if cand["id"][2:].split('-')[0].isdigit() else 0
    box_sz = 30.0 if gen_num >= 4 else 25.0
    box = (box_sz, box_sz, box_sz)
    exhaustiveness = cand.get("exhaustiveness", 48 if gen_num >= 11 else (32 if gen_num >= 3 else 16))

    config = (
        f"receptor = {receptor_pdbqt}\n"
        f"ligand   = {ligand_pdbqt}\n"
        f"center_x = {center[0]}\ncenter_y = {center[1]}\ncenter_z = {center[2]}\n"
        f"size_x = {box[0]}\nsize_y = {box[1]}\nsize_z = {box[2]}\n"
        f"out = {out_pdbqt}\nexhaustiveness = {exhaustiveness}\nnum_modes = 9\nenergy_range = 3\n"
    )
    config_path.write_text(config, encoding="utf-8")

    result = subprocess.run(
        [str(VINA_EXE), "--config", str(config_path)],
        capture_output=True, text=True, timeout=300,
    )
    output = result.stdout + result.stderr

    modes = _parse_vina(output)
    if not modes:
        return {"error": "도킹 결과 없음", "id": cid, "vina_log": output[-300:]}

    best_dg = modes[0]["affinity_kcal_mol"]
    best_ki = _dg_to_ki(best_dg)
    return {
        "id": cid,
        "best_dG": best_dg,
        "best_Ki_uM": best_ki,
        "all_modes": modes,
        "docking_pass": best_ki < KI_TARGET_UM,
    }


def _parse_vina(stdout: str) -> list:
    modes, in_table = [], False
    for line in stdout.splitlines():
        if "affinity" in line.lower():
            in_table = True
            continue
        if in_table and "---" in line:
            continue
        if in_table and line.strip():
            parts = line.split()
            if len(parts) >= 4:
                try:
                    modes.append({
                        "mode": int(parts[0]),
                        "affinity_kcal_mol": float(parts[1]),
                        "rmsd_lb": float(parts[2]),
                        "rmsd_ub": float(parts[3]),
                    })
                except (ValueError, IndexError):
                    pass
    return modes


def _dg_to_ki(dg: float) -> float:
    return round(math.exp(dg / (0.001987 * 298.0)) * 1e6, 4)


# ─────────────────────────────────────────────
# 4. RDKit 계산
# ─────────────────────────────────────────────
def calc_rdkit(smiles: str) -> dict:
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors
        from rdkit.Chem import QED as QEDmod
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"error": "파싱 실패"}
        cooh = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
        nh2  = Chem.MolFromSmarts("[NH2]")
        return {
            "mw":    round(Descriptors.MolWt(mol), 2),
            "logp":  round(Descriptors.MolLogP(mol), 3),
            "qed":   round(QEDmod.qed(mol), 4),
            "hbd":   rdMolDescriptors.CalcNumHBD(mol),
            "hba":   rdMolDescriptors.CalcNumHBA(mol),
            "tpsa":  round(Descriptors.TPSA(mol), 2),
            "rotb":  rdMolDescriptors.CalcNumRotatableBonds(mol),
            "ro5":   sum([Descriptors.MolWt(mol) > 500,
                          Descriptors.MolLogP(mol) > 5,
                          rdMolDescriptors.CalcNumHBD(mol) > 5,
                          rdMolDescriptors.CalcNumHBA(mol) > 10]),
            "lipinski_pass": sum([Descriptors.MolWt(mol) > 500,
                                  Descriptors.MolLogP(mol) > 5,
                                  rdMolDescriptors.CalcNumHBD(mol) > 5,
                                  rdMolDescriptors.CalcNumHBA(mol) > 10]) == 0,
            "has_cooh": mol.HasSubstructMatch(cooh),
            "has_nh2":  mol.HasSubstructMatch(nh2),
            "qed_warn": 0.48 <= QEDmod.qed(mol) <= 0.52,
        }
    except Exception as e:
        return {"error": str(e)}


# ─────────────────────────────────────────────
# 5. ADMET 필터 (Cancer-Gate 기준)
# ─────────────────────────────────────────────
_admet_model = None

def get_admet_model():
    global _admet_model
    if _admet_model is None:
        from admet_ai import ADMETModel
        _admet_model = ADMETModel()
    return _admet_model


def calc_admet_and_filter(smiles: str, has_cooh: bool) -> dict:
    try:
        model = get_admet_model()
        preds = model.predict(smiles=smiles)
        raw = preds.iloc[0].to_dict() if hasattr(preds, "iloc") else preds

        def g(k, alts=None):
            v = raw.get(k)
            if v is None and alts:
                for a in alts:
                    v = raw.get(a)
                    if v is not None:
                        break
            return round(float(v), 4) if v is not None else None

        herg = g("hERG", ["hERG_at_10uM"])
        ames = g("AMES")
        dili = g("DILI")
        hia  = g("HIA_Hou")
        bbb  = g("BBB_Martins")
        cyp  = g("CYP3A4_Substrate_CarbonMangels")

        # COOH 보정 (F% +20%p)
        hia_corrected = min(round(hia + 0.20, 4), 1.0) if (has_cooh and hia is not None) else hia
        cooh_cal_applied = has_cooh and hia is not None

        # Cancer-Gate 필터 적용
        dili_pass = dili is None or dili < DILI_THRESHOLD
        herg_pass = herg is None or herg < HERG_THRESHOLD
        gate_pass = dili_pass and herg_pass

        return {
            "herg": herg,  "ames": ames,  "dili": dili,
            "hia": hia,    "hia_corrected": hia_corrected,
            "bbb": bbb,    "cyp3a4": cyp,
            "cooh_correction_applied": cooh_cal_applied,
            # Cancer-Gate 필터
            "dili_pass":  dili_pass,
            "herg_pass":  herg_pass,
            "gate_pass":  gate_pass,
            "filter_detail": {
                f"DILI {dili:.3f} {'<' if dili_pass else '>='} {DILI_THRESHOLD}" if dili is not None else "DILI N/A": dili_pass,
                f"hERG {herg:.3f} {'<' if herg_pass else '>='} {HERG_THRESHOLD}" if herg is not None else "hERG N/A": herg_pass,
            },
        }
    except Exception as e:
        return {"error": str(e), "gate_pass": False}


# ─────────────────────────────────────────────
# 6. 최종 점수 (Cancer-Gate 기준)
# ─────────────────────────────────────────────
def score_candidate(rdkit_r: dict, admet_r: dict, dock_r: dict) -> dict:
    score = 60.0
    flags = []

    # RDKit (20점)
    if rdkit_r.get("lipinski_pass"):
        score += 10
    qed = rdkit_r.get("qed", 0)
    score += qed * 10

    # Docking (35점) — Ki < 1μM 목표
    dg = dock_r.get("best_dG", -5.0)
    ki = dock_r.get("best_Ki_uM", 999)
    if ki < 0.1:
        score += 35
    elif ki < 1.0:
        score += 28
    elif ki < 10:
        score += 18
    elif ki < 100:
        score += 10
    else:
        score += 3
    if ki >= KI_TARGET_UM:
        flags.append(f"Ki {ki:.1f}μM ≥ 1μM 목표 미달")

    # ADMET (25점) + Cancer-Gate 필터
    if not admet_r.get("dili_pass", True):
        score -= 20
        flags.append(f"DILI {admet_r.get('dili', 'N/A')} ≥ 0.4 FAIL")
    if not admet_r.get("herg_pass", True):
        score -= 15
        flags.append(f"hERG {admet_r.get('herg', 'N/A')} ≥ 0.3 FAIL")
    if admet_r.get("ames", 0) and admet_r["ames"] > 0.5:
        score -= 10
        flags.append(f"AMES {admet_r['ames']:.3f}")
    hia = admet_r.get("hia_corrected") or admet_r.get("hia", 0)
    if hia and hia >= 0.7:
        score += 5
    if admet_r.get("bbb", 0) and admet_r["bbb"] >= 0.5:
        score += 3

    score = max(0, min(100, round(score, 1)))

    # Cancer-Gate 최종 판정
    gate_pass  = admet_r.get("gate_pass", False)
    ki_pass    = ki < KI_TARGET_UM
    if gate_pass and ki_pass and score >= 70:
        verdict = "✅ GO"
    elif gate_pass and (ki_pass or score >= 55):
        verdict = "⚡ CONDITIONAL"
    elif not gate_pass:
        verdict = "❌ FILTERED"
    else:
        verdict = "❌ NO-GO"

    return {
        "score": score,
        "verdict": verdict,
        "gate_pass": gate_pass,
        "ki_pass": ki_pass,
        "flags": flags,
    }


# ─────────────────────────────────────────────
# 7. 메인 파이프라인
# ─────────────────────────────────────────────
def run_cancer_gate():
    print(f"\n{'='*65}")
    print(f"  Cancer-Gate Pipeline v1.0")
    print(f"  타깃: GRM4 (mGluR4)  |  PDB: {PDB_ID}")
    print(f"  필터: DILI < {DILI_THRESHOLD}  |  hERG < {HERG_THRESHOLD}  |  Ki < {KI_TARGET_UM}μM")
    print(f"{'='*65}\n")

    DOCK_DIR.mkdir(parents=True, exist_ok=True)

    # 수용체 준비
    print("[SETUP] 수용체 준비 (7E9H 체인 R = mGlu4)")
    receptor_pdbqt, center = download_and_prepare_receptor()

    results = []
    for i, cand in enumerate(GEN0_CANDIDATES, 1):
        cid = cand["id"]
        smiles = cand["smiles"]
        print(f"\n{'─'*65}")
        print(f"  [{i}/5] {cid} | {cand['name']}")
        print(f"  SMILES : {smiles}")
        print(f"  Scaffold: {cand['scaffold']} | 핸들: {cand['adc_handle']}")
        print(f"{'─'*65}")

        # RDKit
        rdkit_r = calc_rdkit(smiles)
        if "error" in rdkit_r:
            print(f"  [RDKit 오류] {rdkit_r['error']}")
            continue
        print(f"  RDKit  MW={rdkit_r['mw']}  LogP={rdkit_r['logp']}  "
              f"QED={rdkit_r['qed']}  Ro5={'PASS' if rdkit_r['lipinski_pass'] else 'FAIL'}  "
              f"COOH={'Y' if rdkit_r['has_cooh'] else 'N'}  NH2={'Y' if rdkit_r['has_nh2'] else 'N'}")

        # ADMET + Cancer-Gate 필터
        print(f"  ADMET  계산 중...", end=" ", flush=True)
        admet_r = calc_admet_and_filter(smiles, rdkit_r["has_cooh"])
        if "error" in admet_r:
            print(f"오류: {admet_r['error']}")
        else:
            cal_note = "(COOH +20%p 보정)" if admet_r.get("cooh_correction_applied") else ""
            print(f"완료 {cal_note}")
            print(f"  ADMET  hERG={admet_r.get('herg','N/A')}  "
                  f"DILI={admet_r.get('dili','N/A')}  "
                  f"AMES={admet_r.get('ames','N/A')}  "
                  f"HIA={admet_r.get('hia_corrected', admet_r.get('hia','N/A'))}  "
                  f"BBB={admet_r.get('bbb','N/A')}")
            gate = "PASS ✓" if admet_r.get("gate_pass") else "FAIL ✗"
            print(f"  Cancer-Gate 필터: DILI={'PASS' if admet_r.get('dili_pass') else 'FAIL'}  "
                  f"hERG={'PASS' if admet_r.get('herg_pass') else 'FAIL'}  → {gate}")

        # Docking
        _exh = 32 if int(cand["id"][2:].split('-')[0]) >= 3 else 16
        print(f"  Docking 실행 중 (exhaustiveness={_exh})...", end=" ", flush=True)
        dock_r = run_vina_for(cand, receptor_pdbqt, center)
        if "error" in dock_r:
            print(f"오류: {dock_r['error']}")
            dock_r = {"best_dG": -5.0, "best_Ki_uM": 999.9, "docking_pass": False}
        else:
            ki_note = "✅ Ki 목표 달성" if dock_r["docking_pass"] else "⚠ Ki 목표 미달"
            print(f"완료")
            print(f"  Docking ΔG={dock_r['best_dG']} kcal/mol  "
                  f"Ki={dock_r['best_Ki_uM']}μM  {ki_note}")

        # 점수 산출
        scored = score_candidate(rdkit_r, admet_r, dock_r)
        print(f"  종합   점수={scored['score']}/100  판정: {scored['verdict']}")
        if scored["flags"]:
            for fl in scored["flags"]:
                print(f"         ⚠ {fl}")

        results.append({
            "candidate": cand,
            "rdkit":  rdkit_r,
            "admet":  admet_r,
            "docking": dock_r,
            "score":  scored,
        })

    _print_ranking(results, gen="Gen0")

    # JSON 저장
    report_path = BASE_DIR / "cancer_gate_report.json"
    report_path.write_text(
        json.dumps(results, ensure_ascii=False, indent=2, default=str),
        encoding="utf-8"
    )
    print(f"\n  리포트: {report_path}\n")
    return results


def run_cancer_gate_gen1():
    """Gen1 후보 5종 스크리닝 (chain R 수용체 + 수정된 그리드 박스)"""
    print(f"\n{'='*65}")
    print(f"  Cancer-Gate Pipeline v1.0  —  Gen1 스크리닝")
    print(f"  타깃: GRM4 (mGluR4)  |  PDB: {PDB_ID} Chain R")
    print(f"  필터: DILI < {DILI_THRESHOLD}  |  hERG < {HERG_THRESHOLD}  |  Ki < {KI_TARGET_UM}uM")
    print(f"  그리드: TM3-7 Calpha centroid (211.7, 188.1, 202.9) sz=25A")
    print(f"{'='*65}\n")

    DOCK_DIR.mkdir(parents=True, exist_ok=True)

    print("[SETUP] 수용체 준비 (7E9H 체인 R = mGlu4)")
    receptor_pdbqt, center = download_and_prepare_receptor()

    results = []
    for i, cand in enumerate(GEN1_CANDIDATES, 1):
        cid = cand["id"]
        smiles = cand["smiles"]
        print(f"\n{'─'*65}")
        print(f"  [{i}/5] {cid} | {cand['name']}")
        print(f"  SMILES : {smiles}")
        print(f"  Scaffold: {cand['scaffold']} | 핸들: {cand['adc_handle']}")
        print(f"{'─'*65}")

        rdkit_r = calc_rdkit(smiles)
        if "error" in rdkit_r:
            print(f"  [RDKit 오류] {rdkit_r['error']}")
            continue
        print(f"  RDKit  MW={rdkit_r['mw']}  LogP={rdkit_r['logp']}  "
              f"QED={rdkit_r['qed']}  Ro5={'PASS' if rdkit_r['lipinski_pass'] else 'FAIL'}  "
              f"COOH={'Y' if rdkit_r['has_cooh'] else 'N'}")

        print(f"  ADMET  계산 중...", end=" ", flush=True)
        admet_r = calc_admet_and_filter(smiles, rdkit_r["has_cooh"])
        if "error" in admet_r:
            print(f"오류: {admet_r['error']}")
        else:
            cal_note = "(COOH +20%p 보정)" if admet_r.get("cooh_correction_applied") else ""
            print(f"완료 {cal_note}")
            print(f"  ADMET  hERG={admet_r.get('herg','N/A')}  "
                  f"DILI={admet_r.get('dili','N/A')}  "
                  f"AMES={admet_r.get('ames','N/A')}  "
                  f"HIA={admet_r.get('hia_corrected', admet_r.get('hia','N/A'))}  "
                  f"BBB={admet_r.get('bbb','N/A')}")
            gate = "PASS ✓" if admet_r.get("gate_pass") else "FAIL ✗"
            print(f"  Cancer-Gate: DILI={'PASS' if admet_r.get('dili_pass') else 'FAIL'}  "
                  f"hERG={'PASS' if admet_r.get('herg_pass') else 'FAIL'}  -> {gate}")

        _exh = 32 if int(cand["id"][2:].split('-')[0]) >= 3 else 16
        print(f"  Docking 실행 중 (exhaustiveness={_exh})...", end=" ", flush=True)
        dock_r = run_vina_for(cand, receptor_pdbqt, center)
        if "error" in dock_r:
            print(f"오류: {dock_r['error']}")
            dock_r = {"best_dG": -5.0, "best_Ki_uM": 999.9, "docking_pass": False}
        else:
            ki_note = "Ki 목표 달성" if dock_r["docking_pass"] else "Ki 목표 미달"
            print(f"완료")
            print(f"  Docking dG={dock_r['best_dG']} kcal/mol  "
                  f"Ki={dock_r['best_Ki_uM']}uM  {ki_note}")

        scored = score_candidate(rdkit_r, admet_r, dock_r)
        print(f"  종합   점수={scored['score']}/100  판정: {scored['verdict']}")
        if scored["flags"]:
            for fl in scored["flags"]:
                print(f"         * {fl}")

        results.append({
            "candidate": cand,
            "rdkit":  rdkit_r,
            "admet":  admet_r,
            "docking": dock_r,
            "score":  scored,
        })

    _print_ranking(results, gen="Gen1")

    report_path = BASE_DIR / "cancer_gate_gen1_report.json"
    report_path.write_text(
        json.dumps(results, ensure_ascii=False, indent=2, default=str),
        encoding="utf-8"
    )
    print(f"\n  리포트: {report_path}\n")
    return results


def _print_ranking(results: list, gen: str = "Gen0"):
    print(f"\n\n{'='*65}")
    print(f"  Cancer-Gate {gen} 스크리닝 결과 랭킹")
    print(f"  필터 기준: DILI < {DILI_THRESHOLD}  hERG < {HERG_THRESHOLD}  Ki < {KI_TARGET_UM}μM")
    print(f"{'='*65}")
    print(f"  {'순위':<4} {'ID':<8} {'점수':>5} {'판정':<14} "
          f"{'ΔG':>8} {'Ki(μM)':>9} {'DILI':>6} {'hERG':>6} {'QED':>6}")
    print(f"  {'-'*70}")

    ranked = sorted(results, key=lambda r: r["score"]["score"], reverse=True)
    for rank, r in enumerate(ranked, 1):
        cid   = r["candidate"]["id"]
        sc    = r["score"]["score"]
        verd  = r["score"]["verdict"]
        dg    = r["docking"].get("best_dG", "N/A")
        ki    = r["docking"].get("best_Ki_uM", "N/A")
        dili  = r["admet"].get("dili", "N/A")
        herg  = r["admet"].get("herg", "N/A")
        qed   = r["rdkit"].get("qed", "N/A")
        dg_s  = f"{dg:.3f}" if isinstance(dg, float) else str(dg)
        ki_s  = f"{ki:.4f}" if isinstance(ki, float) else str(ki)
        dili_s = f"{dili:.3f}" if isinstance(dili, float) else str(dili)
        herg_s = f"{herg:.3f}" if isinstance(herg, float) else str(herg)
        qed_s  = f"{qed:.4f}" if isinstance(qed, float) else str(qed)
        print(f"  {rank:<4} {cid:<8} {sc:>5} {verd:<14} "
              f"{dg_s:>8} {ki_s:>9} {dili_s:>6} {herg_s:>6} {qed_s:>6}")

    print(f"\n  GO 후보:")
    go_list = [r for r in ranked if "GO" in r["score"]["verdict"] and "NO" not in r["score"]["verdict"]]
    if go_list:
        for r in go_list:
            print(f"    {r['candidate']['id']} — {r['candidate']['name']}")
            print(f"      설계 근거: {r['candidate'].get('rationale', r['candidate'].get('name', ''))}")
    else:
        print("    없음 (Gen1 최적화 필요)")
    print(f"{'='*65}\n")


GEN2_CANDIDATES = [
    {
        "id": "CG2-001",
        "smiles": "OC(=O)C1CCN(Cc2ccc(F)c(F)c2)CC1",
        "name": "1-(3,4-difluorobenzyl)piperidine-4-carboxylic acid",
        "scaffold": "Piperidine-4-COOH + 3,4-diF-benzyl",
        "adc_handle": "COOH",
        "rationale": "CG1-001 + 3,4-diF: 추가 F 소수성 접촉 → TM 소수성 포켓 증진 / DILI 저위험 유지",
    },
    {
        "id": "CG2-002",
        "smiles": "OC(=O)C1CCC(Cc2ccc(F)c(F)c2)CC1",
        "name": "4-(3,4-difluorobenzyl)cyclohexane-1-carboxylic acid",
        "scaffold": "Cyclohexane-COOH + 3,4-diF-benzyl (N-free)",
        "adc_handle": "COOH",
        "rationale": "CG1-004 + 3,4-diF: N-free로 hERG 최소 / 두 F → van der Waals + polar 접촉 동시",
    },
    {
        "id": "CG2-003",
        "smiles": "OC(=O)C1CCN(CC(=O)Nc2ccc(F)cc2)CC1",
        "name": "2-((4-(4-carboxypiperidine-1-yl)acetamido))-4-fluorobenzene",
        "scaffold": "Piperidine-4-COOH + N-CH2-amide-4F-phenyl",
        "adc_handle": "COOH",
        "rationale": "amide NH H-bond donor: TM Ser/Asn 접촉 기대 / COOH+amide 이중 H-bond 전략",
    },
    {
        "id": "CG2-004",
        "smiles": "OC(=O)C1CCN(Cc2ccc(F)c(F)c2)C(C)C1",
        "name": "1-(3,4-difluorobenzyl)-2-methylpiperidine-4-carboxylic acid",
        "scaffold": "2-methyl-piperidine-4-COOH + 3,4-diF-benzyl",
        "adc_handle": "COOH",
        "rationale": "CG2-001 + 2-methyl: piperidine methylation → TM 소수성 채널 shape 최적화 + Fsp3 증가",
    },
    {
        "id": "CG2-005",
        "smiles": "OC(=O)C1CCN(CCc2ccc(F)c(F)c2)CC1",
        "name": "1-(3,4-difluorophenethyl)piperidine-4-carboxylic acid",
        "scaffold": "Piperidine-4-COOH + 3,4-diF-phenethyl (CH2CH2 링커)",
        "adc_handle": "COOH",
        "rationale": "benzyl → phenethyl: CH2 연장 → TM 포켓 심부 탐색 / 더 큰 van der Waals surface",
    },
]


GEN3_CANDIDATES = [
    {
        "id": "CG3-001",
        "smiles": "OC(=O)C1CCN(CCc2ccc(F)c(F)c2OC)CC1",
        "name": "1-(4,5-difluoro-2-methoxyphenethyl)piperidine-4-carboxylic acid",
        "scaffold": "Piperidine-4-COOH + 4,5-diF-2-OMe-phenethyl",
        "adc_handle": "COOH",
        "rationale": "CG2-005 + ortho-OMe: H-bond acceptor 추가 → TM Ser 접촉 / MW=299 포켓 체적 확장",
    },
    {
        "id": "CG3-002",
        "smiles": "OC(=O)C1CCN(CCc2ccc(F)c(F)c2C2CC2)CC1",
        "name": "1-(2-cyclopropyl-4,5-difluorophenethyl)piperidine-4-carboxylic acid",
        "scaffold": "Piperidine-4-COOH + cyclopropyl-diF-phenethyl",
        "adc_handle": "COOH",
        "rationale": "cyclopropyl 입체 확장 → TM5/TM6 소수성 포켓 형상 최적화 / Fsp3=0.588",
    },
    {
        "id": "CG3-003",
        "smiles": "OC(=O)C1CCN(CCCc2ccc(F)c(F)c2OC)CC1",
        "name": "1-(3-(4,5-difluoro-2-methoxyphenyl)propyl)piperidine-4-carboxylic acid",
        "scaffold": "Piperidine-4-COOH + propyl-diF-OMe-phenyl",
        "adc_handle": "COOH",
        "rationale": "링커 CH2 연장(ethyl→propyl) + OMe: 포켓 심부 탐색 / MW=313 최대 체적",
    },
    {
        "id": "CG3-004",
        "smiles": "OC(=O)C1COCCN1CCc2ccc(F)c(F)c2OC",
        "name": "4-(4,5-difluoro-2-methoxyphenethyl)morpholine-3-carboxylic acid",
        "scaffold": "Morpholine-3-COOH + diF-OMe-phenethyl",
        "adc_handle": "COOH",
        "rationale": "morpholine O: hERG 위험 최소화 + H-bond acceptor / OMe: TM 포켓 친화도 증가",
    },
    {
        "id": "CG3-005",
        "smiles": "OC(=O)C12CCN(CCc3ccc(F)c(F)c3)CC1CCC2",
        "name": "2'-azaspiro[cyclobutane-1,4'-piperidine]-1'-(3,4-difluorophenethyl)-5'-carboxylic acid",
        "scaffold": "Spiro[cyclobutane-piperidine]-COOH + 3,4-diF-phenethyl",
        "adc_handle": "COOH",
        "rationale": "spiro 링 제약: TM 포켓 3D 형상 lock / cyclobutane 체적 + rigid scaffold → 엔트로피 이득",
    },
]

# Gen3 PAM 포켓 좌표 (TM5-7 Cα centroid, chain R 7E9H)
GEN3_GRID_CENTER = (207.216, 184.017, 204.370)

GEN4_CANDIDATES = [
    {
        "id": "CG4-001",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)CC1CCN1CCOCC1",
        "name": "1-(3,4-difluorophenyl)-3-(2-morpholinoethyl)cyclopropane-1-carboxylic acid",
        "scaffold": "VU0155041 core + morpholine tail",
        "adc_handle": "COOH",
        "rationale": "cyclopropane geometry -> W773 pi-stacking / morpholine: H-bond acceptor + sp3 / MW=311",
    },
    {
        "id": "CG4-002",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)CC1OC",
        "name": "1-(3,4-difluorophenyl)-2-methoxycyclopropane-1-carboxylic acid",
        "scaffold": "VU0155041 core + 2-OMe cyclopropane",
        "adc_handle": "COOH",
        "rationale": "OMe at C2: H-bond acceptor within TM cavity / compact fit / MW=228",
    },
    {
        "id": "CG4-003",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)CC1CN1CCCCC1",
        "name": "1-(3,4-difluorophenyl)-2-(piperidinylmethyl)cyclopropane-1-carboxylic acid",
        "scaffold": "VU0155041 core + piperidinyl-CH2",
        "adc_handle": "COOH",
        "rationale": "piperidine N H-bond donor/acceptor: TM Ser contact / hydrophobic extension MW=295",
    },
    {
        "id": "CG4-004",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)CC1c1ccncc1",
        "name": "1-(3,4-difluorophenyl)-2-(pyridin-4-yl)cyclopropane-1-carboxylic acid",
        "scaffold": "VU0155041 core + pyridyl (dual pi-stack)",
        "adc_handle": "COOH",
        "rationale": "diF-phenyl + pyridyl dual aromatic: W773 pi-stack + N H-bond acceptor / MW=275",
    },
    {
        "id": "CG4-005",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)CC1(C)C",
        "name": "1-(3,4-difluorophenyl)-2,2-dimethylcyclopropane-1-carboxylic acid",
        "scaffold": "VU0155041 core + gem-dimethyl",
        "adc_handle": "COOH",
        "rationale": "gem-dimethyl: hydrophobic fill of TM pocket / rigid conformation lock / MW=226",
    },
]
# Gen4: key PAM binding residue centroid (W773^6.48 포함, chain R 7E9H)
GEN4_GRID_CENTER = (207.213, 186.578, 202.506)

# ─────────────────────────────────────────────
# Gen5 후보: CG4-003 최적화 (cis/trans + H-bond 강화)
# ─────────────────────────────────────────────
GEN5_CANDIDATES = [
    {
        "id": "CG5-001",
        "smiles": "OC(=O)[C@@]1(c2ccc(F)c(F)c2)C[C@H]1CN1CCCCC1",
        "name": "(1S,3R)-1-(3,4-difluorophenyl)-3-(piperidinylmethyl)cyclopropane-1-COOH (cis)",
        "scaffold": "CG4-003 cis 입체이성질체",
        "adc_handle": "COOH",
        "rationale": "cis: aryl/CH2N 동일면 → TM 포켓 심부 tight fit / 결합 기하학 최적화",
    },
    {
        "id": "CG5-002",
        "smiles": "OC(=O)[C@@]1(c2ccc(F)c(F)c2)C[C@@H]1CN1CCCCC1",
        "name": "(1S,3S)-1-(3,4-difluorophenyl)-3-(piperidinylmethyl)cyclopropane-1-COOH (trans)",
        "scaffold": "CG4-003 trans 입체이성질체",
        "adc_handle": "COOH",
        "rationale": "trans: aryl/CH2N 반대면 → 포켓 입구 탐색 / 다른 결합 mode 비교",
    },
    {
        "id": "CG5-003",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)CC1C[N+]1(CCCCC1)[O-]",
        "name": "1-(3,4-difluorophenyl)-2-((1-oxido-1lambda5-piperidin-1-yl)methyl)cyclopropane-1-COOH",
        "scaffold": "CG4-003 piperidine N-oxide",
        "adc_handle": "COOH",
        "rationale": "N-oxide: N+→O- 쌍극자 → TM Ser/Thr H-bond 극대화 / TPSA+20 → hERG 저감",
    },
    {
        "id": "CG5-004",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)CC1CN1CCNCC1",
        "name": "1-(3,4-difluorophenyl)-2-(piperazin-1-ylmethyl)cyclopropane-1-COOH",
        "scaffold": "CG4-003 piperidine → piperazine (NH 추가)",
        "adc_handle": "COOH",
        "rationale": "piperazine NH: H-bond donor → TM Glu/Asp 직접 접촉 / HBD=2 (NH+COOH)",
    },
    {
        "id": "CG5-005",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)CC1CNC(C)=O",
        "name": "N-((1-(3,4-difluorophenyl)-1-(carboxyl)cyclopropan-2-yl)methyl)acetamide",
        "scaffold": "CG4-003 piperidine → acetylamide (NHCO H-bond pair)",
        "adc_handle": "COOH",
        "rationale": "NH + C=O 이중 H-bond: TM Ser+Asn 동시 접촉 가능 / 링 제거로 entropy 이득",
    },
]
GEN5_GRID_CENTER = GEN4_GRID_CENTER   # 동일 포켓 좌표 유지

# ─────────────────────────────────────────────
# Gen6 후보: CG4-003 라세미 베이스 / Option A+C 동시 적용
#   Option A: cyclopropane C3 소수성 확장 → TM5 Val/Leu 클램프 접촉
#   Option C: COOH → tetrazole bioisostere → 이온 접촉 개선
# ─────────────────────────────────────────────
GEN6_CANDIDATES = [
    {
        "id": "CG6-001",
        "smiles": "C1(c2ccc(F)c(F)c2)(c3nnn[nH]3)CC1CN4CCCCC4",
        "name": "1-(3,4-difluorophenyl)-1-(1H-tetrazol-5-yl)-2-(piperidinylmethyl)cyclopropane",
        "scaffold": "CG4-003 COOH→tetrazole (Option C 기준점)",
        "adc_handle": "tetrazole",
        "rationale": "순수 bioisostere 치환: pKa≈4.9 유지, 4개 N-H-bond acceptor로 이온 접촉 개선",
    },
    {
        "id": "CG6-002",
        "smiles": "C1(c2ccc(F)c(F)c2)(c3nnn[nH]3)CC1(C)CN4CCCCC4",
        "name": "1-(3,4-difluorophenyl)-1-(tetrazol-5-yl)-3-methyl-2-(piperidinylmethyl)cyclopropane",
        "scaffold": "CG4-003 COOH→tetrazole + C3-methyl",
        "adc_handle": "tetrazole",
        "rationale": "Option A+C: C3-methyl → TM5 Val/Leu 소수성 클램프 / Fsp3=0.588 유지",
    },
    {
        "id": "CG6-003",
        "smiles": "C1(c2ccc(F)c(F)c2)(c3nnn[nH]3)CC1(OC)CN4CCCCC4",
        "name": "1-(3,4-difluorophenyl)-1-(tetrazol-5-yl)-3-methoxy-2-(piperidinylmethyl)cyclopropane",
        "scaffold": "CG4-003 COOH→tetrazole + C3-OMe",
        "adc_handle": "tetrazole",
        "rationale": "Option A+C: C3-OMe H-bond acceptor → TM Ser/Thr 추가 극성 접촉 / HBA=5",
    },
    {
        "id": "CG6-004",
        "smiles": "C1(c2ccc(F)c(F)c2)(c3nnn[nH]3)CC1(C2CC2)CN4CCCCC4",
        "name": "1-(3,4-difluorophenyl)-1-(tetrazol-5-yl)-3-cyclopropyl-2-(piperidinylmethyl)cyclopropane",
        "scaffold": "CG4-003 COOH→tetrazole + C3-cyclopropyl",
        "adc_handle": "tetrazole",
        "rationale": "Option A+C: cyclopropyl Walsh orbital → 추가 π-stacking 기여 / Fsp3=0.632 최고",
    },
    {
        "id": "CG6-005",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)CC1(OC)CN3CCCCC3",
        "name": "1-(3,4-difluorophenyl)-1-carboxy-3-methoxy-2-(piperidinylmethyl)cyclopropane",
        "scaffold": "CG4-003 COOH 유지 + C3-OMe (대조군)",
        "adc_handle": "COOH",
        "rationale": "Option A only: COOH 기준 C3-OMe 효과 단독 검증 / tetrazole vs COOH 비교",
    },
]
GEN6_GRID_CENTER = GEN4_GRID_CENTER   # 동일 PAM 포켓 좌표 유지

# ─────────────────────────────────────────────
# Gen7 후보: CG4-003 베이스 / COOH 유지 / tetrazole&C3수식 블랙리스트
#   탐색1: cyclopropane C2 브릿지 수식 (C2-F, C2-Me) — 미탐색 축
#   탐색2: piperidine 링 내부 치환 (4-Me, 4-F, 3,3-diF)
#   exhaustiveness=32 / box=30 (버그 수정 후 자동 적용)
# ─────────────────────────────────────────────
GEN7_CANDIDATES = [
    {
        "id": "CG7-001",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)C(F)C1CN1CCCCC1",
        "name": "1-(3,4-difluorophenyl)-2-fluoro-3-(piperidinylmethyl)cyclopropane-1-COOH",
        "scaffold": "CG4-003 C2-F 수식 (C2 브릿지 탐색)",
        "adc_handle": "COOH",
        "rationale": "C2-H→C2-F: 대사 안정화 + cyclopropane 기하학 미조정 → TM pocket 밀착도 변화",
    },
    {
        "id": "CG7-002",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)C(C)C1CN1CCCCC1",
        "name": "1-(3,4-difluorophenyl)-2-methyl-3-(piperidinylmethyl)cyclopropane-1-COOH",
        "scaffold": "CG4-003 C2-Me 수식 (C2 브릿지 탐색)",
        "adc_handle": "COOH",
        "rationale": "C2-H→C2-Me: TM5 방향 소수성 클램프 탐색 / C2-F 직접 비교기준",
    },
    {
        "id": "CG7-003",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)CC1CN1CCC(C)CC1",
        "name": "1-(3,4-difluorophenyl)-2-((4-methylpiperidin-1-yl)methyl)cyclopropane-1-COOH",
        "scaffold": "CG4-003 pip-4-methyl",
        "adc_handle": "COOH",
        "rationale": "4-Me pip: TM7 Val775/Phe817 소수성 포켓 심층 접촉 / 링 equatorial 고정",
    },
    {
        "id": "CG7-004",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)CC1CN1CCC(F)CC1",
        "name": "1-(3,4-difluorophenyl)-2-((4-fluoropiperidin-1-yl)methyl)cyclopropane-1-COOH",
        "scaffold": "CG4-003 pip-4-fluoro",
        "adc_handle": "COOH",
        "rationale": "4-F pip: gauche 효과 → pip 링 구조 고정 + C-F 전자 효과로 N 기하학 변경",
    },
    {
        "id": "CG7-005",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)CC1CN1CCC(F)(F)CC1",
        "name": "1-(3,4-difluorophenyl)-2-((3,3-difluoropiperidin-1-yl)methyl)cyclopropane-1-COOH",
        "scaffold": "CG4-003 pip-3,3-gem-diF",
        "adc_handle": "COOH",
        "rationale": "3,3-diF pip: gem-diF 대사 차단 + 링 평면성↑ → N lone pair 방향 최적화",
    },
]
GEN7_GRID_CENTER = GEN4_GRID_CENTER

# ─────────────────────────────────────────────
# Gen8 후보: Gen7 최우수 기능 조합 synergy test
#   C2-브릿지 수식 (Me/F) × pip 불소화 (4,4-diF / 4-F) 조합 매트릭스
#   CG8-005: C2-gem-diF (CF2 cyclopropane) 탐색
# ─────────────────────────────────────────────
GEN8_CANDIDATES = [
    {
        "id": "CG8-001",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)C(C)C1CN1CCC(F)(F)CC1",
        "name": "1-(3,4-diF-Ph)-2-methyl-3-((4,4-diF-piperidin-1-yl)methyl)cyclopropane-1-COOH",
        "scaffold": "CG7-002 + CG7-005 조합 (C2-Me × 4,4-diF pip)",
        "adc_handle": "COOH",
        "rationale": "synergy: C2-Me 소수성 클램프 × 4,4-diF pip 링 고정 → additive ΔG 기대",
    },
    {
        "id": "CG8-002",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)C(F)C1CN1CCC(F)(F)CC1",
        "name": "1-(3,4-diF-Ph)-2-fluoro-3-((4,4-diF-piperidin-1-yl)methyl)cyclopropane-1-COOH",
        "scaffold": "CG7-001 + CG7-005 조합 (C2-F × 4,4-diF pip)",
        "adc_handle": "COOH",
        "rationale": "synergy: C2-F 기하학 미조정 × 4,4-diF pip → fluorine-rich variant",
    },
    {
        "id": "CG8-003",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)C(C)C1CN1CCC(F)CC1",
        "name": "1-(3,4-diF-Ph)-2-methyl-3-((4-F-piperidin-1-yl)methyl)cyclopropane-1-COOH",
        "scaffold": "CG7-002 + CG7-004 조합 (C2-Me × 4-F pip)",
        "adc_handle": "COOH",
        "rationale": "synergy: C2-Me × 4-F pip gauche 효과 / MW 최소화 조합",
    },
    {
        "id": "CG8-004",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)C(F)C1CN1CCC(F)CC1",
        "name": "1-(3,4-diF-Ph)-2-fluoro-3-((4-F-piperidin-1-yl)methyl)cyclopropane-1-COOH",
        "scaffold": "CG7-001 + CG7-004 조합 (C2-F × 4-F pip)",
        "adc_handle": "COOH",
        "rationale": "synergy: C2-F × 4-F pip / 극성 낮은 양면 F 도입",
    },
    {
        "id": "CG8-005",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1CCCCC1",
        "name": "1-(3,4-diF-Ph)-2,2-diF-3-(piperidinylmethyl)cyclopropane-1-COOH",
        "scaffold": "C2-gem-diF cyclopropane (CF₂ bridge)",
        "adc_handle": "COOH",
        "rationale": "C2에 gem-diF: Walsh orbital 전자 효과 최대화 + 포켓 기하학 탐색",
    },
]
GEN8_GRID_CENTER = GEN4_GRID_CENTER

# ─────────────────────────────────────────────
# Gen9: CG8-005 (CF₂, Ki=1.51μM) × CG8-001 (4,4-diF pip, Ki=1.97μM) 최대 조합
# ─────────────────────────────────────────────
GEN9_CANDIDATES = [
    {
        "id": "CG9-001",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1CCC(F)(F)CC1",
        "name": "1-(3,4-diF-Ph)-2,2-diF-3-((4,4-diF-pip-1-yl)methyl)cyclopropane-1-COOH",
        "scaffold": "CG8-005 + CG8-001 핵심 조합 (CF₂ × 4,4-diF pip)",
        "adc_handle": "COOH",
        "rationale": "CF₂ Walsh 전자효과 + 4,4-diF pip 링 고정 → additive ΔG 최대화",
    },
    {
        "id": "CG9-002",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1CCC(F)CC1",
        "name": "1-(3,4-diF-Ph)-2,2-diF-3-((4-F-pip-1-yl)methyl)cyclopropane-1-COOH",
        "scaffold": "CF₂-cyclopropane + 4-F pip",
        "adc_handle": "COOH",
        "rationale": "CF₂ × 4-F pip: 단일 pip-F로 부분 최적화 / MW 최소화",
    },
    {
        "id": "CG9-003",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)C(CC)C1CN1CCC(F)(F)CC1",
        "name": "1-(3,4-diF-Ph)-2-ethyl-3-((4,4-diF-pip-1-yl)methyl)cyclopropane-1-COOH",
        "scaffold": "C2-Et + 4,4-diF pip (CG8-001 C2-Me→Et)",
        "adc_handle": "COOH",
        "rationale": "C2-Me를 Et로 확장: TM5 소수성 포켓 심층 접촉 + 4,4-diF pip synergy",
    },
    {
        "id": "CG9-004",
        "smiles": "OC(=O)C1(c2cc(F)c(F)c(F)c2)C(F)(F)C1CN1CCCCC1",
        "name": "1-(2,3,4-triF-Ph)-2,2-diF-3-(pip-1-ylmethyl)cyclopropane-1-COOH",
        "scaffold": "2,3,4-triF-Ph + CF₂-cyclopropane",
        "adc_handle": "COOH",
        "rationale": "phenyl 2-위치 추가 F: ortho-F → C1 근방 전자 밀도 변화 + W773 stacking 기여",
    },
    {
        "id": "CG9-005",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1CC(F)(F)CCC1",
        "name": "1-(3,4-diF-Ph)-2,2-diF-3-((3,3-diF-pip-1-yl)methyl)cyclopropane-1-COOH",
        "scaffold": "CF₂-cyclopropane + 3,3-diF pip (N-adjacent)",
        "adc_handle": "COOH",
        "rationale": "3,3-diF pip: N 인접 gem-diF → N lone pair 방향 제약 + 4,4-diF와 다른 기하학",
    },
]
GEN9_GRID_CENTER = GEN4_GRID_CENTER

# ─────────────────────────────────────────────
# GEN10: Ki<0.5μM 목표 / CG9-005 scaffold 기반
# Base: OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1CC(F)(F)CCC1 (Ki=0.725μM)
# A1: 3-CF3 phenyl (meta-trifluoromethyl)
# A2: 3-F,4-Cl mixed halogen phenyl
# A3: 3,4,5-triF phenyl (add F at 5-position)
# B1: C3-Me (quaternary C3 on cyclopropane, steric TM pocket fill)
# C1: 2-Me-3,3-diF piperidine (alpha-methyl conformational restriction)
# ─────────────────────────────────────────────
GEN10_CANDIDATES = [
    {
        "id": "CG10-001",
        "smiles": "OC(=O)C1(c2cccc(C(F)(F)F)c2)C(F)(F)C1CN1CC(F)(F)CCC1",
        "name": "A1: 3-CF3-Ph",
        "scaffold": "CF2-cyclopropane + 3,3-diF-pip",
        "adc_handle": "COOH",
    },
    {
        "id": "CG10-002",
        "smiles": "OC(=O)C1(c2ccc(Cl)c(F)c2)C(F)(F)C1CN1CC(F)(F)CCC1",
        "name": "A2: 3-F,4-Cl-Ph",
        "scaffold": "CF2-cyclopropane + 3,3-diF-pip",
        "adc_handle": "COOH",
    },
    {
        "id": "CG10-003",
        "smiles": "OC(=O)C1(c2cc(F)c(F)c(F)c2)C(F)(F)C1CN1CC(F)(F)CCC1",
        "name": "A3: 3,4,5-triF-Ph",
        "scaffold": "CF2-cyclopropane + 3,3-diF-pip",
        "adc_handle": "COOH",
    },
    {
        "id": "CG10-004",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1(C)CN1CC(F)(F)CCC1",
        "name": "B1: C3-Me cyclopropane",
        "scaffold": "CF2-cyclopropane + C3-Me + 3,3-diF-pip",
        "adc_handle": "COOH",
    },
    {
        "id": "CG10-005",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1C(C)C(F)(F)CCC1",
        "name": "C1: 2-Me-3,3-diF-pip",
        "scaffold": "CF2-cyclopropane + 2-Me-3,3-diF-pip",
        "adc_handle": "COOH",
    },
]
GEN10_GRID_CENTER = GEN4_GRID_CENTER

# ─────────────────────────────────────────────
# GEN11: Ki<0.2μM 목표 / CG10-005 scaffold 고정
# Base: OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1C(C)C(F)(F)CCC1
#        (Ki=0.345μM, MD PASS 1.69Å, ΔG=-8.811)
# pip α-methyl(2-Me) 위치 변형:
#   CG11-001: 2-Me → 2-Et (크기 확장, TM6 Leu805 접촉면)
#   CG11-002: 2-Me → 2,2-diMe (gem-diMe, 구조 고정)
#   CG11-003: 2-Me → 2-F (전자 인출, 친전자 포켓)
#   CG11-004: 2-Me 유지 + 피페리딘 5위치 F 추가 (추가 fluorination)
#   CG11-005: 2-Me 유지 + cyclopropane C3에 CN 추가 (극성 핸들)
# exhaustiveness=48 (재현성 강화)
# ─────────────────────────────────────────────
GEN11_CANDIDATES = [
    {
        "id": "CG11-001",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1C(CC)C(F)(F)CCC1",
        "name": "2-Et-3,3-diF-pip (alpha-Et)",
        "scaffold": "CF2-cyclopropane + 2-Et-3,3-diF-pip",
        "adc_handle": "COOH",
        "exhaustiveness": 48,
    },
    {
        "id": "CG11-002",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1C(C)(C)C(F)(F)CCC1",
        "name": "2,2-diMe-3,3-diF-pip (gem-diMe)",
        "scaffold": "CF2-cyclopropane + 2,2-diMe-3,3-diF-pip",
        "adc_handle": "COOH",
        "exhaustiveness": 48,
    },
    {
        "id": "CG11-003",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1C(F)C(F)(F)CCC1",
        "name": "2-F-3,3-diF-pip (alpha-F)",
        "scaffold": "CF2-cyclopropane + 2-F-3,3-diF-pip",
        "adc_handle": "COOH",
        "exhaustiveness": 48,
    },
    {
        "id": "CG11-004",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1C(C)C(F)(F)CC(F)C1",
        "name": "2-Me-3,3-diF-5F-pip (additional 5-F)",
        "scaffold": "CF2-cyclopropane + 2-Me-3,3-diF-5F-pip",
        "adc_handle": "COOH",
        "exhaustiveness": 48,
    },
    {
        "id": "CG11-005",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1(CC#N)CN1C(C)C(F)(F)CCC1",
        "name": "C3-CH2CN + 2-Me-3,3-diF-pip (polar nitrile handle)",
        "scaffold": "CF2-cyclopropane + C3-CH2CN + 2-Me-3,3-diF-pip",
        "adc_handle": "COOH",
        "exhaustiveness": 48,
    },
]
GEN11_GRID_CENTER = GEN4_GRID_CENTER

# ─────────────────────────────────────────────
# GEN12: CG11-002 gem-diMe scaffold hERG 감소 전략
# Base ref: CG11-002 Ki=0.251μM, hERG=0.344 (0.044 초과)
# 전략: LogP 감소 (F 제거 / 극성 추가) → hERG < 0.30 목표
#   CG12-001: gem-diMe + 3-monoF pip  (3,3-diF→3-F, LogP 3.76)
#   CG12-002: gem-diMe + plain pip    (pip F 전부 제거, LogP 3.82)
#   CG12-003: gem-diMe + 4,4-diF pip  (diF 위치 3→4, LogP 4.06)
#   CG12-004: CG10-005 base + 4-OH pip (2-Me+3,3-diF+OH, LogP 2.64)
#   CG12-005: gem-diMe + CF-cyclopropane (CF2→CF, LogP 3.76)
# exhaustiveness=48
# ─────────────────────────────────────────────
GEN12_CANDIDATES = [
    {
        "id": "CG12-001",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1C(C)(C)C(F)CCC1",
        "name": "gem-diMe + 3-monoF pip",
        "scaffold": "CF2-cyclopropane + gem-diMe-3F-pip",
        "adc_handle": "COOH",
        "exhaustiveness": 48,
    },
    {
        "id": "CG12-002",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1C(C)(C)CCCC1",
        "name": "gem-diMe + plain pip (no F)",
        "scaffold": "CF2-cyclopropane + gem-diMe-pip",
        "adc_handle": "COOH",
        "exhaustiveness": 48,
    },
    {
        "id": "CG12-003",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1C(C)(C)CC(F)(F)CC1",
        "name": "gem-diMe + 4,4-diF pip",
        "scaffold": "CF2-cyclopropane + gem-diMe-4,4-diF-pip",
        "adc_handle": "COOH",
        "exhaustiveness": 48,
    },
    {
        "id": "CG12-004",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1C(C)C(F)(F)C(O)CC1",
        "name": "2-Me + 3,3-diF + 4-OH pip (polar)",
        "scaffold": "CF2-cyclopropane + 2-Me-3,3-diF-4OH-pip",
        "adc_handle": "COOH",
        "exhaustiveness": 48,
    },
    {
        "id": "CG12-005",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)C(F)C1CN1C(C)(C)C(F)(F)CCC1",
        "name": "gem-diMe + 3,3-diF + CF-cyclopropane",
        "scaffold": "CF-cyclopropane + gem-diMe-3,3-diF-pip",
        "adc_handle": "COOH",
        "exhaustiveness": 48,
    },
]
GEN12_GRID_CENTER = GEN4_GRID_CENTER

# ─────────────────────────────────────────────
# GEN13: Track A (OH위치 변경) + Track B (4-위치 치환체 비교)
# 참조: CG10-005 Ki=0.345μM hERG=0.259 / CG12-004 Ki=0.580μM hERG=0.154 (4-OH)
# Track A: CG12-004 scaffold (2-Me+3,3-diF+CF2) + 5-OH / 6-OH 이동
#   CG13-A1: 5-OH — LogP 2.64 (CG12-004과 동일)
#   CG13-A2: 6-OH — LogP 2.99 (알파-OH, N 인접)
# Track B: CG10-005 base (2-Me+3,3-diF+CF2) + 4-위치 치환체
#   CG13-B1: 4-F  — LogP 3.62 (친유성 유지, H-bond donor 없음)
#   CG13-B2: 4-OMe — LogP 3.30 (극성 있으나 donor 없음)
# exhaustiveness=48
# ─────────────────────────────────────────────
GEN13_CANDIDATES = [
    {
        "id": "CG13-A1",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1C(C)C(F)(F)CC(O)C1",
        "name": "Track A: 5-OH (C4→C5 shift)",
        "scaffold": "CF2-cyclopropane + 2-Me-3,3-diF-5OH-pip",
        "adc_handle": "COOH",
        "exhaustiveness": 48,
    },
    {
        "id": "CG13-A2",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1C(C)C(F)(F)CCC1O",
        "name": "Track A: 6-OH (alpha-OH, C4→C6 shift)",
        "scaffold": "CF2-cyclopropane + 2-Me-3,3-diF-6OH-pip",
        "adc_handle": "COOH",
        "exhaustiveness": 48,
    },
    {
        "id": "CG13-B1",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1C(C)C(F)(F)C(F)CC1",
        "name": "Track B: 4-F (F bioisostere, LogP 3.62)",
        "scaffold": "CF2-cyclopropane + 2-Me-3,3,4-triF-pip",
        "adc_handle": "COOH",
        "exhaustiveness": 48,
    },
    {
        "id": "CG13-B2",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1C(C)C(F)(F)C(OC)CC1",
        "name": "Track B: 4-OMe (methoxy, LogP 3.30)",
        "scaffold": "CF2-cyclopropane + 2-Me-3,3-diF-4OMe-pip",
        "adc_handle": "COOH",
        "exhaustiveness": 48,
    },
]
GEN13_GRID_CENTER = GEN4_GRID_CENTER

# ─────────────────────────────────────────────
# GEN14: DILI 억제 2종 전략 검증
# Gen13 DILI FAIL 원인 분석: 3,3-diF pip 전자밀도 과잉 + C4 이외 OH → DILI↑
# 전략 A: 3-monoF+5-OH pip (3,3-diF→3-F, OH위치 조정)  LogP 2.35
# 전략 B: CHF-cyclopropane (CF2→CHF) + 5-OH pip          LogP 2.35
# 대조군:  CHF-cyclopropane + CG10-005 pip (OH 없음)      LogP 3.37
# exhaustiveness=48
# ─────────────────────────────────────────────
GEN14_CANDIDATES = [
    {
        "id": "CG14-001",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1C(C)C(F)CC(O)C1",
        "name": "3-monoF+5-OH pip (CF2-cyclopropane)",
        "scaffold": "CF2-cyclopropane + 2-Me-3F-5OH-pip",
        "adc_handle": "COOH",
        "exhaustiveness": 48,
    },
    {
        "id": "CG14-002",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)[C@@H](F)C1CN1C(C)C(F)(F)CC(O)C1",
        "name": "CHF-cyclopropane + 3,3-diF+5-OH pip",
        "scaffold": "CHF-cyclopropane + 2-Me-3,3-diF-5OH-pip",
        "adc_handle": "COOH",
        "exhaustiveness": 48,
    },
    {
        "id": "CG14-003",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)[C@@H](F)C1CN1C(C)C(F)(F)CCC1",
        "name": "CHF-cyclopropane + CG10-005 pip (대조군)",
        "scaffold": "CHF-cyclopropane + 2-Me-3,3-diF-pip",
        "adc_handle": "COOH",
        "exhaustiveness": 48,
    },
]
GEN14_GRID_CENTER = GEN4_GRID_CENTER

# ─────────────────────────────────────────────
# GEN15: Ki<0.400μM + DILI<0.400 + hERG<0.300 동시 달성
# 기준: CG13-B1 Ki=0.372 DILI=0.404(FAIL) — DILI만 0.004 초과
# 001: 2-Me 제거 + 3,3-diF+4-F pip → DILI 미세 완화, LogP 3.23
# 002: 3,3-gem-diF → 3,4-vicinal-diF pip → 전자밀도 분산, LogP 3.32
# 003: Ph F위치 3,4→3,5 이성질체 + pip 4-F → 전자분포 변화, LogP 3.62
# 004: pip 4-F → 4-Cl 바이오이소스터 → DILI 패턴 변화, LogP 3.89
# exhaustiveness=48
# ─────────────────────────────────────────────
GEN15_CANDIDATES = [
    {
        "id": "CG15-001",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1CC(F)(F)C(F)CC1",
        "name": "N-deMe+3,3-diF+4-F pip",
        "rationale": "2-Me 제거로 입체/전자 완화 → DILI 미세조정, Ki 유지 기대",
        "scaffold": "CF2-cyclopropane",
        "adc_handle": "COOH",
        "exhaustiveness": 48,
    },
    {
        "id": "CG15-002",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1C(C)C(F)C(F)CC1",
        "name": "3,4-diF pip (gem→vicinal 재배치)",
        "rationale": "3,3-gem-diF→3,4-vicinal-diF: 전자밀도 분산으로 DILI↓ 기대",
        "scaffold": "CF2-cyclopropane",
        "adc_handle": "COOH",
        "exhaustiveness": 48,
    },
    {
        "id": "CG15-003",
        "smiles": "OC(=O)C1(c2cc(F)c(F)cc2)C(F)(F)C1CN1C(C)C(F)(F)C(F)CC1",
        "name": "Ph-3,5-diF + pip-2-Me-3,3-diF-4-F",
        "rationale": "페닐 F위치 3,4→3,5 이성질체: W773 접촉 유지하며 전자분포 변화",
        "scaffold": "CF2-cyclopropane",
        "adc_handle": "COOH",
        "exhaustiveness": 48,
    },
    {
        "id": "CG15-004",
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1C(C)C(F)(F)C(Cl)CC1",
        "name": "4-F→4-Cl 바이오이소스터",
        "rationale": "F→Cl 바이오이소스터: DILI 패턴 변화 가능, LogP 3.89",
        "scaffold": "CF2-cyclopropane",
        "adc_handle": "COOH",
        "exhaustiveness": 48,
    },
]
GEN15_GRID_CENTER = GEN4_GRID_CENTER


def run_cancer_gate_gen(gen_name: str, candidates: list, grid_center=None):
    """공용 Cancer-Gate 스크리닝 실행기."""
    print(f"\n{'='*65}")
    print(f"  Cancer-Gate Pipeline v1.0  —  {gen_name} 스크리닝")
    print(f"  타깃: GRM4 (mGluR4)  |  PDB: {PDB_ID} Chain R")
    print(f"  필터: DILI < {DILI_THRESHOLD}  |  hERG < {HERG_THRESHOLD}  |  Ki < {KI_TARGET_UM}uM")
    if grid_center:
        print(f"  그리드: {grid_center} (PAM TM5-7 pocket) sz=25A  exhaustiveness=32")
    print(f"{'='*65}\n")

    DOCK_DIR.mkdir(parents=True, exist_ok=True)

    print(f"[SETUP] 수용체 준비 (7E9H 체인 R = mGlu4)")
    receptor_pdbqt, center = download_and_prepare_receptor()
    if grid_center:
        center = grid_center
        print(f"  [OVERRIDE] 그리드 중심: {center}")

    results = []
    n = len(candidates)
    for i, cand in enumerate(candidates, 1):
        cid = cand["id"]
        smiles = cand["smiles"]
        print(f"\n{'─'*65}")
        print(f"  [{i}/{n}] {cid} | {cand['name']}")
        print(f"  SMILES : {smiles}")
        print(f"  Scaffold: {cand['scaffold']} | 핸들: {cand['adc_handle']}")
        print(f"{'─'*65}")

        rdkit_r = calc_rdkit(smiles)
        if "error" in rdkit_r:
            print(f"  [RDKit 오류] {rdkit_r['error']}")
            continue
        print(f"  RDKit  MW={rdkit_r['mw']}  LogP={rdkit_r['logp']}  "
              f"QED={rdkit_r['qed']}  Ro5={'PASS' if rdkit_r['lipinski_pass'] else 'FAIL'}  "
              f"COOH={'Y' if rdkit_r['has_cooh'] else 'N'}")

        print(f"  ADMET  계산 중...", end=" ", flush=True)
        admet_r = calc_admet_and_filter(smiles, rdkit_r["has_cooh"])
        if "error" in admet_r:
            print(f"오류: {admet_r['error']}")
        else:
            cal_note = "(COOH +20%p 보정)" if admet_r.get("cooh_correction_applied") else ""
            print(f"완료 {cal_note}")
            print(f"  ADMET  hERG={admet_r.get('herg','N/A')}  "
                  f"DILI={admet_r.get('dili','N/A')}  "
                  f"AMES={admet_r.get('ames','N/A')}  "
                  f"HIA={admet_r.get('hia_corrected', admet_r.get('hia','N/A'))}  "
                  f"BBB={admet_r.get('bbb','N/A')}")
            gate = "PASS ✓" if admet_r.get("gate_pass") else "FAIL ✗"
            print(f"  Cancer-Gate: DILI={'PASS' if admet_r.get('dili_pass') else 'FAIL'}  "
                  f"hERG={'PASS' if admet_r.get('herg_pass') else 'FAIL'}  -> {gate}")

        _exh = 32 if int(cand["id"][2:].split('-')[0]) >= 3 else 16
        print(f"  Docking 실행 중 (exhaustiveness={_exh})...", end=" ", flush=True)
        dock_r = run_vina_for(cand, receptor_pdbqt, center)
        if "error" in dock_r:
            print(f"오류: {dock_r['error']}")
            dock_r = {"best_dG": -5.0, "best_Ki_uM": 999.9, "docking_pass": False}
        else:
            ki_note = "Ki 목표 달성" if dock_r["docking_pass"] else "Ki 목표 미달"
            print(f"완료")
            print(f"  Docking dG={dock_r['best_dG']} kcal/mol  "
                  f"Ki={dock_r['best_Ki_uM']}uM  {ki_note}")

        scored = score_candidate(rdkit_r, admet_r, dock_r)
        print(f"  종합   점수={scored['score']}/100  판정: {scored['verdict']}")
        if scored["flags"]:
            for fl in scored["flags"]:
                print(f"         * {fl}")

        results.append({
            "candidate": cand,
            "rdkit":  rdkit_r,
            "admet":  admet_r,
            "docking": dock_r,
            "score":  scored,
        })

    _print_ranking(results, gen=gen_name)

    report_path = BASE_DIR / f"cancer_gate_{gen_name.lower().replace(' ', '_')}_report.json"
    report_path.write_text(
        json.dumps(results, ensure_ascii=False, indent=2, default=str),
        encoding="utf-8"
    )
    print(f"\n  리포트: {report_path}\n")
    return results


if __name__ == "__main__":
    import sys as _sys
    arg = _sys.argv[1] if len(_sys.argv) > 1 else "gen0"
    if arg == "gen15":
        run_cancer_gate_gen("Gen15", GEN15_CANDIDATES, grid_center=GEN15_GRID_CENTER)
    elif arg == "gen14":
        run_cancer_gate_gen("Gen14", GEN14_CANDIDATES, grid_center=GEN14_GRID_CENTER)
    elif arg == "gen13":
        run_cancer_gate_gen("Gen13", GEN13_CANDIDATES, grid_center=GEN13_GRID_CENTER)
    elif arg == "gen12":
        run_cancer_gate_gen("Gen12", GEN12_CANDIDATES, grid_center=GEN12_GRID_CENTER)
    elif arg == "gen11":
        run_cancer_gate_gen("Gen11", GEN11_CANDIDATES, grid_center=GEN11_GRID_CENTER)
    elif arg == "gen10":
        run_cancer_gate_gen("Gen10", GEN10_CANDIDATES, grid_center=GEN10_GRID_CENTER)
    elif arg == "gen9":
        run_cancer_gate_gen("Gen9", GEN9_CANDIDATES, grid_center=GEN9_GRID_CENTER)
    elif arg == "gen8":
        run_cancer_gate_gen("Gen8", GEN8_CANDIDATES, grid_center=GEN8_GRID_CENTER)
    elif arg == "gen7":
        run_cancer_gate_gen("Gen7", GEN7_CANDIDATES, grid_center=GEN7_GRID_CENTER)
    elif arg == "gen6":
        run_cancer_gate_gen("Gen6", GEN6_CANDIDATES, grid_center=GEN6_GRID_CENTER)
    elif arg == "gen5":
        run_cancer_gate_gen("Gen5", GEN5_CANDIDATES, grid_center=GEN5_GRID_CENTER)
    elif arg == "gen4":
        run_cancer_gate_gen("Gen4", GEN4_CANDIDATES, grid_center=GEN4_GRID_CENTER)
    elif arg == "gen3":
        run_cancer_gate_gen("Gen3", GEN3_CANDIDATES, grid_center=GEN3_GRID_CENTER)
    elif arg == "gen2":
        run_cancer_gate_gen("Gen2", GEN2_CANDIDATES)
    elif arg == "gen1":
        run_cancer_gate_gen("Gen1", GEN1_CANDIDATES)
    else:
        run_cancer_gate()
