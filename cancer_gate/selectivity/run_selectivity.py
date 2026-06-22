"""
Selectivity Counter-Screening: GO candidates vs mGluR5 (5CGC)
=============================================================
목적: mGluR4 선택적 PAM 확인 — mGluR5에 약하게 결합해야 함
선택성 지수 = Ki(mGluR5) / Ki(mGluR4) >= 10 목표

receptor: C:\sigma_protocol\cancer_gate\selectivity\receptor_5CGC.pdbqt
ligands:  GO 후보 (CG10-005, CG12-004, CG14-001)
grid:     5CGC 51D NAM pocket centroid (-23.938, 16.228, 43.562)
"""

import os, json, subprocess, math, re

VINA     = r"C:\sigma_protocol\vina\vina.exe"
RECEPTOR = r"C:\sigma_protocol\cancer_gate\selectivity\receptor_5CGC.pdbqt"
OUT_DIR  = r"C:\sigma_protocol\cancer_gate\selectivity"
DOCK_DIR = r"C:\sigma_protocol\cancer_gate\docking\7E9H"  # ligand pdbqt source

# 5CGC NAM pocket (51D centroid)
GRID = {"x": -23.938, "y": 16.228, "z": 43.562}
BOX  = 25.0

# mGluR4 결과 (이전 스크리닝)
MGR4 = {
    "CG10-005": {"dG": -8.811, "Ki_uM": 0.345},
    "CG12-004": {"dG": -7.948, "Ki_uM": 0.580},
    "CG14-001": {"dG": -7.393, "Ki_uM": 0.982},
}

GO_CANDIDATES = [
    {"id": "CG10-005", "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1C(C)C(F)(F)CCC1"},
    {"id": "CG12-004", "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1C(C)C(F)(F)CC(O)C1"},
    {"id": "CG14-001", "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1C(C)C(F)CC(O)C1"},
]

RT = 0.592  # kcal/mol @ 298K

def dG_to_Ki(dG):
    return math.exp(dG / RT) * 1e6  # uM

def parse_vina_output(text):
    for line in text.splitlines():
        m = re.match(r'\s+1\s+([-\d.]+)', line)
        if m:
            return float(m.group(1))
    return None

def run_docking(cand_id):
    ligand_src = os.path.join(DOCK_DIR, f"ligand_{cand_id}.pdbqt")
    if not os.path.exists(ligand_src):
        return None, None, f"ligand not found: {ligand_src}"

    out_pdbqt = os.path.join(OUT_DIR, f"out_{cand_id}_5CGC.pdbqt")
    config    = os.path.join(OUT_DIR, f"config_{cand_id}_5CGC.txt")

    with open(config, "w") as f:
        f.write(f"receptor = {RECEPTOR}\n")
        f.write(f"ligand   = {ligand_src}\n")
        f.write(f"center_x = {GRID['x']}\n")
        f.write(f"center_y = {GRID['y']}\n")
        f.write(f"center_z = {GRID['z']}\n")
        f.write(f"size_x = {BOX}\nsize_y = {BOX}\nsize_z = {BOX}\n")
        f.write("exhaustiveness = 32\nnum_modes = 9\nenergy_range = 3\n")

    result = subprocess.run(
        [VINA, "--config", config, "--out", out_pdbqt],
        capture_output=True, text=True
    )
    output = result.stdout + result.stderr
    dG = parse_vina_output(output)
    if dG is None:
        return None, None, output[-500:]
    Ki = dG_to_Ki(dG)
    return dG, Ki, None

def main():
    print("Selectivity Counter-Screening: GO candidates vs mGluR5 (5CGC)")
    print(f"Grid: ({GRID['x']}, {GRID['y']}, {GRID['z']})  Box={BOX}A")
    print(f"exhaustiveness=32\n")

    results = []
    for cand in GO_CANDIDATES:
        cid = cand["id"]
        print(f"  {cid} ...", end=" ", flush=True)
        dG5, Ki5, err = run_docking(cid)

        if err:
            print(f"ERROR: {err}")
            results.append({"id": cid, "error": err})
            continue

        r4 = MGR4[cid]
        sel_idx = Ki5 / r4["Ki_uM"] if r4["Ki_uM"] > 0 else 0
        status  = "OK" if sel_idx >= 10 else ("COND" if sel_idx >= 5 else "FAIL")

        print(f"dG={dG5:.3f}  Ki(5CGC)={Ki5:.3f}uM  SI={sel_idx:.1f}x  [{status}]")
        results.append({
            "id":           cid,
            "dG_mGluR5":    dG5,
            "Ki_mGluR5_uM": Ki5,
            "dG_mGluR4":    r4["dG"],
            "Ki_mGluR4_uM": r4["Ki_uM"],
            "selectivity_index": sel_idx,
            "status":       status,
        })

    # 저장
    json_path = os.path.join(OUT_DIR, "selectivity_results.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(results, f, ensure_ascii=False, indent=2)

    # 요약
    print("\n" + "="*65)
    print(f"  {'ID':<12} {'Ki4(uM)':>9} {'Ki5(uM)':>9} {'SI':>7}  판정")
    print("="*65)
    for r in results:
        if "error" in r: continue
        print(f"  {r['id']:<12} {r['Ki_mGluR4_uM']:>9.3f} {r['Ki_mGluR5_uM']:>9.3f} {r['selectivity_index']:>7.1f}x  {r['status']}")
    print(f"\n결과: {json_path}")

if __name__ == "__main__":
    main()
