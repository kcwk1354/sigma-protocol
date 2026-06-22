"""
Cancer-Gate MD Simulation
방법: OpenFF Sage 2.2.0 + Gasteiger charges + 진공 MD + Kabsch RMSD
"""

import os, json, time
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

OUT_DIR = r"C:\sigma_protocol\cancer_gate\md"
os.makedirs(OUT_DIR, exist_ok=True)

CANDIDATES = {
    "CG10-005": {
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1C(C)C(F)(F)CCC1",
        "ki_uM": 0.345, "dili": 0.349, "herg": 0.259, "note": "Ki최고, hERG 경계"
    },
    "CG12-004": {
        "smiles": "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1C(C)C(F)(F)CC(O)C1",
        "ki_uM": 0.580, "dili": 0.312, "herg": 0.154, "note": "hERG최우수, DILI OK"
    },
}

FAST_STEPS = 500
FULL_STEPS = 5000
RMSD_FAST_PASS = 2.5
RMSD_FULL_PASS = 3.0


def kabsch_rmsd(P_nm, Q_nm):
    """Kabsch RMSD in Angstrom. P, Q: Nx3 float arrays in nm."""
    p = (P_nm - P_nm.mean(axis=0)) * 10.0
    q = (Q_nm - Q_nm.mean(axis=0)) * 10.0
    H = p.T @ q
    U, _, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    return float(np.sqrt(np.mean(np.sum((p @ R.T - q) ** 2, axis=1))))


def run_md(compound_id, smiles, n_steps):
    from openff.toolkit.topology import Molecule, Topology
    from openff.toolkit.typing.engines.smirnoff import ForceField
    from openmm import app, unit, LangevinMiddleIntegrator, Platform

    mol_rd = Chem.MolFromSmiles(smiles)
    mol_rd = Chem.AddHs(mol_rd)
    AllChem.EmbedMolecule(mol_rd, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(mol_rd, mmffVariant="MMFF94s")

    sdf_path = os.path.join(OUT_DIR, f"{compound_id}_ligand.sdf")
    w = Chem.SDWriter(sdf_path); w.write(mol_rd); w.close()

    mol_off = Molecule.from_file(sdf_path, allow_undefined_stereo=True)
    mol_off.assign_partial_charges("gasteiger")

    top_off = Topology.from_molecules([mol_off])
    ff_off  = ForceField("openff-2.2.0.offxml")
    interchange = ff_off.create_interchange(top_off, charge_from_molecules=[mol_off])

    omm_sys   = interchange.to_openmm(combine_nonbonded_forces=True)
    omm_top   = interchange.to_openmm_topology()
    positions = np.array(interchange.positions.magnitude)  # nm, plain numpy

    integrator = LangevinMiddleIntegrator(
        300 * unit.kelvin, 1.0 / unit.picosecond, 2.0 * unit.femtosecond
    )
    sim = app.Simulation(omm_top, omm_sys, integrator, Platform.getPlatformByName("CPU"))
    sim.context.setPositions(positions * unit.nanometer)
    sim.minimizeEnergy(maxIterations=500)

    def get_pos_nm():
        state = sim.context.getState(getPositions=True)
        pos_q = state.getPositions(asNumpy=True)
        # strip OpenMM Quantity -> plain numpy in nm
        return np.array(pos_q.value_in_unit(unit.nanometer))

    ref_pos = get_pos_nm()
    heavy   = [a.index for a in omm_top.atoms() if a.element.symbol != 'H']

    t0 = time.time()
    sim.step(n_steps)
    final_pos = get_pos_nm()
    elapsed   = time.time() - t0

    rmsd = kabsch_rmsd(final_pos[heavy], ref_pos[heavy])
    return rmsd, elapsed


def main():
    print("Cancer-Gate MD Simulation")
    print("방법: OpenFF Sage 2.2.0 / 진공 / CPU / Kabsch-aligned 중원자 RMSD")
    print()

    all_results = []
    winner = None

    for cid, info in CANDIDATES.items():
        print(f"\n{'='*55}")
        print(f"  {cid}  ({info['note']})")
        print(f"{'='*55}")
        result = {"compound_id": cid, **info}

        try:
            print(f"  [FAST {FAST_STEPS} steps] ...", end=" ", flush=True)
            rmsd_fast, t_fast = run_md(cid, info["smiles"], FAST_STEPS)
            print(f"RMSD={rmsd_fast:.2f}A  ({t_fast:.0f}s)")

            fast_pass = rmsd_fast < RMSD_FAST_PASS
            result["rmsd_fast_A"] = rmsd_fast
            result["fast_pass"]   = fast_pass

            if fast_pass:
                print(f"  [FULL {FULL_STEPS} steps] ...", end=" ", flush=True)
                rmsd_full, t_full = run_md(cid, info["smiles"], FULL_STEPS)
                print(f"RMSD={rmsd_full:.2f}A  ({t_full:.0f}s)")
                result["rmsd_full_A"] = rmsd_full
                result["stable"]      = rmsd_full < RMSD_FULL_PASS
            else:
                result["rmsd_full_A"] = None
                result["stable"]      = False

            print(f"  -> {'STABLE' if result['stable'] else 'UNSTABLE'}")

        except Exception as e:
            import traceback
            print(f"\n[ERROR] {e}")
            traceback.print_exc()
            result["error"]  = str(e)
            result["stable"] = False

        all_results.append(result)

    json_path = os.path.join(OUT_DIR, "cancer_gate_md_results.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(all_results, f, ensure_ascii=False, indent=2)

    valid = [r for r in all_results if "error" not in r]
    report_path = os.path.join(OUT_DIR, "cancer_gate_md_report.txt")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("Cancer-Gate MD 최종 판정\n" + "="*50 + "\n\n")
        f.write("방법: OpenFF Sage 2.2.0 / 진공 / Kabsch-aligned 중원자 RMSD\n\n")
        for r in valid:
            full_v = r.get("rmsd_full_A")
            full_s = f"{full_v:.2f}" if full_v is not None else "skip"
            f.write(f"{r['compound_id']}\n")
            f.write(f"  Ki={r['ki_uM']}uM  DILI={r['dili']}  hERG={r['herg']}\n")
            f.write(f"  FAST={r.get('rmsd_fast_A',0):.2f}A  FULL={full_s}A\n")
            f.write(f"  안정성: {'STABLE' if r['stable'] else 'UNSTABLE'}\n\n")

        stable = [r for r in valid if r.get("stable")]
        reason = ""
        if len(stable) == 2:
            winner = min(stable, key=lambda r: r.get("rmsd_full_A") or 99)
            reason = "RMSD 더 낮음"
        elif len(stable) == 1:
            winner = stable[0]; reason = "유일한 STABLE 후보"
        elif valid:
            winner = min(valid, key=lambda r: r["ki_uM"]); reason = "둘 다 UNSTABLE -> Ki 낮은 쪽"

        if winner:
            f.write(f"★ 최종 선택: {winner['compound_id']}\n")
            f.write(f"  근거: {reason}\n  {winner['note']}\n  Ki={winner['ki_uM']}uM\n")

    print("\n" + "="*55 + "\n  MD 완료\n" + "="*55)
    for r in valid:
        full = r.get("rmsd_full_A")
        full_str = f"  FULL={full:.2f}A" if full is not None else ""
        print(f"  {r['compound_id']}: FAST={r.get('rmsd_fast_A',0):.2f}A{full_str}  {'STABLE' if r.get('stable') else 'UNSTABLE'}")
    if winner:
        print(f"\n★ 최종 선택: {winner['compound_id']}  ({reason})")
    print(f"\n결과: {json_path}\n판정: {report_path}")


if __name__ == "__main__":
    main()
