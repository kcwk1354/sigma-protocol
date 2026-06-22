"""CG10-005 MD — OpenFF Sage 2.2.0 + Langevin 300K (GAFF2 antechamber 불필요)"""
import json
import numpy as np

SMILES = "OC(=O)C1(c2ccc(F)c(F)c2)C(F)(F)C1CN1C(C)C(F)(F)CCC1"
STEPS_FAST = 500
STEPS_FULL = 5000
OUTPUT = "C:/sigma_protocol/md_results/CG10_005_md.json"


def run_md_sage(smiles: str, steps: int) -> dict:
    try:
        from openff.toolkit import Molecule, ForceField as OFFForceField
        from openmmforcefields.generators import SMIRNOFFTemplateGenerator
        from openmm.app import ForceField, Simulation
        from openmm import LangevinMiddleIntegrator, Platform
        from openmm.unit import kelvin, picosecond, femtoseconds, nanometer, kilojoules_per_mole

        mol = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
        mol.generate_conformers(n_conformers=1)
        mol.assign_partial_charges("gasteiger")

        smirnoff = SMIRNOFFTemplateGenerator(
            molecules=[mol],
            forcefield="openff-2.2.0.offxml",
        )
        forcefield = ForceField()
        forcefield.registerTemplateGenerator(smirnoff.generator)

        topology = mol.to_topology().to_openmm()
        positions = mol.conformers[0].to_openmm()

        system = forcefield.createSystem(topology)

        integrator = LangevinMiddleIntegrator(300 * kelvin, 1 / picosecond, 2 * femtoseconds)
        platform = Platform.getPlatformByName("CPU")
        simulation = Simulation(topology, system, integrator, platform)
        simulation.context.setPositions(positions)
        simulation.minimizeEnergy(maxIterations=200)

        state0 = simulation.context.getState(getPositions=True)
        pos0 = np.array([[v.x, v.y, v.z]
                         for v in state0.getPositions(asNumpy=False)])

        simulation.context.setVelocitiesToTemperature(300 * kelvin)
        simulation.step(steps)

        state1 = simulation.context.getState(getPositions=True, getEnergy=True)
        pos1 = np.array([[v.x, v.y, v.z]
                         for v in state1.getPositions(asNumpy=False)])
        pot_energy = state1.getPotentialEnergy().value_in_unit(kilojoules_per_mole)

        # heavy-atom indices (exclude H)
        from rdkit import Chem as _Chem
        _rdmol = _Chem.MolFromSmiles(smiles)
        _rdmol = _Chem.AddHs(_rdmol)
        heavy_idx = [i for i in range(_rdmol.GetNumAtoms())
                     if _rdmol.GetAtomWithIdx(i).GetAtomicNum() != 1]

        # Kabsch-aligned RMSD (removes translation + rotation)
        def kabsch_rmsd(p, q):
            p = p - p.mean(axis=0)
            q = q - q.mean(axis=0)
            H = p.T @ q
            U, S, Vt = np.linalg.svd(H)
            d = np.linalg.det(Vt.T @ U.T)
            R = Vt.T @ np.diag([1, 1, d]) @ U.T
            p_rot = p @ R.T
            diff = p_rot - q
            return float(np.sqrt(np.mean(np.sum(diff**2, axis=1))))

        p0h = pos0[heavy_idx]
        p1h = pos1[heavy_idx]
        rmsd_heavy = kabsch_rmsd(p1h, p0h)
        rmsd_all = kabsch_rmsd(pos1, pos0)

        return {
            "steps": steps,
            "rmsd_all_nm": round(rmsd_all, 4),
            "rmsd_all_angstrom": round(rmsd_all * 10, 4),
            "rmsd_heavy_nm": round(rmsd_heavy, 4),
            "rmsd_heavy_angstrom": round(rmsd_heavy * 10, 4),
            "rmsd_nm": round(rmsd_heavy, 4),
            "rmsd_angstrom": round(rmsd_heavy * 10, 4),
            "potential_energy_kJ_mol": round(pot_energy, 2),
            "md_pass": rmsd_heavy * 10 < 2.0,  # Kabsch heavy-atom RMSD < 2.0 Å
            "status": "OK",
            "forcefield": "OpenFF Sage 2.2.0 + Langevin 300K (Kabsch-aligned)",
        }

    except Exception as e:
        import traceback
        return {
            "steps": steps,
            "status": "ERROR",
            "error": str(e),
            "traceback": traceback.format_exc()[-800:],
            "md_pass": False,
        }


if __name__ == "__main__":
    print("=== CG10-005 MD (OpenFF Sage 2.2.0 / 300K Langevin / vacuum) ===")
    print(f"SMILES: {SMILES}")
    print()

    print(f"FAST ({STEPS_FAST} steps)...")
    result_fast = run_md_sage(SMILES, STEPS_FAST)
    if result_fast["status"] == "OK":
        print(f"  RMSD: {result_fast['rmsd_angstrom']} A  |  PASS: {result_fast['md_pass']}")
    else:
        print(f"  ERROR: {result_fast['error']}")
        print(result_fast.get("traceback", ""))

    # 2.0 Å에 매우 근접하면 FULL 실행 (진공 시뮬레이션 노이즈 고려)
    run_full = result_fast.get("md_pass") or (
        result_fast.get("status") == "OK" and result_fast.get("rmsd_heavy_angstrom", 99) < 2.5
    )
    if run_full:
        print(f"\nFULL ({STEPS_FULL} steps)...")
        result_full = run_md_sage(SMILES, STEPS_FULL)
        if result_full["status"] == "OK":
            print(f"  RMSD (heavy): {result_full['rmsd_heavy_angstrom']} A  all: {result_full['rmsd_all_angstrom']} A")
            print(f"  Energy: {result_full['potential_energy_kJ_mol']} kJ/mol")
            print(f"  PASS: {result_full['md_pass']}")
        else:
            print(f"  ERROR: {result_full['error']}")
    else:
        result_full = {"md_pass": False, "note": "FAST RMSD too high - FULL skipped"}
        print("\nFAST RMSD too high - FULL skipped")

    output = {
        "candidate": "CG10-005",
        "smiles": SMILES,
        "vina_dG_best": -8.811,
        "vina_Ki_uM_best": 0.345,
        "forcefield_note": "OpenFF Sage 2.2.0 (SMIRNOFF), Gasteiger charges, vacuum, 300K",
        "md_fast": result_fast,
        "md_full": result_full,
        "final_md_pass": result_full.get("md_pass", False),
    }

    with open(OUTPUT, "w", encoding="utf-8") as f:
        json.dump(output, f, ensure_ascii=False, indent=2)

    status = "PASS" if output["final_md_pass"] else ("FAST_ONLY_PASS" if result_fast.get("md_pass") else "FAIL")
    print(f"\nSaved: {OUTPUT}")
    print(f"MD result: {status}")
