"""
Epsilon-Gate MD — OpenMM 진공 분자동역학
=========================================
도킹 통과 분자의 구조 안정성 검증

실행 시간 (CPU 기준):
  10ps  = ~30초/분자   ← 기본값 (이 환경)
  100ps = ~5분/분자    ← 코랩 CPU
  1ns   = ~50분/분자   ← 코랩 GPU 권장

타임스케일 변경:
  N_STEPS 상수 하나만 수정하면 됨
  10ps  → N_STEPS = 20_000
  100ps → N_STEPS = 200_000
  1ns   → N_STEPS = 2_000_000

설치:
  pip install rdkit openmm pandas

사용:
  python epsilon_gate_md.py
  또는 도킹 결과 CSV 지정:
  python epsilon_gate_md.py --csv gen2_docking_results.csv
"""

import os, sys, tempfile, argparse
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdmolfiles import MolToPDBBlock
import openmm as mm
from openmm import app, unit

# ──────────────────────────────────────────
# 타임스케일 설정 (여기만 바꾸면 됨)
# ──────────────────────────────────────────
TIMESTEP_FS  = 0.5           # 타임스텝 (fs)
N_STEPS      = 20_000        # 10ps  ← 이 숫자만 바꾸면 스케일 변경
#N_STEPS     = 200_000       # 100ps (코랩 CPU)
#N_STEPS     = 2_000_000     # 1ns   (코랩 GPU)
REPORT_EVERY = 2_000         # 로그 출력 간격
TEMPERATURE  = 300           # K

def ps_total():
    return N_STEPS * TIMESTEP_FS / 1000

# ──────────────────────────────────────────
# MMFF → OpenMM System 변환
# ──────────────────────────────────────────
def build_system(mol):
    """RDKit MMFF94 파라미터를 OpenMM System으로 변환"""
    mmff_props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94')
    if mmff_props is None:
        raise ValueError("MMFF 파라미터 생성 실패 (분자 확인 필요)")

    system = mm.System()
    n = mol.GetNumAtoms()

    # 파티클 질량
    for i in range(n):
        system.addParticle(mol.GetAtomWithIdx(i).GetMass() * unit.dalton)

    # Bond stretching
    bf = mm.HarmonicBondForce()
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        p = mmff_props.GetMMFFBondStretchParams(mol, i, j)
        kb = (p[1] if p else 300) * 4.184 * 2 * 100  # kcal/Å² → kJ/nm²
        r0 = (p[2] if p else 1.5) / 10                # Å → nm
        bf.addBond(i, j, r0*unit.nanometer,
                   kb*unit.kilojoule_per_mole/unit.nanometer**2)
    system.addForce(bf)

    # Angle bending
    af = mm.HarmonicAngleForce()
    for atom in mol.GetAtoms():
        j = atom.GetIdx()
        nbs = [a.GetIdx() for a in atom.GetNeighbors()]
        for ii in range(len(nbs)):
            for kk in range(ii+1, len(nbs)):
                i, k = nbs[ii], nbs[kk]
                p = mmff_props.GetMMFFAngleBendParams(mol, i, j, k)
                ka = (p[1] if p else 100) * 4.184 * 2
                th = np.radians(p[2] if p else 109.5)
                af.addAngle(i, j, k, th,
                            ka*unit.kilojoule_per_mole/unit.radian**2)
    system.addForce(af)

    # Torsions
    tf = mm.PeriodicTorsionForce()
    for bond in mol.GetBonds():
        j, k = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        jnb = [a.GetIdx() for a in mol.GetAtomWithIdx(j).GetNeighbors() if a.GetIdx()!=k]
        knb = [a.GetIdx() for a in mol.GetAtomWithIdx(k).GetNeighbors() if a.GetIdx()!=j]
        for i in jnb:
            for l in knb:
                p = mmff_props.GetMMFFTorsionParams(mol, i, j, k, l)
                if p:
                    for per in [1, 2, 3]:
                        v = abs(p[per-1]) * 4.184 / 2
                        if v > 0.01:
                            phase = 0.0 if per % 2 == 1 else np.pi
                            tf.addTorsion(i, j, k, l, per, phase,
                                         v*unit.kilojoule_per_mole)
    system.addForce(tf)

    # Nonbonded (LJ, 간소화)
    nbf = mm.NonbondedForce()
    nbf.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
    vdw = {'C':(0.34,0.36),'N':(0.31,0.30),'O':(0.30,0.21),
           'H':(0.24,0.03),'F':(0.29,0.25),'S':(0.36,0.43)}
    for i in range(n):
        sym = mol.GetAtomWithIdx(i).GetSymbol()
        r, e = vdw.get(sym, (0.34, 0.36))
        nbf.addParticle(0.0,
                        r/10*unit.nanometer,
                        e*4.184*unit.kilojoule_per_mole)

    bonded = set()
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        nbf.addException(i, j, 0.0, 1.0, 0.0)
        bonded.add((min(i,j), max(i,j)))
    for atom in mol.GetAtoms():
        j = atom.GetIdx()
        nbs = [a.GetIdx() for a in atom.GetNeighbors()]
        for ii in range(len(nbs)):
            for kk in range(ii+1, len(nbs)):
                i, k = nbs[ii], nbs[kk]
                pair = (min(i,k), max(i,k))
                if pair not in bonded:
                    nbf.addException(i, k, 0.0, 1.0, 0.0)
                    bonded.add(pair)
    system.addForce(nbf)
    return system

# ──────────────────────────────────────────
# 단일 분자 MD 실행
# ──────────────────────────────────────────
def run_md(mol_id, smiles, verbose=True):
    """
    한 분자의 MD 실행
    반환: {'rmsd_final', 'rmsd_mean', 'rmsd_max', 'stable', 'traj'}
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)

    pdb_block = MolToPDBBlock(mol)
    with tempfile.NamedTemporaryFile(suffix='.pdb', mode='w', delete=False) as f:
        f.write(pdb_block); pdb_path = f.name

    try:
        pdb = app.PDBFile(pdb_path)
        system = build_system(mol)
        integrator = mm.LangevinMiddleIntegrator(
            TEMPERATURE*unit.kelvin,
            10/unit.picosecond,
            TIMESTEP_FS*0.001*unit.picoseconds
        )
        sim = app.Simulation(pdb.topology, system, integrator,
                             mm.Platform.getPlatformByName('CPU'))
        sim.context.setPositions(pdb.positions)

        # 에너지 최소화
        sim.minimizeEnergy(maxIterations=1000)
        state0 = sim.context.getState(getPositions=True)
        pos0 = np.array(state0.getPositions().value_in_unit(unit.angstrom))

        # 점진적 가열 (100→200→300K)
        for T in [100, 200, TEMPERATURE]:
            integrator.setTemperature(T*unit.kelvin)
            sim.step(500)

        # MD 실행 + 트래킹
        traj = []
        for step_i in range(N_STEPS // REPORT_EVERY):
            sim.step(REPORT_EVERY)
            st = sim.context.getState(getPositions=True, getEnergy=True)
            pos = np.array(st.getPositions().value_in_unit(unit.angstrom))
            rmsd = float(np.sqrt(np.mean(np.sum((pos-pos0)**2, axis=1))))
            pe   = st.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
            ke   = st.getKineticEnergy().value_in_unit(unit.kilocalories_per_mole)
            t_ps = (step_i+1) * REPORT_EVERY * TIMESTEP_FS / 1000
            traj.append({'t_ps': round(t_ps,2), 'RMSD': round(rmsd,3),
                         'PE': round(pe,2), 'KE': round(ke,2)})
            if verbose:
                print(f"    {t_ps:.1f}ps  RMSD={rmsd:.2f}Å  "
                      f"PE={pe:.1f}  KE={ke:.1f} kcal/mol")

        rmsds = [r['RMSD'] for r in traj]
        # 안정 기준: 마지막 절반 RMSD 평균 < 3.0 Å
        late_rmsd = rmsds[len(rmsds)//2:]
        stable = np.mean(late_rmsd) < 3.0

        return {
            'rmsd_final': rmsds[-1],
            'rmsd_mean':  round(float(np.mean(late_rmsd)), 3),
            'rmsd_max':   round(float(max(rmsds)), 3),
            'stable':     stable,
            'traj':       traj,
        }
    finally:
        os.unlink(pdb_path)

# ──────────────────────────────────────────
# 안정성 판정
# ──────────────────────────────────────────
def grade_stability(result):
    if result is None:
        return "측정불가"
    m = result['rmsd_mean']
    if m < 1.5:   return "⭐⭐⭐ 매우 안정"
    elif m < 2.5: return "⭐⭐  안정"
    elif m < 3.5: return "⭐   보통"
    else:         return "△   불안정"

# ──────────────────────────────────────────
# 메인
# ──────────────────────────────────────────
def main(csv_path="gen2_docking_results.csv", top_n=3):
    print("\n" + "="*65)
    print(f"  Epsilon-Gate MD — OpenMM 진공 MD ({ps_total():.0f}ps)")
    print(f"  {N_STEPS:,} 스텝 × {TIMESTEP_FS}fs  |  {TEMPERATURE}K")
    print(f"  ※ N_STEPS 변경 → 100ps(200,000) / 1ns(2,000,000)")
    print("="*65)

    # 도킹 결과에서 상위 분자 선택
    if os.path.exists(csv_path):
        df = pd.read_csv(csv_path).sort_values('Docking_dG').head(top_n)
    else:
        # 도킹 결과 없으면 gen2_library에서 QED 상위
        lib = os.path.join(os.path.dirname(csv_path), 'gen2_library.csv')
        df = pd.read_csv(lib)
        df = df[df['STATUS']=='PASS'].nlargest(top_n, 'QED')

    results = []
    for _, row in df.iterrows():
        mol_id = row['ID']
        smiles = row['SMILES']
        dg = row.get('Docking_dG', '—')
        print(f"\n  [{mol_id}]")
        print(f"  도킹 ΔG: {dg}  →  {ps_total():.0f}ps MD 실행")
        print("  " + "-"*55)

        md = run_md(mol_id, smiles, verbose=True)
        grade = grade_stability(md)

        print(f"\n  → {grade}  |  RMSD(평균)={md['rmsd_mean']:.2f}Å  "
              f"최대={md['rmsd_max']:.2f}Å")

        results.append({
            'ID':          mol_id,
            'Docking_dG':  dg,
            'RMSD_mean':   md['rmsd_mean'],
            'RMSD_max':    md['rmsd_max'],
            'RMSD_final':  md['rmsd_final'],
            'Stable':      md['stable'],
            'Grade':       grade,
            'MD_ps':       ps_total(),
        })

    # 최종 결과
    res_df = pd.DataFrame(results)
    print("\n" + "="*65)
    print("  🏆 종합 순위 (도킹 × MD 안정성)")
    print("="*65)
    print(f"  {'순위':<4} {'ID':<35} {'ΔG':>7}  {'RMSD':>6}  {'판정'}")
    print("  " + "-"*65)
    for rank, (_, r) in enumerate(res_df.iterrows(), 1):
        print(f"  {rank:<4} {r['ID']:<35} {str(r['Docking_dG']):>7}  "
              f"{r['RMSD_mean']:>5.2f}Å  {r['Grade']}")

    out = os.path.join(os.path.dirname(csv_path), 'gen2_md_results.csv')
    res_df.to_csv(out, index=False)
    print(f"\n  결과 저장: {out}")

    best = res_df.sort_values('RMSD_mean').iloc[0]
    print(f"\n  📌 Gen3 힌트: {best['ID']}")
    print(f"     도킹 {best['Docking_dG']} kcal/mol  |  MD RMSD {best['RMSD_mean']:.2f}Å")
    print(f"\n  ⚠️  진공 MD → 용매 효과 미포함")
    print(f"     실제 검증: 코랩 + 명시적 수분자(TIP3P) + N_STEPS=2,000,000\n")
    return res_df

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--csv', default='gen2_docking_results.csv')
    parser.add_argument('--top', type=int, default=3)
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    csv = os.path.join(script_dir, args.csv)
    main(csv_path=csv, top_n=args.top)
