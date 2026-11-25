"""对多帧做 MM/GBSA，给出 ΔG 的帧间波动。"""
from pathlib import Path
import sys

import numpy as np
import openmm
from openmm import unit
from openmm.app import ForceField, PDBFile, NoCutoff

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mmgbsa_obc2 import (KCAL, ligand_parameters, split_protein,
                         add_ligand_to_system, make_ligand_only_system,
                         protein_gb_radii)

ROOT = Path(__file__).resolve().parent.parent
PLATFORM = openmm.Platform.getPlatformByName("Reference")


def parts(system, positions):
    for i, f in enumerate(system.getForces()):
        f.setForceGroup(i)
    ctx = openmm.Context(system, openmm.VerletIntegrator(0.001), PLATFORM)
    ctx.setPositions(positions)
    out = {}
    for i, f in enumerate(system.getForces()):
        out[f.getName()] = ctx.getState(getEnergy=True, groups={i}).getPotentialEnergy(
        ).value_in_unit(unit.kilojoule_per_mole) / KCAL
    del ctx
    return out


def one_frame(frame, ff, rad_tab, sca_tab):
    prot = PDBFile(str(split_protein(frame)))
    system = ff.createSystem(prot.topology, nonbondedMethod=NoCutoff, constraints=None)
    nb = next(f for f in system.getForces() if isinstance(f, openmm.NonbondedForce))
    gb = next(f for f in system.getForces() if isinstance(f, openmm.CustomGBForce))
    lig = ligand_parameters(frame, rad_tab, sca_tab)
    add_ligand_to_system(system, nb, gb, lig)

    pos_prot = np.array(prot.positions.value_in_unit(unit.nanometer))
    c = parts(system, np.concatenate([pos_prot, lig[0]]))

    rec = ff.createSystem(prot.topology, nonbondedMethod=NoCutoff, constraints=None)
    r = parts(rec, pos_prot)

    ls, lp = make_ligand_only_system(lig)
    l = parts(ls, lp)

    dn = c["NonbondedForce"] - r["NonbondedForce"] - l.get("NonbondedForce", 0)
    dgb = c["CustomGBForce"] - r["CustomGBForce"] - l.get("CustomGBForce", 0)
    return dn, dgb, dn + dgb


def main():
    frames = sys.argv[1:] or ["800", "850", "901", "950", "1001"]
    ff = ForceField("amber14-all.xml", "implicit/obc2.xml")
    prot0 = PDBFile(str(split_protein(frames[0])))
    sys0 = ff.createSystem(prot0.topology, nonbondedMethod=NoCutoff, constraints=None)
    gb0 = next(f for f in sys0.getForces() if isinstance(f, openmm.CustomGBForce))
    rad_tab, sca_tab = protein_gb_radii(prot0.topology, gb0)

    rows = []
    print(f"{'帧':>6} {'ΔE(非键)':>12} {'ΔG(GB/SA)':>12} {'ΔG':>10}")
    for fr in frames:
        dn, dgb, dg = one_frame(fr, ff, rad_tab, sca_tab)
        rows.append((fr, dn, dgb, dg))
        print(f"{fr:>6} {dn:12.2f} {dgb:12.2f} {dg:10.2f}")

    arr = np.array([[r[1], r[2], r[3]] for r in rows])
    print("-" * 46)
    print(f"{'均值':>6} {arr[:, 0].mean():12.2f} {arr[:, 1].mean():12.2f} {arr[:, 2].mean():10.2f}")
    print(f"{'SD':>6} {arr[:, 0].std(ddof=1):12.2f} {arr[:, 1].std(ddof=1):12.2f} {arr[:, 2].std(ddof=1):10.2f}")


if __name__ == "__main__":
    main()
