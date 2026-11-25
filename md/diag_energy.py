"""诊断：把 MM/GBSA 的总能量拆成各个力的分量，定位 +649 kcal/mol 的来源。"""
from pathlib import Path
import sys

import numpy as np
import openmm
from openmm import unit
from openmm.app import ForceField, PDBFile, NoCutoff

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mmgbsa_obc2 import (KCAL, ligand_parameters, split_protein,
                         add_ligand_to_system, make_ligand_only_system)

ROOT = Path(__file__).resolve().parent.parent
MD = ROOT / "md"


def components(system, positions, platform, label):
    for i, f in enumerate(system.getForces()):
        f.setForceGroup(i)
    ctx = openmm.Context(system, openmm.VerletIntegrator(0.001), platform)
    ctx.setPositions(positions)
    out = {}
    for i, f in enumerate(system.getForces()):
        e = ctx.getState(getEnergy=True, groups={i}).getPotentialEnergy()
        out[f.getName()] = e.value_in_unit(unit.kilojoule_per_mole) / KCAL
    del ctx
    print(f"--- {label} (kcal/mol) ---")
    for k, v in out.items():
        print(f"    {k:28s} {v:12.2f}")
    print(f"    {'TOTAL':28s} {sum(out.values()):12.2f}")
    return out


def main():
    frame = sys.argv[1] if len(sys.argv) > 1 else "901"
    prot_pdb = split_protein(frame)
    prot = PDBFile(str(prot_pdb))
    lig = ligand_parameters(frame)

    ff = ForceField("amber14-all.xml", "implicit/obc2.xml")
    system = ff.createSystem(prot.topology, nonbondedMethod=NoCutoff, constraints=None)
    nb = next(f for f in system.getForces() if isinstance(f, openmm.NonbondedForce))
    gb = next(f for f in system.getForces() if isinstance(f, openmm.CustomGBForce))
    add_ligand_to_system(system, nb, gb, lig)

    pos_prot = np.array(prot.positions.value_in_unit(unit.nanometer))
    pos_complex = np.concatenate([pos_prot, lig[0]])

    platform = openmm.Platform.getPlatformByName("Reference")
    c = components(system, pos_complex, platform, "复合物")

    rec_system = ff.createSystem(prot.topology, nonbondedMethod=NoCutoff, constraints=None)
    r = components(rec_system, pos_prot, platform, "受体")

    lig_system, lig_pos = make_ligand_only_system(lig)
    l = components(lig_system, lig_pos, platform, "游离配体")

    print("\n=== 差值 (复合物 - 受体 - 配体) ===")
    for key in set(c) | set(r) | set(l):
        print(f"    {key:28s} {c.get(key, 0) - r.get(key, 0) - l.get(key, 0):12.2f}")
    print(f"    {'ΔG':28s} {sum(c.values()) - sum(r.values()) - sum(l.values()):12.2f}")


if __name__ == "__main__":
    main()
