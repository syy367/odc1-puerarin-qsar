"""单帧 MM/GBSA（OBC2, igb=5）——不依赖 OpenFF 的实现。

为什么不用 OpenFF：openff-forcefields / openff-interchange / openff-units 在 PyPI 上
要么缺失要么全部被 yank，Python 3.13 + Windows 装不出可用的 SMIRNOFF 栈。
但 MM/GBSA 的分解式

    ΔG = ΔE_vdW + ΔE_ele + ΔG_polar(GB) + ΔG_nonpolar(SA)

中每一项都是"复合物 − 受体 − 配体"的差，配体的键合项（bond/angle/torsion）
在复合物和游离配体中完全相同、自动抵消。因此配体只需要非键参数：
电荷、LJ（sigma/epsilon）、GB 半径与 scale。

做法：
  * 蛋白走 OpenMM 自带 amber14-all + implicit/obc2.xml（模板来自 PDB 残基名）；
  * 配体按元素/杂化给 GAFF 量级的 LJ、mbondi2 半径、OBC2 的 scale；
  * 电荷用 RDKit Gasteiger（MD 阶段是 ORCA/B3LYP-D3 + RESP，这里有偏差，需在论文写明）。

用法：python md/mmgbsa_obc2.py 901
"""
from pathlib import Path
import sys
import time

import numpy as np
import openmm
from openmm import unit, Platform
from openmm.app import ForceField, PDBFile, NoCutoff
from openmm.app.internal.customgbforces import GBSAOBC2Force
from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = Path(__file__).resolve().parent.parent
MD = ROOT / "md"

KCAL = 4.184                      # kJ per kcal
STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HIE", "HID", "HIP", "CYX",
}

# GAFF 量级的元素级 LJ 参数：(sigma / nm, epsilon / kJ·mol⁻¹)
LJ = {
    "C.ar": (0.1908, 0.086 * KCAL),     # 芳香碳 ca
    "C.sp2": (0.1908, 0.086 * KCAL),    # 羰基/烯碳 c2
    "C.sp3": (0.1908, 0.1094 * KCAL),   # 饱和碳 c3
    "O.oh": (0.1721, 0.2104 * KCAL),    # 羟基氧 oh
    "O.os": (0.1684, 0.1700 * KCAL),    # 醚/糖环氧 os
    "O.o": (0.1661, 0.2100 * KCAL),     # 羰基氧 o
    "H.pol": (0.0000, 0.0000),          # 连 O/N 的极性氢 HO
    "H.np": (0.1487, 0.0157 * KCAL),    # 非极性氢 HC
}

# mbondi2 半径（OBC2 标准）与 OBC2 scale
GB_RADIUS = {"H": 0.12, "C": 0.17, "N": 0.155, "O": 0.15, "S": 0.18,
             "P": 0.185, "F": 0.15, "CL": 0.17, "BR": 0.185, "I": 0.198}


def classify(atom: Chem.Atom) -> str:
    """把 RDKit 原子映射到上面的 LJ 类别。"""
    sym = atom.GetSymbol()
    if sym == "C":
        if atom.GetIsAromatic():
            return "C.ar"
        hyb = atom.GetHybridization()
        if hyb == Chem.HybridizationType.SP2:
            return "C.sp2"
        return "C.sp3"
    if sym == "O":
        neigh = list(atom.GetNeighbors())
        if len(neigh) == 2:
            return "O.os"
        bond = atom.GetBonds()[0]
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            return "O.o"
        return "O.oh"
    if sym == "H":
        neigh = list(atom.GetNeighbors())
        if neigh and neigh[0].GetSymbol() in ("O", "N", "S"):
            return "H.pol"
        return "H.np"
    if sym == "N":
        return "C.sp3"      # 配体里没有 N，兜底用不到
    raise ValueError(f"未覆盖的元素 {sym}")


def protein_gb_radii(topology, gb):
    """从蛋白自带的 OBC2 参数里统计每个元素的 radius / scale 众数。

    自己猜 mbondi2 半径会和 OpenMM 蛋白侧那套（OBC2 专用半径）对不上，
    配体原子半径偏大或偏小都会让 GB 自能在结合前后剧变（实测差 640 kcal/mol）。
    """
    from collections import Counter
    radii, scales = {}, {}
    for i, atom in enumerate(topology.atoms()):
        p = gb.getParticleParameters(i)
        sym = atom.element.symbol.upper()
        radii.setdefault(sym, []).append(round(p[1], 5))
        scales.setdefault(sym, []).append(round(p[2], 3))
    return ({k: Counter(v).most_common(1)[0][0] for k, v in radii.items()},
            {k: Counter(v).most_common(1)[0][0] for k, v in scales.items()})


def ligand_parameters(frame: str, radius_table=None, scale_table=None):
    """返回配体的 (坐标 nm, 质量, 电荷 e, sigma nm, epsilon kJ, GB 半径 nm, scale)。"""
    sdf = (MD / "frames" / f"lig_{frame}_fix.sdf").read_text()
    mol = Chem.MolFromMolBlock(sdf, removeHs=False, sanitize=True)
    if mol is None:
        raise SystemExit("配体 SDF 解析失败")
    AllChem.ComputeGasteigerCharges(mol)
    q = np.array([float(a.GetProp("_GasteigerCharge")) for a in mol.GetAtoms()])
    q = np.nan_to_num(q, nan=0.0, posinf=0.0, neginf=0.0)
    q -= q.mean()

    conf = mol.GetConformer()
    pos = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
                     conf.GetAtomPosition(i).z] for i in range(mol.GetNumAtoms())])
    pos = pos / 10.0                                # Å -> nm

    mass, sigma, eps, gbr, scale = [], [], [], [], []
    radius_table = radius_table or {}
    scale_table = scale_table or {}
    for atom in mol.GetAtoms():
        mass.append(atom.GetMass())
        s, e = LJ[classify(atom)]
        sigma.append(s)
        eps.append(e)
        sym = atom.GetSymbol().upper()
        # 优先沿用蛋白侧同一元素在 OBC2 下的取值，保证 GB 参数口径一致
        gbr.append(radius_table.get(sym, GB_RADIUS[sym]))
        if sym == "H":
            heavy = atom.GetNeighbors()[0].GetSymbol().upper()
            gbr[-1] = radius_table.get("H", gbr[-1])
            scale.append(scale_table.get("H", 0.85 if heavy in ("O", "N") else 0.8))
        else:
            scale.append(scale_table.get(sym, 0.8))
    return (pos, np.array(mass), q, np.array(sigma), np.array(eps),
            np.array(gbr), np.array(scale), mol)


def split_protein(frame: str) -> Path:
    src = MD / "frames" / f"frame_{frame}.pdb"
    keep = [ln for ln in src.read_text().splitlines()
            if ln.startswith(("ATOM", "HETATM")) and ln[17:20].strip() in STANDARD_AA]
    out = MD / "frames" / f"protein_{frame}.pdb"
    out.write_text("\n".join(keep) + "\nEND\n")
    return out


def add_ligand_to_system(system, nb, gb, lig_params):
    pos, mass, q, sigma, eps, gbr, scale, _ = lig_params
    n0 = system.getNumParticles()
    for i in range(len(mass)):
        system.addParticle(float(mass[i]))
        nb.addParticle(float(q[i]), float(sigma[i]), float(eps[i]))
        gb.addParticle([float(q[i]), float(gbr[i]), float(scale[i])])
    # 配体内部对全部设为 0 例外：键合项本来就没加，内部非键项在复合物与游离态
    # 完全一致必然抵消，但成键原子对按正常 LJ 算会产生 +900 kcal/mol 的虚假排斥，
    # 数值上会淹没真实信号，这里直接置零（两边同步置零，抵消关系不变）。
    for i in range(n0, n0 + len(mass)):
        for j in range(i + 1, n0 + len(mass)):
            nb.addException(i, j, 0.0, 1.0, 0.0)


def make_ligand_only_system(lig_params):
    pos, mass, q, sigma, eps, gbr, scale, _ = lig_params
    system = openmm.System()
    for m in mass:
        system.addParticle(float(m))
    nb = openmm.NonbondedForce()
    nb.setNonbondedMethod(openmm.NonbondedForce.NoCutoff)
    nb.setReactionFieldDielectric(1.0)
    gb = GBSAOBC2Force(soluteDielectric=1.0, solventDielectric=78.5, SA="ACE")
    for i in range(len(mass)):
        nb.addParticle(float(q[i]), float(sigma[i]), float(eps[i]))
        gb.addParticle([float(q[i]), float(gbr[i]), float(scale[i])])
    for i in range(len(mass)):
        for j in range(i + 1, len(mass)):
            nb.addException(i, j, 0.0, 1.0, 0.0)
    gb.finalize()
    gb.setNonbondedMethod(openmm.CustomNonbondedForce.NoCutoff)
    system.addForce(nb)
    system.addForce(gb)
    return system, pos


def energy(system, positions, platform):
    ctx = openmm.Context(system, openmm.VerletIntegrator(0.001), platform)
    ctx.setPositions(positions)
    e = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    del ctx
    return e


def main():
    frame = sys.argv[1] if len(sys.argv) > 1 else "901"
    t0 = time.time()

    prot_pdb = split_protein(frame)
    prot = PDBFile(str(prot_pdb))
    n_prot = prot.topology.getNumAtoms()
    print(f"[1] 蛋白 {n_prot} 原子 / {prot.topology.getNumResidues()} 残基")

    ff = ForceField("amber14-all.xml", "implicit/obc2.xml")
    system = ff.createSystem(prot.topology, nonbondedMethod=NoCutoff, constraints=None)
    nb = next(f for f in system.getForces() if isinstance(f, openmm.NonbondedForce))
    gb = next(f for f in system.getForces() if isinstance(f, openmm.CustomGBForce))

    rad_tab, sca_tab = protein_gb_radii(prot.topology, gb)
    print(f"[2] 蛋白侧 OBC2 元素半径：{ {k: round(v, 4) for k, v in rad_tab.items()} }")
    lig = ligand_parameters(frame, rad_tab, sca_tab)
    print(f"    配体 {len(lig[1])} 原子，净电荷 {lig[2].sum():+.3f} e")
    print(f"[3] 蛋白 system：{system.getNumParticles()} 粒子，"
          f"GB 力场 {gb.getNumParticles()} 粒子 ({time.time() - t0:.0f} s)")

    add_ligand_to_system(system, nb, gb, lig)
    print(f"[4] 加入配体后 {system.getNumParticles()} 粒子")

    pos_prot = np.array(prot.positions.value_in_unit(unit.nanometer))
    pos_complex = np.concatenate([pos_prot, lig[0]])

    lig_system, lig_pos = make_ligand_only_system(lig)
    print(f"[5] 游离配体 system：{lig_system.getNumParticles()} 粒子")

    platform = Platform.getPlatformByName("Reference")

    def parts(sys_, pos, label):
        for i, f in enumerate(sys_.getForces()):
            f.setForceGroup(i)
        ctx = openmm.Context(sys_, openmm.VerletIntegrator(0.001), platform)
        ctx.setPositions(pos)
        res = {}
        for i, f in enumerate(sys_.getForces()):
            key = f.getName()
            res[key] = res.get(key, 0.0) + ctx.getState(
                getEnergy=True, groups={i}).getPotentialEnergy().value_in_unit(
                    unit.kilojoule_per_mole) / KCAL
        del ctx
        print(f"    {label:10s} 非键 {res.get('NonbondedForce', 0):10.2f}   "
              f"GB/SA {res.get('CustomGBForce', 0):10.2f}   "
              f"键合 {sum(v for k, v in res.items() if 'Bond' in k or 'Angle' in k or 'Torsion' in k):9.2f}")
        return res

    c = parts(system, pos_complex, "复合物")
    rec_system = ff.createSystem(prot.topology, nonbondedMethod=NoCutoff, constraints=None)
    r = parts(rec_system, pos_prot, "受体")
    l = parts(lig_system, lig_pos, "游离配体")

    dn = c["NonbondedForce"] - r["NonbondedForce"] - l.get("NonbondedForce", 0)
    dgb = c["CustomGBForce"] - r["CustomGBForce"] - l.get("CustomGBForce", 0)
    print()
    print(f"[6] ΔE(非键) = {dn:8.2f} kcal/mol        （vdW + 静电）")
    print(f"    ΔG(GB/SA) = {dgb:8.2f} kcal/mol        （极性去溶剂化 + 非极性）")
    print(f"    ΔG(MM/GBSA) = {dn + dgb:8.2f} kcal/mol   单帧，无熵校正，总耗时 {time.time() - t0:.0f} s")

    print(f"    总耗时 {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
