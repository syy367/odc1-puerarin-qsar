"""单帧 MM/GBSA 验证：OpenMM(amber14 + GBn2) + OpenFF Sage 配体参数。

思路：MD 轨迹给的是 GROMACS/AMBER14SB 的坐标，这里不重跑 MD，
只把某一帧的坐标拿来做气相分子力学能 + GB 隐式溶剂能，按
    ΔG_bind ≈ E(复合物) - E(受体) - E(配体)
做能量分解。

已知偏差（必须写进论文方法）：
  * 配体电荷用 Gasteiger（MD 阶段是 ORCA/B3LYP-D3 单点 + RESP）；
  * 未做熵校正、未加非极性 SASA 表面项（GBn2 自带极性部分）；
  * 单帧，无构象平均。

用法：python md/mmgbsa.py 901
"""
from pathlib import Path
import sys
import time

import numpy as np
from openmm import app, unit
from openmmforcefields.generators import SystemGenerator
from openff.toolkit import Molecule as OffMolecule
from openff.units import unit as off_unit
from rdkit import Chem
from rdkit.Chem import AllChem

# 本机装不到 openff 依赖的 constraint（只在 conda-forge 发布），而它只用于给 RDKit
# 分子指派双键立体化学；MM/GBSA 只看键级与坐标，不依赖 E/Z 标记，这里直接跳过。
from openff.toolkit.utils import rdkit_wrapper as _rdkit_wrapper

_rdkit_wrapper.RDKitToolkitWrapper._assign_rdmol_bonds_stereo = lambda self, *a, **k: None

ROOT = Path(__file__).resolve().parent.parent
MD = ROOT / "md"
OFFXML = MD / "forcefields" / "openff_unconstrained-2.2.1.offxml"

STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HIE", "HID", "HIP", "CYX",
}


def split_frame(frame: str):
    """把一帧拆成蛋白部分（去掉配体与离子）与配体部分。"""
    src = MD / "frames" / f"frame_{frame}.pdb"
    prot_lines, lig_lines = [], []
    for line in src.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        resname = line[17:20].strip()
        if resname == "MOL":
            lig_lines.append(line)
        elif resname in STANDARD_AA:
            prot_lines.append(line)
    prot = MD / "frames" / f"protein_{frame}.pdb"
    prot.write_text("\n".join(prot_lines) + "\nEND\n")
    return prot, lig_lines


def merge_topologies(prot_top, lig_top):
    """把配体残基追加到蛋白拓扑末尾（OpenMM 没有内置 merge）。

    注意必须连键一起搬过去：只加原子的话 OpenMM 会认为残基内部没有键，
    模板匹配直接失败（"The set of atoms matches NPHE, but the residue has no bonds"）。
    """
    merged = app.Topology()
    box = prot_top.getPeriodicBoxVectors()
    if box is not None:
        merged.setPeriodicBoxVectors(*box)

    atom_map = {}

    def copy_into(src_top, chain_id, rename=None):
        for chain in src_top.chains():
            new_chain = merged.addChain(chain.id if rename is None else chain_id)
            for res in chain.residues():
                new_res = merged.addResidue(res.name if rename is None else rename,
                                            new_chain, res.id)
                for atom in res.atoms():
                    atom_map[atom] = merged.addAtom(atom.name, atom.element, new_res)
        for bond in src_top.bonds():
            merged.addBond(atom_map[bond[0]], atom_map[bond[1]])

    copy_into(prot_top, None)
    copy_into(lig_top, "L", rename="LIG")
    return merged


def build_ligand_molecule(frame: str):
    sdf_text = (MD / "frames" / f"lig_{frame}_fix.sdf").read_text()
    rd_mol = Chem.MolFromMolBlock(sdf_text, removeHs=False, sanitize=True)
    if rd_mol is None:
        raise SystemExit("配体 SDF 解析失败")
    # openff 的 toolkit registry 在本机缺 constraint 包，Gasteiger 走 RDKit 直接算
    AllChem.ComputeGasteigerCharges(rd_mol)
    charges = np.array(
        [float(a.GetProp("_GasteigerCharge")) for a in rd_mol.GetAtoms()]
    )
    charges = np.nan_to_num(charges, nan=0.0, posinf=0.0, neginf=0.0)
    charges -= charges.mean()          # 保证净电荷为整数（中性分子）
    off_mol = OffMolecule.from_rdkit(rd_mol, allow_undefined_stereo=True)
    off_mol.partial_charges = charges * off_unit.elementary_charge
    return off_mol


def main():
    frame = sys.argv[1] if len(sys.argv) > 1 else "901"
    t0 = time.time()

    prot_pdb, _ = split_frame(frame)
    print(f"[1] 蛋白部分写出 {prot_pdb.relative_to(ROOT)} "
          f"({len(prot_pdb.read_text().splitlines()) - 1} 行)")

    off_mol = build_ligand_molecule(frame)
    print(f"[2] 配体 OpenFF 分子：{off_mol.n_atoms} 原子，净电荷 {off_mol.total_charge}")

    generator = SystemGenerator(
        forcefields=["amber14-all.xml", "implicit/gbn2.xml"],
        small_molecule_forcefield=str(OFFXML),
        forcefield_kwargs={
            "constraints": None,
            "rigidWater": False,
            "implicitSolvent": app.GBn2,
        },
        nonperiodic_forcefield_kwargs={"nonbondedMethod": app.NoCutoff},
    )

    prot = app.PDBFile(str(prot_pdb))
    print(f"[3] 蛋白拓扑：{prot.topology.getNumAtoms()} 原子 / "
          f"{prot.topology.getNumResidues()} 残基 / {prot.topology.getNumChains()} 链")

    lig_top_off = off_mol.to_topology()
    lig_top = lig_top_off.to_openmm()
    lig_pos = lig_top_off.get_positions().to("nanometer")

    complex_top = merge_topologies(prot.topology, lig_top)
    n_prot = prot.topology.getNumAtoms()
    print(f"[4] 复合物拓扑：{complex_top.getNumAtoms()} 原子 "
          f"(蛋白 {n_prot} + 配体 {complex_top.getNumAtoms() - n_prot})")

    system_complex = generator.create_system(complex_top, molecules=[off_mol])
    system_rec = generator.create_system(prot.topology)
    system_lig = generator.create_system(lig_top, molecules=[off_mol])
    print(f"[5] 三个 system 建好：复合物 {system_complex.getNumParticles()} 原子 / "
          f"受体 {system_rec.getNumParticles()} / 配体 {system_lig.getNumParticles()}"
          f"   ({time.time() - t0:.0f} s)")

    platform = app.Platform.getPlatformByName("Reference")

    def energy(system, positions):
        ctx = app.Context(system, app.VerletIntegrator(0.001), platform)
        ctx.setPositions(positions)
        return ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
            unit.kilocalorie_per_mole)

    pos_complex = prot.positions + lig_pos
    e_complex = energy(system_complex, pos_complex)
    print(f"    E(complex) = {e_complex:12.1f} kcal/mol")
    e_rec = energy(system_rec, prot.positions)
    print(f"    E(receptor) = {e_rec:12.1f} kcal/mol")
    e_lig = energy(system_lig, lig_pos)
    print(f"    E(ligand)  = {e_lig:12.1f} kcal/mol")

    dg = e_complex - e_rec - e_lig
    print()
    print(f"[6] ΔG(GB) = {dg:.1f} kcal/mol   （单帧、无熵、无 SA 表面项）")
    print(f"    总耗时 {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
