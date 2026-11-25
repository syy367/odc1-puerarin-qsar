"""从 MD 轨迹帧中提取葛根素配体，修正键级并补齐缺失的氢，输出 SDF。

轨迹里的配体只有 36 个原子（21C + 9O + 6 个羟基 H），缺少 14 个非极性氢，
且 PDB 无 CONECT，OpenBabel 直接感知键级会把糖环误判成芳香体系。
这里改用"参考 SMILES + 重原子图匹配"把正确键级搬到轨迹坐标上。

踩过的坑：
  * RDKit 的 C++ 层打不开含中文的路径，必须走 MolFromPDBBlock/MolToMolBlock；
  * 葛根素是 8-C-葡萄糖基而非 6-C，SMILES 里糖的位置写错会导致图匹配失败；
  * AssignBondOrdersFromTemplate 前要先把目标原子设成 NoImplicit，否则隐形氢会造成价态冲突。

用法：python md/prep_ligand.py [帧号，默认 901]
"""
from pathlib import Path
import sys

from rdkit import Chem
from rdkit.Chem import AllChem, rdDetermineBonds
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

ROOT = Path(__file__).resolve().parent.parent
MD = ROOT / "md"

# 葛根素：7,4'-二羟基-8-C-葡萄糖基异黄酮
SMILES = "O=C1C(=COc2c([C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)c(O)ccc12)c1ccc(O)cc1"
EXPECTED_FORMULA = "C21H20O9"


def build(frame: str) -> Chem.Mol:
    pdb = MD / "frames" / f"lig_{frame}.pdb"
    ref = Chem.MolFromSmiles(SMILES)
    if ref is None:
        raise SystemExit("SMILES 解析失败")

    traj = Chem.MolFromPDBBlock(pdb.read_text(), removeHs=False, sanitize=False)
    if traj is None:
        raise SystemExit(f"无法解析 {pdb}")
    heavy = Chem.RemoveHs(traj)
    if heavy.GetNumAtoms() != ref.GetNumAtoms():
        raise SystemExit(f"重原子数不符：轨迹 {heavy.GetNumAtoms()} vs 参考 {ref.GetNumAtoms()}")

    work = Chem.RWMol(heavy)
    rdDetermineBonds.DetermineConnectivity(work, charge=0)
    work = work.GetMol()
    for atom in work.GetAtoms():          # 关键：关掉隐式氢，否则键级指派会冲突
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(0)

    try:
        fixed = AllChem.AssignBondOrdersFromTemplate(ref, work)
    except Exception as exc:  # noqa: BLE001
        raise SystemExit(f"键级指派失败（SMILES 与轨迹拓扑不匹配）：{exc}") from exc

    # 只保留元素重建成干净分子，让 RDKit 按价态重新补氢（原对象带着 NoImplicit 标记）
    rebuilt = Chem.RWMol()
    for atom in fixed.GetAtoms():
        new_atom = Chem.Atom(atom.GetAtomicNum())
        new_atom.SetNoImplicit(False)
        rebuilt.AddAtom(new_atom)
    for bond in fixed.GetBonds():
        rebuilt.AddBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondType())
    mol = rebuilt.GetMol()
    mol.AddConformer(fixed.GetConformer(), assignId=True)
    Chem.SanitizeMol(mol)
    mol_h = Chem.AddHs(mol, addCoords=True)   # addCoords：在轨迹重原子坐标上补氢
    Chem.AssignStereochemistryFrom3D(mol_h)

    formula = CalcMolFormula(mol_h)
    if formula != EXPECTED_FORMULA:
        raise SystemExit(f"补氢后分子式不符：{formula}（应为 {EXPECTED_FORMULA}）")
    return mol_h


def main():
    frame = sys.argv[1] if len(sys.argv) > 1 else "901"
    out = MD / "frames" / f"lig_{frame}_fix.sdf"
    mol = build(frame)
    print(f"[{frame}] {CalcMolFormula(mol)}  原子 {mol.GetNumAtoms()}  "
          f"MW {Descriptors.MolWt(mol):.2f}")
    print(f"      SMILES {Chem.MolToSmiles(Chem.RemoveHs(mol))}")
    print(f"      手性 {Chem.FindMolChiralCenters(mol, includeUnassigned=True)}")
    mol.SetProp("_Name", f"puerarin_frame{frame}")
    out.write_text(Chem.MolToMolBlock(mol))
    print(f"      -> {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
