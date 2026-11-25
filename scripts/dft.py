"""DFT 电子结构计算（PySCF, B3LYP/6-31G*）。

输出 HOMO/LUMO 能量与能隙、偶极矩、Mulliken 电荷，用于解释葛根素与口袋残基
形成氢键/静电作用的电子结构基础。

几何：RDKit ETKDGv3 + MMFF94 优化（单点能计算前不做 DFT 几何优化，以控制耗时）。
用法: python scripts/dft.py [分子名，可省略=全部]
"""
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
RES.mkdir(exist_ok=True)

BASIS = "6-31g*"
XC = "b3lyp"


def embed(smiles, seed=42):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        AllChem.EmbedMolecule(mol, randomSeed=seed, useRandomCoords=True)
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=5000)
    except Exception:
        pass
    conf = mol.GetConformer()
    atoms = []
    for i, a in enumerate(mol.GetAtoms()):
        p = conf.GetAtomPosition(i)
        atoms.append((a.GetSymbol(), (p.x, p.y, p.z)))
    return mol, atoms


def run_dft(atoms, verbose=3):
    from pyscf import dft, gto

    mol = gto.M(atom=atoms, basis=BASIS, charge=0, spin=0, verbose=verbose,
                max_memory=3000)
    mf = dft.RKS(mol)
    mf.xc = XC
    mf = mf.density_fit(auxbasis="weigend")   # 密度拟合，显著降低内存与耗时
    mf.conv_tol = 1e-7
    t0 = time.time()
    mf.kernel()
    dt = time.time() - t0
    occ = mf.mo_occ
    energies = mf.mo_energy
    homo = float(energies[occ > 0][-1])
    lumo = float(energies[occ == 0][0])
    pop, charges = mf.mulliken_pop(verbose=0)
    dip = mf.dip_moment(verbose=0)
    return dict(
        homo_eV=round(homo * 27.2114, 3),
        lumo_eV=round(lumo * 27.2114, 3),
        gap_eV=round((lumo - homo) * 27.2114, 3),
        dipole_D=round(float(np.linalg.norm(dip)), 3),
        n_basis=int(mol.nao_nr()),
        charges=[round(float(c), 4) for c in charges],
        seconds=round(dt, 1),
        e_tot=round(float(mf.e_tot), 6),
    )


def main():
    q = pd.read_csv(ROOT / "data" / "query_molecules.csv")
    target = sys.argv[1] if len(sys.argv) > 1 else None
    if target:
        q = q[q["name"].str.contains(target, case=False, regex=False)]
    out_path = RES / "dft_results.json"
    allres = json.load(open(out_path, encoding="utf-8")) if out_path.exists() else {}
    for _, r in q.iterrows():
        name = r["name"]
        if name in allres and not os.environ.get("DFT_FORCE"):
            print(f"skip {name} (已算)", flush=True)
            continue
        print(f"\n>>> {name}  ({BASIS}/{XC})", flush=True)
        mol, atoms = embed(r["smiles"])
        res = run_dft(atoms)
        res["smiles"] = r["smiles"]
        res["formula"] = Chem.rdMolDescriptors.CalcMolFormula(mol)
        res["elements"] = [a.GetSymbol() for a in mol.GetAtoms()]
        allres[name] = res
        print(f"    HOMO={res['homo_eV']} eV  LUMO={res['lumo_eV']} eV  gap={res['gap_eV']} eV  "
              f"dipole={res['dipole_D']} D  ({res['seconds']} s, {res['n_basis']} basis)", flush=True)
        with open(out_path, "w", encoding="utf-8") as fh:
            json.dump(allres, fh, indent=1, ensure_ascii=False)


if __name__ == "__main__":
    main()
