"""AutoDock Vina 1.2.5 分子对接：ODC1 (1D7K, A 链) vs 葛根素及对照分子。

受体为 apo 口袋（已去除 LLP 共价加合物与水），盒子覆盖 LLP 所在活性位点。
配体 3D 构象由 RDKit ETKDG + MMFF94 生成，pdbqt 由 Meeko 写出。
"""
import json
import subprocess
from pathlib import Path

import pandas as pd
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parent.parent
DOCK = ROOT / "docking"
OUT = DOCK / "poses"
OUT.mkdir(parents=True, exist_ok=True)
VINA = ROOT / "tools" / "bin" / "vina.exe"

CENTER = (6.69, 11.74, 54.76)  # 1D7K A 链 LLP 质心
SIZE = (22.0, 22.0, 22.0)


def make_pdbqt(name, smiles, seed=42):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"非法 SMILES: {name}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        AllChem.EmbedMolecule(mol, randomSeed=seed, useRandomCoords=True)
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
    except Exception:
        pass
    prep = MoleculePreparation()
    setups = prep.prepare(mol)
    pdbqt, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        raise RuntimeError(f"{name} pdbqt 生成失败: {err}")
    safe = "".join(c if c.isalnum() else "_" for c in name)[:30]
    path = OUT / f"{safe}.pdbqt"
    path.write_text(pdbqt, encoding="utf-8")
    return path, safe


def run_vina(lig_path, safe):
    out_pdbqt = OUT / f"{safe}_out.pdbqt"
    log = OUT / f"{safe}.log"
    cmd = [
        str(VINA),
        "--receptor", str(DOCK / "receptor.pdbqt"),
        "--ligand", str(lig_path),
        "--center_x", str(CENTER[0]), "--center_y", str(CENTER[1]), "--center_z", str(CENTER[2]),
        "--size_x", str(SIZE[0]), "--size_y", str(SIZE[1]), "--size_z", str(SIZE[2]),
        "--exhaustiveness", "32", "--num_modes", "9", "--seed", "42",
        "--out", str(out_pdbqt),
    ]
    r = subprocess.run(cmd, capture_output=True, timeout=1800)
    log.write_text((r.stdout + r.stderr).decode("utf-8", errors="ignore"), encoding="utf-8")
    scores = []
    if out_pdbqt.exists():
        for line in out_pdbqt.read_text(encoding="utf-8", errors="ignore").splitlines():
            if line.startswith("REMARK VINA RESULT"):
                parts = line.split()
                scores.append(float(parts[3]))
    return scores, r.returncode


def main():
    q = pd.read_csv(ROOT / "data" / "query_molecules.csv")
    rows = []
    for _, r in q.iterrows():
        try:
            path, safe = make_pdbqt(r["name"], r["smiles"])
            scores, code = run_vina(path, safe)
            best = scores[0] if scores else None
            print(f"{r['name'][:34]:36s} best={best} modes={len(scores)} rc={code}")
            rows.append(dict(
                name=r["name"], role=r["role"], smiles=r["smiles"],
                best_affinity_kcal=best,
                all_modes=json.dumps(scores),
                n_modes=len(scores),
            ))
        except Exception as e:
            print(f"{r['name']}: 失败 {type(e).__name__} {e}")
            rows.append(dict(name=r["name"], role=r["role"], smiles=r["smiles"],
                             best_affinity_kcal=None, all_modes="[]", n_modes=0))
    df = pd.DataFrame(rows).sort_values("best_affinity_kcal")
    df.to_csv(ROOT / "results" / "docking_results.csv", index=False)
    print("\n=== 对接打分 (kcal/mol) ===")
    print(df[["name", "best_affinity_kcal", "n_modes", "role"]].to_string(index=False))


if __name__ == "__main__":
    main()
