"""对接结果分析：配体效率、活性口袋占据、关键残基接触。

对比基准：1D7K 中 LLP（PLP-DFMO 共价加合物）所界定的活性位点残基。
"""
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parent.parent
DOCK = ROOT / "docking"
POSES = DOCK / "poses"


def parse_pdbqt_heavy(path):
    """读取 pdbqt 第一模型的配体重原子坐标。"""
    xyz, names = [], []
    with open(path, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("MODEL"):
                if xyz:
                    break
                continue
            if line.startswith(("ATOM", "HETATM")):
                el = line[77:79].strip()
                if el.upper() == "HD" or el.upper() == "H":
                    continue
                xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                names.append(el)
    return np.array(xyz)


def pocket_residues(cut=4.0):
    """1D7K 中距 LLP 任一原子 < cut 的蛋白残基。"""
    llp, prot = [], {}
    for line in open(ROOT / "data" / "protein" / "1D7K.pdb", encoding="utf-8", errors="ignore"):
        if line.startswith("HETATM") and line[17:20].strip() == "LLP" and line[21] == "A":
            llp.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        elif line.startswith("ATOM") and line[21] == "A":
            key = (line[17:20].strip(), int(line[22:26]))
            prot.setdefault(key, []).append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    llp = np.array(llp)
    out = {}
    for key, coords in prot.items():
        d = np.linalg.norm(np.array(coords)[:, None, :] - llp[None, :, :], axis=2).min()
        if d < cut:
            out[key] = d
    return out


def main():
    pres = pocket_residues(4.0)
    print(f"活性位点残基（LLP 4Å 内）: {len(pres)} 个")
    print("  " + ", ".join(f"{r}{n}" for r, n in sorted(pres, key=lambda x: x[1])))

    rows = []
    for f in sorted(POSES.glob("*_out.pdbqt")):
        xyz = parse_pdbqt_heavy(f)
        if len(xyz) == 0:
            continue
        # 与口袋残基的最近距离
        close = 0
        for (res, num), _ in pres.items():
            pass
        # 重算残基-配体接触
        contacts = {}
        for line in open(ROOT / "data" / "protein" / "1D7K.pdb", encoding="utf-8", errors="ignore"):
            if line.startswith("ATOM") and line[21] == "A":
                key = f"{line[17:20].strip()}{line[22:26].strip()}"
                c = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                d = np.linalg.norm(xyz - np.array(c), axis=1).min()
                if d < 4.0:
                    contacts[key] = round(min(d, contacts.get(key, 99)), 2)

        # 打分
        scores = [float(l.split()[3]) for l in open(POSES / f.name, encoding="utf-8", errors="ignore")
                  if l.startswith("REMARK VINA RESULT")]
        best = min(scores) if scores else None
        nha = len(xyz)
        rows.append(dict(
            pose=f.stem.replace("_out", ""),
            n_heavy=nha,
            best_affinity=best,
            ligand_efficiency=round(-best / nha, 3) if best else None,
            n_pocket_contacts=len(contacts),
            contacts=json.dumps(contacts, ensure_ascii=False),
        ))
    df = pd.DataFrame(rows)
    # 恢复可读名称
    q = pd.read_csv(ROOT / "data" / "query_molecules.csv")
    q["key"] = q["name"].map(lambda s: "".join(c if c.isalnum() else "_" for c in s)[:30])
    df = df.merge(q[["key", "name", "role"]], left_on="pose", right_on="key", how="left")
    df = df.drop(columns=["key"]).sort_values("ligand_efficiency", ascending=False)
    df.to_csv(ROOT / "results" / "docking_analysis.csv", index=False)
    print("\n=== 对接分析 ===")
    print(df[["name", "best_affinity", "n_heavy", "ligand_efficiency", "n_pocket_contacts"]]
          .to_string(index=False, float_format=lambda v: f"{v:.2f}"))
    with open(ROOT / "results" / "pocket_residues.json", "w", encoding="utf-8") as fh:
        json.dump({f"{r}{n}": round(d, 2) for (r, n), d in pres.items()}, fh, indent=1)


if __name__ == "__main__":
    main()
