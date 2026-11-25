"""对接构象的氢键分析：配体极性原子 vs 口袋残基极性原子（距离 <= 3.6 Å）。"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
POSES = ROOT / "docking" / "poses"
PDB = ROOT / "data" / "protein" / "1D7K.pdb"
HB_CUT = 3.6

# 蛋白极性原子（N/O 侧链及主链）
prot_polar = []
for line in open(PDB, encoding="utf-8", errors="ignore"):
    if line.startswith("ATOM") and line[21] == "A":
        el = line[76:78].strip().upper()
        if el in ("N", "O"):
            prot_polar.append((
                f"{line[17:20].strip()}{int(line[22:26])}",
                line[12:16].strip(),
                el,
                np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
            ))


def pose_polar(path):
    out = []
    for line in open(path, encoding="utf-8", errors="ignore"):
        if line.startswith("MODEL") and out:
            break
        if line.startswith(("ATOM", "HETATM")):
            t = line[77:79].strip().upper()
            if t.startswith("O") or t.startswith("N"):
                out.append((t, np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])))
    return out


rows = []
detail = {}
for f in sorted(POSES.glob("*_out.pdbqt")):
    lig = pose_polar(f)
    hbs = []
    for t, xyz in lig:
        for res, aname, el, pxyz in prot_polar:
            d = float(np.linalg.norm(xyz - pxyz))
            if d <= HB_CUT:
                hbs.append((res, aname, t, round(d, 2)))
    key = f.stem.replace("_out", "")
    detail[key] = hbs
    rows.append(dict(pose=key, n_hbond=len(hbs),
                     hbond_residues=";".join(sorted({h[0] for h in hbs}))))

q = pd.read_csv(ROOT / "data" / "query_molecules.csv")
q["key"] = q["name"].map(lambda s: "".join(c if c.isalnum() else "_" for c in s)[:30])
df = pd.DataFrame(rows).merge(q[["key", "name"]], left_on="pose", right_on="key", how="left").drop(columns=["key"])
df.to_csv(ROOT / "results" / "hbond_summary.csv", index=False)
with open(ROOT / "results" / "hbond_detail.json", "w", encoding="utf-8") as fh:
    json.dump(detail, fh, indent=1, ensure_ascii=False)

print(df[["name", "n_hbond", "hbond_residues"]].to_string(index=False))
print("\n--- 葛根素氢键明细 ---")
for k, v in detail.items():
    if "Puerarin" in k:
        for res, aname, t, d in sorted(v, key=lambda x: x[3]):
            print(f"  {res:10s} {aname:5s} ... {t:4s}O/N  {d} Å")
