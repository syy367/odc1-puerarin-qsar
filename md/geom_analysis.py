"""60–100 ns 平衡段的几何证据：接触占有率、氢键占有率、距离/角度统计。

判据
  * 接触：残基任一原子到配体任一重原子的最小距离 < 4.0 Å
  * 氢键：重原子 D···A < 3.5 Å，H···A < 2.5 Å，D–H···A 角 > 120°
  * 占有率：满足条件的帧数 / 总帧数

轨迹已去水、配体保留 6 个羟基氢，配体原子名无区分性（全是元素符号），
因此按几何就近原则确定每个氢的母原子。

用法：python md/geom_analysis.py
"""
from pathlib import Path
import re
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent.parent
MD = ROOT / "md"
FRAMES = MD / "equil_frames"
OUT = MD / "results"

CONTACT_CUT = 4.0          # Å
DA_CUT = 3.5               # D···A
HA_CUT = 2.5               # H···A
ANG_CUT = 120.0            # D-H···A

STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}


def parse(path: Path):
    """返回 (蛋白记录 list[(resname,resseq,atomname,elem,xyz)], 配体重原子 xyz, 配体H xyz, 配体H母原子索引)。"""
    prot = []
    lig_xyz, lig_el, lig_name = [], [], []
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        resname = line[17:20].strip()
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        elem = line[76:78].strip() or re.sub(r"[^A-Za-z]", "", line[12:16])[:1]
        if resname == "MOL":
            lig_xyz.append(xyz)
            lig_el.append(elem.upper())
            lig_name.append(line[12:16].strip())
        elif resname in STANDARD_AA:
            prot.append((resname, line[22:26].strip(), line[12:16].strip(),
                         elem.upper(), xyz))
    return prot, np.array(lig_xyz), lig_el, lig_name


def analyse_frame(path: Path):
    prot, lig_xyz, lig_el, lig_name = parse(path)
    lig_el = np.array(lig_el)
    heavy_mask = lig_el != "H"
    lig_heavy = lig_xyz[heavy_mask]
    lig_heavy_el = lig_el[heavy_mask]
    lig_h = lig_xyz[~heavy_mask]

    p_xyz = np.array([a[4] for a in prot])
    p_el = np.array([a[3] for a in prot])
    p_res = np.array([f"{a[0]}{a[1]}" for a in prot])
    p_name = np.array([a[2] for a in prot])

    tree_lig = cKDTree(lig_heavy)
    d_all, idx_all = tree_lig.query(p_xyz)

    # 1) 残基级接触（<4.0 Å）
    res_min = defaultdict(lambda: 9e9)
    for d, r in zip(d_all, p_res):
        if d < res_min[r]:
            res_min[r] = d
    contacts = {r: d for r, d in res_min.items() if d < CONTACT_CUT}

    # 2) 氢键
    hbonds = []          # (残基, 供体/受体原子, 配体原子, D···A, H···A, 角度)
    nolig = lig_heavy_el != "C"

    # 2a 配体作供体：配体 O-H ··· 蛋白受体
    if len(lig_h):
        tree_h = cKDTree(lig_h)
        for i, h in enumerate(lig_h):
            _, j = tree_h.query(h, k=2)                     # 自身 + 最近重原子?? 见下
            dj = np.linalg.norm(lig_heavy - h, axis=1)
            donor = int(np.argmin(dj))
            if lig_heavy_el[donor] != "O":
                continue
            acc_mask = np.isin(p_el, ("O", "N", "S"))
            if not acc_mask.any():
                continue
            acc_xyz = p_xyz[acc_mask]
            dha, k = cKDTree(acc_xyz).query(h)
            if dha > HA_CUT:
                continue
            a = acc_xyz[k]
            dd = np.linalg.norm(lig_heavy[donor] - a)
            if dd > DA_CUT:
                continue
            v1, v2 = lig_heavy[donor] - h, a - h
            ang = np.degrees(np.arccos(np.clip(v1 @ v2 / (np.linalg.norm(v1) * np.linalg.norm(v2)),
                                               -1, 1)))
            if ang < ANG_CUT:
                continue
            gi = np.where(acc_mask)[0][k]
            hbonds.append((p_res[gi], p_name[gi], f"O{i + 1}", dd, dha, ang, "ligand-donor"))

    # 2b 配体作受体：蛋白 N/O-H ··· 配体 O
    prot_h = np.where(p_el == "H")[0]
    if len(prot_h):
        tree_all = cKDTree(p_xyz)
        acc = lig_heavy_el != "C"
        for ih in prot_h:
            h = p_xyz[ih]
            kd = p_name[ih][1:] if p_name[ih].startswith("H") else ""
            # 母重原子：同残基内最近的重原子
            same = np.where((p_res == p_res[ih]) & (p_el != "H"))[0]
            if not len(same):
                continue
            dd_vec = np.linalg.norm(p_xyz[same] - h, axis=1)
            donor_i = same[int(np.argmin(dd_vec))]
            if p_el[donor_i] not in ("N", "O"):
                continue
            if not acc.any():
                continue
            dh = np.linalg.norm(lig_heavy[acc] - h, axis=1)
            k = int(np.argmin(dh))
            if dh[k] > HA_CUT:
                continue
            a = lig_heavy[acc][k]
            dd = np.linalg.norm(p_xyz[donor_i] - a)
            if dd > DA_CUT:
                continue
            v1, v2 = p_xyz[donor_i] - h, a - h
            ang = np.degrees(np.arccos(np.clip(v1 @ v2 / (np.linalg.norm(v1) * np.linalg.norm(v2)),
                                               -1, 1)))
            if ang < ANG_CUT:
                continue
            hbonds.append((p_res[donor_i], p_name[donor_i],
                           f"{lig_heavy_el[acc][k]}{np.where(acc)[0][k] + 1}",
                           dd, dh[k], ang, "protein-donor"))

    # 配体重原子转成 1-based 编号用于记录
    return contacts, hbonds, len(lig_heavy), res_min


def main():
    files = sorted(FRAMES.glob("f*.pdb"), key=lambda p: int(p.stem[1:]))
    if not files:
        raise SystemExit("没有找到平衡段帧")
    print(f"帧数 {len(files)}（{files[0].stem} … {files[-1].stem}）")

    occ = defaultdict(int)
    dist_sum = defaultdict(float)
    hb_occ = defaultdict(int)
    hb_stats = defaultdict(list)
    lig_atom_of = {}

    for n, f in enumerate(files, 1):
        contacts, hbonds, nheavy, _ = analyse_frame(f)
        for r, d in contacts.items():
            occ[r] += 1
            dist_sum[r] += d
        for r, atom, latom, dd, ha, ang, kind in hbonds:
            key = (r, atom, kind)
            hb_occ[key] += 1
            hb_stats[key].append((dd, ha, ang))
            lig_atom_of[key] = latom
        if n % 10 == 0:
            print(f"  已处理 {n}/{len(files)}")

    nf = len(files)
    rows = [{"residue": r, "contact_occupancy_pct": 100.0 * c / nf,
             "mean_min_dist_A": dist_sum[r] / c, "n_frames": c}
            for r, c in occ.items()]
    df = pd.DataFrame(rows).sort_values(["contact_occupancy_pct", "mean_min_dist_A"],
                                        ascending=[False, True])
    OUT.mkdir(exist_ok=True)
    df.to_csv(OUT / "geom_contact_occupancy.csv", index=False, float_format="%.2f")
    print(f"\n接触残基 {len(df)} 个，占有率达 50% 的 {int((df.contact_occupancy_pct >= 50).sum())} 个")

    hrows = []
    for (r, atom, kind), c in sorted(hb_occ.items(), key=lambda x: -x[1]):
        st = np.array(hb_stats[(r, atom, kind)])
        hrows.append({"residue": r, "protein_atom": atom, "type": kind,
                      "hbond_occupancy_pct": 100.0 * c / nf,
                      "mean_DA_A": st[:, 0].mean(), "mean_HA_A": st[:, 1].mean(),
                      "mean_angle_deg": st[:, 2].mean(),
                      "ligand_atom": lig_atom_of[(r, atom, kind)]})
    hdf = pd.DataFrame(hrows)
    if len(hdf):
        hdf.to_csv(OUT / "geom_hbonds.csv", index=False, float_format="%.2f")
        print(f"氢键条目 {len(hdf)} 个（占有率≥20% 的 {int((hdf.hbond_occupancy_pct >= 20).sum())} 个）")
    else:
        print("未检出满足判据的氢键")

    print("\n=== 接触占有率前 20 ===")
    print(df.head(20).to_string(index=False))
    if len(hdf):
        print("\n=== 氢键 ===")
        print(hdf.head(15).to_string(index=False))


if __name__ == "__main__":
    main()
