"""平衡段几何证据出图：接触占有率 + 关键残基距离时间序列。

出图规范：单图、Times New Roman、全英文标签。
"""
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent.parent
MD = ROOT / "md"
OUT = MD / "results"

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 10
plt.rcParams["axes.linewidth"] = 0.8
plt.rcParams["xtick.direction"] = "out"
plt.rcParams["ytick.direction"] = "out"

HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "PHE", "TRP", "MET", "PRO", "GLY", "CYS"}
POLAR = {"SER", "THR", "ASN", "GLN", "TYR"}
CHARGED = {"ASP", "GLU", "LYS", "ARG", "HIS"}


def residue_class(name: str) -> str:
    aa = name[:3]
    if aa in CHARGED:
        return "charged"
    if aa in POLAR:
        return "polar"
    return "hydrophobic"


def plot_occupancy(df: pd.DataFrame):
    df = df.sort_values("contact_occupancy_pct", ascending=True)
    colors = {"hydrophobic": "#4C72B0", "polar": "#55A868", "charged": "#C44E52"}
    cs = [colors[residue_class(r)] for r in df.residue]

    fig, ax = plt.subplots(figsize=(3.5, 4.6))
    ax.barh(df.residue, df.contact_occupancy_pct, color=cs, height=0.72,
            edgecolor="none")
    ax.axvline(50, ls="--", lw=0.8, color="#888888")
    ax.axvline(70, ls=":", lw=0.8, color="#888888")
    ax.set_xlim(0, 105)
    ax.set_xlabel("Contact occupancy (%)", fontsize=10)
    ax.set_ylabel("Residue", fontsize=10)
    ax.tick_params(labelsize=8)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.text(51, len(df) - 0.45, "50%", fontsize=7, color="#666666", va="bottom")
    ax.text(71, len(df) - 0.45, "70%", fontsize=7, color="#666666", va="bottom")
    ax.set_ylim(-0.7, len(df) + 0.2)

    handles = [plt.Rectangle((0, 0), 1, 1, color=colors[k]) for k in colors]
    ax.legend(handles, ["hydrophobic", "polar", "charged"], fontsize=7,
              frameon=False, loc="lower right", bbox_to_anchor=(1.02, -0.02))
    fig.tight_layout()
    fig.savefig(OUT / "fig_contact_occupancy.png", dpi=600)
    plt.close(fig)
    print("写出 fig_contact_occupancy.png")


def frame_min_distances(targets):
    files = sorted((MD / "equil_frames").glob("f*.pdb"), key=lambda p: int(p.stem[1:]))
    series = {t: [] for t in targets}
    times = []
    for f in files:
        lines = f.read_text().splitlines()
        lig = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                        for l in lines if l[17:20].strip() == "MOL"])
        tree = cKDTree(lig)
        acc = {t: 9e9 for t in targets}
        for l in lines:
            if not l.startswith("ATOM"):
                continue
            rn = l[17:20].strip()
            if rn == "MOL":
                continue
            key = f"{rn}{l[22:26].strip()}"
            if key not in acc:
                continue
            d, _ = tree.query((float(l[30:38]), float(l[38:46]), float(l[46:54])))
            if d < acc[key]:
                acc[key] = d
        times.append((int(f.stem[1:]) - 1) * 0.1)
        for t in targets:
            series[t].append(acc[t])
    return np.array(times), {k: np.array(v) for k, v in series.items()}


def plot_timeseries():
    targets = ["ILE251", "LEU246", "PHE208", "CYS202", "LYS245", "ASP204"]
    times, series = frame_min_distances(targets)
    fig, ax = plt.subplots(figsize=(3.4, 2.4))
    colors = ["#4C72B0", "#55A868", "#C44E52", "#8172B2", "#937860", "#DA8BC3"]
    for (k, v), c in zip(series.items(), colors):
        ax.plot(times, v, lw=0.9, color=c, label=k.capitalize().replace("Ile", "Ile"))
    ax.axhline(4.0, ls="--", lw=0.8, color="#888888")
    ax.set_xlabel("Time (ns)", fontsize=9)
    ax.set_ylabel("Min. distance to ligand (Å)", fontsize=9)
    ax.set_ylim(1.5, 6.0)
    ax.tick_params(labelsize=8)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.legend(fontsize=6.5, frameon=False, ncol=2, loc="upper right")
    fig.tight_layout()
    fig.savefig(OUT / "fig_contact_timeseries.png", dpi=600)
    plt.close(fig)
    print("写出 fig_contact_timeseries.png")


if __name__ == "__main__":
    df = pd.read_csv(OUT / "geom_contact_occupancy.csv")
    plot_occupancy(df)
    plot_timeseries()
