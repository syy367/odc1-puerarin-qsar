"""DFT 前线轨道能级图（单图，Times New Roman）。"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 7
plt.rcParams["axes.titlesize"] = 7.5
plt.rcParams["axes.labelsize"] = 7
plt.rcParams["xtick.labelsize"] = 6.5
plt.rcParams["ytick.labelsize"] = 6.5
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["figure.dpi"] = 300

ROOT = Path(__file__).resolve().parent.parent
dft = json.load(open(ROOT / "results" / "dft_results.json", encoding="utf-8"))


def find(sub):
    for k in dft:
        if sub.lower() in k.lower():
            return k
    return None


PLAN = [("Puerarin", "Puerarin"), ("Quercetin", "Quercetin"), ("Genistein", "Genistein"),
        ("Daidzein", "Daidzein"), ("Resveratrol", "Resveratrol"),
        ("Eflornithine", "DFMO"), ("aminooxy", "APA")]
items = [(find(k), lab) for k, lab in PLAN]
items = [(k, lab) for k, lab in items if k]
print("绘图分子:", [lab for _, lab in items])

fig, ax = plt.subplots(figsize=(3.346, 2.087))
for i, (k, lab) in enumerate(items):
    v = dft[k]
    c = "#E45756" if lab == "Puerarin" else ("#4C78A8" if lab in ("DFMO", "APA") else "#8C8C8C")
    ax.plot([i - 0.30, i + 0.30], [v["homo_eV"]] * 2, color=c, lw=3.4, solid_capstyle="butt")
    ax.plot([i - 0.30, i + 0.30], [v["lumo_eV"]] * 2, color=c, lw=1.5, ls=(0, (3, 2)))
    ax.annotate("", xy=(i, v["lumo_eV"]), xytext=(i, v["homo_eV"]),
                arrowprops=dict(arrowstyle="<->", color=c, lw=0.8, shrinkA=0, shrinkB=0))
    ax.text(i - 0.36, (v["homo_eV"] + v["lumo_eV"]) / 2, f"{v['gap_eV']:.2f}",
            fontsize=5.8, va="center", ha="right", color=c)

ax.set_xticks(range(len(items)))
ax.set_xticklabels([lab for _, lab in items], rotation=32, ha="right", fontsize=6)
ax.set_ylabel("Orbital energy (eV)")
ax.set_title("Frontier molecular orbitals (B3LYP/6-31G*)")
lo = min(v["homo_eV"] for v in dft.values()) - 0.8
hi = max(v["lumo_eV"] for v in dft.values()) + 0.9
ax.set_ylim(lo, hi)
ax.set_xlim(-0.55, len(items) - 0.15)
ax.axhline(0, color="#cccccc", lw=0.6, ls=":")
ax.legend(handles=[Line2D([], [], color="k", lw=3.2, label="HOMO"),
                   Line2D([], [], color="k", lw=1.5, ls=(0, (3, 2)), label="LUMO")],
          frameon=False, loc="upper left", fontsize=9, ncol=2)
fig.tight_layout()
fig.savefig(ROOT / "figures" / "journal" / "fig11_dft_orbitals.png", dpi=300, bbox_inches="tight", facecolor="white")
print("fig11 已输出")
