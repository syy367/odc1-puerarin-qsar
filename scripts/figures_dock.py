"""对接结果与分子结构图（单图输出，图内文字 Times New Roman）。"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Draw

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parent.parent
FIG = ROOT / "figures" / "journal"
FIG.mkdir(exist_ok=True)

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 7
plt.rcParams["axes.titlesize"] = 7.5
plt.rcParams["axes.labelsize"] = 7
plt.rcParams["xtick.labelsize"] = 6.5
plt.rcParams["ytick.labelsize"] = 6.5
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["figure.dpi"] = 300

COL_TARGET = "#E45756"
COL_REF = "#4C78A8"
COL_NP = "#B4B4B4"


def short(n):
    return {
        "Puerarin": "Puerarin",
        "Eflornithine (DFMO)": "DFMO",
        "Genistein": "Genistein",
        "Daidzein": "Daidzein",
        "Quercetin": "Quercetin",
        "Resveratrol": "Resveratrol",
        "APA (3-aminooxy-1-aminopropane)": "APA",
    }.get(n, n)


def color(n):
    if n == "Puerarin":
        return COL_TARGET
    if "DFMO" in n or "APA" in n:
        return COL_REF
    return COL_NP


def save(fig, name):
    fig.tight_layout()
    fig.savefig(FIG / name, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


da = pd.read_csv(ROOT / "results" / "docking_analysis.csv")
da["label"] = da["name"].map(short)

# 图 7：对接打分
d = da.sort_values("best_affinity")
fig, ax = plt.subplots(figsize=(3.346, 2.323))
ax.barh(range(len(d)), d.best_affinity, color=[color(n) for n in d.name], height=0.66)
ax.set_yticks(range(len(d)))
ax.set_yticklabels(d.label)
ax.set_xlabel("Vina docking score (kcal/mol)")
ax.set_xlim(min(d.best_affinity) - 1.2, 0)
for i, v in enumerate(d.best_affinity):
    ax.text(v - 0.12, i, f"{v:.2f}", va="center", ha="right", fontsize=8)
ax.set_title("Docking score against ODC1 (1D7K)")
save(fig, "fig7_docking_score.png")

# 图 8：配体效率
d = da.sort_values("ligand_efficiency")
fig, ax = plt.subplots(figsize=(3.346, 2.323))
ax.barh(range(len(d)), d.ligand_efficiency, color=[color(n) for n in d.name], height=0.66)
ax.set_yticks(range(len(d)))
ax.set_yticklabels(d.label)
ax.set_xlabel("Ligand efficiency (kcal/mol per heavy atom)")
for i, v in enumerate(d.ligand_efficiency):
    ax.text(v + 0.008, i, f"{v:.2f}", va="center", fontsize=8)
ax.set_title("Ligand efficiency against ODC1")
save(fig, "fig8_ligand_efficiency.png")

# 图 9：QSAR 预测 vs 配体效率
pr = pd.read_csv(ROOT / "results" / "predictions.csv")
pr["label"] = pr["name"].map(short).fillna(pr["name"])
m = pr.dropna(subset=["pIC50_ensemble"]).merge(da[["name", "ligand_efficiency"]], on="name", how="inner")
fig, ax = plt.subplots(figsize=(3.346, 2.638))
for _, r in m.iterrows():
    ax.scatter(r.ligand_efficiency, r.pIC50_ensemble, s=70, color=color(r["name"]),
               edgecolor="white", linewidth=0.6, zorder=3)
    off = {"Daidzein": (-42, -1), "Genistein": (9, -11), "Quercetin": (9, 5),
           "Resveratrol": (9, -11), "Puerarin": (9, 3), "DFMO": (9, 3), "APA": (-20, 7)}
    ax.annotate(r["label"], (r["ligand_efficiency"], r["pIC50_ensemble"]),
                textcoords="offset points", xytext=off.get(r["label"], (6, 5)), fontsize=8)
ax.set_xlim(0.19, 0.72); ax.set_ylim(2.0, 7.3)
ax.axvspan(0.19, 0.72, color="#F7F7F7", zorder=0)
ax.set_xlabel("Ligand efficiency (docking)")
ax.set_ylabel("Predicted pIC$_{50}$ (QSAR)")
ax.set_title("Structure-based vs ligand-based ranking")
save(fig, "fig9_le_vs_pIC50.png")

# 图 10：两个关键分子的二维结构
q = pd.read_csv(ROOT / "data" / "query_molecules.csv")
mols = []
for nm in ["Puerarin", "Eflornithine (DFMO)"]:
    smi = q.loc[q.name == nm, "smiles"].iloc[0]
    mols.append(Chem.MolFromSmiles(smi))
img = Draw.MolsToGridImage(mols, molsPerRow=2, subImgSize=(500, 380),
                           legends=["Puerarin (C21H20O9)", "Eflornithine / DFMO (C6H12F2N2O2)"])
img.save(FIG / "fig10_structures.png")
print("图已输出: fig7, fig8, fig9, fig10")
