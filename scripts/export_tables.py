"""补齐投稿表格所需的三份数据：
  table_s2_descriptors.csv   保留描述符的定义/类型/重要性/相关系数
  table_s3_docking_detail.csv 每个分子 9 个对接构象的打分
  table_s5_sensitivity.csv   仅用 Tier A+B（纯 IC50）重训的敏感性分析
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.cross_decomposition import PLSRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler

import qsar as Q

ROOT = Path(__file__).resolve().parent.parent
DATA, RES, DOCK = ROOT / "data", ROOT / "results", ROOT / "docking"

DEF = {
    "FpDensityMorgan2": ("Morgan fingerprint (radius 2) bit density", "2D topology"),
    "BalabanJ": ("Balaban topological connectivity index", "2D topology"),
    "PEOE_VSA3": ("Gasteiger charge-weighted VdW surface area, bin 3", "Surface / charge"),
    "SlogP_VSA3": ("Crippen logP-weighted VdW surface area, bin 3", "Surface / hydrophobicity"),
    "BCUT2D_MWHI": ("BCUT high eigenvalue (mass-weighted)", "2D topology"),
    "VSA_EState8": ("VSA-corrected atom-type E-state, bin 8", "Surface / electronic"),
    "SMR_VSA4": ("Crippen molar refractivity-weighted VdW surface area, bin 4", "Surface / polarizability"),
    "VSA_EState4": ("VSA-corrected atom-type E-state, bin 4", "Surface / electronic"),
}


def main():
    df = pd.read_csv(DATA / "dataset_odc1.csv")
    y = df["pIC50"].values
    X_all = Q.descriptors(df["smiles"].tolist())
    X_f, feats = Q.clean_features(X_all, y)
    scaler = StandardScaler().fit(X_f)
    Xs = pd.DataFrame(scaler.transform(X_f), columns=feats, index=X_f.index)

    sel = json.load(open(RES / "qsar_summary.json", encoding="utf-8"))["selected_descriptors"]
    tr, te = Q.kennard_stone(Xs.values, n_test=12)
    rf = RandomForestRegressor(n_estimators=800, random_state=42, n_jobs=-1)
    rf.fit(Xs[sel].iloc[tr], y[tr])
    imp = pd.Series(rf.feature_importances_, index=sel).sort_values(ascending=False)

    rows = []
    for name in imp.index:
        r = float(np.corrcoef(Xs[name], y)[0, 1])
        d, t = DEF.get(name, ("—", "—"))
        rows.append(dict(descriptor=name, definition=d, type=t,
                         rf_importance=round(float(imp[name]), 4),
                         pearson_r_pIC50=round(r, 3)))
    s2 = pd.DataFrame(rows)
    s2.to_csv(RES / "table_s2_descriptors.csv", index=False)
    print("Table S2 ->"); print(s2.to_string(index=False))

    # ---------- Table S3: 9 个构象打分 ----------
    q = pd.read_csv(DATA / "query_molecules.csv")
    q["key"] = q["name"].map(lambda s: "".join(c if c.isalnum() else "_" for c in s)[:30])
    rows = []
    for _, r in q.iterrows():
        f = DOCK / "poses" / f"{r['key']}_out.pdbqt"
        scores = []
        if f.exists():
            for line in open(f, encoding="utf-8", errors="ignore"):
                if line.startswith("REMARK VINA RESULT"):
                    scores.append(float(line.split()[3]))
        while len(scores) < 9:
            scores.append(np.nan)
        rows.append({"compound": r["name"], "role": r["role"],
                     **{f"mode{i+1}": scores[i] for i in range(9)}})
    s3 = pd.DataFrame(rows)
    s3["best"] = s3[[f"mode{i+1}" for i in range(9)]].min(axis=1)
    s3.to_csv(RES / "table_s3_docking_detail.csv", index=False)
    print("\nTable S3 ->", s3[["compound", "best"]].to_string(index=False))

    # ---------- Table S5: 仅 Tier A+B 重训 ----------
    sub = df[df.tier.isin(["A", "B"])].reset_index(drop=True)
    ys = sub["pIC50"].values
    Xs_raw = Q.descriptors(sub["smiles"].tolist())
    Xs_f, feats_s = Q.clean_features(Xs_raw, ys)
    sc = StandardScaler().fit(Xs_f)
    Xss = pd.DataFrame(sc.transform(Xs_f), columns=feats_s)
    trs, tes = Q.kennard_stone(Xss.values, n_test=max(5, int(round(0.2 * len(sub)))))

    rf0 = RandomForestRegressor(n_estimators=500, random_state=42, n_jobs=-1).fit(Xss.iloc[trs], ys[trs])
    imps = pd.Series(rf0.feature_importances_, index=feats_s).sort_values(ascending=False)
    best = (-9e9, None)
    for k in [4, 6, 8, 10, 12]:
        if k > len(feats_s):
            continue
        q2 = Q.cv_eval(PLSRegression(n_components=min(3, k)), Xss[imps.index[:k]].iloc[trs], ys[trs])[0]
        if q2 > best[0]:
            best = (q2, k)
    kbest = best[1]
    Xsel_s = Xss[imps.index[:kbest]]

    rows = []
    for name, mdl in Q.build_models(kbest).items():
        q2, rc, _ = Q.cv_eval(mdl, Xsel_s.iloc[trs].values, ys[trs])
        mdl.fit(Xsel_s.iloc[trs].values, ys[trs])
        ptr = mdl.predict(Xsel_s.iloc[trs].values)
        pte = mdl.predict(Xsel_s.iloc[tes].values)
        mtr, mte = Q.metrics(ys[trs], ptr), Q.metrics(ys[tes], pte)
        rows.append(dict(model=name, n_compounds=len(sub), n_train=len(trs), n_test=len(tes),
                         n_descriptors=kbest, Q2_CV=round(q2, 3), RMSE_CV=round(rc, 3),
                         R2_train=round(mtr["R2"], 3), R2_test=round(mte["R2"], 3),
                         RMSE_test=round(mte["RMSE"], 3)))
    s5 = pd.DataFrame(rows)
    s5.to_csv(RES / "table_s5_sensitivity.csv", index=False)
    print("\nTable S5（仅 Tier A+B，n=%d）->" % len(sub))
    print(s5.to_string(index=False))


if __name__ == "__main__":
    main()
