"""ODC1 抑制剂 2D-QSAR：描述符 -> 筛选 -> 划分 -> PLS/RF/XGB -> 验证 -> 预测。

全部路径基于脚本自身目录推导。
"""
import json
import math
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from scipy.spatial.distance import cdist, pdist
from sklearn.cross_decomposition import PLSRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.preprocessing import StandardScaler
from xgboost import XGBRegressor

warnings.filterwarnings("ignore")
RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parent.parent
DATA, RES, FIG = ROOT / "data", ROOT / "results", ROOT / "figures" / "journal"
for d in (RES, FIG):
    d.mkdir(exist_ok=True)

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["font.size"] = 7
plt.rcParams["axes.titlesize"] = 7.5
plt.rcParams["axes.labelsize"] = 7
plt.rcParams["xtick.labelsize"] = 6.5
plt.rcParams["ytick.labelsize"] = 6.5
plt.rcParams["legend.fontsize"] = 6.5
plt.rcParams["legend.handlelength"] = 1.2
plt.rcParams["figure.dpi"] = 300
SEED = 42


# ---------------- 描述符 ----------------
CALC = MoleculeDescriptors.MolecularDescriptorCalculator([d[0] for d in Descriptors.descList])


def descriptors(smiles_list):
    names = [d[0] for d in Descriptors.descList]
    rows = []
    for s in smiles_list:
        m = Chem.MolFromSmiles(s)
        vals = CALC.CalcDescriptors(m) if m is not None else [np.nan] * len(names)
        rows.append(vals)
    X = pd.DataFrame(rows, columns=names)
    X = X.replace([np.inf, -np.inf], np.nan)
    return X


def clean_features(X, y):
    """去全空/常数、去高相关冗余（保留与 y 相关性高者）。"""
    X = X.dropna(axis=1, thresh=int(0.95 * len(X)))
    X = X.fillna(X.median(numeric_only=True))
    nun = X.nunique()
    X = X.loc[:, nun > 1]
    corr_y = X.apply(lambda c: abs(np.corrcoef(c, y)[0, 1]) if c.std() > 0 else 0)
    order = corr_y.sort_values(ascending=False).index.tolist()
    keep, kept = [], []
    cm = X.corr().abs()
    for c in order:
        if all(cm.loc[c, k] <= 0.90 for k in keep):
            keep.append(c)
    return X[keep], keep


def kennard_stone(X, n_test):
    """KS 划分：先取最近质心的点，再贪心取最远点，最后 n_test 个作测试集。"""
    Xa = np.asarray(X, dtype=float)
    D = cdist(Xa, Xa)
    center = Xa.mean(0)
    first = int(np.argmin(((Xa - center) ** 2).sum(1)))
    sel = [first]
    while len(sel) < n_test:
        dmin = D[:, sel].min(1)
        dmin[sel] = -1
        sel.append(int(np.argmax(dmin)))
    test = sel[-n_test:]
    train = [i for i in range(len(Xa)) if i not in test]
    return np.array(train), np.array(test)


# ---------------- 建模 ----------------
def metrics(y, yhat):
    return dict(
        R2=r2_score(y, yhat),
        RMSE=math.sqrt(mean_squared_error(y, yhat)),
        MAE=mean_absolute_error(y, yhat),
    )


def build_models(nfeat):
    return {
        "PLS": PLSRegression(n_components=min(4, nfeat)),
        "RF": RandomForestRegressor(n_estimators=500, max_features="sqrt", min_samples_leaf=2,
                                    random_state=SEED, n_jobs=-1),
        "XGB": XGBRegressor(n_estimators=400, max_depth=3, learning_rate=0.05, subsample=0.8,
                            colsample_bytree=0.8, reg_lambda=1.0, random_state=SEED, n_jobs=2),
    }


def cv_eval(model, X, y, folds=5):
    kf = KFold(n_splits=folds, shuffle=True, random_state=SEED)
    pred = cross_val_predict(model, X, y, cv=kf)
    sc = metrics(y, pred)
    return sc["R2"], sc["RMSE"], pred


def main():
    df = pd.read_csv(DATA / "dataset_odc1.csv")
    X_all = descriptors(df["smiles"].tolist())
    y = df["pIC50"].values

    # ---- 防泄漏的顺序：基础清洗不接触 y -> KS 划分 -> 仅用训练集 y 去冗余 ----
    # (1) 基础清洗：只依据 X 本身（全空、常数），不看活性
    X0 = X_all.dropna(axis=1, thresh=int(0.95 * len(X_all)))
    X0 = X0.fillna(X0.median(numeric_only=True))
    X0 = X0.loc[:, X0.nunique() > 1]
    print(f"候选描述符 {X_all.shape[1]} -> 基础清洗后 {X0.shape[1]}")

    # (2) Kennard-Stone 划分：仅依赖描述符几何，不使用 y
    X0s = StandardScaler().fit_transform(X0.values)
    tr, te = kennard_stone(X0s, n_test=max(10, int(round(0.2 * len(df)))))
    print(f"KS 划分: 训练 {len(tr)} / 测试 {len(te)}")

    # (3) 去冗余（|r|>0.90，保留与活性相关更高者）只在训练集上计算，避免特征选择泄漏
    _, feats = clean_features(X0.iloc[tr], y[tr])
    print(f"训练集去冗余后候选描述符 {len(feats)} 个")

    # (4) 标准化参数只用训练集拟合
    scaler = StandardScaler().fit(X0[feats].iloc[tr])
    Xs = pd.DataFrame(scaler.transform(X0[feats]), columns=feats, index=X0.index)
    X_f = X0[feats]

    # 用 RF 重要性做特征选择，k 由训练集内部 5 折 CV 决定
    rf0 = RandomForestRegressor(n_estimators=500, random_state=SEED, n_jobs=-1).fit(Xs.iloc[tr], y[tr])
    imp = pd.Series(rf0.feature_importances_, index=feats).sort_values(ascending=False)

    best_k, best_q2, k_records = None, -9e9, []
    for k in [6, 8, 10, 12, 15, 20, 25]:
        if k > len(feats):
            continue
        Xk = Xs[imp.index[:k].tolist()].iloc[tr]
        m = PLSRegression(n_components=min(4, k))
        q2, rmse, _ = cv_eval(m, Xk, y[tr])
        k_records.append(dict(k=k, Q2=q2, RMSEcv=rmse))
        if q2 > best_q2:
            best_q2, best_k = q2, k
    print("特征数寻优:", {r["k"]: round(r["Q2"], 3) for r in k_records})
    print(f"选定 {best_k} 个描述符, 训练集 Q2={best_q2:.3f}")

    sel = imp.index[:best_k].tolist()
    Xsel = Xs[sel]
    Xtr, Xte = Xsel.iloc[tr], Xsel.iloc[te]
    ytr, yte = y[tr], y[te]

    rows, preds_store = [], {}
    for name, model in build_models(best_k).items():
        q2, rmsecv, cvpred = cv_eval(model, Xtr.values, ytr)
        model.fit(Xtr.values, ytr)
        yhat_tr = model.predict(Xtr.values)
        yhat_te = model.predict(Xte.values)
        trm, tem = metrics(ytr, yhat_tr), metrics(yte, yhat_te)
        rows.append(dict(model=name, Q2_CV=q2, RMSE_CV=rmsecv,
                         R2_train=trm["R2"], RMSE_train=trm["RMSE"],
                         R2_test=tem["R2"], RMSE_test=tem["RMSE"], MAE_test=tem["MAE"]))
        preds_store[name] = dict(model=model, yhat_tr=yhat_tr, yhat_te=yhat_te, cvpred=cvpred)
        print(f"{name}: Q2={q2:.3f} R2train={trm['R2']:.3f} R2test={tem['R2']:.3f} RMSEtest={tem['RMSE']:.3f}")

    perf = pd.DataFrame(rows)
    perf.to_csv(RES / "model_performance.csv", index=False)

    # 主模型预先按训练集内部 5 折 CV 的 Q² 确定；测试集只用于报告，不参与模型选择
    best_name = perf.sort_values("Q2_CV", ascending=False).iloc[0]["model"]
    print("主模型:", best_name, "(按训练集 CV Q² 选定)")
    B = preds_store[best_name]

    # ---------------- 适用性域 (William plot) ----------------
    Xd = np.hstack([np.ones((len(Xsel), 1)), Xsel.values])
    XtX_inv = np.linalg.pinv(Xd[tr].T @ Xd[tr])
    h_all = np.einsum("ij,jk,ik->i", Xd, XtX_inv, Xd)
    h_star = 3 * (best_k + 1) / len(tr)
    # Williams 图纵轴改用交叉验证残差（拟合残差会被 RF 压小，失去判别意义）
    cv_resid_raw = y[tr] - B["cvpred"]
    resid = cv_resid_raw / (np.std(cv_resid_raw) + 1e-9)
    ad = pd.DataFrame(dict(idx=df["idx"].values[tr], leverage=h_all[tr], std_resid=resid))

    # ---------------- Y-randomization ----------------
    rng = np.random.default_rng(SEED)
    yr_r2 = []
    for _ in range(30):
        yp = rng.permutation(ytr)
        m = build_models(best_k)[best_name]
        q2p, _, _ = cv_eval(m, Xtr.values, yp)
        yr_r2.append(q2p)
    yr_r2 = np.array(yr_r2)
    print(f"Y-randomization: 均值 Q2={yr_r2.mean():.3f}, 95%分位={np.percentile(yr_r2,95):.3f}, 真实 Q2={rows[0]['Q2_CV']:.3f}")

    # ---------------- 外部验证 DFMO + 查询分子 ----------------
    ext = pd.read_csv(DATA / "external_validation.csv")
    q = pd.read_csv(DATA / "query_molecules.csv")
    pred_rows = []
    for _, r in pd.concat([ext.assign(name="Eflornithine (DFMO, external)"), q], ignore_index=True).iterrows():
        Xq = descriptors([r["smiles"]])
        Xq = Xq.reindex(columns=feats).fillna(X0[feats].iloc[tr].median(numeric_only=True))
        Xq_s = pd.DataFrame(scaler.transform(Xq), columns=feats)[sel]
        Xd_q = np.hstack([np.ones((1, 1)), Xq_s.values])
        h_q = float(np.squeeze(Xd_q @ XtX_inv @ Xd_q.T))
        yq = {n: float(s["model"].predict(Xq_s.values)[0]) for n, s in preds_store.items()}
        ensemble = float(np.mean(list(yq.values())))
        pred_rows.append(dict(
            name=r.get("name"), role=r.get("role", "external validation"),
            pubchem_cid=r.get("pubchem_cid"), smiles=r["smiles"],
            pIC50_PLS=yq["PLS"], pIC50_RF=yq["RF"], pIC50_XGB=yq["XGB"],
            pIC50_ensemble=ensemble,
            IC50_uM_ensemble=10 ** (6 - ensemble),
            leverage=h_q, in_AD=bool(h_q <= h_star),
            observed_pIC50=r.get("pIC50", np.nan),
        ))
    pred = pd.DataFrame(pred_rows)
    pred.to_csv(RES / "predictions.csv", index=False)
    print("\n=== 预测结果 ===")
    print(pred[["name", "pIC50_PLS", "pIC50_RF", "pIC50_XGB", "pIC50_ensemble",
                "IC50_uM_ensemble", "leverage", "in_AD", "observed_pIC50"]]
          .to_string(index=False, float_format=lambda v: f"{v:.2f}"))

    summary = dict(
        n_total=int(len(df)), n_train=int(len(tr)), n_test=int(len(te)),
        n_descriptors_raw=int(X_all.shape[1]), n_descriptors_selected=int(best_k),
        n_descriptors_after_cleanup=int(X0.shape[1]), n_descriptors_after_dedup=int(len(feats)),
        selected_descriptors=sel,
        tier_counts=df["tier"].value_counts().to_dict(),
        pIC50_min=float(y.min()), pIC50_max=float(y.max()), pIC50_mean=float(y.mean()),
        h_star=float(h_star), best_model=best_name,
        y_randomization=dict(mean=float(yr_r2.mean()), p95=float(np.percentile(yr_r2, 95)),
                             max=float(yr_r2.max())),
        performance=perf.to_dict("records"),
        feature_selection_curve=k_records,
    )
    with open(RES / "qsar_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)

    # ---------------- 图（单图，Times New Roman） ----------------
    def save(fig, name):
        fig.tight_layout()
        fig.savefig(FIG / name, dpi=300, bbox_inches="tight", facecolor="white")
        plt.close(fig)

    # 1 活性分布
    fig, ax = plt.subplots(figsize=(2.205, 1.496))
    ax.hist(y, bins=14, color="#4C78A8", edgecolor="white")
    ax.set_xlabel("pIC$_{50}$"); ax.set_ylabel("Count")
    ax.set_title(f"ODC1 inhibitor dataset (n = {len(y)})")
    save(fig, "fig1_pIC50_distribution.png")

    # 2 实验-预测（主模型）
    fig, ax = plt.subplots(figsize=(2.205, 2.087))
    ax.scatter(ytr, B["yhat_tr"], s=34, c="#4C78A8", label="Training", edgecolor="white", linewidth=0.4)
    ax.scatter(yte, B["yhat_te"], s=44, c="#E45756", marker="^", label="Test", edgecolor="white", linewidth=0.4)
    lo = min(y.min(), B["yhat_tr"].min(), B["yhat_te"].min()) - 0.3
    hi = max(y.max(), B["yhat_tr"].max(), B["yhat_te"].max()) + 0.3
    ax.plot([lo, hi], [lo, hi], "k--", lw=0.9)
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_xlabel("Observed pIC$_{50}$"); ax.set_ylabel("Predicted pIC$_{50}$")
    ax.legend(frameon=False, loc="upper left")
    _pf = perf[perf.model == best_name].iloc[0]
    ax.set_title(f"{best_name}: $Q^2$={_pf.Q2_CV:.2f}, $R^2_{{test}}$={_pf.R2_test:.2f}")
    save(fig, "fig2_observed_vs_predicted.png")

    # 3 Y-randomization
    fig, ax = plt.subplots(figsize=(2.205, 1.693))
    ax.hist(yr_r2, bins=10, color="#B4B4B4", edgecolor="white")
    ax.axvline(summary["y_randomization"]["p95"], color="#E45756", ls="--", lw=0.9,
               label=f"95% of random ({summary['y_randomization']['p95']:.2f})")
    ax.axvline(perf.loc[perf.model == best_name, "Q2_CV"].iloc[0], color="#1D9E75", lw=1.1,
               label=f"Real model ({perf.loc[perf.model == best_name, 'Q2_CV'].iloc[0]:.2f})")
    ax.legend(frameon=False, loc="upper right", fontsize=5.5,
              labels=["Random 95th pct", "Real model"], handlelength=1.0)
    ax.set_xlabel("Q$^2$ of randomized models"); ax.set_ylabel("Frequency")
    ax.set_title("Y-randomization (30 runs)")
    save(fig, "fig3_y_randomization.png")

    # 4 William plot：横轴按最大杠杆动态延伸，确保全部查询分子可入图；纵轴为 CV 残差
    xmax = max(0.95, float(ad.leverage.max()) * 1.12, float(pred.leverage.max()) * 1.12)
    fig, ax = plt.subplots(figsize=(3.346, 2.441))
    ax.axvspan(0, h_star, color="#F0F6F2", zorder=0)
    ax.axvspan(h_star, xmax, color="#FBEEEC", zorder=0)
    ax.axhline(3, color="grey", ls=":", lw=0.8, zorder=1)
    ax.axhline(-3, color="grey", ls=":", lw=0.8, zorder=1)
    ax.axvline(h_star, color="#E45756", ls="--", lw=0.9, zorder=2)
    ax.scatter(ad.leverage, ad.std_resid, s=26, c="#54A24B", edgecolor="white",
               linewidth=0.35, label="Training (CV residual)", zorder=3)
    pq = pred.drop_duplicates(subset=["smiles"])
    ybase = -3.55
    ax.scatter(pq.leverage, np.full(len(pq), ybase), marker="^", s=34, c="#C0392B",
               edgecolor="white", linewidth=0.35, label=f"Query (n={len(pq)}, h only)", zorder=4)
    _pue = pq[pq.name.str.contains("Puerarin", na=False)]
    if len(_pue):
        _x = float(_pue.leverage.iloc[0])
        ax.annotate("Puerarin", xy=(_x, ybase + 0.18), xytext=(_x + 0.12, ybase + 1.15),
                    fontsize=5.6, color="#C0392B", ha="left", va="center",
                    arrowprops=dict(arrowstyle="-", lw=0.5, color="#C0392B",
                                    shrinkA=0, shrinkB=2))
    ax.set_xlim(0, xmax)
    ax.set_ylim(-4.05, 4.35)
    ax.text(h_star * 0.5, 4.0, "Inside AD", ha="center", va="center", fontsize=6, color="#20794A")
    ax.text((h_star + xmax) / 2, 4.0, f"Outside AD (h* = {h_star:.2f})", ha="center", va="center",
            fontsize=6, color="#B23A2E")
    ax.set_xlabel("Leverage ($h$)")
    ax.set_ylabel("Standardized CV residual")
    ax.set_title("Applicability domain")
    ax.legend(frameon=False, fontsize=5.4, loc="center right", handlelength=1.0,
              borderaxespad=0.3, labelspacing=0.4)
    save(fig, "fig4_william_plot.png")

    # 5 描述符重要性
    fig, ax = plt.subplots(figsize=(3.346, 2.598))
    top = imp.head(best_k).sort_values()
    ax.barh(range(len(top)), top.values, color="#4C78A8")
    ax.set_yticks(range(len(top))); ax.set_yticklabels(top.index, fontsize=8)
    ax.set_xlabel("RF importance")
    save(fig, "fig5_descriptor_importance.png")

    # 6 预测对比
    fig, ax = plt.subplots(figsize=(3.346, 2.047))
    SHORT = {"APA (3-aminooxy-1-aminopropane)": "APA", "Eflornithine (DFMO)": "DFMO",
             "Eflornithine (DFMO, external)": "DFMO"}
    d = pred.drop_duplicates(subset=["smiles"]).sort_values("pIC50_ensemble")
    colors = ["#E45756" if not v else "#4C78A8" for v in d.in_AD]
    ax.barh(range(len(d)), d.pIC50_ensemble, color=colors, height=0.68)
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels([SHORT.get(n, n) for n in d.name], fontsize=6.5)
    ax.set_xlabel("Predicted pIC$_{50}$")
    ax.set_xlim(0, 7.6)
    for i, v in enumerate(d.pIC50_ensemble):
        ax.text(v + 0.06, i, f"{v:.2f}", va="center", fontsize=6)
    save(fig, "fig6_query_predictions.png")
    print("\n图已输出至 figures/")


if __name__ == "__main__":
    main()
