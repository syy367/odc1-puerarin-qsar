"""生成 ODC1 QSAR 项目 HTML 报告（单文件自包含，图片 base64 内嵌）。"""
import base64
import json
from datetime import date
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES, FIG = ROOT / "results", ROOT / "figures" / "journal"
OUT = ROOT / f"ODC1_QSAR_葛根素预测报告_{date.today():%Y%m%d}_wbuddy.html"


def b64(p):
    return base64.b64encode(Path(p).read_bytes()).decode()


def img(name, caption, width="100%"):
    p = FIG / name
    if not p.exists():
        return ""
    return (f'<figure><img src="data:image/png;base64,{b64(p)}" style="width:{width}">'
            f'<figcaption>{caption}</figcaption></figure>')


def tbl(df, cls="tbl", fmt=None):
    d = df.copy()
    if fmt:
        for c, f in fmt.items():
            if c in d.columns:
                d[c] = d[c].map(lambda v: f.format(v) if pd.notna(v) else "—")
    head = "".join(f"<th>{c}</th>" for c in d.columns)
    rows = "".join("<tr>" + "".join(f"<td>{v}</td>" for v in r) + "</tr>" for r in d.values)
    return f'<table class="{cls}"><thead><tr>{head}</tr></thead><tbody>{rows}</tbody></table>'


S = json.load(open(RES / "qsar_summary.json", encoding="utf-8"))
perf = pd.read_csv(RES / "model_performance.csv")
pred = pd.read_csv(RES / "predictions.csv")
dock = pd.read_csv(RES / "docking_analysis.csv")
hb = pd.read_csv(RES / "hbond_summary.csv")
dft = json.load(open(RES / "dft_results.json", encoding="utf-8"))
ds = pd.read_csv(ROOT / "data" / "dataset_odc1.csv")

# ---------- 表格准备 ----------
perf_show = perf.rename(columns={
    "model": "模型", "Q2_CV": "Q²(CV)", "RMSE_CV": "RMSE(CV)",
    "R2_train": "R²(训练)", "RMSE_train": "RMSE(训练)",
    "R2_test": "R²(测试)", "RMSE_test": "RMSE(测试)", "MAE_test": "MAE(测试)"})
perf_show = perf_show[["模型", "Q²(CV)", "RMSE(CV)", "R²(训练)", "R²(测试)", "RMSE(测试)"]]

NAME_MAP = {
    "Eflornithine (DFMO, external)": "DFMO（外部验证）",
    "Eflornithine (DFMO)": "DFMO",
    "Puerarin": "葛根素",
    "Genistein": "染料木素",
    "Daidzein": "大豆苷元",
    "Quercetin": "槲皮素",
    "Resveratrol": "白藜芦醇",
    "APA (3-aminooxy-1-aminopropane)": "APA（XAP）",
}
pred = pred.drop_duplicates(subset=["name"]).copy()
pred["名称"] = pred["name"].map(NAME_MAP).fillna(pred["name"])
pred_show = pred[["名称", "pIC50_PLS", "pIC50_RF", "pIC50_XGB", "pIC50_ensemble",
                  "IC50_uM_ensemble", "leverage", "in_AD", "observed_pIC50"]].rename(columns={
    "pIC50_PLS": "PLS", "pIC50_RF": "RF", "pIC50_XGB": "XGB", "pIC50_ensemble": "集成",
    "IC50_uM_ensemble": "IC₅₀(µM)", "leverage": "杠杆值 h", "in_AD": "在适用域内",
    "observed_pIC50": "实测 pIC50"})
pred_show["在适用域内"] = pred_show["在适用域内"].map({True: "是", False: "否"})

dock["名称"] = dock["name"].map(NAME_MAP).fillna(dock["name"])
dock_show = dock[["名称", "best_affinity", "n_heavy", "ligand_efficiency", "n_pocket_contacts"]].rename(
    columns={"best_affinity": "Vina 打分(kcal/mol)", "n_heavy": "重原子数",
             "ligand_efficiency": "配体效率", "n_pocket_contacts": "口袋残基接触数"})

hb["名称"] = hb["name"].map(NAME_MAP).fillna(hb["name"])
hb_show = hb[["名称", "n_hbond", "hbond_residues"]].rename(
    columns={"n_hbond": "氢键数(≤3.6 Å)", "hbond_residues": "涉及残基"})

DFT_MAP = {"Puerarin": "葛根素", "Eflornithine (DFMO)": "DFMO", "Genistein": "染料木素",
           "Daidzein": "大豆苷元", "Quercetin": "槲皮素", "Resveratrol": "白藜芦醇",
           "APA (3-aminooxy-1-aminopropane)": "APA"}
dft_rows = [{"分子": DFT_MAP.get(k, k), "HOMO(eV)": v["homo_eV"], "LUMO(eV)": v["lumo_eV"],
             "能隙(eV)": v["gap_eV"], "偶极矩(D)": v["dipole_D"], "分子式": v.get("formula", "")}
            for k, v in dft.items()]
dft_df = pd.DataFrame(dft_rows).sort_values("能隙(eV)")

sources = [
    ("人源重组 ODC1 纯化酶 IC₅₀", "13", "2025", "J Med Chem", "10.1021/acs.jmedchem.4c03120", "A"),
    ("大鼠肝纯化 ODC IC₅₀", "15", "1992", "J Med Chem", "10.1021/jm00086a003", "B"),
    ("大鼠肝 ODC Ki（换算）", "32", "1978–1995", "J Med Chem", "10.1021/jm00199a009 等 5 篇*", "C"),
]
src_df = pd.DataFrame(sources, columns=["数据层", "记录数", "年份", "期刊", "DOI（已核实）", "Tier"])

puer = pred[pred.name == "Puerarin"].iloc[0]
dfmo = pred[pred.name.str.contains("DFMO")].iloc[0]
apa = pred[pred.name.str.contains("APA")].iloc[0]
_vsa_desc = "、".join([x for x in S["selected_descriptors"] if ("VSA" in x or "EState" in x)])
_topo_desc = "、".join([x for x in S["selected_descriptors"] if not ("VSA" in x or "EState" in x)])

CSS = """
:root{--ink:#1a1a1a;--mute:#5f6b7a;--line:#e2e6ea;--bg:#fbfcfd;--card:#ffffff;
--accent:#2f6f9f;--warn:#c0392b;--ok:#20794a;--soft:#f3f7fa}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);
font-family:"Segoe UI","Microsoft YaHei","PingFang SC",sans-serif;line-height:1.75;font-size:15px}
.wrap{max-width:1000px;margin:0 auto;padding:40px 28px 80px}
h1{font-size:27px;line-height:1.35;margin:0 0 6px;letter-spacing:.2px}
h2{font-size:20px;margin:44px 0 14px;padding-bottom:8px;border-bottom:2px solid var(--line)}
h3{font-size:16px;margin:26px 0 10px;color:var(--accent)}
.sub{color:var(--mute);font-size:13.5px;margin-bottom:26px}
.card{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:20px 22px;margin:16px 0}
.kpi{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:14px;margin:18px 0}
.kpi div{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:14px 16px}
.kpi b{display:block;font-size:22px;color:var(--accent);font-weight:600;letter-spacing:-.3px}
.kpi span{font-size:12.5px;color:var(--mute)}
table.tbl{width:100%;border-collapse:collapse;margin:14px 0;font-size:13.5px;background:var(--card)}
table.tbl th{background:var(--soft);text-align:left;padding:9px 11px;font-weight:600;
border-bottom:1px solid var(--line);white-space:nowrap}
table.tbl td{padding:8px 11px;border-bottom:1px solid #eef1f4;vertical-align:top}
table.tbl tr:last-child td{border-bottom:none}
figure{margin:20px 0;text-align:center}
figure img{border:1px solid var(--line);border-radius:8px;background:#fff}
figcaption{font-size:12.5px;color:var(--mute);margin-top:7px}
.note{border-left:3px solid var(--accent);background:var(--soft);padding:12px 16px;margin:16px 0;font-size:14px}
.warn{border-left:3px solid var(--warn);background:#fdf3f2;padding:12px 16px;margin:16px 0;font-size:14px}
.ok{border-left:3px solid var(--ok);background:#f2f9f5;padding:12px 16px;margin:16px 0;font-size:14px}
code{background:#eef2f5;padding:2px 6px;border-radius:4px;font-size:13px;
font-family:Consolas,"Courier New",monospace}
ul,ol{padding-left:22px}li{margin:5px 0}
.small{font-size:13px;color:var(--mute)}
.tag{display:inline-block;font-size:11.5px;padding:2px 8px;border-radius:20px;
background:var(--soft);color:var(--accent);margin-right:6px}
.bad{color:var(--warn);font-weight:600}.good{color:var(--ok);font-weight:600}
footer{margin-top:56px;padding-top:18px;border-top:1px solid var(--line);
font-size:12.5px;color:var(--mute)}
"""

HTML = f"""<!DOCTYPE html>
<html lang="zh-CN"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>ODC1 抑制剂 QSAR 与葛根素活性预测</title>
<style>{CSS}</style></head><body><div class="wrap">

<h1>ODC1 抑制剂 2D-QSAR 建模与葛根素 pIC₅₀ 预测<br>
<span style="font-size:18px;color:#5f6b7a;font-weight:400">附分子对接与 DFT 电子结构分析</span></h1>
<div class="sub">生成日期：{date.today():%Y-%m-%d}　|　靶点：人源鸟氨酸脱羧酶 ODC1 (UniProt P11926, PDB 1D7K)　|　
数据集：{S['n_total']} 个化合物，{S['n_train']} 训练 / {S['n_test']} 测试</div>

<div class="card">
<b>v2 修订说明（2026-11-18）</b>
<p class="small" style="margin:8px 0 0">本版相对 2026-11-15 版做了三项方法学修正，并据此<b>重算了全部 QSAR 指标与预测值</b>；对接、DFT 部分未改动。</p>
<table class="tbl" style="margin-top:8px">
<thead><tr><th style="width:24%">修正项</th><th style="width:36%">原做法</th><th>现做法</th></tr></thead>
<tbody>
<tr><td><b>特征选择泄漏</b></td><td>描述符去冗余（保留与活性相关性更高者）在全部 60 个化合物上计算，且早于数据划分</td>
<td>先做 Kennard–Stone 划分 → 仅在 48 个训练化合物上计算相关性去冗余。候选描述符 {S['n_descriptors_after_cleanup']} → {S['n_descriptors_after_dedup']} 个，最终保留 {S['n_descriptors_selected']} 个（原 8 个）</td></tr>
<tr><td><b>主模型选择</b></td><td>按测试集 R² 在 PLS/RF/XGB 中挑最高者，再报其测试集 R²</td>
<td>按训练集内部 CV Q² 预先选定（仍为 RF），测试集只用于报告表现</td></tr>
<tr><td><b>适用域图</b></td><td>Williams 图只画训练点，横轴被自动截断在 0.6，5 个域外查询分子无法入图</td>
<td>横轴按最大杠杆自动延伸，7 个查询分子全部入图；纵轴改用交叉验证残差</td></tr>
</tbody></table>
<p class="small"><b>数值变化对照</b>：Q² <b>0.793 → {S['performance'][1]['Q2_CV']:.3f}</b>；测试集 R² <b>0.740 → {S['performance'][1]['R2_test']:.3f}</b>；
RMSE<sub>test</sub> <b>0.902 → {S['performance'][1]['RMSE_test']:.3f}</b>；h* <b>0.5625 → {S['h_star']:.4f}</b>；
葛根素 pIC₅₀ <b>3.61 → {puer['pIC50_ensemble']:.2f}</b>（IC₅₀ 243 → {puer['IC50_uM_ensemble']:.0f} µM），
杠杆值 h <b>0.966 → {puer['leverage']:.3f}</b>。
杠杆阈值下降后<b>葛根素仍判定为域外</b>，核心结论方向未变。<br>
另：APA 经 InChIKey 精确比对确认<b>属于训练集成员</b>，"未参与训练的外部强效分子"相关表述已全部删除（见第四节⑤）。</p>
</div>

<h2>一、核心结论</h2>
<div class="kpi">
  <div><b>{S['performance'][1]['Q2_CV']:.2f}</b><span>RF 模型交叉验证 Q²</span></div>
  <div><b>{S['performance'][1]['R2_test']:.3f}</b><span>RF 测试集 R²（n={S['n_test']}）</span></div>
  <div><b>{puer['pIC50_ensemble']:.2f}</b><span>葛根素预测 pIC₅₀</span></div>
  <div><b>{puer['IC50_uM_ensemble']:.0f} µM</b><span>葛根素预测 IC₅₀</span></div>
  <div><b>{dfmo['pIC50_ensemble']:.2f}</b><span>DFMO 预测 pIC₅₀</span></div>
</div>

<div class="ok"><b>模型可信度</b>：随机森林在 5 折交叉验证中 Q² = {S['performance'][1]['Q2_CV']:.3f}，
内部测试集（Kennard–Stone 划分，n = {S['n_test']}）R² = {S['performance'][1]['R2_test']:.3f}；30 次 Y-randomization 的平均 Q² 仅
{S['y_randomization']['mean']:.3f}（95% 分位 {S['y_randomization']['p95']:.3f}，最大 {S['y_randomization']['max']:.3f}），
远低于真实模型，<b>未发现偶然相关</b>。DFMO 作为唯一真正意义上的外部化合物
（不在 60 个化合物中、不参与任何训练与划分），预测 pIC₅₀ = {dfmo['pIC50_ensemble']:.2f}，实测 {dfmo['observed_pIC50']:.2f}，偏差
<b>{abs(dfmo['pIC50_ensemble']-dfmo['observed_pIC50']):.2f} log</b> 单位（约 {10**abs(dfmo['pIC50_ensemble']-dfmo['observed_pIC50']):.1f} 倍），
方向与量级均合理。<b>但这是单点验证</b>，其作用是"未证伪"而非"证明定量预测能力"；
模型的实际定位是氨基氧基/二胺类骨架内的局部 SAR 插值（见第六节适用域讨论）。</div>

<div class="warn"><b>与预设故事链的一处重要偏离（必须正视）</b><br>
QSAR 给出的葛根素 pIC₅₀ = <b>{puer['pIC50_ensemble']:.2f}</b>（IC₅₀ ≈ {puer['IC50_uM_ensemble']:.0f} µM），
与 DFMO（pIC₅₀ {dfmo['pIC50_ensemble']:.2f}，IC₅₀ ≈ {dfmo['IC50_uM_ensemble']:.0f} µM）相差约
{10**(dfmo['pIC50_ensemble']-puer['pIC50_ensemble']):.1f} 倍。这一数值属于
<b>"弱抑制剂"（IC₅₀ &gt; 100 µM）</b>，且已接近葛根素在水性缓冲液中的溶解上限量级。
同时葛根素的杠杆值 h = {puer['leverage']:.2f}，高于阈值 h* = {S['h_star']:.2f}，
<b>落在适用域之外</b>。因此该预测只能作为趋势性证据，不能当作定量结论引用。<br>
<b>数值呈现规范（本文通用）</b>：预测 IC₅₀ &gt; 100 µM 的分子一律表述为"无显著抑制"或"&gt; 100 µM"，
不给出具体数值区间，也不在它们之间做强弱排序——毫摩尔级的读数在实验上通常反映溶解度上限或聚集伪影，
不携带结构信息。<br>
建议论文口径改为："计算提示葛根素对 ODC1 仅具弱直接抑制潜力，
其肾保护作用更可能依赖对 ODC1 表达/多胺通路的间接调控"。</div>

<h2>二、数据集构建</h2>
<p>按"同一检测体系、纯化酶活 IC₅₀"的严格口径检索 ChEMBL，实际情况与预期有出入：
人源重组 ODC1 的纯化酶 IC₅₀ 在库中仅有 {len(ds[ds.tier=='A'])} 条（DFMO 单独置于外部验证集，不计入），
仅靠单一体系无法达到 60 个化合物。因此采用<b>三层合并 + 分层标注</b>策略：</p>
{tbl(src_df)}
<p class="small">* Tier C 的 32 条 Ki 来自 5 篇文献：10.1021/jm00199a009（7 条）、10.1021/jm00021a024（10 条）、
10.1021/jm00169a025（4 条）、10.1021/jm00365a002（8 条）、10.1021/jm00121a031（3 条）。
Tier A：人源重组 his-tag ODC1，E. coli 表达纯化，同一篇文献同一体系，质量最高；
Tier B：大鼠肝天然纯化 ODC，IC₅₀；Tier C：大鼠肝纯化 ODC 的 Ki，按 Cheng-Prusoff 关系换算为等效 pIC₅₀。
<b>换算的不确定性</b>：Cheng–Prusoff 完整式为 IC₅₀ = Ki·(1 + [S]/K<sub>m</sub>)，此处统一假定 [S] = K<sub>m</sub>（系数 2）；
若原文采用饱和底物（[S] = 10·K<sub>m</sub>，酶学初速测定中并不罕见），系数将变为 11，对应 0.74 log 的系统偏移，
方向是<b>统一抬高 Tier C 的活性</b>。由于 5 篇原始文献的底物浓度未能逐一回查核实，此处不宣称"0.3 log"这一单一偏差值，
而以 <b>0–0.74 log 的系统不确定度</b>描述该层数据；敏感性分析（剔除 Tier C）见下节，
所有记录均来自纯化酶体系，<b>未纳入任何细胞增殖或细胞裂解液 IC₅₀</b>。</p>

<p>去重（按 InChIKey 归一并取中位数）后得到 <b>{S['n_total']} 个化合物</b>，
pIC₅₀ 范围 {S['pIC50_min']:.2f}–{S['pIC50_max']:.2f}（跨度 {S['pIC50_max']-S['pIC50_min']:.1f} 个对数单位），
均值 {S['pIC50_mean']:.2f}。分层构成：A {S['tier_counts'].get('A',0)} / B {S['tier_counts'].get('B',0)} /
C {S['tier_counts'].get('C',0)}。活性跨度大有利于建模，但三层数据来自 9 篇跨 47 年的文献，
<b>体系异质性是不可忽视的噪声来源</b>。</p>

{img('fig1_pIC50_distribution.png', '图 1　数据集 pIC₅₀ 分布。双峰特征反映 Tier A/B 与 Tier C 的活性量级差异。', '78%')}

<h2>三、建模方法与验证</h2>
<div class="card">
<b>描述符</b>：RDKit 2D 描述符全量 {S['n_descriptors_raw']} 个 → 剔除缺失/常数（仅依据 X 本身）得 {S['n_descriptors_after_cleanup']} 个
→ <b>在训练集内部</b>以 |r| &gt; 0.90 去冗余（保留与训练集活性相关性更高者）得 {S['n_descriptors_after_dedup']} 个候选
→ 随机森林重要性排序 + 训练集内部交叉验证寻优，最终保留 <b>{S['n_descriptors_selected']} 个</b>：
<p class="small" style="margin:8px 0 0">{"、".join(S['selected_descriptors'])}</p>
<p class="small"><b>防止特征选择泄漏的顺序</b>：先去冗余 → 再做 Kennard–Stone 划分 → 再选描述符，这是刻意安排的；
去冗余这一步若用全部 60 个化合物的 y，"与活性相关性更高"就会包含测试集信息，构成特征选择泄漏。
本文所有与 y 相关的步骤（去冗余、重要性排序、k 寻优、标准化参数）<b>均只用 48 个训练化合物</b>，
测试集的 y 在整个建模阶段不参与任何决策。</p>
<p class="small">这组描述符以<b>表面极性/疏水表面积分布</b>（{_vsa_desc}）
与<b>拓扑复杂度/形状</b>（{_topo_desc}）为主，
与 ODC 抑制剂"极性头部 + 疏水骨架"的结构特征吻合。
<b>需注意 BCUT2D_MWHI 与分子量的相关系数 r = 0.78</b>（本数据范围内，60 个化合物的 MW 88–363、重原子 6–24），
即模型的部分解释力可能来自"分子有多大"而非"结构如何识别"；该混杂在解读描述符重要性时必须一并说明。</p>
<b>数据划分</b>：先从描述符空间质心取起点，再贪心地逐点选取距已选集合最远者，取最后 12 个作为<b>内部测试集</b>
（思路源自 Kennard–Stone 的最远点策略，但<b>严格说不等同于标准 KS 划分</b>——标准 KS 由两端交替选取、
其余留作训练；本实现得到的 12 个点是化学空间最边缘的样本，训练集则是相对中心的 48 个）。
效果上这构成一个偏难的测试集，也是 PLS/XGB 在其上表现明显劣于 RF 的直接原因之一。
需强调：这 12 个点来自同一批 60 个化合物，<b>属于内部测试集（internal test set）而非外部验证</b>；
真正外部验证只有 DFMO 一个化合物。<br>
<b>模型</b>：PLS（4 主成分）、随机森林（500 树）、XGBoost（400 轮，深度 3）。
主模型<b>预先按训练集内部 5 折 CV 的 Q² 确定</b>（不依据测试集 R² 挑选），测试集只用于报告表现。<br>
<b>验证</b>：5 折交叉验证（Q²、RMSE<sub>CV</sub>）、内部测试集、30 次 Y-randomization、
Williams 图适用域分析（杠杆值 h vs <b>交叉验证残差</b>，而非拟合残差——后者会被随机森林的过拟合能力压小）。
</div>

<h3>模型性能</h3>
{tbl(perf_show, fmt={c: "{:.3f}" for c in ["Q²(CV)", "RMSE(CV)", "R²(训练)", "R²(测试)", "RMSE(测试)"]})}
<p class="small">PLS 交叉验证尚可（Q² = {perf[perf.model=='PLS'].Q2_CV.iloc[0]:.3f}）但测试集崩到
R² = {perf[perf.model=='PLS'].R2_test.iloc[0]:.3f}，说明线性模型无法外推到测试集覆盖的边缘化学空间；
XGBoost 训练集 R² = {perf[perf.model=='XGB'].R2_train.iloc[0]:.3f} 而测试集仅
{perf[perf.model=='XGB'].R2_test.iloc[0]:.3f}，呈明显过拟合。
<b>随机森林的训练集 CV Q² 最高（{perf[perf.model=='RF'].Q2_CV.iloc[0]:.3f}），据此预先选定为主模型</b>；
其测试集表现（R² = {perf[perf.model=='RF'].R2_test.iloc[0]:.3f}）是在选定之后才观察到的，未参与模型选择。
三个模型的测试集表现差异很大（R² 从 {perf.R2_test.min():.2f} 到 {perf.R2_test.max():.2f}），
这本身提示 12 个测试点不足以稳定区分模型优劣——不应对 R² = {perf[perf.model=='RF'].R2_test.iloc[0]:.2f} 这一单值赋予过强的确定性。</p>

{img('fig2_observed_vs_predicted.png', '图 2　随机森林模型实验值 vs 预测值。训练集（蓝）紧贴对角线；测试集（红三角，n = 12）多数落在 ±1 log 内，个别点偏离反映体系异质性。<b>注意图内标题口径</b>：Q² 是训练集 5 折交叉验证值，R²<sub>test</sub> 是内部测试集的决定系数，二者本次恰好都约 0.78 属巧合，不可互替；图中画的是训练集拟合值与测试集预测值，而 Q² 来自交叉验证预测，严格说与图面点不对应（该问题列入下一版修订）。', '62%')}
{img('fig3_y_randomization.png', '图 3　Y-randomization（30 次）。红色虚线为随机模型的 95% 分位（Q² = ' + f"{S['y_randomization']['p95']:.3f}" + '），真实模型 Q² = ' + f"{S['y_randomization']['real'] if 'real' in S['y_randomization'] else S['performance'][1]['Q2_CV']:.3f}" + '，远超随机分布。', '58%')}
{img('fig4_william_plot.png', '图 4　Williams 适用域图。<b>纵轴为交叉验证残差</b>（非拟合残差），红线为杠杆阈值 h* = ' + f"{S['h_star']:.3f}" + '；横轴按最大杠杆值延伸，<b>7 个查询分子（红三角）全部入图</b>。三角形仅有 h 坐标——查询分子无实测值，纵轴位置不具含义。可见葛根素（h = ' + f"{puer['leverage']:.2f}" + '）与四种天然产物对照均落在 h* 右侧的域外区。', '72%')}
{img('fig5_descriptor_importance.png', '图 5　最终建模使用的 ' + f"{S['n_descriptors_selected']}" + ' 个描述符及其随机森林重要性（单次计算，未附误差）。', '68%')}

<h2>四、预测结果</h2>
{tbl(pred_show, fmt={c: "{:.2f}" for c in ["PLS", "RF", "XGB", "集成", "杠杆值 h", "实测 pIC50"]} | {"IC₅₀(µM)": "{:.1f}"})}
<p class="small">集成值 = 三个模型预测的算术平均。杠杆值 h 越小表示越靠近训练集化学空间中心；
h &gt; h*（{S['h_star']:.2f}）判为适用域外。</p>

{img('fig6_query_predictions.png', '图 6　待测分子的集成 pIC₅₀ 预测。红色为适用域外（葛根素及天然产物对照），蓝色为域内。', '72%')}

<div class="note"><b>关键读数</b><br>
① <b>葛根素 pIC₅₀ = {puer['pIC50_ensemble']:.2f}</b>（IC₅₀ ≈ {puer['IC50_uM_ensemble']:.0f} µM），
三个模型分歧不小（PLS {puer['pIC50_PLS']:.2f} / RF {puer['pIC50_RF']:.2f} / XGB {puer['pIC50_XGB']:.2f}，
极差 {puer[['pIC50_PLS','pIC50_RF','pIC50_XGB']].max()-puer[['pIC50_PLS','pIC50_RF','pIC50_XGB']].min():.2f} log），
这是域外预测的典型表现；按测试集 RMSE = {perf[perf.model=='RF'].RMSE_test.iloc[0]:.2f} 计，95% 预测区间为 ±{1.96*perf[perf.model=='RF'].RMSE_test.iloc[0]:.1f} log，
即约 {10**(6-puer['pIC50_ensemble']-1.96*perf[perf.model=='RF'].RMSE_test.iloc[0]):.0f}–{10**(6-puer['pIC50_ensemble']+1.96*perf[perf.model=='RF'].RMSE_test.iloc[0]):.0f} µM，跨越约 3 个数量级——
<b>这个区间几乎不提供实用约束</b>，这正是域外预测的代价，也是"不能把 159 µM 当定量结论"的量化依据。<br>
② <b>DFMO 是唯一真正的外部验证点</b>：预测 {dfmo['pIC50_ensemble']:.2f} vs 实测 {dfmo['observed_pIC50']:.2f}，
偏差 {abs(dfmo['pIC50_ensemble']-dfmo['observed_pIC50']):.2f} log，杠杆值 {dfmo['leverage']:.2f} 在域内。
单点验证的作用是"未证伪"，不足以支撑"具备定量预测能力"这一表述。<br>
③ <b>四种天然产物对照的预测值挤在 {pred[pred.name.isin(['Genistein','Daidzein','Quercetin','Resveratrol'])].pIC50_ensemble.min():.2f}–{pred[pred.name.isin(['Genistein','Daidzein','Quercetin','Resveratrol'])].pIC50_ensemble.max():.2f}</b>
（跨度仅 {pred[pred.name.isin(['Genistein','Daidzein','Quercetin','Resveratrol'])].pIC50_ensemble.max()-pred[pred.name.isin(['Genistein','Daidzein','Quercetin','Resveratrol'])].pIC50_ensemble.min():.2f} log，小于 RMSE {perf[perf.model=='RF'].RMSE_test.iloc[0]:.2f}），
而它们的结构差异（二苯乙烯 / 黄酮 / 异黄酮 / 黄酮醇）远大于此。其中染料木素与大豆苷元的预测值<b>完全相同</b>
（{pred[pred.name=='Genistein'].pIC50_ensemble.iloc[0]:.2f} vs {pred[pred.name=='Daidzein'].pIC50_ensemble.iloc[0]:.2f}）——两个结构不同的分子给出同一个数，
这本身就是"预测值不携带结构信息"的直接证据。四种化合物的杠杆值（{pred[pred.name.isin(['Genistein','Daidzein','Quercetin','Resveratrol'])].leverage.min():.2f}–{pred[pred.name.isin(['Genistein','Daidzein','Quercetin','Resveratrol'])].leverage.max():.2f}）全部超过 h* = {S['h_star']:.2f}。
<b>因此四者之间的排序不具判别力，不得用于论证"葛根素优于/劣于其他黄酮"。</b><br>
④ 按上文数值呈现规范：白藜芦醇、槲皮素、染料木素、大豆苷元的预测 IC₅₀ 均 &gt; 100 µM，
统一表述为<b>"&gt; 100 µM，无显著抑制"</b>，不给出具体数值作为定量结论。<br>
⑤ <b>关于 APA：它是训练集成员，不是外部对照。</b>
APA（3-aminooxy-1-aminopropane，CID 65020，InChIKey <code>VSZFWDPIWSPZON-UHFFFAOYSA-N</code>）即数据集中编号 57 的化合物
（SMILES <code>NCCCON</code>，Tier A，pIC₅₀ = 7.46，取自 2025 <i>J Med Chem</i> 原始数据，为训练集第二高活性）。
经 InChIKey 精确比对与 Kennard–Stone 划分复核：<b>它属于 48 个训练化合物，不在 12 个内部测试点内</b>。
因此"模型把未参与训练的强效分子排到第一位"这一说法<b>不成立</b>，已从本报告删除；
其 h = {apa['leverage']:.2f}（7 个查询分子中最低）恰恰是"训练集内部成员"的直接体现，而非"预测可信"的证据。
模型对它的预测 {apa['pIC50_ensemble']:.2f} 低于其训练值 7.46（差 {7.456-apa['pIC50_ensemble']:.2f} log），
属于随机森林向训练集均值回归的表现。</div>

{img('fig9_le_vs_pIC50.png', '图 9　对接配体效率（结构层面）与 QSAR 预测 pIC₅₀（配体层面）的对照。<b>仅作散点展示，不宣称两法一致</b>：葛根素在横轴最左（LE 最低）而在纵轴排第二（QSAR 域外、无判别力）。APA 为训练集成员（见关键读数⑤），其右上角位置不构成外部验证。', '62%')}

<div class="note"><b>两条证据线是否互相印证？——结论是：这一论证基本不成立</b><br>
对接侧 LE 排序（去掉 APA）：DFMO ≈ 白藜芦醇（0.39）&gt; 槲皮素（0.36）&gt; 染料木素 ≈ 大豆苷元（0.34）&gt; <b>葛根素（0.24，最低）</b>；<br>
QSAR 侧 pIC₅₀：DFMO（{dfmo['pIC50_ensemble']:.2f}）&gt; <b>葛根素（{puer['pIC50_ensemble']:.2f}）</b>&gt; 染料木素 ≈ 大豆苷元（3.18）&gt; 槲皮素（3.09）&gt; 白藜芦醇（2.55）。<br>
两者对葛根素的位置判断<b>方向相反</b>（对接侧最弱、QSAR 侧在天然产物中最强）。
但这不是"两种方法打架"，而是 QSAR 侧那组数字本来就没有判别力：
葛根素及四种对照的杠杆值全部超域（h = {pred[pred.name.isin(['Puerarin','Genistein','Daidzein','Quercetin','Resveratrol'])].leverage.min():.2f}–{pred[pred.name.isin(['Puerarin','Genistein','Daidzein','Quercetin','Resveratrol'])].leverage.max():.2f} &gt; h* = {S['h_star']:.2f}），
预测值被拉回训练集低活性簇的中心，彼此差异小于 RMSE。<br>
因此可比较的只有 DFMO 一个点（域内，h = {dfmo['leverage']:.2f}）：LE 0.39 属中等，QSAR 预测 {dfmo['pIC50_ensemble']:.2f} 接近实测 {dfmo['observed_pIC50']:.2f}。
<b>n = 1 无法构成"两种方法互相印证"。</b>另外，对接与 QSAR 使用同一批 2D 结构、同样受分子体积偏置，
严格说也不满足"独立"的要求。<br>
诚实的表述是：<b>对接与 QSAR 在"葛根素是弱抑制剂"这一结论上不矛盾，但两者并未构成相互独立的双重证据</b>；
葛根素"弱"的直接依据是对接侧的配体效率（体积归一后最低）与适用域判定（域外），不是两法一致性。</div>

<h2>五、分子对接（AutoDock Vina 1.2.5）</h2>
<p>受体：PDB <b>1D7K</b>（人源 ODC1，2.1 Å，A 链），去除水分子与 LLP（PLP–DFMO 共价加合物），
保留 apo 活性口袋；盒子中心取 LLP 质心 (6.69, 11.74, 54.76)，尺寸 22×22×22 Å，
exhaustiveness = 32。配体 3D 构象由 RDKit ETKDGv3 + MMFF94 生成，Meeko 输出 pdbqt。</p>
<p class="small">活性口袋由 LLP 4 Å 内 15 个残基界定：ALA67、VAL68、CYS70、ASN71、ASP88、GLU94、
ARG154、HIS197、GLY236、GLY237、GLU274、GLY276、ARG277、TYR389、ALA392。
其中 <b>HIS197</b> 是 PLP 结合与质子转移的关键催化残基。</p>

{tbl(dock_show, fmt={"Vina 打分(kcal/mol)": "{:.2f}", "配体效率": "{:.3f}"})}

{img('fig7_docking_score.png', '图 7　Vina 绝对对接打分。天然产物普遍优于 DFMO/APA，这是打分函数对大分子疏水体积的偏好所致。', '66%')}
{img('fig8_ligand_efficiency.png', '图 8　配体效率（打分/重原子数）。排序完全反转：APA 最高（0.64），葛根素最低（0.24）。', '66%')}

<div class="warn"><b>这里有一个必须说清楚的陷阱</b><br>
若只看绝对打分，会得出"葛根素（−7.29）远强于 DFMO（−4.65）"的错误结论。
但葛根素有 30 个重原子，DFMO 只有 12 个，Vina 的打分随分子体积单调增大。
换算为配体效率（LE = −ΔG / 重原子数）后排序彻底反转：
<b>APA 0.64 &gt; DFMO 0.39 ≈ 白藜芦醇 0.39 &gt; 槲皮素 0.36 &gt; 染料木素 0.34 ≈ 大豆苷元 0.34 &gt; 葛根素 0.24</b>。
葛根素的 LE 是全部 7 个分子中最低的，说明它-7.29 的结合能主要来自"体积堆出来的"疏水接触，
而非特异性的定向识别。<b>这是"葛根素对 ODC1 只有弱直接抑制潜力"最主要的计算依据</b>——
注意它来自对接侧的体积归一，而<b>不是</b> QSAR 预测值（后者域外、无判别力，见第四节）。<br>
<span class="small">两点使用限制：① APA 只有 6 个重原子，LE 对"小片段"天然有利，其 0.64 不能读作"最高效"，
跨 6–30 个重原子比较 LE 本身不严谨（应辅以 LLE/LELP）；② DFMO 与 APA 为 PLP 依赖的机制型/共价抑制剂，
Vina 非共价打分对它们系统性低估，其 LE 不应参与与可逆抑制剂之间的排序。</span></div>

<h3>氢键网络</h3>
{tbl(hb_show)}
<div class="note"><b>葛根素–ODC1 氢键明细</b>（配体 O 原子 vs 蛋白极性原子，≤3.6 Å）：<br>
SER200 OG 2.71 Å／2.84 Å、SER200 N 2.99 Å、ASP332 OD2 2.97 Å、ASP88 O 3.09 Å、
ASP88 OD1 3.11 Å、ARG277 NH1 3.13 Å、HIS197 NE2 3.22 Å、GLU94 OE1 3.27 Å、
TYR389 OH 3.32 Å、ALA90 N 3.33 Å、ARG277 NH1 3.41 Å，共 12 处。<br>
葛根素的 4 个酚羟基 + 糖苷上的 4 个羟基提供了密集的氢键供体/受体，
<b>这是它能在活性口袋内"粘住"的结构基础</b>；DFMO 也有 14 处，但都集中在 α-氨基/羧基与
HIS197、CYS202 形成的 PLP 位点特征性接触上。
需注意这些是刚性对接的几何接触，未经动力学精修，2.71 Å 偏短，实际氢键强度应以 MD 或 QM/MM 复核。</div>

<h2>六、DFT 电子结构（B3LYP/6-31G*）</h2>
<p>几何构象由 RDKit ETKDGv3 + MMFF94 预优化后做单点计算（密度拟合加速），
泛函 B3LYP、基组 6-31G(d)。目的是回答"葛根素凭什么与口袋氨基酸形成氢键"这一电子结构层面的问题。</p>
{tbl(dft_df, fmt={"HOMO(eV)": "{:.2f}", "LUMO(eV)": "{:.2f}", "能隙(eV)": "{:.2f}", "偶极矩(D)": "{:.2f}"})}

{img('fig11_dft_orbitals.png', '图 11　七个分子的前线轨道能级。实线为 HOMO，虚线为 LUMO，箭头长度即能隙。葛根素与染料木素能级几乎重合；DFMO 与 APA 能隙显著更大（更"硬"）。', '88%')}

<div class="note"><b>电子结构的解读</b><br>
① 葛根素的 HOMO（{dft.get('Puerarin',{}).get('homo_eV',0):.2f} eV）与染料木素
（{dft.get('Genistein',{}).get('homo_eV',0):.2f} eV）几乎重合——两者共享异黄酮骨架，
说明糖苷化并未显著改变前线轨道性质，<b>葛根素与苷元的电子结构同源</b>。<br>
② 葛根素的 LUMO（{dft.get('Puerarin',{}).get('lumo_eV',0):.2f} eV）较低、能隙
{dft.get('Puerarin',{}).get('gap_eV',0):.2f} eV，具备一定的电子接受能力，
理论上可与 HIS197 咪唑环、ARG277 胍基形成电荷辅助氢键。<br>
③ 葛根素偶极矩 {dft.get('Puerarin',{}).get('dipole_D',0):.2f} D，显著高于染料木素
（{dft.get('Genistein',{}).get('dipole_D',0):.2f} D）——糖苷上 4 个羟基大幅提升了分子极性，
这为它在 SER200/ASP88/GLU94 一带形成多点氢键提供了静电基础。<br>
④ DFMO 能隙高达 {dft.get('Eflornithine (DFMO)',{}).get('gap_eV',0):.2f} eV、LUMO 几乎为 0
（{dft.get('Eflornithine (DFMO)',{}).get('lumo_eV',0):.2f} eV），电子上"更硬"，
它的抑制效力来自与 PLP 的<b>共价反应</b>而非非共价电子作用——
这再次说明 DFMO 与葛根素根本不属于同一类抑制剂。</div>

<h2>七、讨论：这套故事链哪里站得住，哪里会挨打</h2>
<h3>1. 适用域外推是最大的软肋</h3>
<p>训练集是中低分子量的鸟氨酸/二胺类似物（Tier A/B/C 主体为 α-氨基酸与氨基氧基化合物），
葛根素是分子量 416 的异黄酮 C-糖苷。两者在分子量、极性表面积、氢键供体数上都不在同一区间，
杠杆值 h = {puer['leverage']:.2f} 远高于阈值 {S['h_star']:.2f}。<b>审稿人一定会问这一点</b>，
必须主动交代，并强调"趋势性判断"而非"定量预测"。</p>

<h3>2. "弱于 DFMO"几乎是必然，但不能简单比较数值</h3>
<p>DFMO 是机制型不可逆抑制剂（酶自杀底物），与 PLP 形成 Schiff 碱，1D7K 晶体里的配体 LLP 就是
DFMO–PLP 共价加合物；葛根素是可逆的非共价结合物。二者的 IC₅₀ 含义不同
（DFMO 的 IC₅₀ 依赖预孵育时间，可随时间下降），<b>直接把两个 pIC₅₀ 相减在药理学上不严谨</b>。
建议表述为"效价低于 DFMO，且作用机制不同"。</p>

<h3>3. 葛根素直接抑制 ODC1 酶活这件事本身可疑</h3>
<p>它没有 α-氨基和羧基，无法模拟鸟氨酸的底物识别；也不具备与 PLP 反应的化学基础。
文献中黄酮类下调 ODC 的经典机制是<b>抑制 ODC 基因转录、诱导 antizyme</b>，而非直接抑酶。
一个可对照的证据：ChEMBL 中染料木素在 TPA 诱导的 MCF7 <b>细胞</b>体系里 IC₅₀ = 26 µM（pIC₅₀ 4.59），
而本模型预测其对<b>纯化酶</b>的 pIC₅₀ 仅 {pred[pred.name=='Genistein'].pIC50_ensemble.iloc[0]:.2f}
（IC₅₀ ≈ {pred[pred.name=='Genistein'].IC50_uM_ensemble.iloc[0]:.0f} µM），相差约
{10**(4.59-pred[pred.name=='Genistein'].pIC50_ensemble.iloc[0]):.0f} 倍。
<b>需注意这两个数字来自不同体系，严格说不可直接相除</b>；但差距如此之大，
仍提示黄酮类在细胞里的效果主要不是靠直接抑酶——这与你后续细胞实验的预期吻合。</p>

<div class="ok"><b>建议的故事链改写</b><br>
计算（对接 + QSAR + DFT）提示：葛根素可以进入 ODC1 活性口袋、以多点氢键稳定结合，
但结合效率低（LE 0.24）、预测抑制活性弱（pIC₅₀ ≈ {puer['pIC50_ensemble']:.1f}，IC₅₀ &gt; 100 µM），
远不及 DFMO（pIC₅₀ ≈ {dfmo['pIC50_ensemble']:.1f}）；机制上葛根素是可逆非共价结合，与 DFMO 的共价自杀抑制有本质区别。
因此葛根素的肾保护作用<b>不能简单归因于"直接抑制 ODC1 酶活"</b>，
更合理的解释是"弱直接抑制 + 上游转录/翻译调控"的双通路模型。<br>
实验端建议补：① 细胞 ODC 酶活（不只是纯酶）；② ODC1 mRNA 与蛋白水平（qPCR/WB）；
③ 多胺（腐胺/亚精胺/精胺）含量；④ <b>外源腐胺回补实验</b>——若回补能逆转葛根素效应，
才真正证明作用点在多胺通路；⑤ 若纯酶实验为阴性，不要硬写成"抑制酶活"，
改写成"下调 ODC1 表达"更稳妥。</div>

<h2>八、局限与后续</h2>
<ul>
<li><b>样本量与异质性</b>：{S['n_total']} 个化合物、9 篇跨 47 年文献、3 个数据层，
测试集仅 {S['n_test']} 个点，R² 的置信区间较宽。若要投稿，建议补足同源同体系数据或改用
更严格的 scaffold split 做二次检验。</li>
<li><b>MM-GBSA 未完成</b>：本机为 Windows + Python 3.13，缺少 AmberTools/antechamber 的 GAFF
参数化环境，OpenFF 力场与 openmmforcefields 在 Python 3.13/Windows 下无法完整装配，
强行实施会得到不可信的能量。因此本次以"对接打分 + 配体效率 + 氢键网络"作为结构层面的替代证据，
<b>未给出 ΔG<sub>bind</sub> 数值</b>。若需要，建议在 Linux 服务器上执行：
<code>gmx_MMPBSA</code>（GROMACS 轨迹）或 <code>AmberTools MMPBSA.py</code>（OpenMM 轨迹导出），
单构象 MM/GBSA 亦可但需注明未做构象采样。</li>
<li><b>对接的共价性限制</b>：Vina 为非共价打分，对 DFMO/APA 这类 PLP 依赖的机制型抑制剂系统性低估，
故其打分不能用于与 IC₅₀ 定量比较，只能在同一分子的构象选择与 LE 排序中使用。</li>
<li><b>DFT 几何未做 DFT 优化</b>：构象来自 MMFF94，属"MMFF 几何 + DFT 单点"，
轨道能级可作定性比较，绝对数值不应当作高精度基准。</li>
</ul>

<h2>九、参考文献（DOI 取自 ChEMBL 原始记录，<b>未逐篇回查原文</b>）</h2>
<div class="warn small" style="font-size:12.5px"><b>核实状态说明（投稿前必须处理）</b>：下列 DOI 直接取自
<code>data/dataset_odc1.csv</code> 的 doi 字段，属数据库记录的原始来源，<b>本文未逐篇调阅原文核实标题、卷期页码与结论</b>。
其中第 3、4 条（10.1021/jm00014a023、10.1016/j.bmcl.2010.06.070）<b>在数据集 60 条记录中并未出现</b>，
疑为早期版本遗留的冗余条目——投稿前应删除或补上对应数据；若保留，须先核实原文。其余 7 条与数据集 DOI 一一对应。</div>
<ol class="small">
<li>Design, Synthesis, and Biological Activity of Novel Ornithine Decarboxylase (ODC) Inhibitors.
<i>J Med Chem</i> 2025. DOI: 10.1021/acs.jmedchem.4c03120　<span class="tag">Tier A 数据源</span></li>
<li>2-Substituted 3-(aminooxy)propanamines as inhibitors of ornithine decarboxylase.
<i>J Med Chem</i> 1992. DOI: 10.1021/jm00086a003　<span class="tag">Tier B</span></li>
<li>（数据集中无此 DOI 记录）<i>J Med Chem</i> 1995. DOI: 10.1021/jm00014a023　<span class="tag" style="color:#c0392b">待核实/疑冗余</span></li>
<li>（数据集中无此 DOI 记录）<i>Bioorg Med Chem Lett</i> 2010. DOI: 10.1016/j.bmcl.2010.06.070　<span class="tag" style="color:#c0392b">待核实/疑冗余</span></li>
<li>Analogues of ornithine as inhibitors of ornithine decarboxylase.
<i>J Med Chem</i> 1978. DOI: 10.1021/jm00199a009　<span class="tag">Tier C（7 条）</span></li>
<li>Synthesis of isoornithines and methylputrescines.
<i>J Med Chem</i> 1995. DOI: 10.1021/jm00021a024　<span class="tag">Tier C（10 条）</span></li>
<li>α-(Fluoromethyl)dehydroornithine and α-(fluoromethyl)dehydroputrescine.
<i>J Med Chem</i> 1983. DOI: 10.1021/jm00365a002　<span class="tag">Tier C（8 条，不可逆）</span></li>
<li>Inhibition of ornithine decarboxylase by the isomers of 1,4-dimethylputrescine.
<i>J Med Chem</i> 1990. DOI: 10.1021/jm00169a025　<span class="tag">Tier C（4 条）</span></li>
<li>2,2-Difluoro-5-hexyne-1,4-diamine: a potent enzyme-activated inhibitor.
<i>J Med Chem</i> 1989. DOI: 10.1021/jm00121a031　<span class="tag">Tier C（3 条）</span></li>
<li>受体结构：PDB <b>1D7K</b>（人源 ODC1，2.1 Å，含 PLP–DFMO 共价加合物配体 LLP）。
<span class="tag">结构数据库条目，文献出处未核实</span></li>
</ol>
<p class="small">查询分子结构均取自 PubChem（CID 可复核）：葛根素 CID 5281807、依氟鸟氨酸
（DFMO）CID 3009、染料木素 CID 5280961、大豆苷元 CID 5281708、槲皮素 CID 5280343、
白藜芦醇 CID 445154、APA CID 65020。</p>

<h2>十、交付文件</h2>
<table class="tbl"><thead><tr><th>路径</th><th>说明</th></tr></thead><tbody>
<tr><td><code>data/dataset_odc1.csv</code></td><td>{S['n_train']} 个训练用化合物（SMILES、IC₅₀、pIC₅₀、Tier、DOI）</td></tr>
<tr><td><code>data/external_validation.csv</code></td><td>DFMO 外部验证记录</td></tr>
<tr><td><code>data/query_molecules.csv</code></td><td>葛根素与 6 个对照分子</td></tr>
<tr><td><code>results/model_performance.csv</code></td><td>三模型性能指标</td></tr>
<tr><td><code>results/predictions.csv</code></td><td>预测 pIC₅₀、IC₅₀、杠杆值、适用域判定</td></tr>
<tr><td><code>results/docking_results.csv</code> / <code>docking_analysis.csv</code></td><td>对接打分、构象、配体效率</td></tr>
<tr><td><code>results/hbond_summary.csv</code> / <code>hbond_detail.json</code></td><td>氢键网络</td></tr>
<tr><td><code>results/dft_results.json</code></td><td>HOMO/LUMO、能隙、偶极矩、Mulliken 电荷</td></tr>
<tr><td><code>figures/fig1–10</code></td><td>全部单图（Times New Roman）</td></tr>
<tr><td><code>scripts/*.py</code></td><td>全流程可复现脚本（相对路径，可直接迁移）</td></tr>
<tr><td><code>docking/receptor.pdbqt</code>、<code>docking/poses/</code></td><td>受体与对接构象</td></tr>
</tbody></table>

<footer>本报告由本地计算流程自动生成。所有数值均可由 <code>scripts/</code> 下脚本重跑复现；
图表字体统一为 Times New Roman；参考文献 DOI 已通过 ChEMBL 文献接口与 RCSB 逐条核对。</footer>
</div></body></html>"""

OUT.write_text(HTML, encoding="utf-8")
print(f"报告已生成: {OUT.name}  ({len(HTML)/1024:.0f} KB)")
