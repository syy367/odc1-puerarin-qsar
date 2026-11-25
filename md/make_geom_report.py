"""把平衡段几何证据整理成单文件 HTML 报告（图片 base64 内嵌）。"""
from pathlib import Path
import base64
import io

import pandas as pd
from PIL import Image

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "md" / "results"
TARGET = ROOT / "几何证据分析_20261117_wbuddy.html"


def embed(png: Path, width: int = 1200) -> str:
    im = Image.open(png).convert("RGB")
    if im.width > width:
        im = im.resize((width, round(im.height * width / im.width)), Image.LANCZOS)
    buf = io.BytesIO()
    im.save(buf, "JPEG", quality=90, optimize=True)
    b64 = base64.b64encode(buf.getvalue()).decode()
    return f"data:image/jpeg;base64,{b64}"


def main():
    df = pd.read_csv(OUT / "geom_contact_occupancy.csv")
    hb = pd.read_csv(OUT / "geom_hbonds.csv")

    rows = "\n".join(
        f"<tr><td>{r.residue}</td><td class='num'>{r.contact_occupancy_pct:.1f}</td>"
        f"<td class='num'>{r.mean_min_dist_A:.2f}</td><td class='num'>{r.n_frames}</td></tr>"
        for r in df.itertuples())

    hrows = "\n".join(
        f"<tr><td>{r.residue}</td><td>{r.protein_atom}</td><td>{r.ligand_atom}</td>"
        f"<td class='num'>{r.hbond_occupancy_pct:.1f}</td><td class='num'>{r.mean_DA_A:.2f}</td>"
        f"<td class='num'>{r.mean_angle_deg:.0f}</td></tr>"
        for r in hb.itertuples())

    img_occ = embed(OUT / "fig_contact_occupancy.png")
    img_ts = embed(OUT / "fig_contact_timeseries.png")

    html = f"""<!DOCTYPE html>
<html lang="zh-CN"><head><meta charset="utf-8">
<title>葛根素–ODC1 平衡段几何证据（2026-11-17）</title>
<style>
  :root{{--ink:#1f2328;--muted:#57606a;--line:#d8dee4;--soft:#f6f8fa;
        --bad:#a40e26;--badbg:#ffebe9;--warn:#9a6700;--warnbg:#fff8c5;--ok:#1a7f37;--okbg:#dafbe1;}}
  *{{box-sizing:border-box}}
  body{{margin:0;padding:40px 24px 64px;background:#fff;color:var(--ink);
       font-family:"Microsoft YaHei","PingFang SC","Segoe UI",system-ui,sans-serif;
       font-size:15px;line-height:1.75;}}
  .wrap{{max-width:960px;margin:0 auto}}
  h1{{font-size:23px;margin:0 0 6px}}
  .sub{{color:var(--muted);font-size:13px;margin-bottom:26px}}
  h2{{font-size:18px;margin:32px 0 12px;padding-bottom:6px;border-bottom:1px solid var(--line)}}
  h3{{font-size:15px;margin:18px 0 8px}}
  table{{border-collapse:collapse;width:100%;margin:12px 0;font-size:13.5px}}
  th,td{{border:1px solid var(--line);padding:6px 10px;text-align:left;vertical-align:top}}
  th{{background:var(--soft);font-weight:600}}
  td.num,th.num{{text-align:right;font-variant-numeric:tabular-nums}}
  code{{background:var(--soft);padding:1px 5px;border-radius:4px;font-size:12.5px;
       font-family:Consolas,"Courier New",monospace}}
  .box{{border-left:3px solid var(--line);background:var(--soft);padding:12px 16px;margin:14px 0;border-radius:0 6px 6px 0}}
  .box.bad{{border-left-color:var(--bad);background:var(--badbg)}}
  .box.warn{{border-left-color:var(--warn);background:var(--warnbg)}}
  .box.ok{{border-left-color:var(--ok);background:var(--okbg)}}
  ul{{margin:8px 0 8px 20px;padding:0}} li{{margin:5px 0}}
  .fig{{margin:20px 0 8px;text-align:center}}
  .fig img{{max-width:100%;border:1px solid var(--line);border-radius:4px}}
  .cap{{font-size:12.5px;color:var(--muted);margin-top:6px;text-align:left}}
  .foot{{margin-top:34px;padding-top:14px;border-top:1px solid var(--line);color:var(--muted);font-size:12.5px}}
</style></head><body><div class="wrap">

<h1>葛根素–ODC1 复合物平衡段几何证据</h1>
<div class="sub">60–100 ns，每 1 ns 抽取 1 帧共 41 帧 · 轨迹 <code>movie_noWater.pdb</code> · 配体葛根素 C<sub>21</sub>H<sub>20</sub>O<sub>9</sub></div>

<div class="box ok">
<b>结论：配体在 60–100 ns 内稳定占据一个高度疏水的表面口袋，17 个残基接触占有率全部 ≥ 50%（其中 11 个 100%），
全过程中只检出一条占有率 22% 的瞬时氢键。</b><br>
即：结合靠范德华与疏水堆积驱动，没有稳定的氢键网络。这与"弱直接抑制"（IC<sub>50</sub> ≈ 243 µM）以及多酚类易非特异结合的告诫完全自洽，
但与对接预测的极性氢键模型（Ser200 / Asp88 / Arg277 / His197 / Glu94 / Tyr389）在构象上不一致。
</div>

<h2>1. 方法与判据</h2>
<table>
<tbody>
<tr><th style="width:26%">采样</th><td>60–100 ns 平衡段，每 1 ns 一帧，共 41 帧（<code>md/equil_frames/f601–f1001.pdb</code>）</td></tr>
<tr><th>接触判据</th><td>残基任一原子到配体任一重原子的最小距离 &lt; 4.0 Å；占有率 = 满足帧数 / 41</td></tr>
<tr><th>氢键判据</th><td>D···A &lt; 3.5 Å，H···A &lt; 2.5 Å，D–H···A 角 &gt; 120°（双向扫描：配体作供体、配体作受体）</td></tr>
<tr><th>配体原子指认</th><td>轨迹中配体原子名无区分性（均为元素符号），按几何就近原则确定每个氢的母重原子</td></tr>
</tbody>
</table>

<h2>2. 接触占有率</h2>
<table>
<thead><tr><th>残基</th><th class="num">占有率 / %</th><th class="num">平均最小距离 / Å</th><th class="num">命中帧数</th></tr></thead>
<tbody>
{rows}
</tbody>
</table>

<div class="fig">
<img src="{img_occ}" alt="contact occupancy">
<div class="cap"><b>图 1</b> 各残基在 60–100 ns 内的接触占有率（4 Å 判据）。颜色区分侧链性质：蓝 = 疏水，绿 = 极性，红 = 带电。
虚线为 50% 与 70% 参考线。</div>
</div>

<div class="fig">
<img src="{img_ts}" alt="contact time series">
<div class="cap"><b>图 2</b> 六个代表性残基到配体的最小距离随时间变化。水平虚线为 4 Å 接触判据。
所有曲线全程在 1.7–4.8 Å 之间波动，无解离事件。</div>
</div>

<h2>3. 氢键</h2>
<table>
<thead><tr><th>残基</th><th>蛋白原子</th><th>配体原子</th><th class="num">占有率 / %</th><th class="num">D···A / Å</th><th class="num">角度 / °</th></tr></thead>
<tbody>
{hrows}
</tbody>
</table>
<p>仅 Gly236 主链羰基氧与配体羟基形成氢键，占有率 22%，且 D···A 距离偏长（3.22 Å）。这是典型的<b>瞬时、非特异性</b>接触。
与 <code>hbond.xvg</code> 实测的 0.18 ± 0.40 个氢键一致。</p>

<h2>4. 化学解读</h2>
<ul>
<li><b>疏水占绝对主导</b>：17 个接触残基中 14 个为疏水（Val198、Pro205、Phe208、Ile234、Gly236、Gly237、Phe238、Pro239、Val244、Leu246、Ile251、Val254、Ile255、Cys202），
仅 1 个极性（Thr203）、2 个带电（Lys245、Asp204）。</li>
<li><b>极性残基只是"贴边"</b>：Lys245 与 Asp204 虽然带电，但占有率与平均距离都偏低（97.6% / 2.80 Å、58.5% / 3.52 Å），且没有形成氢键——它们的极性头朝溶剂。</li>
<li><b>这解释了热力学结果</b>：ΔE 非键 −21.8 kcal/mol 是范德华堆积的典型量级；6 个羟基被埋进疏水环境却没有氢键受体，
GB 去溶剂化惩罚因此高达 +22.2 kcal/mol，净结合自由能接近零。</li>
<li><b>与酶学数据自洽</b>：IC<sub>50</sub> ≈ 243 µM 对应 ΔG ≈ −4.9 kcal/mol，本就属弱结合；多酚的非特异性疏水吸附是已知现象。</li>
</ul>

<div class="box bad">
<b>与对接结果的根本冲突（必须在投稿前解决）</b>
<ul style="margin-bottom:0">
<li>对接（Vina，本项目 <code>results/hbond_detail.json</code>）给出的葛根素氢键残基是
Tyr389 / Asp88 / Ala90 / Glu94 / His197 / Ser200 / Arg277 / Asp332；</li>
<li>而这 8 个残基在 MD 轨迹中距配体 <b>5.3–21.6 Å</b>，全部远离界面；</li>
<li>MD 的实际口袋与本方案的 17 个残基完全落在另一处（Phe196–Val198–Thr203–Asp204…—Ile255 一带），
两者相距 7–20 Å，是两个不同位点。</li>
</ul>
也就是说：<b>对接与 MD 描述的不是同一个结合模式</b>。论文中不能同时用两套残基，须先查明 MD 起始结构的配体来源。
</div>

<h2>5. 对丙氨酸扫描的直接影响</h2>
<p>按本方案的三级判据框架，候选残基应当从这 17 个接触残基中挑选，而不是从别处引用。具体分级建议：</p>
<table>
<thead><tr><th style="width:24%">层级</th><th>残基</th><th style="width:26%">依据</th></tr></thead>
<tbody>
<tr><td>首选（100% 占有 + 距离 &lt; 3.1 Å）</td><td>Thr203、Pro239、Ile251、Pro205、Leu246、Phe208、Phe238、Cys202、Ile255、Val198</td><td>全程稳定接触，突变最可能显著改变 ΔΔG</td></tr>
<tr><td>次选（占有 85–100%）</td><td>Gly236、Lys245、Gly237、Val244</td><td>含唯一氢键来源 Gly236 主链</td></tr>
<tr><td>边缘（占 50–76%）</td><td>Val254、Asp204、Ile234</td><td>接触时断时续，适合作为灰区样本</td></tr>
<tr><td>负对照</td><td>任一距配体 &gt; 10 Å 的表面残基（如 Glu94、Ala90）</td><td>与本方案要求一致</td></tr>
</tbody>
</table>
<p>用户原方案列的 His335 / Asn385 / Ser241 / Cys202 中，只有 <b>Cys202</b> 落在界面内（100% / 2.98 Å）；
His335（14.5–20.8 Å）、Asn385（16.3–19.9 Å）、Ser241（7.4–9.2 Å）完全够不到配体，作为扫描对象无意义。</p>

<h2>6. 局限</h2>
<ul>
<li>接触与氢键均为<b>几何判据</b>，不含能量项，不能替代 ΔΔG。</li>
<li>轨迹已去水，<b>无法做水桥分析</b>——而疏水口袋中的水分子恰恰可能是极性羟基的替代伴侣，这是本分析最大的信息缺口。</li>
<li>轨迹为<b>单体 390 残基、无 PLP</b>，结论只能限定为"单体模型下的表面结合"。</li>
<li>未计算埋藏面积（SASA 变化），可后续补充以支持"疏水埋藏"的判断。</li>
</ul>

<div class="foot">
数据文件：<code>md/results/geom_contact_occupancy.csv</code>、<code>md/results/geom_hbonds.csv</code>、
<code>md/results/fig_contact_occupancy.png</code>、<code>md/results/fig_contact_timeseries.png</code>。
分析脚本：<code>md/geom_analysis.py</code>（判据与统计）、<code>md/fig_geom.py</code>（出图，Times New Roman）。
</div>

</div></body></html>"""

    TARGET.write_text(html, encoding="utf-8")
    print(f"写出 {TARGET.name}  ({TARGET.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
