# ODC1 Inhibitor Virtual Screening — Puerarin

本仓库收录 **ODC1（鸟氨酸脱羧酶 1）抑制剂虚拟筛选**的完整可复现流程，围绕一个核心问题：*葛根素（Puerarin）能否直接结合 ODC1？*

项目整合了四条相互独立的证据链——**2D-QSAR 建模、分子对接、DFT 电子结构计算、分子动力学几何分析**——共同回答上述问题。最终结论为：葛根素与 ODC1 存在**弱而松散的直接结合**（以疏水作用为主、缺乏特异性氢键识别），因此论文口径为"弱直接抑制 + 上游转录调控"双通路，而非强效直接抑制。

> 数据集、脚本与结论均为本项目原始产出；只提交代码、数据、结果与图，不提交虚拟环境、外部工具、汇报文档与备份文件。

[English version](README.md)

---

## 🚀 Features（功能）

- **2D-QSAR 建模**：从 RDKit 2D 描述符出发，经描述符筛选与训练集划分，构建 PLS / Random Forest / XGBoost 回归模型；含 Y-randomization、外部验证（DFMO）、applicability domain（William 图）与敏感性分析。
- **分子对接**：AutoDock Vina 1.2.5，受体为 ODC1（PDB 1D7K, A 链）活性位点；配体 3D 构象由 RDKit ETKDG + MMFF94 生成，pdbqt 由 Meeko 写出；计算配体效率（LE）以校正打分随分子体积增大的偏差。
- **DFT 电子结构**：PySCF，B3LYP/6-31G* 单点能，输出 HOMO/LUMO 能隙、偶极矩、Mulliken 电荷，用于解释结合的电子结构基础。
- **分子动力学几何证据**：100 ns 轨迹抽帧，统计接触残基占有率、氢键、RMSD 与 SASA，判定结合稳定性与作用力类型。
- **投稿级图表**：全部图片 Times New Roman 字体，含正文图与补充图（journal 版）。

---

## 📂 Project Structure（目录结构）

```
├── scripts/                    # 数据处理与建模脚本
│   ├── build_dataset.py        #   从 raw 数据构建分层 QSAR 数据集（A/B/C 三层 + DFMO 外部验证）
│   ├── qsar.py                 #   2D-QSAR 全流程（描述符 → 模型 → 验证 → 预测）
│   ├── dock.py                 #   AutoDock Vina 对接
│   ├── analyze_docking.py      #   对接结果分析（配体效率、残基接触）
│   ├── dft.py                  #   PySCF 单点能计算
│   ├── figures_*.py            #   各步骤出图脚本
│   ├── export_tables.py        #   投稿表格导出
│   ├── hbond.py / make_panel.py / make_report.py / make_word.py / fetch_wheels.py
├── md/                         # 分子动力学分析脚本与结果
│   ├── prep_ligand.py          #   配体准备
│   ├── run_frames.py           #   轨迹抽帧
│   ├── geom_analysis.py        #   几何证据：接触占有率、氢键
│   ├── mmgbsa.py / mmgbsa_obc2.py / diag_energy.py / fig_geom.py / make_geom_report.py
│   ├── frames/                 #   抽帧结构（pdb / sdf）
│   ├── forcefields/            #   OpenFF 力场 xml
│   └── results/                #   几何与能量结果（csv / png）
├── data/                       # 数据
│   ├── dataset_odc1.csv        #   建模数据集（61 个化合物，不含 DFMO）
│   ├── external_validation.csv #   DFMO 外部验证集
│   ├── query_molecules.csv     #   待预测分子（葛根素及对照）
│   ├── protein/1D7K.pdb        #   受体结构
│   └── raw/                    #   ChEMBL / BindingDB 原始 json
├── docking/                    # 对接输入与输出
│   ├── receptor*.pdbqt         #   受体
│   └── poses/                  #   各分子 pdbqt 与 log
├── results/                    # 各步骤结果（csv / json）
├── figures/                    # 图片（journal/ 为投稿尺寸版）
├── requirements.txt            # Python 依赖版本锁定
├── README.md                   # 英文说明
└── README_cn.md                # 中文说明（本文件）
```

> 注：`tools/`（vina.exe、xtb.exe、离线 wheel）与 `.venv/` 未纳入仓库，需按下方说明自行准备。

---

## 🛠️ Installation（安装）

### 1. 克隆仓库

```bash
git clone <your-repo-url>.git
cd <repo-name>
```

### 2. 创建并激活环境

```bash
conda create -n odc1 python=3.13 -y
conda activate odc1
```

### 3. 安装 Python 依赖

```bash
pip install -r requirements.txt
```

### 4. 准备外部工具（未随仓库提交）

| 工具 | 用途 | 放置位置 |
|---|---|---|
| AutoDock Vina 1.2.5 | 分子对接 | `tools/bin/vina.exe` |
| xtb 6.6.1 | 构象/半经验计算（可选） | `tools/xtb/xtb-6.6.1/bin/xtb.exe` |

Windows 下的离线 wheel 可从 PyPI 或对应项目主页下载后本地安装。所有脚本内的路径均基于脚本自身目录推导，**无需修改绝对路径即可在任意机器运行**。

---

## ▶️ Usage（使用）

按顺序运行即可复现完整流程（在项目根目录执行）：

```bash
# 1) 由 raw 数据构建数据集
python scripts/build_dataset.py

# 2) 2D-QSAR：建模 + 验证 + 对查询分子预测
python scripts/qsar.py

# 3) 分子对接（需先放好 vina.exe）
python scripts/dock.py
python scripts/analyze_docking.py

# 4) DFT 单点能
python scripts/dft.py

# 5) 出图与投稿表格
python scripts/figures_dock.py
python scripts/figures_dft.py
python scripts/export_tables.py

# 6) 分子动力学几何证据（需有轨迹抽帧结构）
python md/geom_analysis.py
```

核心结果见 `results/`（`model_performance.csv`、`predictions.csv`、`docking_results.csv`、`qsar_summary.json` 等）。

---

## 📊 Example Input Format（输入格式）

查询分子文件 `data/query_molecules.csv`（逗号分隔）：

```csv
name,pubchem_cid,smiles,inchikey,role
Puerarin,5281807,O=c1c(-c2ccc(O)cc2)coc2c([C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)c(O)ccc12,HKEAFJYKMMKDOR-VPRICQMDSA-N,目标分子
Eflornithine (DFMO),3009,NCCCC(N)(C(=O)O)C(F)F,VLCYCQAOQCDTCN-UHFFFAOYSA-N,阳性对照/金标准
Genistein,5280961,O=c1c(-c2ccc(O)cc2)coc2cc(O)cc(O)c12,TZBJGXHYKVUXJN-UHFFFAOYSA-N,天然产物对照（葛根素苷元）
```

---

## 📖 Data Sources & Tools（数据与工具）

- **数据来源**：ChEMBL、BindingDB 及原始文献；数据集分层与来源 DOI 见 `scripts/build_dataset.py` 顶部注释与 `data/dataset_odc1.csv` 的 `doi` 列。
- **化学信息学与建模**：RDKit（https://www.rdkit.org）、scikit-learn（https://scikit-learn.org）、XGBoost（https://xgboost.readthedocs.io）。
- **分子对接**：AutoDock Vina（https://vina.scripps.edu）、Meeko（https://github.com/forlilab/Meeko）。
- **量子化学**：PySCF（https://pyscf.org）。
- **分子动力学**：OpenMM（https://openmm.org）、OpenFF（https://openforcefield.org）。

---

## 👤 Credits（致谢）

原始数据集、脚本与结果由本项目独立整理与计算完成。
