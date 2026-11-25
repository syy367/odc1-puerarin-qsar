# ODC1 Inhibitor Virtual Screening — Puerarin

This repository contains a fully reproducible workflow for the **virtual screening of ODC1 (ornithine decarboxylase 1) inhibitors**, built around a single question: *Can puerarin bind directly to ODC1?*

Four independent lines of evidence — **2D-QSAR modeling, molecular docking, DFT electronic-structure calculation, and molecular-dynamics geometric analysis** — are integrated to address this question. The overall conclusion is that puerarin exhibits **weak, loose direct binding** to ODC1 (dominated by hydrophobic contacts, with no specific hydrogen-bond network). The reported interpretation is therefore a dual-pathway model of "weak direct inhibition plus upstream transcriptional regulation," rather than potent direct inhibition.

> The dataset, scripts, and conclusions are original outputs of this project. Only code, data, results, and figures are committed; virtual environments, external tools, reports, and backup files are excluded.

[Chinese version](README_cn.md)

---

## 🚀 Features

- **2D-QSAR modeling**: RDKit 2D descriptors are filtered and split into training/test sets to build PLS, Random Forest, and XGBoost regression models, with Y-randomization, external validation (DFMO), applicability-domain assessment (William's plot), and sensitivity analysis.
- **Molecular docking**: AutoDock Vina 1.2.5 against the ODC1 active site (PDB 1D7K, chain A). Ligand 3D conformers are generated with RDKit ETKDG + MMFF94 and converted to pdbqt with Meeko. Ligand efficiency (LE) is computed to correct for the size-dependent inflation of docking scores.
- **DFT electronic structure**: PySCF single-point energies at B3LYP/6-31G*, reporting HOMO/LUMO gap, dipole moment, and Mulliken charges to rationalize the electronic basis of binding.
- **MD geometric evidence**: Frame extraction from a 100 ns trajectory to quantify contact-residue occupancy, hydrogen bonds, RMSD, and SASA, determining binding stability and interaction type.
- **Publication-ready figures**: All plots use Times New Roman; both main and supplementary figures (journal versions) are included.

---

## 📂 Project Structure

```
├── scripts/                    # Data processing and modeling scripts
│   ├── build_dataset.py        #   Build the tiered QSAR dataset from raw data (tiers A/B/C + DFMO external set)
│   ├── qsar.py                 #   End-to-end 2D-QSAR (descriptors → model → validation → prediction)
│   ├── dock.py                 #   AutoDock Vina docking
│   ├── analyze_docking.py      #   Docking analysis (ligand efficiency, residue contacts)
│   ├── dft.py                  #   PySCF single-point energies
│   ├── figures_*.py            #   Plotting scripts for each step
│   ├── export_tables.py        #   Export publication tables
│   ├── hbond.py / make_panel.py / make_report.py / make_word.py / fetch_wheels.py
├── md/                         # MD analysis scripts and results
│   ├── prep_ligand.py          #   Ligand preparation
│   ├── run_frames.py           #   Trajectory frame extraction
│   ├── geom_analysis.py        #   Geometric evidence: contact occupancy, hydrogen bonds
│   ├── mmgbsa.py / mmgbsa_obc2.py / diag_energy.py / fig_geom.py / make_geom_report.py
│   ├── frames/                 #   Extracted structures (pdb / sdf)
│   ├── forcefields/            #   OpenFF force-field xml
│   └── results/                #   Geometric and energy results (csv / png)
├── data/                       # Data
│   ├── dataset_odc1.csv        #   Modeling dataset (61 compounds, DFMO excluded)
│   ├── external_validation.csv #   DFMO external validation set
│   ├── query_molecules.csv     #   Query molecules (puerarin and controls)
│   ├── protein/1D7K.pdb        #   Receptor structure
│   └── raw/                    #   Raw ChEMBL / BindingDB json
├── docking/                    # Docking inputs and outputs
│   ├── receptor*.pdbqt         #   Receptor
│   └── poses/                  #   Per-molecule pdbqt and logs
├── results/                    # Step-wise results (csv / json)
├── figures/                    # Figures (journal/ holds publication-size versions)
├── requirements.txt            # Pinned Python dependencies
├── README.md                   # This file (English)
└── README_cn.md                # Chinese version
```

> Note: `tools/` (vina.exe, xtb.exe, offline wheels) and `.venv/` are not included in the repository and must be prepared as described below.

---

## 🛠️ Installation

### 1. Clone the repository

```bash
git clone <your-repo-url>.git
cd <repo-name>
```

### 2. Create and activate the environment

```bash
conda create -n odc1 python=3.13 -y
conda activate odc1
```

### 3. Install Python dependencies

```bash
pip install -r requirements.txt
```

### 4. Prepare external tools (not committed)

| Tool | Purpose | Location |
|---|---|---|
| AutoDock Vina 1.2.5 | Molecular docking | `tools/bin/vina.exe` |
| xtb 6.6.1 | Conformer / semi-empirical calculation (optional) | `tools/xtb/xtb-6.6.1/bin/xtb.exe` |

Offline wheels for Windows can be downloaded from PyPI or the respective project homepages and installed locally. All paths in the scripts are resolved relative to each script's own directory, so **no absolute paths need to be edited to run them on another machine**.

---

## ▶️ Usage

Run the steps below in order from the project root to reproduce the full workflow:

```bash
# 1) Build the dataset from raw data
python scripts/build_dataset.py

# 2) 2D-QSAR: modeling + validation + prediction on query molecules
python scripts/qsar.py

# 3) Molecular docking (requires vina.exe)
python scripts/dock.py
python scripts/analyze_docking.py

# 4) DFT single-point energies
python scripts/dft.py

# 5) Figures and publication tables
python scripts/figures_dock.py
python scripts/figures_dft.py
python scripts/export_tables.py

# 6) MD geometric evidence (requires extracted trajectory frames)
python md/geom_analysis.py
```

Key results are in `results/` (`model_performance.csv`, `predictions.csv`, `docking_results.csv`, `qsar_summary.json`, etc.).

---

## 📊 Example Input Format

Query molecules file `data/query_molecules.csv` (comma-separated):

```csv
name,pubchem_cid,smiles,inchikey,role
Puerarin,5281807,O=c1c(-c2ccc(O)cc2)coc2c([C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)c(O)ccc12,HKEAFJYKMMKDOR-VPRICQMDSA-N,target molecule
Eflornithine (DFMO),3009,NCCCC(N)(C(=O)O)C(F)F,VLCYCQAOQCDTCN-UHFFFAOYSA-N,positive control / gold standard
Genistein,5280961,O=c1c(-c2ccc(O)cc2)coc2cc(O)cc(O)c12,TZBJGXHYKVUXJN-UHFFFAOYSA-N,natural-product control (puerarin aglycone)
```

---

## 📖 Data Sources & Tools

- **Data sources**: ChEMBL, BindingDB, and original literature. Dataset tiers and source DOIs are documented in the header comments of `scripts/build_dataset.py` and the `doi` column of `data/dataset_odc1.csv`.
- **Cheminformatics and modeling**: RDKit (https://www.rdkit.org), scikit-learn (https://scikit-learn.org), XGBoost (https://xgboost.readthedocs.io).
- **Molecular docking**: AutoDock Vina (https://vina.scripps.edu), Meeko (https://github.com/forlilab/Meeko).
- **Quantum chemistry**: PySCF (https://pyscf.org).
- **Molecular dynamics**: OpenMM (https://openmm.org), OpenFF (https://openforcefield.org).

---

## 👤 Credits

The original dataset, scripts, and results were independently compiled and computed in this project.
