"""构建 ODC1 抑制剂 pIC50 数据集。

来源：ChEMBL（已核对 DOI，见 data/sources.md）
分层：
  A = 人源重组 ODC1 纯化酶 IC50（2025 J Med Chem, 10.1021/acs.jmedchem.4c03120）
  B = 大鼠肝纯化 ODC IC50（1992/1995 J Med Chem, 2010 BMCL）
  C = 大鼠肝纯化 ODC Ki（1978-1995 J Med Chem），按 Cheng-Prusoff（[S]=Km）换算 IC50 = 2*Ki

DFMO 被单独拎出作为外部验证化合物，不进训练/测试集。

输出（相对脚本目录）：
  ../data/dataset_odc1.csv      训练/测试用（不含 DFMO）
  ../data/external_validation.csv  DFMO 外部验证记录
  ../data/query_molecules.csv   查询分子（葛根素及对照）
"""
import json
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"

# 文献 DOI（已通过 ChEMBL document API 核对）
DOI = {
    "CHEMBL6103633": "10.1021/acs.jmedchem.4c03120",
    "CHEMBL1125431": "10.1021/jm00086a003",
    "CHEMBL1128469": "10.1021/jm00014a023",
    "CHEMBL1208752": "10.1016/j.bmcl.2010.06.070",
    "CHEMBL3243972": "10.1021/jm00199a009",
    "CHEMBL1128697": "10.1021/jm00021a024",
    "CHEMBL1122386": "10.1021/jm00365a002",
    "CHEMBL1125165": "10.1021/jm00169a025",
    "CHEMBL1124320": "10.1021/jm00121a031",
}

_lfc = rdMolStandardize.LargestFragmentChooser()


def norm(smi):
    """去盐/取最大片段并标准化，返回 (canonical_smiles, inchikey)。"""
    m = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
    if m is None:
        return None, None
    try:
        m = _lfc.choose(m)
    except Exception:
        pass
    return Chem.MolToSmiles(m), Chem.MolToInchiKey(m)


def load_raw(fname, tier, endpoint):
    with open(RAW / fname, encoding="utf-8") as fh:
        acts = json.load(fh)["activities"]
    rows = []
    for a in acts:
        desc = (a.get("assay_description") or "").lower()
        if tier == "A" and "his-6 tagged human odc" not in desc:
            continue
        if tier in ("B", "C") and "ornithine decarboxylase" not in desc:
            continue
        val = a.get("standard_value")
        smi = a.get("canonical_smiles")
        if val in (None, "") or not smi:
            continue
        try:
            val = float(val)
        except (TypeError, ValueError):
            continue
        if val <= 0:
            continue
        rows.append(
            {
                "chembl_id": a.get("molecule_chembl_id"),
                "pref_name": a.get("molecule_pref_name"),
                "smiles_raw": smi,
                "value_nM": val,
                "endpoint": endpoint,
                "tier": tier,
                "species": "Homo sapiens" if tier == "A" else "Rattus norvegicus",
                "assay_desc": a.get("assay_description"),
                "doc": a.get("document_chembl_id"),
                "year": a.get("document_year"),
                "journal": a.get("document_journal"),
            }
        )
    return rows


def main():
    rows = []
    rows += load_raw("human_ic50.json", "A", "IC50")
    rows += load_raw("rat_ic50.json", "B", "IC50")
    rows += load_raw("rat_ki.json", "C", "Ki")
    df = pd.DataFrame(rows)
    print(f"原始记录: {len(df)}  (A={sum(df.tier=='A')}, B={sum(df.tier=='B')}, C={sum(df.tier=='C')})")

    # pIC50：IC50 直接换算；Ki 按 IC50 = 2*Ki 换算
    df["pIC50"] = df.apply(
        lambda r: 9 - __import__("math").log10(r["value_nM"] * (2 if r["endpoint"] == "Ki" else 1)),
        axis=1,
    )

    std = df["smiles_raw"].map(norm)
    df["smiles"] = [s for s, _ in std]
    df["inchikey"] = [k for _, k in std]
    df = df.dropna(subset=["smiles"])
    df["doi"] = df["doc"].map(DOI)

    # 同一结构多条记录：按 tier 优先（A>B>C）再取中位数
    df["tier_rank"] = df["tier"].map({"A": 0, "B": 1, "C": 2})
    agg = (
        df.sort_values(["inchikey", "tier_rank", "value_nM"])
        .groupby("inchikey", as_index=False)
        .agg(
            chembl_id=("chembl_id", "first"),
            pref_name=("pref_name", "first"),
            smiles=("smiles", "first"),
            value_nM=("value_nM", "median"),
            pIC50=("pIC50", "median"),
            endpoint=("endpoint", "first"),
            tier=("tier", "first"),
            species=("species", "first"),
            assay_desc=("assay_desc", "first"),
            doi=("doi", "first"),
            year=("year", "first"),
            journal=("journal", "first"),
            n_records=("value_nM", "size"),
        )
    )
    print(f"去重后化合物: {len(agg)}")

    # DFMO 拎出做外部验证
    dfmo = agg[agg["chembl_id"] == "CHEMBL830"]
    if len(dfmo):
        dfmo.to_csv(ROOT / "data" / "external_validation.csv", index=False)
        print("DFMO 外部验证:", dfmo[["chembl_id", "value_nM", "pIC50", "tier"]].to_dict("records"))
    train = agg[agg["chembl_id"] != "CHEMBL830"].copy()
    train = train.sort_values("pIC50").reset_index(drop=True)
    train.insert(0, "idx", range(1, len(train) + 1))
    train.to_csv(ROOT / "data" / "dataset_odc1.csv", index=False)
    print(f"训练集候选: {len(train)}")
    print("pIC50 范围 %.2f - %.2f, 均值 %.2f, 标准差 %.2f"
          % (train.pIC50.min(), train.pIC50.max(), train.pIC50.mean(), train.pIC50.std()))
    print("tier 分布:", train.tier.value_counts().to_dict())

    # 查询分子（结构来自 PubChem，CID 可核）
    queries = [
        ("Puerarin", "5281807", "C1=CC(=CC=C1C2=COC3=C(C2=O)C=CC(=C3[C@H]4[C@@H]([C@H]([C@@H]([C@H](O4)CO)O)O)O)O)O", "目标分子"),
        ("Eflornithine (DFMO)", "3009", "C(CC(C(F)F)(C(=O)O)N)CN", "阳性对照/金标准"),
        ("Genistein", "5280961", "C1=CC(=CC=C1C2=COC3=CC(=CC(=C3C2=O)O)O)O", "天然产物对照（葛根素苷元）"),
        ("Daidzein", "5281708", "C1=CC(=CC=C1C2=COC3=C(C2=O)C=CC(=C3)O)O", "天然产物对照（异黄酮）"),
        ("Quercetin", "5280343", "C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O)O", "天然产物对照（黄酮醇）"),
        ("Resveratrol", "445154", "C1=CC(=CC=C1/C=C/C2=CC(=CC(=C2)O)O)O", "天然产物对照（芪类）"),
        ("APA (3-aminooxy-1-aminopropane)", "65020", "C(CN)CON", "ODC 强效抑制剂（共晶配体 XAP）"),
    ]
    qrows = []
    for name, cid, smi, role in queries:
        s, k = norm(smi)
        qrows.append({"name": name, "pubchem_cid": cid, "smiles": s, "inchikey": k, "role": role})
    pd.DataFrame(qrows).to_csv(ROOT / "data" / "query_molecules.csv", index=False)
    print(f"\n查询分子 {len(qrows)} 个已写入 data/query_molecules.csv")


if __name__ == "__main__":
    main()
