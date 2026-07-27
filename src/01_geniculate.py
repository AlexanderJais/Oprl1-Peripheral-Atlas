"""Geniculate ganglion: opioid receptor and ligand levels, recomputed from GEO.

Two independent datasets, deliberately different in sensitivity:
  GSE102443  96 neurons, full-length SMART-seq, FPKM. High per-cell sensitivity,
             so this carries the expression-level claims.
  GSE135801  454 Phox2b+ neurons, 3' scRNA-seq counts. Lower sensitivity; used
             to check the receptor ordering holds on a second platform.

Levels are summarised as pseudobulk means, never as cross-platform detection
percentages. Output feeds src/04_synthesis.py.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import atlas_common as ac

FPKM_FILE = ("GSE102443_GEO-ID_Dvoryanchikov_2017_Datatable_FPKM.txt.gz")
ZUKER_XLSX = "gse135801/GSM4037432_GG_scRNAseq_Phox2b_expressed454.xlsx"


def load_dvoryanchikov():
    """96 geniculate neurons, transcript FPKM collapsed to gene level."""
    df = pd.read_csv(ac.DATA / FPKM_FILE, sep="\t", low_memory=False)
    cell_cols = [c for c in df.columns
                 if c not in ("Transcript", "Chromosome", "Begin", "End",
                              "Strand", "Gene")]
    # Cell labels carry plate coordinates; keep the trailing sample ID.
    names = [c.split("]]")[-1].strip() if "]]" in c else c.strip()
                 for c in cell_cols]
    mat = df[cell_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    mat.columns = names
    mat["Gene"] = df["Gene"].astype(str).values
    genes = mat.groupby("Gene").sum()
    print(f"  GSE102443: {genes.shape[1]} cells x {genes.shape[0]:,} genes (FPKM)")
    return genes


def load_zuker():
    """454 Phox2b+ geniculate neurons, raw counts -> CPM."""
    x = pd.read_excel(ac.DATA / ZUKER_XLSX, index_col=0)
    x = x.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    # Orient genes as rows.
    if x.shape[0] < x.shape[1]:
        x = x.T
    lib = x.sum(axis=0)
    lib[lib == 0] = 1.0
    cpm = x.div(lib, axis=1) * 1e6
    print(f"  GSE135801: {cpm.shape[1]} cells x {cpm.shape[0]:,} genes (CPM)")
    return cpm


def summarise(mat, label, unit):
    """Pseudobulk mean and detection for one matrix, plus sanity checks."""
    mean = mat.mean(axis=1)
    det = (mat > 0).mean(axis=1) * 100.0

    missing = [g for g in ac.SANITY_GENES if g not in mat.index]
    print(f"  {label} sanity: " + ", ".join(
        f"{g} {mean[g]:.1f}" for g in ac.SANITY_GENES if g in mat.index)
        + (f"  [absent: {missing}]" if missing else ""))

    rows = []
    for g in ac.OPIOID_GENES:
        rows.append({
            "dataset": label,
            "tissue": "geniculate",
            "assay": ac.ASSAY_LABEL[label],
            "gene": g,
            "unit": unit,
            "mean_level": round(float(mean.get(g, 0.0)), 4),
            "pct_detected": round(float(det.get(g, 0.0)), 2),
            "n_cells": int(mat.shape[1]),
            "in_matrix": bool(g in mat.index),
        })
    return pd.DataFrame(rows), mean


def main() -> int:
    ac.RES.mkdir(parents=True, exist_ok=True)

    dv = load_dvoryanchikov()
    zk = load_zuker()

    dv_tbl, dv_mean = summarise(dv, "GSE102443", "FPKM")
    zk_tbl, zk_mean = summarise(zk, "GSE135801", "CPM")
    tbl = pd.concat([dv_tbl, zk_tbl], ignore_index=True)
    ac.save_table(tbl, "geniculate_opioid_levels.csv")

    print("\n  Opioid gene levels, geniculate:")
    print(tbl.pivot(index="gene", columns="dataset",
                    values="mean_level").reindex(ac.OPIOID_GENES).to_string())

    # Receptor ordering, the statistic that survives the platform difference.
    ranks = []
    for label, mean in (("GSE102443", dv_mean), ("GSE135801", zk_mean)):
        r = ac.receptor_rank(mean)
        r.insert(0, "dataset", label)
        r.insert(0, "tissue", "geniculate")
        ranks.append(r)
    rank_tbl = pd.concat(ranks, ignore_index=True)
    ac.save_table(rank_tbl, "geniculate_receptor_rank.csv")
    print("\n  Receptor rank within each dataset:")
    print(rank_tbl.to_string(index=False))

    # Receptor against its own ligand, in the same cells.
    lr = []
    for label, mean in (("GSE102443", dv_mean), ("GSE135801", zk_mean)):
        lr.append({
            "tissue": "geniculate", "dataset": label,
            "Oprl1": round(float(mean.get("Oprl1", 0.0)), 4),
            "Pnoc": round(float(mean.get("Pnoc", 0.0)), 4),
            "Oprl1_over_Pnoc": ac.ligand_receptor_ratio(mean),
            "Pnoc_in_matrix": bool("Pnoc" in mean.index),
            "Oprl1_transcriptome_percentile":
                round(ac.transcriptome_percentile(mean, "Oprl1"), 1),
        })
    lr_tbl = pd.DataFrame(lr)
    ac.save_table(lr_tbl, "geniculate_ligand_receptor.csv")
    print("\n  Receptor vs ligand in the same cells:")
    print(lr_tbl.to_string(index=False))

    # Per-cell Oprl1 for the deep dataset, split by gustatory/somatosensory as
    # the source project does (Phox2b FPKM > 5).
    if "Phox2b" in dv.index and "Oprl1" in dv.index:
        gust = dv.loc["Phox2b"] > 5
        per_cell = pd.DataFrame({
            "cell": dv.columns,
            "division": np.where(gust, "gustatory (Phox2b+)",
                                 "somatosensory (Phox2b-)"),
            "Oprl1_FPKM": dv.loc["Oprl1"].values,
            "Phox2b_FPKM": dv.loc["Phox2b"].values,
            "Penk_FPKM": dv.loc["Penk"].values if "Penk" in dv.index else np.nan,
        })
        ac.save_table(per_cell, "geniculate_per_cell_GSE102443.csv")
        print("\n  Oprl1 FPKM by division (GSE102443):")
        print(per_cell.groupby("division")["Oprl1_FPKM"]
              .agg(["size", "mean", "median"]).round(2).to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
