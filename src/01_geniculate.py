"""Geniculate ganglion: opioid receptor and ligand levels, recomputed from GEO.

Two independent datasets, deliberately different in sensitivity:
  GSE102443  96 neurons, full-length SMART-seq, FPKM. High per-cell sensitivity,
             so this carries the expression-level claims.
  GSE135801  454 Phox2b+ neurons, 3' scRNA-seq counts. Lower sensitivity; used
             to check the receptor ordering holds on a second platform.

The subject is Oprl1: its level, its detection rate, and which division of the
ganglion carries it. The other opioid genes are reported for context.

Levels are summarised as pseudobulk means, never as cross-platform detection
percentages. Each matrix passes `check_markers` before any opioid number is read
out of it, and every opioid gene carries an `in_matrix` flag taken from the
matrix rather than assumed: GSE102443 does not quantify Pnoc at all, so that
gene is recorded as unmeasured rather than as a zero.

Output feeds src/04_synthesis.py.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import atlas_common as ac
import atlas_style as st

N_BOOT = 2000


def load_dvoryanchikov():
    """96 geniculate neurons, transcript FPKM collapsed to gene level."""
    df = pd.read_csv(ac.DATA / ac.GENICULATE_FPKM, sep="\t", low_memory=False)
    cell_cols = [c for c in df.columns
                 if c not in ("Transcript", "Chromosome", "Begin", "End",
                              "Strand", "Gene")]
    # Cell labels carry plate coordinates; keep the trailing sample ID.
    names = [c.split("]]")[-1].strip() if "]]" in c else c.strip()
             for c in cell_cols]
    if len(set(names)) != len(names):
        raise ValueError("cell labels are not unique after stripping plate "
                         "coordinates; gene means would double-count cells")
    mat = df[cell_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    mat.columns = names
    mat["Gene"] = df["Gene"].astype(str).values
    genes = mat.groupby("Gene").sum()
    print(f"  GSE102443: {genes.shape[1]} cells x {genes.shape[0]:,} genes (FPKM)")
    return genes


def load_zuker():
    """454 Phox2b+ geniculate neurons, raw counts -> CPM."""
    x = pd.read_excel(ac.DATA / ac.GENICULATE_XLSX, index_col=0)
    x = x.apply(pd.to_numeric, errors="coerce").fillna(0.0)

    # Orient genes as rows. Decided by which axis actually carries gene symbols,
    # not by which axis is longer: the shape heuristic silently transposes any
    # matrix with more cells than genes.
    probe = set(ac.SANITY_GENES) | set(ac.OPIOID_GENES)
    in_rows = len(probe & set(map(str, x.index)))
    in_cols = len(probe & set(map(str, x.columns)))
    if in_cols > in_rows:
        x = x.T
    elif in_rows == 0:
        raise ValueError("neither axis of the GSE135801 matrix carries "
                         "recognisable gene symbols")

    # CPM assumes raw counts. Say so out loud if the file is already normalised:
    # ranks and ratios would survive it, but the absolute CPM would not.
    sample = x.to_numpy()[:2000]
    if not np.allclose(sample, np.round(sample)):
        print("  [warn] GSE135801 values are not integers; the file may already "
              "be normalised, in which case absolute CPM is not meaningful "
              "(ranks and within-sample ratios are unaffected)")

    lib = x.sum(axis=0)
    lib[lib == 0] = 1.0
    cpm = x.div(lib, axis=1) * 1e6
    print(f"  GSE135801: {cpm.shape[1]} cells x {cpm.shape[0]:,} genes (CPM)")
    return cpm


def summarise(mat, label, unit):
    """Pseudobulk mean and detection for one matrix, plus the opioid table."""
    mean = mat.mean(axis=1)
    det = (mat > 0).mean(axis=1) * 100.0

    rows = []
    for g in ac.OPIOID_GENES:
        here = g in mat.index
        rows.append({
            "dataset": label,
            "tissue": "geniculate",
            "assay": ac.ASSAY_LABEL[label],
            "gene": g,
            "unit": unit,
            # A gene absent from the annotation gets a null level, not 0.0:
            # "never measured" and "measured at zero" are different evidence.
            "mean_level": round(float(mean[g]), 4) if here else np.nan,
            "pct_detected": round(float(det[g]), 2) if here else np.nan,
            "n_cells": int(mat.shape[1]),
            "n_genes_quantified": int(mat.shape[0]),
            "in_matrix": bool(here),
            "transcriptome_percentile":
                round(ac.transcriptome_percentile(mean, g), 1) if here else np.nan,
        })
    return pd.DataFrame(rows), mean


def main() -> int:
    ac.RES.mkdir(parents=True, exist_ok=True)

    dv = load_dvoryanchikov()
    zk = load_zuker()

    # Marker gate: raises before any opioid number is read if the matrix is not
    # what it claims to be.
    checks = pd.concat([ac.check_markers(dv.mean(axis=1), "GSE102443"),
                        ac.check_markers(zk.mean(axis=1), "GSE135801")],
                       ignore_index=True)
    ac.save_table(checks, "geniculate_marker_checks.csv")

    dv_tbl, dv_mean = summarise(dv, "GSE102443", "FPKM")
    zk_tbl, zk_mean = summarise(zk, "GSE135801", "CPM")
    tbl = pd.concat([dv_tbl, zk_tbl], ignore_index=True)
    ac.save_table(tbl, "geniculate_opioid_levels.csv")

    print("\n  Opioid gene levels, geniculate:")
    print(tbl.pivot(index="gene", columns="dataset",
                    values="mean_level").reindex(ac.OPIOID_GENES).to_string())
    print("\n  Transcriptome percentile within each dataset:")
    print(tbl.pivot(index="gene", columns="dataset",
                    values="transcriptome_percentile")
          .reindex(ac.OPIOID_GENES).to_string())

    # Receptor ordering, the statistic that survives the platform difference,
    # with a cell bootstrap so a narrow ordering is not read as a wide one.
    ranks, support = [], []
    for label, mean, mat in (("GSE102443", dv_mean, dv), ("GSE135801", zk_mean, zk)):
        r = ac.receptor_rank(mean)
        r.insert(0, "dataset", label)
        r.insert(0, "tissue", "geniculate")
        ranks.append(r)

        per_cell = mat.loc[[g for g in ac.RECEPTORS if g in mat.index]].T
        b = ac.bootstrap_receptor_support(per_cell, n_boot=N_BOOT)
        b.update({"dataset": label, "tissue": "geniculate"})
        support.append(b)

    rank_tbl = pd.concat(ranks, ignore_index=True)
    ac.save_table(rank_tbl, "geniculate_receptor_rank.csv")
    print("\n  Receptor rank within each dataset:")
    print(rank_tbl.to_string(index=False))

    sup = pd.DataFrame(support)[["tissue", "dataset", "top_gene", "runner_up",
                                 "margin", "margin_lo", "margin_hi", "support",
                                 "n_cells", "n_boot"]]
    ac.save_table(sup, "geniculate_rank_support.csv")
    print("\n  Bootstrap support for the top receptor:")
    print(sup.round(3).to_string(index=False))

    # Per-cell Oprl1 for the deep dataset, split by gustatory/somatosensory as
    # the source project does (Phox2b FPKM > 5).
    per_cell = None
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

    figures(dv, zk, per_cell)
    return 0


def figures(dv, zk, per_cell):
    """Oprl1 in the geniculate ganglion: mean expression, and per-neuron spread.

    One message per panel, each read directly off the axis in the unit the
    dataset was measured in.
    """
    st.set_theme()
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.6))

    # (a) The four opioid receptors, deep full-length data.
    ax = axes[0]
    vals = [float(dv.loc[g].mean()) if g in dv.index else np.nan
            for g in ac.RECEPTORS]
    order = np.argsort(-np.array(vals))
    st.expression_bars(ax, [vals[i] for i in order],
                       [ac.RECEPTORS[i] for i in order],
                       "Mean expression (FPKM)")
    st.panel_letter(ax, "a")
    ax.set_title("GSE102443 — 96 geniculate neurons", fontsize=11, pad=12)

    # (b) Every neuron, both divisions of the ganglion.
    ax = axes[1]
    if per_cell is not None:
        groups = ["Gustatory", "Somato-\nsensory"]
        by = {"Gustatory":
              per_cell.loc[per_cell.division.str.startswith("gustatory"),
                           "Oprl1_FPKM"].values,
              "Somato-\nsensory":
              per_cell.loc[per_cell.division.str.startswith("somato"),
                           "Oprl1_FPKM"].values}
        st.violin_points(ax, groups, by, "Oprl1 FPKM")
    st.panel_letter(ax, "b")
    ax.set_title("Oprl1 in individual neurons", fontsize=11, pad=12)

    # (c) The same four receptors on an independent platform.
    ax = axes[2]
    vals = [float(zk.loc[g].mean()) if g in zk.index else np.nan
            for g in ac.RECEPTORS]
    order = np.argsort(-np.array(vals))
    st.expression_bars(ax, [vals[i] for i in order],
                       [ac.RECEPTORS[i] for i in order],
                       "Mean expression (CPM)")
    st.panel_letter(ax, "c")
    ax.set_title("GSE135801 — 454 geniculate neurons", fontsize=11, pad=12)

    fig.tight_layout()
    st.save(fig, "figure1_geniculate_oprl1")


if __name__ == "__main__":
    raise SystemExit(main())
