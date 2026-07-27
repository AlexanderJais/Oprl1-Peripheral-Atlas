"""Nodose and jugular ganglia: where Oprl1 sits in the NodoMap atlas.

Reads the integrated atlas published with NodoMap (Cheng et al.), 106,436 cells
in 52 author clusters, and asks the only question this repo is about: which
cells carry Oprl1, and how strongly.

The figures deliberately use the same idiom as the sibling PNOC-Nodose project
(`src/nodomap_style.py`): a FeaturePlot on the published UMAP, ranked bars of
percent-expressing across neuronal clusters, and a dot plot of the opioid panel
across all 52 clusters.

Two differences from the sibling project's own analysis, both required for
comparability with the geniculate and NTS data here:

  * levels are pseudobulk mean CPM computed from raw counts, so the geniculate
    means and these are the same kind of quantity (they are still not
    interchangeable in magnitude, which is why only ranks cross tissues);
  * the five constituent datasets are summarised separately as well as pooled,
    because the in-house snRNA-seq inflates long-intron genes such as Oprm1 and
    would otherwise distort the receptor ordering.

Preparation type is read from the atlas's own `suspension_type` field rather
than hardcoded, and checked against the registry in atlas_common.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import anndata as ad
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import atlas_common as ac
import atlas_style as st

NEURON_CLASSES = ["Nodose ganglion neuron", "Jugular ganglion neuron"]

# Cluster order as published by the NodoMap authors, matching the sibling
# project's `nodomap_style.CLUSTER_ORDER`.
CLUSTER_ORDER = (
    ["EC1", "EC2", "EC3"]
    + ["FB1", "FB2", "FB3", "FB4"]
    + ["HC1", "HC2"]
    + ["GC1", "GC2", "GC3"]
    + [f"MGC{i}" for i in range(1, 7)]
    + [f"SGC{i}" for i in range(1, 9)]
    + [f"JGN{i}" for i in range(1, 6)]
    + [f"NGN{i}" for i in range(1, 22)]
)
NEURON_CLUSTERS = [c for c in CLUSTER_ORDER if c.startswith(("NGN", "JGN"))]

N_BOOT = 2000


def load_counts(adata, path, genes):
    """Raw counts for a few genes, as a cells x genes frame.

    Returns (frame, absent). A symbol the atlas does not carry is reported as
    absent rather than silently becoming a column of zeros: "not in this
    annotation" and "measured at zero" are different evidence.
    """
    sym = adata.raw.var["feature_name"].astype(str).values
    by_symbol = {}
    for j, s in enumerate(sym):
        by_symbol.setdefault(s, []).append(j)

    cols = {g: by_symbol[g] for g in genes if g in by_symbol}
    absent = [g for g in genes if g not in by_symbol]
    dupes = {g: len(v) for g, v in cols.items() if len(v) > 1}
    if dupes:
        print(f"  [note] summed duplicate annotation rows for {dupes}")
    if absent:
        print(f"  [warn] absent from the atlas annotation: {absent}")

    with h5py.File(path, "r") as fh:
        mat = ac.csr_gene_columns(fh["raw/X"], cols, adata.n_obs)
    return pd.DataFrame(mat, columns=list(cols), index=adata.obs_names), absent


def pseudobulk_cpm(counts, total_counts, mask):
    """Mean per-cell CPM over `mask`, matching how the geniculate CPM is built."""
    sub = counts.loc[mask]
    lib = total_counts[mask].astype(float)
    lib[lib == 0] = 1.0
    return (sub.div(lib, axis=0) * 1e6).mean()


def library_size(obs, counts):
    """The authors' recorded library size, or the observed sum if it is missing."""
    if "nCount_RNA" in obs.columns:
        total = obs["nCount_RNA"].values.astype(float)
        if np.isfinite(total).all():
            return total
        print("  [warn] nCount_RNA has non-finite entries; falling back to row sums")
    else:
        print("  [warn] nCount_RNA absent; falling back to the sum over loaded genes")
    return counts.sum(axis=1).values.astype(float)


def opioid_rows(counts, total, mask, tissue, label, prep, absent):
    rows = []
    m = pseudobulk_cpm(counts, total, mask)
    det = (counts.loc[mask] > 0).mean() * 100.0
    for g in ac.OPIOID_GENES:
        here = g not in absent
        rows.append({
            "dataset": label, "tissue": tissue,
            "assay": ac.ASSAY_LABEL["NodoMap"], "prep": prep,
            "gene": g, "unit": "CPM",
            "mean_level": round(float(m[g]), 4) if here else np.nan,
            "pct_detected": round(float(det[g]), 2) if here else np.nan,
            "n_cells": int(mask.sum()),
            "in_matrix": bool(here),
        })
    return rows, m


def main() -> int:
    ac.RES.mkdir(parents=True, exist_ok=True)
    if not ac.NODOSE_H5AD.exists():
        raise SystemExit(
            f"NodoMap atlas not found at {ac.NODOSE_H5AD}.\n"
            "Set NODOSE_ROOT, or fetch it from CZ CELLxGENE:\n"
            "  curl -fsSL -o nodomap_integrated.h5ad \\\n"
            "    https://datasets.cellxgene.cziscience.com/"
            "a503329b-7bac-4a9f-be64-03cefc978452.h5ad")

    adata = ad.read_h5ad(ac.NODOSE_H5AD, backed="r")
    obs = adata.obs
    counts, absent = load_counts(adata, ac.NODOSE_H5AD,
                                 ac.OPIOID_GENES + ac.SANITY_GENES)
    total = library_size(obs, counts)
    print(f"  NodoMap: {adata.n_obs:,} cells, median library {np.median(total):,.0f}")

    # Preparation type from the atlas itself, then checked against the registry.
    prep_of = (obs.groupby("dataset", observed=True)["suspension_type"]
               .agg(lambda s: s.mode().iloc[0]).to_dict())
    prep_of = {k: ("nuclear" if v == "nucleus" else "whole cell")
               for k, v in prep_of.items()}
    for ds, prep in prep_of.items():
        declared = ac.DATASETS.get(f"NodoMap:{ds}", {}).get("prep")
        if declared and declared != prep:
            raise ac.SanityCheckError(
                f"NodoMap:{ds} is declared '{declared}' but the atlas says '{prep}'")
    print(f"  preparation from suspension_type: {prep_of}")

    is_neuron = obs["cell_class"].isin(NEURON_CLASSES).values
    is_nodose = (obs["cell_class"] == "Nodose ganglion neuron").values
    is_jugular = (obs["cell_class"] == "Jugular ganglion neuron").values

    ac.save_table(ac.check_markers(pseudobulk_cpm(counts, total, is_neuron),
                                   "NodoMap (neurons)"),
                  "nodose_marker_checks.csv")

    # ------------------------------------------------ levels, pooled and split
    rows, means = [], {}
    groups = [("nodose", "NodoMap", is_nodose, "mixed (pooled)"),
              ("jugular", "NodoMap", is_jugular, "mixed (pooled)"),
              ("nodose+jugular", "NodoMap", is_neuron, "mixed (pooled)")]
    for ds in obs["dataset"].cat.categories:
        groups.append(("nodose+jugular", f"NodoMap:{ds}",
                       is_neuron & (obs["dataset"] == ds).values, prep_of[ds]))

    for tissue, label, mask, prep in groups:
        if mask.sum() == 0:
            continue
        r, m = opioid_rows(counts, total, mask, tissue, label, prep, absent)
        rows += r
        means[(tissue, label)] = m
    tbl = pd.DataFrame(rows)
    ac.save_table(tbl, "nodose_opioid_levels.csv")

    print("\n  Opioid gene levels, mean CPM:")
    print(tbl[tbl.dataset == "NodoMap"].pivot(index="gene", columns="tissue",
          values="mean_level").reindex(ac.OPIOID_GENES).to_string())

    # ------------------------------------------------------- receptor ordering
    ranks, support = [], []
    for (tissue, label), m in means.items():
        r = ac.receptor_rank(m)
        r.insert(0, "dataset", label)
        r.insert(0, "tissue", tissue)
        ranks.append(r)
    rank_tbl = pd.concat(ranks, ignore_index=True)
    ac.save_table(rank_tbl, "nodose_receptor_rank.csv")

    for ds in obs["dataset"].cat.categories:
        mask = is_neuron & (obs["dataset"] == ds).values
        lib = total[mask].astype(float)
        lib[lib == 0] = 1.0
        per_cell = counts.loc[mask, [g for g in ac.RECEPTORS
                                     if g in counts.columns]].div(lib, axis=0) * 1e6
        b = ac.bootstrap_receptor_support(per_cell, n_boot=N_BOOT)
        b.update({"dataset": f"NodoMap:{ds}", "tissue": "nodose+jugular"})
        support.append(b)
    sup = pd.DataFrame(support)[["tissue", "dataset", "top_gene", "runner_up",
                                 "margin", "margin_lo", "margin_hi", "support",
                                 "n_cells", "n_boot"]]
    ac.save_table(sup, "nodose_rank_support.csv")
    print("\n  Bootstrap support for the top receptor:")
    print(sup.round(3).to_string(index=False))

    # ------------------------------------------------------- Oprl1 per cluster
    cluster = obs["author_cell_type"].astype(str).values
    lib_all = total.astype(float).copy()
    lib_all[lib_all == 0] = 1.0
    cpm_all = counts.div(lib_all, axis=0) * 1e6
    present = [c for c in CLUSTER_ORDER if c in set(cluster)]
    per_cluster = pd.DataFrame({
        "cluster": present,
        "n_cells": [int((cluster == c).sum()) for c in present],
        "cell_class": [obs["cell_class"][cluster == c].mode().iloc[0] for c in present],
        "is_neuron": [c.startswith(("NGN", "JGN")) for c in present],
        "Oprl1_CPM": [round(float(cpm_all.loc[cluster == c, "Oprl1"].mean()), 3)
                      for c in present],
        "Oprl1_pct": [round(float((counts.loc[cluster == c, "Oprl1"] > 0).mean() * 100), 2)
                      for c in present],
    })
    ac.save_table(per_cluster, "nodose_oprl1_by_cluster.csv")
    print("\n  Top clusters by Oprl1 detection:")
    print(per_cluster.sort_values("Oprl1_pct", ascending=False)
          .head(10).to_string(index=False))

    # --------------------------------------- neuronal enrichment, whole cell only
    # Pooling nuclear with whole-cell data here would contradict the stratification
    # the rest of this project rests on, so the in-house nuclear set is excluded.
    wc = obs["suspension_type"].values == "cell"
    m_neu = pseudobulk_cpm(counts, total, is_neuron & wc)
    m_non = pseudobulk_cpm(counts, total, (~is_neuron) & wc)
    enr = pd.DataFrame({
        "gene": ac.OPIOID_GENES,
        "neuron_CPM": [round(float(m_neu[g]), 4) if g not in absent else np.nan
                       for g in ac.OPIOID_GENES],
        "non_neuron_CPM": [round(float(m_non[g]), 4) if g not in absent else np.nan
                           for g in ac.OPIOID_GENES],
        "prep": "whole cell only",
        "n_neurons": int((is_neuron & wc).sum()),
        "n_non_neurons": int(((~is_neuron) & wc).sum()),
    })
    enr["log2_enrichment"] = np.log2((enr.neuron_CPM + 0.01) /
                                     (enr.non_neuron_CPM + 0.01)).round(3)
    ac.save_table(enr, "nodose_neuron_enrichment.csv")
    print("\n  Neuronal enrichment, whole-cell datasets only (mean CPM):")
    print(enr[["gene", "neuron_CPM", "non_neuron_CPM", "log2_enrichment"]]
          .to_string(index=False))

    figures(adata, obs, counts, cpm_all, cluster, per_cluster, tbl)
    return 0


def figures(adata, obs, counts, cpm_all, cluster, per_cluster, tbl):
    st.set_theme()
    umap = np.asarray(adata.obsm["X_umap"])

    # ---- Figure 2: Oprl1 across the atlas ---------------------------------
    fig = plt.figure(figsize=(15.5, 5.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.15, 0.95], wspace=0.32)

    ax = fig.add_subplot(gs[0])
    cmap, vmax = st.feature_plot(ax, umap, cpm_all["Oprl1"].values,
                                 "Oprl1", point_size=0.7)
    st.add_colorbar(fig, ax, cmap, vmax, label="Oprl1 (CPM)")
    ax.set_title("(A) Oprl1 on the published NodoMap UMAP\n"
                 f"{adata.n_obs:,} cells", fontsize=10)

    ax = fig.add_subplot(gs[1])
    neu = per_cluster[per_cluster.is_neuron]
    st.ranked_bars(ax, neu.Oprl1_pct, neu.cluster,
                   "% of cells expressing Oprl1",
                   annotate=[f"{v:.1f} CPM" for v in neu.Oprl1_CPM], fontsize=7,
                   threshold=30, threshold_label="NodoMap 30% positivity rule")
    ax.set_title("(B) Oprl1 across the 26 neuronal clusters\n"
                 "nodose (NGN) and jugular (JGN)", fontsize=10)

    ax = fig.add_subplot(gs[2])
    ds = tbl[(tbl.gene == "Oprl1") & tbl.dataset.str.contains(":")]
    st.ranked_bars(ax, ds.mean_level, [d.replace("NodoMap:", "") for d in ds.dataset],
                   "Oprl1 (mean CPM)",
                   annotate=[f"{p}" for p in ds.prep], fontsize=8)
    ax.set_title("(C) Oprl1 by constituent dataset", fontsize=10)

    st.save(fig, "figure2_nodose_oprl1")

    # ---- Figure 3: the opioid panel across all 52 clusters ----------------
    genes = [g for g in ac.OPIOID_GENES if g in counts.columns]
    order = per_cluster.cluster.tolist()
    idx = pd.Categorical(cluster, categories=order, ordered=True)
    pct = st.percent_expressing(counts[genes], idx).loc[order]
    lvl = st.mean_expression(cpm_all[genes], idx).loc[order]

    # Colour on a log scale: Oprm1 and Penk reach two orders of magnitude above
    # the rest, and on a linear scale they flatten every other row to one shade.
    fig, ax = plt.subplots(figsize=(17.5, 3.6))
    sc = st.dot_plot(ax, pct, np.log1p(lvl))
    # Separate non-neuronal from neuronal clusters, as the sibling project does.
    first_neuron = next(i for i, c in enumerate(order) if c.startswith(("JGN", "NGN")))
    ax.axvline(first_neuron - 0.5, color="black", lw=0.8, ls="--")
    ax.set_title("Opioid genes across all 52 NodoMap clusters "
                 "(dot size = % of cells expressing)", fontsize=11)
    cb = fig.colorbar(sc, ax=ax, fraction=0.015, pad=0.16, shrink=0.9)
    cb.set_label("log(1 + mean CPM)", size=7)
    cb.ax.tick_params(labelsize=6)
    st.dot_size_legend(ax, values=(1, 5, 10, 20, 30))
    st.save(fig, "figure3_nodose_opioid_dotplot")


if __name__ == "__main__":
    raise SystemExit(main())
