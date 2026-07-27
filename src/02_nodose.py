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
    # Whole-cell-only groups carry the section 1 comparison; the pooled rows mix
    # in the nuclear dataset and are reported for completeness only.
    wc = obs["suspension_type"].values == "cell"
    groups = [("nodose", "NodoMap:whole-cell", is_nodose & wc, "whole cell"),
              ("jugular", "NodoMap:whole-cell", is_jugular & wc, "whole cell"),
              ("nodose", "NodoMap", is_nodose, "mixed (pooled)"),
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

    # Jugular is the only neural-crest-derived ganglion in this atlas, so it
    # tests whether the receptor ordering holds outside the placodal series.
    for label, mask in (("nodose", is_nodose & wc), ("jugular", is_jugular & wc)):
        lib = total[mask].astype(float)
        lib[lib == 0] = 1.0
        pc = counts.loc[mask, [g for g in ac.RECEPTORS
                               if g in counts.columns]].div(lib, axis=0) * 1e6
        b = ac.bootstrap_receptor_support(pc, n_boot=N_BOOT)
        b.update({"dataset": "NodoMap:whole-cell", "tissue": label})
        support.append(b)

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
    """Mean expression, read straight off the axis. No ratios, no ranks."""
    st.set_theme()
    whole_cell = obs["suspension_type"].values == "cell"
    is_nodose = (obs["cell_class"] == "Nodose ganglion neuron").values & whole_cell

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.8),
                             gridspec_kw={"width_ratios": [1.0, 2.1, 1.0]})

    # (a) The four opioid receptors in nodose neurons.
    ax = axes[0]
    vals = [float(cpm_all.loc[is_nodose, g].mean()) if g in cpm_all.columns else np.nan
            for g in ac.RECEPTORS]
    order = np.argsort(-np.array(vals))
    st.expression_bars(ax, [vals[i] for i in order],
                       [ac.RECEPTORS[i] for i in order],
                       "Mean expression (CPM)")
    st.panel_letter(ax, "a", dx=-0.28)
    ax.set_title(f"Nodose neurons\n{is_nodose.sum():,} cells, whole-cell datasets",
                 fontsize=11, pad=10)

    # (b) Oprl1 in every nodose cluster.
    ax = axes[1]
    ngn = per_cluster[per_cluster.cluster.str.startswith("NGN")] \
        .sort_values("Oprl1_CPM", ascending=False)
    st.expression_bars(ax, ngn.Oprl1_CPM.values, ngn.cluster.values,
                       "Oprl1 (mean CPM)", italic=False, rotation=90, fontsize=11)
    ax.axhline(float(cpm_all.loc[is_nodose, "Oprl1"].mean()), color="black",
               lw=1.2, ls="--")
    ax.text(len(ngn) - 0.4, float(cpm_all.loc[is_nodose, "Oprl1"].mean()),
            " ganglion mean", fontsize=9, va="bottom", ha="right")
    st.panel_letter(ax, "b", dx=-0.07)
    ax.set_title("Oprl1 is expressed in every nodose cluster", fontsize=11, pad=10)

    # (c) Oprl1 in each constituent whole-cell dataset.
    ax = axes[2]
    ds = tbl[(tbl.gene == "Oprl1") & tbl.dataset.str.contains(":")
             & (tbl.prep == "whole cell")]
    ds = ds.sort_values("mean_level", ascending=False)
    st.expression_bars(ax, ds.mean_level.values,
                       [d.replace("NodoMap:", "") for d in ds.dataset],
                       "Oprl1 (mean CPM)", italic=False, fontsize=11)
    st.panel_letter(ax, "c", dx=-0.28)
    ax.set_title("Reproduces across all four\nwhole-cell datasets", fontsize=11, pad=10)

    fig.tight_layout()
    st.save(fig, "figure2_nodose_oprl1")

    # ---- Supplementary: the full opioid panel, neurons against everything else
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.6), sharey=True)
    for ax, (mask, title) in zip(axes, (
            (is_nodose, f"Nodose neurons (n = {is_nodose.sum():,})"),
            ((~obs["cell_class"].isin(NEURON_CLASSES).values) & whole_cell,
             "Non-neuronal cells"))):
        vals = [float(cpm_all.loc[mask, g].mean()) if g in cpm_all.columns else np.nan
                for g in ac.OPIOID_GENES]
        colours = [st.BAR_BLUE if g in ac.RECEPTORS else st.BAR_GREY
                   for g in ac.OPIOID_GENES]
        st.expression_bars(ax, vals, ac.OPIOID_GENES, "Mean expression (CPM)",
                           colors=colours, fontsize=11)
        ax.set_title(title, fontsize=11, pad=10)
        ax.set_ylim(0, 17)
    st.panel_letter(axes[0], "a", dx=-0.16)
    st.panel_letter(axes[1], "b", dx=-0.10)
    fig.suptitle("Opioid receptors (blue) are neuronal; the peptides (grey) are not",
                 fontsize=12, y=1.02)
    fig.tight_layout()
    st.save(fig, "figureS2_nodose_opioid_panel")


if __name__ == "__main__":
    raise SystemExit(main())
