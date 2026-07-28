"""Where do the Oprl1 neurons sit? An unbiased answer.

Oprl1 is not a marker of the Glp1r/Cckar populations, so the obvious next
question is what it IS a marker of. Two analyses, neither of which starts from a
candidate gene:

  1. The NodoMap authors annotated every nodose cluster by organ projection,
     fibre type, sensor type and sodium-channel class. Oprl1 mean expression is
     computed in each, so the answer is stated in the atlas's own vocabulary.

  2. A transcriptome-wide correlation. Pseudobulk mean CPM is computed for all
     54,640 genes in each of the 21 nodose clusters, and every gene is
     correlated with Oprl1 across those clusters. Nothing is pre-selected; the
     partners come out of the data.

Cluster-level pseudobulk is the right resolution for this. Per-cell correlation
between two sparsely detected transcripts is dominated by capture depth;
averaging within a cluster removes dropout before the correlation is taken.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import anndata as ad
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

import atlas_common as ac
import atlas_style as st

_nodose = __import__("02_nodose")
library_size = _nodose.library_size

ANNOTATIONS = [
    ("organ_projection", "Organ projection"),
    ("fibre_type", "Fibre type"),
    ("sensor_type", "Sensor type"),
    ("sodium_channel_type", "Sodium channel"),
]
MIN_CELLS = 30
MIN_CPM = 1.0          # a gene must reach this in at least MIN_CLUSTERS clusters
MIN_CLUSTERS = 3
CACHE = ac.DATA / "nodomap_cluster_pseudobulk.npz"


def cluster_pseudobulk(path, codes, n_codes, n_genes, lib, force=False):
    """Mean per-cell CPM for every gene in every cluster, in one streamed pass."""
    if CACHE.exists() and not force:
        z = np.load(CACHE)
        if z["mat"].shape == (n_codes, n_genes):
            print(f"  [cache] {CACHE.name}")
            return z["mat"]

    print(f"  streaming the full matrix for {n_genes:,} genes x {n_codes} clusters ...")
    sums = np.zeros((n_codes, n_genes), dtype=np.float64)
    counts = np.bincount(codes[codes >= 0], minlength=n_codes).astype(float)
    with h5py.File(path, "r") as fh:
        grp = fh["raw/X"]
        indptr = grp["indptr"][:]
        data_ds, idx_ds = grp["data"], grp["indices"]
        n_cells = len(indptr) - 1
        block = 20_000
        for start in range(0, n_cells, block):
            stop = min(start + block, n_cells)
            lo, hi = int(indptr[start]), int(indptr[stop])
            m = sp.csr_matrix(
                (data_ds[lo:hi], idx_ds[lo:hi], indptr[start:stop + 1] - lo),
                shape=(stop - start, n_genes))
            # Per-cell CPM first, then the mean over cells: the same quantity
            # every other table in this project reports.
            m = sp.diags(1e6 / lib[start:stop]) @ m
            blk = codes[start:stop]
            for c in np.unique(blk):
                if c < 0:
                    continue
                sums[c] += np.asarray(m[blk == c].sum(axis=0)).ravel()
            print(f"    {stop:,}/{n_cells:,} cells", flush=True)

    mat = sums / np.maximum(counts, 1)[:, None]
    np.savez_compressed(CACHE, mat=mat)
    print(f"  [cache] wrote {CACHE.name}")
    return mat


def benjamini_hochberg(p):
    p = np.asarray(p, dtype=float)
    n = p.size
    order = np.argsort(p)
    q = np.empty(n)
    q[order] = np.minimum.accumulate((p[order] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def main() -> int:
    ac.RES.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(ac.NODOSE_H5AD, backed="r")
    obs = adata.obs
    symbols = adata.raw.var["feature_name"].astype(str).values

    whole_cell = obs["suspension_type"].values == "cell"
    is_nodose = (obs["cell_class"] == "Nodose ganglion neuron").values & whole_cell
    cluster = obs["author_cell_type"].astype(str).values

    # Library size for CPM, from a small gene load (Oprl1 only is enough).
    counts, _ = _nodose.load_counts(adata, ac.NODOSE_H5AD, ["Oprl1"])
    lib = library_size(obs, counts).astype(float)
    lib[lib == 0] = 1.0
    oprl1_cpm = (counts["Oprl1"].values / lib) * 1e6
    detected = counts["Oprl1"].values > 0

    # ---------------------------------------- 1. the atlas's own annotations
    rows = []
    for col, pretty in ANNOTATIONS:
        vals = obs[col].astype(str).values
        for level in pd.unique(vals[is_nodose]):
            m = is_nodose & (vals == level)
            if m.sum() < MIN_CELLS or level in ("NA", "nan"):
                continue
            rows.append({
                "annotation": pretty, "level": level, "n_cells": int(m.sum()),
                "Oprl1_CPM": round(float(oprl1_cpm[m].mean()), 3),
                "Oprl1_pct": round(float(detected[m].mean() * 100), 2),
            })
    ann = pd.DataFrame(rows)
    ac.save_table(ann, "nodose_oprl1_by_annotation.csv")
    print("\n  Oprl1 by the atlas's own annotations (nodose neurons, whole cell):")
    for _, pretty in ANNOTATIONS:
        sub = ann[ann.annotation == pretty].sort_values("Oprl1_CPM", ascending=False)
        print(f"\n  {pretty}:")
        print(sub[["level", "n_cells", "Oprl1_CPM", "Oprl1_pct"]].to_string(index=False))

    # Every annotation above is constant within a cluster, so a cell-level test
    # would treat 26,047 correlated observations as independent. The honest unit
    # is the cluster; the per-cell means are kept only as a description.
    cl_rows = []
    for c in pd.unique(cluster[is_nodose]):
        m = is_nodose & (cluster == c)
        if m.sum() < MIN_CELLS:
            continue
        row = {"cluster": c, "n_cells": int(m.sum()),
               "Oprl1_CPM": round(float(oprl1_cpm[m].mean()), 3),
               "Oprl1_pct": round(float(detected[m].mean() * 100), 2)}
        for col, pretty in ANNOTATIONS:
            row[col] = str(pd.Series(obs[col].astype(str).values[m]).mode().iloc[0])
        cl_rows.append(row)
    per_cluster = pd.DataFrame(cl_rows).sort_values("Oprl1_CPM", ascending=False)
    ac.save_table(per_cluster, "nodose_oprl1_by_cluster_annotated.csv")

    print("\n  Cluster-level test (n = %d clusters, not cells):" % len(per_cluster))
    kw_rows = []
    for col, pretty in ANNOTATIONS:
        groups = [g["Oprl1_CPM"].values for _, g in per_cluster.groupby(col)
                  if len(g) >= 2]
        if len(groups) < 2:
            continue
        h, p = stats.kruskal(*groups)
        kw_rows.append({"annotation": pretty, "n_groups": len(groups),
                        "n_clusters": len(per_cluster), "kruskal_H": round(float(h), 3),
                        "p_value": float(p)})
        print(f"    {pretty:18s} Kruskal-Wallis H = {h:6.3f}  p = {p:.4f}"
              f"  ({len(groups)} groups)")
    ac.save_table(pd.DataFrame(kw_rows), "nodose_oprl1_annotation_tests.csv")

    # ------------------------------------ 2. transcriptome-wide correlation
    ngn = sorted({c for c in cluster[is_nodose] if c.startswith("NGN")},
                 key=lambda s: int(s[3:]))
    keep = [c for c in ngn if (is_nodose & (cluster == c)).sum() >= MIN_CELLS]
    code_of = {c: i for i, c in enumerate(keep)}
    codes = np.full(adata.n_obs, -1, dtype=np.int64)
    for c, i in code_of.items():
        codes[is_nodose & (cluster == c)] = i

    mat = cluster_pseudobulk(ac.NODOSE_H5AD, codes, len(keep), adata.n_vars, lib)
    pb = pd.DataFrame(mat, index=keep, columns=symbols)
    pb = pb.T.groupby(level=0).sum().T          # sum duplicate symbols

    expressed = (pb >= MIN_CPM).sum(axis=0) >= MIN_CLUSTERS
    pb = pb.loc[:, expressed]
    print(f"\n  {pb.shape[1]:,} genes expressed in >= {MIN_CLUSTERS} of "
          f"{len(keep)} nodose clusters; correlating each with Oprl1")

    # One gene against all others. spearmanr() on the whole frame would build a
    # 16k x 16k matrix; ranking once and correlating a single vector is the same
    # statistic in a fraction of the time and memory.
    R = pb.rank(axis=0).values.astype(float)
    names = np.asarray(pb.columns)
    t = R[:, int(np.nonzero(names == "Oprl1")[0][0])]
    Rc = R - R.mean(axis=0)
    tc = t - t.mean()
    with np.errstate(invalid="ignore", divide="ignore"):
        rho = ((Rc * tc[:, None]).sum(axis=0)
               / np.sqrt((Rc ** 2).sum(axis=0) * (tc ** 2).sum()))
        n = R.shape[0]
        tstat = rho * np.sqrt((n - 2) / np.clip(1 - rho ** 2, 1e-12, None))
    pval = 2 * stats.t.sf(np.abs(tstat), n - 2)

    keep_g = np.isfinite(rho) & (names != "Oprl1")
    corr = pd.DataFrame({
        "gene": names[keep_g],
        "spearman_rho": np.round(rho[keep_g], 4),
        "p_value": pval[keep_g],
    })
    corr["q_value"] = np.round(benjamini_hochberg(corr.p_value.values), 4)
    corr["mean_CPM_across_clusters"] = [round(float(pb[g].mean()), 2)
                                        for g in corr.gene]
    corr = corr.sort_values("spearman_rho", ascending=False).reset_index(drop=True)
    ac.save_table(corr, "nodose_oprl1_gene_correlations.csv")

    sig = corr[corr.q_value < 0.05]
    print(f"\n  {len(sig):,} genes correlate with Oprl1 at FDR < 5% "
          f"(n = {len(keep)} clusters)")
    print("\n  Top 25 positive correlates:")
    print(corr.head(25).to_string(index=False))
    print("\n  Top 10 negative correlates:")
    print(corr.tail(10).iloc[::-1].to_string(index=False))

    # Where the named satiation receptors fall in that ranking, for contrast.
    named = corr[corr.gene.isin(ac.SATIATION_GENES)]
    if len(named):
        print("\n  For contrast, the satiation receptors in the same ranking:")
        for _, r in named.iterrows():
            pos = int(corr.index[corr.gene == r.gene][0]) + 1
            print(f"    {r.gene:8s} rho = {r.spearman_rho:+.3f}  q = {r.q_value:.3f}"
                  f"   rank {pos:,} of {len(corr):,}")

    figures(corr, per_cluster, np.asarray(adata.obsm["X_umap"]), oprl1_cpm,
            is_nodose, obs["fibre_type"].astype(str).values)
    return 0


def figures(corr, per_cluster, umap, oprl1_cpm, is_nodose, fibre_percell):
    """Deliberately mixed panel types: an embedding, distributions, a scatter.

    The finding is that several independent annotations agree, and five bar
    charts side by side make that harder to see, not easier.
    """
    st.set_theme()
    fig = plt.figure(figsize=(16.5, 9.0))
    gs = fig.add_gridspec(2, 3, hspace=0.36, wspace=0.30,
                          height_ratios=[1.05, 0.95])

    fibre_order = ["Myelinated", "Lightly myelinated", "Unmyelinated"]
    fibre_cols = {"Myelinated": "#08306B", "Lightly myelinated": "#6BAED6",
                  "Unmyelinated": "#F0A202"}

    # (a) the ganglion, coloured by fibre type
    ax = fig.add_subplot(gs[0, 0])
    lab = np.where(is_nodose, fibre_percell, "")
    st.dim_plot(ax, umap, lab, fibre_order, fibre_cols,
                background=~is_nodose, point_size=1.1)
    ax.legend(loc="lower left", fontsize=7.5, markerscale=5, handletextpad=0.2)
    st.panel_letter(ax, "a", dx=-0.02, dy=1.02)
    ax.set_title("Nodose neurons by fibre type", fontsize=11)

    # (b) Oprl1 on the same embedding
    ax = fig.add_subplot(gs[0, 1])
    v = np.where(is_nodose, oprl1_cpm, np.nan)
    ax.scatter(umap[~is_nodose, 0], umap[~is_nodose, 1], s=0.35, c="#EDEDED",
               linewidths=0, rasterized=True)
    idx = np.argsort(np.nan_to_num(v))
    idx = idx[is_nodose[idx]]
    sc = ax.scatter(umap[idx, 0], umap[idx, 1], c=v[idx], cmap="viridis", s=1.4,
                    vmin=0, vmax=float(np.nanpercentile(v[is_nodose], 99)),
                    linewidths=0, rasterized=True)
    ax.set_aspect("equal")
    ax.axis("off")
    cb = fig.colorbar(sc, ax=ax, fraction=0.03, pad=0.01, shrink=0.6)
    cb.set_label("Oprl1 (CPM)", size=8)
    cb.ax.tick_params(labelsize=7)
    st.panel_letter(ax, "b", dx=-0.02, dy=1.02)
    ax.set_title("Oprl1", fontsize=11, fontstyle="italic")

    # (c) Oprl1 by organ projection, one point per cluster
    ax = fig.add_subplot(gs[0, 2])
    organ_order = ["Pancreas", "Duodenum", "Heart", "Gut", "Broad projection",
                   "Jejunum/Ileum"]
    organ_order = [o for o in organ_order if o in set(per_cluster.organ_projection)]
    st.strip_by_group(ax, organ_order,
                      {o: per_cluster.loc[per_cluster.organ_projection == o,
                                          "Oprl1_CPM"] for o in organ_order},
                      "Oprl1 (mean CPM per cluster)")
    ax.set_xticklabels(organ_order, rotation=35, ha="right", fontsize=9)
    st.panel_letter(ax, "c", dx=-0.24)
    ax.set_title("Organ projection", fontsize=11)

    # (d) fibre type, one point per cluster
    ax = fig.add_subplot(gs[1, 0])
    st.strip_by_group(ax, fibre_order,
                      {f: per_cluster.loc[per_cluster.fibre_type == f, "Oprl1_CPM"]
                       for f in fibre_order}, "Oprl1 (mean CPM per cluster)")
    ax.set_xticklabels(["Myelin-\nated", "Lightly\nmyelinated", "Unmyelin-\nated"],
                       fontsize=9)
    st.panel_letter(ax, "d", dx=-0.24)
    ax.set_title("Fibre type", fontsize=11)

    # (e) sodium channel class, one point per cluster
    ax = fig.add_subplot(gs[1, 1])
    nav_order = [n for n in ["Nav1.1", "Nav1.1/Nav1.8", "Nav1.8"]
                 if n in set(per_cluster.sodium_channel_type)]
    st.strip_by_group(ax, nav_order,
                      {n: per_cluster.loc[per_cluster.sodium_channel_type == n,
                                          "Oprl1_CPM"] for n in nav_order},
                      "Oprl1 (mean CPM per cluster)")
    st.panel_letter(ax, "e", dx=-0.24)
    ax.set_title("Sodium channel class", fontsize=11)

    # (f) the unbiased scan: every expressed gene, correlation against abundance
    ax = fig.add_subplot(gs[1, 2])
    ax.scatter(corr.spearman_rho, corr.mean_CPM_across_clusters, s=4,
               c="#D8D8D8", linewidths=0, rasterized=True)
    top = corr.head(8)
    ax.scatter(top.spearman_rho, top.mean_CPM_across_clusters, s=32,
               c=st.BAR_BLUE, edgecolors="black", linewidths=0.5, zorder=3)
    for _, r in top.iterrows():
        ax.annotate(r.gene, (r.spearman_rho, r.mean_CPM_across_clusters),
                    fontsize=7, fontstyle="italic", xytext=(-6, 5),
                    textcoords="offset points", ha="right")
    named = corr[corr.gene.isin(["Glp1r", "Cckar", "Cckbr"])]
    ax.scatter(named.spearman_rho, named.mean_CPM_across_clusters, s=44,
               c=st.HIGHLIGHT, edgecolors="black", linewidths=0.6, zorder=4)
    for _, r in named.iterrows():
        ax.annotate(r.gene, (r.spearman_rho, r.mean_CPM_across_clusters),
                    fontsize=8, fontstyle="italic", color=st.HIGHLIGHT,
                    xytext=(6, -3), textcoords="offset points")
    ax.axvline(0, color="black", lw=0.8, ls="--")
    ax.set_yscale("log")
    ax.set_xlabel("Spearman rho with Oprl1 across clusters", fontsize=10)
    ax.set_ylabel("mean CPM across clusters", fontsize=10)
    st.panel_letter(ax, "f", dx=-0.24)
    ax.set_title(f"All {len(corr):,} expressed genes", fontsize=11)

    st.save(fig, "figureS5_oprl1_annotations")


if __name__ == "__main__":
    raise SystemExit(main())
