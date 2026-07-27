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

    figures(ann, corr, pb, keep)
    return 0


def figures(ann, corr, pb, clusters):
    st.set_theme()
    fig, axes = plt.subplots(1, 5, figsize=(21.0, 5.0))

    for k, (_, pretty) in enumerate(ANNOTATIONS):
        ax = axes[k]
        sub = ann[ann.annotation == pretty].sort_values("Oprl1_CPM", ascending=False)
        colours = [st.BAR_BLUE] * len(sub)
        st.expression_bars(ax, sub.Oprl1_CPM.values,
                           [f"{lv}" for lv in sub.level],
                           "Oprl1 (mean CPM)",
                           colors=colours, italic=False, annotate=True,
                           fontsize=11)
        for i, n in enumerate(sub.n_cells):
            ax.text(i, 0.4, f"n={n:,}", ha="center", va="bottom", fontsize=8,
                    color="white", rotation=90)
        st.panel_letter(ax, "abcd"[k], dx=-0.26)
        ax.set_title(pretty, fontsize=12, pad=10)

    # (d) the transcriptome-wide correlates, which is the unbiased answer.
    ax = axes[4]
    top = corr.head(15)
    st.expression_bars(ax, top.spearman_rho.values, top.gene.values,
                       "Spearman rho with Oprl1\nacross nodose clusters",
                       colors=[st.HIGHLIGHT] * len(top), annotate=False,
                       fontsize=11)
    ax.set_ylim(0, 1.05)
    st.panel_letter(ax, "e", dx=-0.30)
    ax.set_title(f"Top transcriptome-wide correlates\n"
                 f"(n = {len(clusters)} clusters, {len(corr):,} genes tested)",
                 fontsize=12, pad=10)

    fig.tight_layout()
    st.save(fig, "figure3_oprl1_localisation")


if __name__ == "__main__":
    raise SystemExit(main())
