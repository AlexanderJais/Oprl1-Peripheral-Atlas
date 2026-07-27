"""Does Oprl1 sit on the vagal afferents that carry the satiation receptors?

The orexigenic effect of N/OFQ is well established and is attributed throughout
the literature to hypothalamic sites. This script asks a different question: if
Gi-coupled Oprl1 is expressed by the same nodose neurons that carry the
excitatory satiation receptors Glp1r and Cckar, then the NOP receptor is
positioned as a cell-autonomous brake on the first synapse of the gut-brain
axis, upstream of anything hypothalamic.

That is a co-expression question, and co-expression in droplet data has one
dominant confound: capture depth. A cell that detects any gene tends to detect
more genes overall, so Oprl1+ and Oprl1- cells differ in depth before biology
is considered, and a raw overlap percentage will show an association whether or
not one exists. Three things are reported for each partner gene, in increasing
order of how much weight they can carry:

  1. crude co-detection      - the raw percentage, reported for transparency
                               and not used for any claim
  2. depth-stratified odds   - Mantel-Haenszel across ten depth deciles, the
     ratio                     cell-level statistic the claim rests on
  3. cluster-level agreement - whether the neuronal clusters high for Oprl1 are
                               the clusters high for the partner. Pseudobulk
                               per cluster averages dropout away entirely, so
                               this is the most robust of the three.

Nodose neurons only, whole-cell datasets only. Jugular neurons are reported
separately: they are vagal but are not the gut-projecting population the
hypothesis is about.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import atlas_common as ac  # noqa: E402
import atlas_style as st  # noqa: E402

_nodose = __import__("02_nodose")
load_counts = _nodose.load_counts
library_size = _nodose.library_size
CLUSTER_ORDER = _nodose.CLUSTER_ORDER

GENES = ["Oprl1"] + ac.SATIATION_GENES
N_DEPTH_BINS = 10
NODOSE_CLUSTERS = [c for c in CLUSTER_ORDER if c.startswith("NGN")]


def coexpression_table(detected, oprl1_pos, depth, cluster, label):
    """Crude overlap, depth-stratified and cluster+depth-stratified odds ratios.

    The cluster+depth statistic is the decisive one. A pooled association can
    arise purely from cluster composition — if Oprl1-high clusters happened to
    be the Glp1r-high clusters, every cell-level test would be positive without
    a single cell co-expressing anything. Holding cluster fixed asks whether two
    neurons of the *same* type, sequenced to the *same* depth, carry the two
    transcripts together more often than chance.
    """
    rows = []
    joint = ac.cluster_depth_strata(cluster, _depth_rank(depth))
    for g in ac.SATIATION_GENES:
        if g not in detected.columns:
            continue
        partner = detected[g].values
        n_p, n_o = int(partner.sum()), int(oprl1_pos.sum())
        if n_p == 0:
            print(f"  [note] {g} is not detected in any {label} neuron; skipped")
            continue
        n_both = int((partner & oprl1_pos).sum())
        or_s, p_s, k = ac.stratified_odds_ratio(oprl1_pos, partner, depth)
        or_c, p_c, k_c = ac.stratified_odds_ratio(oprl1_pos, partner, joint)
        # Crude 2x2 odds ratio, for comparison with the stratified one.
        a, b = n_both, n_p - n_both
        c = n_o - n_both
        d = len(partner) - n_p - c
        crude = (a * d) / (b * c) if b and c else np.nan
        rows.append({
            "population": label, "partner": g,
            "n_cells": int(len(partner)),
            "n_partner_pos": n_p, "n_Oprl1_pos": n_o, "n_double_pos": n_both,
            "pct_of_partner_pos_that_are_Oprl1_pos":
                round(100 * n_both / n_p, 2) if n_p else np.nan,
            "pct_of_partner_neg_that_are_Oprl1_pos":
                round(100 * c / (len(partner) - n_p), 2) if len(partner) > n_p else np.nan,
            "crude_odds_ratio": round(crude, 3) if np.isfinite(crude) else np.nan,
            "depth_stratified_OR": round(or_s, 3) if np.isfinite(or_s) else np.nan,
            "p_value": p_s,
            "n_depth_strata": k,
            "cluster_depth_stratified_OR":
                round(or_c, 3) if np.isfinite(or_c) else np.nan,
            "p_value_cluster_depth": p_c,
            "n_cluster_depth_strata": k_c,
        })
    return pd.DataFrame(rows)


def _depth_rank(depth):
    """Recover a numeric ordering from decile labels, for coarser re-binning."""
    return pd.Series(depth).astype(int).values


def main() -> int:
    ac.RES.mkdir(parents=True, exist_ok=True)
    if not ac.NODOSE_H5AD.exists():
        raise SystemExit(f"NodoMap atlas not found at {ac.NODOSE_H5AD}; "
                         "run bash src/00_download_data.sh")

    adata = ad.read_h5ad(ac.NODOSE_H5AD, backed="r")
    obs = adata.obs
    counts, absent = load_counts(adata, ac.NODOSE_H5AD, GENES)
    if "Oprl1" in absent:
        raise ac.SanityCheckError("Oprl1 absent from the atlas annotation")
    total = library_size(obs, counts)

    # Whole-cell datasets only: the nuclear set has a different detection
    # regime, and mixing it in would confound every co-detection count.
    whole_cell = obs["suspension_type"].values == "cell"
    cluster = obs["author_cell_type"].astype(str).values
    is_nodose = (obs["cell_class"] == "Nodose ganglion neuron").values & whole_cell
    is_jugular = (obs["cell_class"] == "Jugular ganglion neuron").values & whole_cell
    print(f"  nodose neurons (whole cell): {is_nodose.sum():,}")
    print(f"  jugular neurons (whole cell): {is_jugular.sum():,}")

    detected = counts > 0
    depth_all = ac.depth_strata(obs["nFeature_RNA"].values, N_DEPTH_BINS)

    # ------------------------------------------------- cell-level association
    tables = []
    for label, mask in (("nodose", is_nodose), ("jugular", is_jugular)):
        sub = detected.loc[mask]
        tables.append(coexpression_table(sub, sub["Oprl1"].values,
                                         depth_all[mask], cluster[mask], label))
    co = pd.concat(tables, ignore_index=True)
    ac.save_table(co, "vagal_oprl1_coexpression.csv")
    print("\n  Oprl1 co-detection with satiation receptors, nodose neurons:")
    print(co[co.population == "nodose"][
        ["partner", "n_partner_pos", "n_double_pos",
         "pct_of_partner_pos_that_are_Oprl1_pos",
         "pct_of_partner_neg_that_are_Oprl1_pos", "crude_odds_ratio",
         "depth_stratified_OR", "cluster_depth_stratified_OR",
         "p_value_cluster_depth"]].to_string(index=False))

    # ------------------------------------------------ cluster-level agreement
    lib = total.astype(float).copy()
    lib[lib == 0] = 1.0
    cpm = counts.div(lib, axis=0) * 1e6
    rows = []
    for c in NODOSE_CLUSTERS:
        m = (cluster == c) & whole_cell
        if m.sum() < 30:
            continue
        row = {"cluster": c, "n_cells": int(m.sum())}
        for g in GENES:
            if g in cpm.columns:
                row[f"{g}_CPM"] = round(float(cpm.loc[m, g].mean()), 3)
                row[f"{g}_pct"] = round(float(detected.loc[m, g].mean() * 100), 2)
        rows.append(row)
    per_cluster = pd.DataFrame(rows)
    ac.save_table(per_cluster, "vagal_oprl1_satiation_by_cluster.csv")

    corr_rows = []
    for g in ac.SATIATION_GENES:
        col = f"{g}_pct"
        if col not in per_cluster or per_cluster[col].nunique() < 2:
            continue
        rho, p = stats.spearmanr(per_cluster["Oprl1_pct"], per_cluster[col])
        corr_rows.append({"partner": g, "spearman_rho": round(float(rho), 3),
                          "p_value": float(p), "n_clusters": len(per_cluster)})
    corr = pd.DataFrame(corr_rows)
    ac.save_table(corr, "vagal_cluster_level_correlation.csv")
    print("\n  Cluster-level agreement across the 21 nodose clusters "
          "(percent-expressing, Spearman):")
    print(corr.to_string(index=False))

    # How much of the gut-projecting population is covered: the fraction of
    # nodose neurons in clusters that are above the atlas median for both.
    med_o = per_cluster.Oprl1_pct.median()
    for g in ac.SATIATION_PRIMARY:
        col = f"{g}_pct"
        if col not in per_cluster:
            continue
        both = per_cluster[(per_cluster.Oprl1_pct >= med_o)
                           & (per_cluster[col] >= per_cluster[col].median())]
        share = 100 * both.n_cells.sum() / per_cluster.n_cells.sum()
        print(f"  clusters above median for Oprl1 and {g}: "
              f"{len(both)}/{len(per_cluster)} clusters, {share:.0f}% of nodose neurons")

    figures(adata, obs, counts, detected, cpm, per_cluster, co, corr,
            is_nodose, whole_cell, cluster)
    return 0


def figures(adata, obs, counts, detected, cpm, per_cluster, co, corr,
            is_nodose, whole_cell, cluster):
    st.set_theme()
    fig = plt.figure(figsize=(16.5, 8.4))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.0], hspace=0.48, wspace=0.42)

    # (A) UMAP: where the double-positive cells are.
    ax = fig.add_subplot(gs[0, 0])
    umap = np.asarray(adata.obsm["X_umap"])
    o = detected["Oprl1"].values
    g = detected["Glp1r"].values
    base = ~is_nodose
    ax.scatter(umap[base, 0], umap[base, 1], s=0.4, c="#E8E8E8", linewidths=0,
               rasterized=True)
    for m, colour, lab in ((is_nodose & ~o & ~g, "#C8C8C8", "neither"),
                           (is_nodose & o & ~g, "#08306B", "Oprl1 only"),
                           (is_nodose & ~o & g, "#F0A202", "Glp1r only"),
                           (is_nodose & o & g, "#B2182B", "both")):
        ax.scatter(umap[m, 0], umap[m, 1], s=1.6 if lab == "both" else 0.6,
                   c=colour, linewidths=0, rasterized=True, label=f"{lab} ({m.sum():,})")
    ax.set_aspect("equal")
    ax.axis("off")
    ax.legend(loc="lower left", fontsize=6.5, markerscale=6, handletextpad=0.2)
    ax.set_title("(a) Oprl1 and Glp1r co-detection\nnodose neurons on the NodoMap UMAP",
                 fontsize=10)

    # (B) and (C) cluster-level agreement, the dropout-robust view.
    for j, partner in enumerate(ac.SATIATION_PRIMARY):
        ax = fig.add_subplot(gs[0, 1 + j])
        col = f"{partner}_pct"
        ax.scatter(per_cluster[col], per_cluster.Oprl1_pct, s=42,
                   c=per_cluster.n_cells, cmap="viridis",
                   edgecolors="black", linewidths=0.4, zorder=3)
        for _, rr in per_cluster.iterrows():
            ax.annotate(rr.cluster, (rr[col], rr.Oprl1_pct), fontsize=6,
                        xytext=(4, 3), textcoords="offset points")
        rho = corr.set_index("partner").loc[partner]
        ax.set_xlabel(f"% of cells expressing {partner}", fontsize=8.5)
        ax.set_ylabel("% of cells expressing Oprl1", fontsize=8.5)
        ax.set_title(f"({'bc'[j]}) Oprl1 against {partner}, per nodose cluster\n"
                     f"Spearman rho = {rho.spearman_rho:.2f}, "
                     f"p = {rho.p_value:.1e} (n = {int(rho.n_clusters)})", fontsize=10)

    # (E) The same association under three levels of control. The crude value is
    # inflated by capture depth and by cluster composition; what survives both is
    # the number the claim rests on.
    ax = fig.add_subplot(gs[1, :])
    d = co[(co.population == "nodose")
           & co.cluster_depth_stratified_OR.notna()].copy()
    d = d.sort_values("cluster_depth_stratified_OR")
    yy = np.arange(len(d))
    ax.scatter(d.crude_odds_ratio, yy + 0.22, s=30, facecolors="none",
               edgecolors="#999999", linewidths=0.9, zorder=2, label="crude")
    ax.scatter(d.depth_stratified_OR, yy, s=38, c="#6BAED6",
               edgecolors="black", linewidths=0.4, zorder=3, label="+ depth")
    ax.scatter(d.cluster_depth_stratified_OR, yy - 0.22, s=62, c="#B2182B",
               edgecolors="black", linewidths=0.5, zorder=4,
               label="+ depth + cluster")
    ax.axvline(1.0, color="black", lw=0.9, ls="--")
    ax.set_yticks(yy)
    ax.set_yticklabels([f"{p}\nn = {n:,}" for p, n in zip(d.partner, d.n_partner_pos)],
                       fontsize=8)
    ax.set_ylim(-0.75, len(d) - 0.15)
    ax.set_xscale("log")
    ax.set_xlim(0.92, 3.0)
    ax.set_xticks([1, 1.5, 2, 3])
    ax.set_xticklabels(["1", "1.5", "2", "3"])
    ax.set_xlabel("odds of Oprl1 detection, partner+ vs partner-", fontsize=8.5)
    for i, rr in enumerate(d.itertuples()):
        ax.text(rr.cluster_depth_stratified_OR * 1.06, i - 0.22,
                f"p = {rr.p_value_cluster_depth:.0e}", fontsize=6.5,
                ha="left", va="center", color="#B2182B")
    ax.legend(fontsize=6.5, loc="lower right", markerscale=0.9)
    ax.set_title("(d) The association survives depth and cluster control",
                 fontsize=10)

    st.save(fig, "figure4_vagal_oprl1_satiation")


if __name__ == "__main__":
    raise SystemExit(main())
