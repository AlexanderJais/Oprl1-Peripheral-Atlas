"""Mechanotransduction and Gi effector genes against Oprl1.

Section 5 of the README places Oprl1 on Nav1.1 nodose neurons. If N/OFQ reduces
firing in those cells, three conditions hold on the same neurons, and each is
checkable in this atlas:

  1. the mechanotransducer is present      Piezo2
  2. the Gi effector genes are present     Kcnj3/6/9 (GIRK), Gnai, Gnao,
                                           Cacna1b, which a Gi-coupled receptor
                                           signals through
  3. the nociceptor programme is absent    Trpv1, Trpa1, Scn10a

Condition 3 carries as much weight as the other two. Oprl1 correlating with
everything would make the first two uninformative, so the nociceptor module
serves as the internal negative control and should run in the opposite
direction.

This script computes nothing new from the matrix: it reads the transcriptome-wide
correlation produced by 06_oprl1_localisation.py and the cluster pseudobulk it
cached, then adds one per-cell test for Piezo2.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import atlas_common as ac
import atlas_style as st

_nodose = __import__("02_nodose")
_loc = __import__("06_oprl1_localisation")

MODULES = {
    "Mechanotransduction\n& large-soma identity":
        ["Piezo2", "Ntrk3", "Nefh", "Pvalb", "Adgrg6"],
    "Gi effector genes\n(what a Gi receptor signals through)":
        ["Kcnj3", "Kcnj6", "Kcnj9", "Gnai1", "Gnai2", "Gnao1", "Cacna1b"],
    "Nociceptor programme\n(negative control)":
        ["Trpv1", "Trpa1", "Scn10a", "Scn9a"],
}
MODULE_COLORS = ["#08306B", "#1B7837", "#B2182B"]
MIN_CELLS = 30


def main() -> int:
    ac.RES.mkdir(parents=True, exist_ok=True)
    corr = pd.read_csv(ac.RES / "nodose_oprl1_gene_correlations.csv")
    lookup = corr.set_index("gene")

    rows = []
    for module, genes in MODULES.items():
        for g in genes:
            if g not in lookup.index:
                print(f"  [note] {g} is not expressed above threshold in the "
                      "nodose clusters; no correlation to report")
                continue
            r = lookup.loc[g]
            rows.append({"module": module.replace("\n", " "), "gene": g,
                         "spearman_rho": float(r.spearman_rho),
                         "q_value": float(r.q_value),
                         "mean_CPM": float(r.mean_CPM_across_clusters)})
    mod = pd.DataFrame(rows)
    ac.save_table(mod, "vagal_gene_module_correlations.csv")
    print("\n  Correlation of each module gene with Oprl1 across nodose clusters:")
    for module in mod.module.unique():
        sub = mod[mod.module == module].sort_values("spearman_rho", ascending=False)
        print(f"\n  {module}:")
        print(sub[["gene", "spearman_rho", "q_value", "mean_CPM"]].to_string(index=False))

    # Do the three modules differ? The nociceptor module is the control.
    groups = [mod.loc[mod.module == m, "spearman_rho"].values
              for m in mod.module.unique()]
    if len(groups) == 3 and all(len(g) >= 2 for g in groups):
        h, p = stats.kruskal(*groups)
        print(f"\n  Kruskal-Wallis across the three modules: H = {h:.3f}, p = {p:.4f}")

    # ------------------------------------------------- per-cell Piezo2 test
    adata = ad.read_h5ad(ac.NODOSE_H5AD, backed="r")
    obs = adata.obs
    counts, absent = _nodose.load_counts(adata, ac.NODOSE_H5AD,
                                         ["Oprl1", "Piezo2", "Trpv1"])
    whole_cell = obs["suspension_type"].values == "cell"
    is_nodose = (obs["cell_class"] == "Nodose ganglion neuron").values & whole_cell
    cluster = obs["author_cell_type"].astype(str).values
    nf = obs["nFeature_RNA"].values
    depth = np.full(len(nf), "", dtype=object)
    depth[is_nodose] = ac.depth_strata(nf[is_nodose], 10)
    joint = np.full(len(nf), "", dtype=object)
    joint[is_nodose] = ac.cluster_depth_strata(cluster[is_nodose], nf[is_nodose])

    det = (counts > 0)
    prows = []
    for partner in ("Piezo2", "Trpv1"):
        if partner in absent:
            continue
        sub = det.loc[is_nodose]
        p_pos = sub[partner].values
        o_pos = sub["Oprl1"].values
        or_d, pd_, kd = ac.stratified_odds_ratio(o_pos, p_pos, depth[is_nodose])
        or_c, pc_, kc = ac.stratified_odds_ratio(o_pos, p_pos, joint[is_nodose])
        prows.append({
            "partner": partner, "n_partner_pos": int(p_pos.sum()),
            "pct_of_partner_pos_that_are_Oprl1_pos":
                round(100 * float((p_pos & o_pos).sum()) / max(p_pos.sum(), 1), 2),
            "pct_of_partner_neg_that_are_Oprl1_pos":
                round(100 * float((~p_pos & o_pos).sum()) / max((~p_pos).sum(), 1), 2),
            "depth_stratified_OR": round(or_d, 3), "p_depth": pd_,
            "cluster_depth_stratified_OR": round(or_c, 3), "p_cluster_depth": pc_,
            "n_cluster_depth_strata": kc,
        })
    piezo = pd.DataFrame(prows)
    ac.save_table(piezo, "vagal_piezo2_oprl1_coexpression.csv")
    print("\n  Per-cell co-detection, nodose neurons:")
    print(piezo.to_string(index=False))

    figures(mod, piezo)
    return 0


def figures(mod, piezo):
    st.set_theme()
    pb, clusters, symbols = _cluster_frame()

    fig = plt.figure(figsize=(15.5, 8.6))
    gs = fig.add_gridspec(2, 3, hspace=0.42, wspace=0.34)

    # (a, b) the two scatters that carry the prediction, cluster by cluster
    for k, partner in enumerate(["Piezo2", "Kcnj3"]):
        ax = fig.add_subplot(gs[0, k])
        if partner not in pb.columns:
            continue
        x, y = pb[partner].values, pb["Oprl1"].values
        ax.scatter(x, y, s=46, c=st.BAR_BLUE, edgecolors="black", linewidths=0.6,
                   zorder=3)
        for c, xi, yi in zip(clusters, x, y):
            ax.annotate(c, (xi, yi), fontsize=6, xytext=(4, 3),
                        textcoords="offset points")
        rho, p = stats.spearmanr(x, y)
        ax.set_xlabel(f"{partner} (mean CPM per cluster)", fontsize=10)
        ax.set_ylabel("Oprl1 (mean CPM per cluster)", fontsize=10)
        st.panel_letter(ax, "ab"[k], dx=-0.22)
        ax.set_title(f"rho = {rho:.2f}, p = {p:.1e}", fontsize=11)

    # (c) Trpv1, the negative control, on the same axes
    ax = fig.add_subplot(gs[0, 2])
    if "Trpv1" in pb.columns:
        x, y = pb["Trpv1"].values, pb["Oprl1"].values
        ax.scatter(x, y, s=46, c=st.HIGHLIGHT, edgecolors="black", linewidths=0.6,
                   zorder=3)
        for c, xi, yi in zip(clusters, x, y):
            ax.annotate(c, (xi, yi), fontsize=6, xytext=(4, 3),
                        textcoords="offset points")
        rho, p = stats.spearmanr(x, y)
        ax.set_xlabel("Trpv1 (mean CPM per cluster)", fontsize=10)
        ax.set_ylabel("Oprl1 (mean CPM per cluster)", fontsize=10)
        ax.set_title(f"rho = {rho:.2f}, p = {p:.1e}", fontsize=11)
    st.panel_letter(ax, "c", dx=-0.22)

    # (d) the three modules side by side, one point per gene
    ax = fig.add_subplot(gs[1, :2])
    modules = list(mod.module.unique())
    rng = np.random.default_rng(0)
    for i, (m, colour) in enumerate(zip(modules, MODULE_COLORS)):
        sub = mod[mod.module == m]
        ax.scatter(i + rng.uniform(-0.14, 0.14, len(sub)), sub.spearman_rho,
                   s=70, c=colour, edgecolors="black", linewidths=0.6, zorder=3)
        # Push labels apart where the points nearly coincide.
        ys = sub.sort_values("spearman_rho", ascending=False)
        placed = []
        for _, r in ys.iterrows():
            y = r.spearman_rho
            while any(abs(y - q) < 0.075 for q in placed):
                y -= 0.075
            placed.append(y)
            ax.annotate(r.gene, (i + 0.20, y), fontsize=7.5,
                        fontstyle="italic", va="center")
        ax.hlines(sub.spearman_rho.median(), i - 0.3, i + 0.3, color="black",
                  lw=2.4, zorder=4)
    ax.axhline(0, color="black", lw=0.9, ls="--")
    ax.set_xticks(range(len(modules)))
    ax.set_xticklabels([m.replace(" (", "\n(").replace(" & ", "\n& ")
                        for m in modules], fontsize=9)
    ax.set_xlim(-0.5, len(modules) - 0.3)
    ax.set_ylim(-1.05, 1.05)
    ax.set_ylabel("Spearman rho with Oprl1\nacross the 21 nodose clusters",
                  fontsize=10)
    st.panel_letter(ax, "d", dx=-0.10)
    ax.set_title("Oprl1 tracks the GIRK channels and runs against the nociceptor "
                 "programme\nand the G-protein subunits",
                 fontsize=11)

    # (e) the per-cell Piezo2 result
    ax = fig.add_subplot(gs[1, 2])
    if len(piezo):
        labels, vals = [], []
        for _, r in piezo.iterrows():
            labels += [f"{r.partner}+", f"{r.partner}−"]
            vals += [r.pct_of_partner_pos_that_are_Oprl1_pos,
                     r.pct_of_partner_neg_that_are_Oprl1_pos]
        colours = [st.BAR_BLUE if i % 2 == 0 else st.BAR_GREY
                   for i in range(len(vals))]
        st.expression_bars(ax, vals, labels, "% of neurons expressing Oprl1",
                           colors=colours, italic=False, rotation=0, fontsize=10)
        for i, (_, r) in enumerate(piezo.iterrows()):
            ax.text(i * 2 + 0.5, max(vals) * 1.04,
                    f"OR {r.cluster_depth_stratified_OR:.2f}\n"
                    f"p = {r.p_cluster_depth:.0e}", ha="center", fontsize=7.5)
    st.panel_letter(ax, "e", dx=-0.28)
    ax.set_title("Per cell, holding cluster\nand depth fixed", fontsize=11)

    st.save(fig, "figure7_transduction_effector_genes")


def _cluster_frame():
    """Rebuild the cluster x gene pseudobulk frame from the cache 06 wrote."""
    adata = ad.read_h5ad(ac.NODOSE_H5AD, backed="r")
    obs = adata.obs
    symbols = adata.raw.var["feature_name"].astype(str).values
    whole_cell = obs["suspension_type"].values == "cell"
    is_nodose = (obs["cell_class"] == "Nodose ganglion neuron").values & whole_cell
    cluster = obs["author_cell_type"].astype(str).values
    ngn = sorted({c for c in cluster[is_nodose] if c.startswith("NGN")},
                 key=lambda s: int(s[3:]))
    keep = [c for c in ngn if (is_nodose & (cluster == c)).sum() >= MIN_CELLS]
    z = np.load(_loc.CACHE)
    pb = pd.DataFrame(z["mat"], index=keep, columns=symbols)
    return pb.T.groupby(level=0).sum().T, keep, symbols


if __name__ == "__main__":
    raise SystemExit(main())
