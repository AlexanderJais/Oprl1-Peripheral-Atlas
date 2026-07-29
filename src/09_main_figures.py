"""Consolidated main figures.

Figure 1 is the claim: the four opioid receptors in every ganglion measured, one
panel per ganglion, each in its own unit. Figure 3 is the sodium-channel
gradient in the two ganglia where it was tested.

Both read tables the numbered pipeline wrote, so nothing is recomputed here.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import atlas_common as ac
import atlas_style as st


NCOL = 4

# Panels carry identifiers only. Every statement about what the figure shows,
# including the margins, the bootstrap support and the preparation of each
# dataset, belongs to the legend in paper/figure_legends.md.


def figure1():
    """The four opioid receptors in each peripheral neuronal population.

    Built to the Cell Press two-column width so nothing is rescaled at
    submission. Panels are blocked by division of the peripheral nervous
    system; the division names are the only text in the figure that is not a
    tissue name, a sample size, a gene symbol or an axis value.
    """
    d = pd.read_csv(ac.RES / "receptor_levels_by_ganglion.csv")
    groups = list(dict.fromkeys(d.group))
    block = {g: int(np.ceil((d.group == g).sum() / NCOL)) for g in groups}

    # A thin empty row between blocks. Without it the division rule lands on the
    # rotated gene labels of the row above, which hang below their axes.
    SPACER = 0.36
    heights, start = [], {}
    for i, g in enumerate(groups):
        if i:
            heights.append(SPACER)
        start[g] = len(heights)
        heights += [1.0] * block[g]

    st.set_theme()
    fig = plt.figure(figsize=(st.W_2COL, 1.50 * sum(heights) + 0.30))
    gs = fig.add_gridspec(len(heights), NCOL, height_ratios=heights,
                          hspace=0.95, wspace=0.62,
                          left=0.062, right=0.995, top=0.955, bottom=0.035)

    for g in groups:
        sub = d[d.group == g]
        for k, r in enumerate(sub.itertuples()):
            ax = fig.add_subplot(gs[start[g] + k // NCOL, k % NCOL])
            vals = [getattr(r, gene) for gene in ac.RECEPTORS]
            order = np.argsort(-np.array(vals))
            genes = [ac.RECEPTORS[i] for i in order]
            st.expression_bars(
                ax, [vals[i] for i in order], genes, r.unit,
                colors=[st.BAR_BLUE if gene == "Oprl1" else st.BAR_GREY
                        for gene in genes],
                fontsize=st.FS_TICK, linewidth=0.5, rotation=45)
            st.panel_letter(ax, r.panel, dx=-0.40, dy=1.30)
            tissue = r.tissue[0].upper() + r.tissue[1:]
            ax.set_title(f"{tissue}\nn = {r.n:,}", fontsize=st.FS_NOTE, pad=3)
            if k == 0:
                y = ax.get_position().y1 + 0.031
                fig.text(0.004, y + 0.003, g[0].upper() + g[1:],
                         fontsize=st.FS_LABEL, fontweight="bold", ha="left",
                         va="bottom")
                fig.lines.append(plt.Line2D(
                    [0.004, 0.996], [y] * 2, transform=fig.transFigure,
                    color="black", lw=0.5))

    st.save(fig, "figure1_oprl1_across_ganglia")


def figure3():
    ann = pd.read_csv(ac.RES / "nodose_oprl1_by_cluster_annotated.csv")
    drg = pd.read_csv(ac.RES / "drg_oprl1_by_subtype.csv")
    st.set_theme()
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.8),
                             gridspec_kw={"width_ratios": [0.9, 1.0, 1.2]})

    ax = axes[0]
    nav = [n for n in ["Nav1.1", "Nav1.1/Nav1.8", "Nav1.8"]
           if n in set(ann.sodium_channel_type)]
    st.strip_by_group(ax, nav,
                      {n: ann.loc[ann.sodium_channel_type == n, "Oprl1_CPM"]
                       for n in nav}, "Oprl1 (mean CPM per cluster)")
    st.panel_letter(ax, "a", dx=-0.26)
    ax.set_title("Nodose, 21 clusters\np = 0.00095", fontsize=11, pad=8)

    ax = axes[1]
    x, y = drg.Scn1a.values, drg.Oprl1.values
    ax.scatter(x, y, s=70, c=[st.HIGHLIGHT if s == "NF2" else st.BAR_BLUE
                              for s in drg.subtype],
               edgecolors="black", linewidths=0.7, zorder=3)
    for _, r in drg.iterrows():
        ax.annotate(r.subtype, (r.Scn1a, r.Oprl1), fontsize=8, xytext=(5, 3),
                    textcoords="offset points")
    rho, p = stats.spearmanr(x, y)
    ax.set_xlabel("Scn1a (mean CPM per subtype)", fontsize=11)
    ax.set_ylabel("Oprl1 (mean CPM per subtype)", fontsize=11)
    st.panel_letter(ax, "b", dx=-0.24)
    ax.set_title(f"Dorsal root ganglion, 9 subtypes\nrho = {rho:.2f}, p = {p:.4f}",
                 fontsize=11, pad=8)

    ax = axes[2]
    s = drg.sort_values("Oprl1", ascending=False)
    st.expression_bars(ax, s.Oprl1.values, s.subtype.values, "Oprl1 (mean CPM)",
                       colors=[st.HIGHLIGHT if t == "NF2" else st.BAR_BLUE
                               for t in s.subtype],
                       italic=False, rotation=90, fontsize=11)
    st.panel_letter(ax, "c", dx=-0.14)
    ax.set_title("NF2, the proprioceptors, rank first for\nOprl1, Pvalb, Runx3 and Scn1a",
                 fontsize=11, pad=8)

    fig.tight_layout()
    st.save(fig, "figure5_sodium_channel_gradient")


if __name__ == "__main__":
    figure1()
    figure3()
