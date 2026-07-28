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


def figure1():
    """The four opioid receptors in every peripheral population measured.

    Panels are blocked by division of the peripheral nervous system, in the
    order sensory, sympathetic, parasympathetic and enteric, so that the
    comparison between divisions is read down the figure rather than assembled
    from the caption.
    """
    d = pd.read_csv(ac.RES / "receptor_levels_by_ganglion.csv")
    groups = list(dict.fromkeys(d.group))
    block = {g: int(np.ceil((d.group == g).sum() / NCOL)) for g in groups}
    nrow = sum(block.values())

    st.set_theme()
    fig = plt.figure(figsize=(2.9 * NCOL, 4.6 * nrow))
    gs = fig.add_gridspec(nrow, NCOL, hspace=0.85, wspace=0.42,
                          left=0.06, right=0.99, top=0.915, bottom=0.03)

    row0 = 0
    for g in groups:
        sub = d[d.group == g]
        for k, r in enumerate(sub.itertuples()):
            ax = fig.add_subplot(gs[row0 + k // NCOL, k % NCOL])
            vals = [getattr(r, gene) for gene in ac.RECEPTORS]
            order = np.argsort(-np.array(vals))
            genes = [ac.RECEPTORS[i] for i in order]
            st.expression_bars(
                ax, [vals[i] for i in order], genes,
                f"Mean expression ({r.unit})",
                colors=[st.BAR_BLUE if gene == "Oprl1" else st.BAR_GREY
                        for gene in genes],
                fontsize=12)
            st.panel_letter(ax, r.panel, dx=-0.32)
            ax.set_title(f"{r.tissue}\n{r.dataset}, n = {r.n:,}", fontsize=11,
                         pad=8)
            if k == 0:
                # The rule clears the two-line panel title, so the group name
                # never sits beside the first panel's tissue name.
                y = ax.get_position().y1 + 0.031
                fig.text(0.005, y + 0.005, g.capitalize(), fontsize=14,
                         fontweight="bold", ha="left", va="bottom")
                fig.lines.append(plt.Line2D(
                    [0.005, 0.995], [y] * 2, transform=fig.transFigure,
                    color="black", lw=1.0))
        row0 += block[g]

    fig.suptitle("The four opioid receptors across the peripheral nervous system",
                 fontsize=15, y=0.995)
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
    st.save(fig, "figure4_sodium_channel_gradient")


if __name__ == "__main__":
    figure1()
    figure3()
