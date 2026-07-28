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


def figure1():
    d = pd.read_csv(ac.RES / "receptor_levels_by_ganglion.csv")
    st.set_theme()
    fig, axes = plt.subplots(1, len(d), figsize=(3.0 * len(d), 4.4))
    for ax, r in zip(np.atleast_1d(axes), d.itertuples()):
        vals = [getattr(r, g) for g in ac.RECEPTORS]
        order = np.argsort(-np.array(vals))
        genes = [ac.RECEPTORS[i] for i in order]
        st.expression_bars(
            ax, [vals[i] for i in order], genes,
            f"Mean expression ({r.unit})",
            colors=[st.BAR_BLUE if g == "Oprl1" else st.BAR_GREY for g in genes],
            fontsize=12)
        st.panel_letter(ax, r.panel, dx=-0.30)
        ax.set_title(f"{r.tissue}\n{r.dataset}, n = {r.n:,}", fontsize=11, pad=8)
        ax.text(0.5, -0.30, f"{r.margin:.2f}× over runner-up\nsupport {r.support:.2f}",
                transform=ax.transAxes, ha="center", va="top", fontsize=8.5,
                color="#444444")
    fig.suptitle("Oprl1 is the highest-expressed opioid receptor in every ganglion measured",
                 fontsize=13, y=1.04)
    fig.tight_layout()
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
    st.save(fig, "figure3_sodium_channel_gradient")


if __name__ == "__main__":
    figure1()
    figure3()
