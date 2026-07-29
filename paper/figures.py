"""Manuscript figures, built to Cell Press specification.

Everything under paper/figures/ is a figure that appears in the paper, named as
the journal will number it. The figures under figures/ in the repository root
are exploratory output from the analysis pipeline and are not manuscript
figures.

Run from the repository root: python3 paper/figures.py
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import atlas_common as ac
import atlas_style as st

OUT = Path(__file__).resolve().parent / "figures"
NCOL = 4

# Preparation is hatched, not coloured, so colour stays free to carry the gene.
NUC_HATCH = "////"

# Panels carry identifiers only: panel letter, tissue, sample size, gene symbol,
# axis value, unit, and the division name. Margins, bootstrap support, dataset
# accessions and preparation live in paper/figure_legends.md.


def figure1():
    """The four opioid receptor genes in 19 peripheral neuronal populations."""
    d = pd.read_csv(ac.RES / "receptor_levels_by_ganglion.csv")
    groups = list(dict.fromkeys(d.group))
    rows = {g: int(np.ceil((d.group == g).sum() / NCOL)) for g in groups}

    # Laid out in inches rather than by gridspec spacing. A spacer row cannot
    # give a small gap between blocks, because hspace applies on both sides of
    # it, so a block gap is never less than twice a row gap however thin the
    # spacer is. Each block gets its own gridspec at an explicit position.
    PANEL_H = 0.95     # axes height
    ABOVE = 0.42       # two-line tissue title and the panel letter above it
    BELOW = 0.30       # gene labels at 45 degrees
    WITHIN = ABOVE + BELOW
    HEAD = ABOVE       # top of the axes to the rule
    BETWEEN = BELOW + 0.22 + ABOVE      # the rule and the division name between
    TOP_PAD = ABOVE + 0.16
    BOT_PAD = BELOW + 0.06

    n_rows = sum(rows.values())
    height = (TOP_PAD + BOT_PAD + n_rows * PANEL_H
              + sum(r - 1 for r in rows.values()) * WITHIN
              + (len(groups) - 1) * BETWEEN)

    st.set_theme()
    fig = plt.figure(figsize=(st.W_2COL, height))
    left, right = 0.075, 0.995
    y = height - TOP_PAD

    for gi, g in enumerate(groups):
        nr = rows[g]
        block_h = nr * PANEL_H + (nr - 1) * WITHIN
        gs = fig.add_gridspec(nr, NCOL, left=left, right=right,
                              top=y / height, bottom=(y - block_h) / height,
                              hspace=WITHIN / PANEL_H, wspace=0.52)
        fig.text(0.004, (y + HEAD + 0.02) / height, g[0].upper() + g[1:],
                 fontsize=st.FS_LABEL, fontweight="bold", ha="left", va="bottom")
        fig.lines.append(plt.Line2D([0.004, 0.996], [(y + HEAD) / height] * 2,
                                    transform=fig.transFigure, color="black",
                                    lw=0.5))
        for k, r in enumerate(d[d.group == g].itertuples()):
            ax = fig.add_subplot(gs[k // NCOL, k % NCOL])
            vals = [getattr(r, gene) for gene in ac.RECEPTORS]
            order = np.argsort(-np.array(vals))
            genes = [ac.RECEPTORS[i] for i in order]
            st.expression_bars(
                ax, [vals[i] for i in order], genes, r.unit,
                colors=[st.BAR_BLUE if gene == "Oprl1" else st.BAR_GREY
                        for gene in genes],
                fontsize=st.FS_TICK, linewidth=0.5, rotation=45)
            st.panel_letter(ax, r.panel, dx=-0.40, dy=1.34)
            tissue = r.tissue[0].upper() + r.tissue[1:]
            ax.set_title(f"{tissue}\nn = {r.n:,}", fontsize=st.FS_NOTE, pad=3)
        y -= block_h + BETWEEN
    return fig


def _it(*symbols):
    """Gene symbols italicised inside an axis label."""
    return " / ".join(rf"$\it{{{s}}}$" for s in symbols)


def _prep_key(fig, x, y):
    """The whole-cell against nuclear key, at a figure position.

    Drawn on the canvas rather than inside a panel. Inside panel A it sat on
    top of the nuclear Oprm1 bar, which is the tallest thing in the figure, and
    no in-axes anchor clears it.
    """
    solid = plt.Rectangle((0, 0), 1, 1, facecolor="white", edgecolor="black",
                          linewidth=0.5)
    hatched = plt.Rectangle((0, 0), 1, 1, facecolor="white", edgecolor="black",
                            linewidth=0.5, hatch=NUC_HATCH)
    fig.legend([solid, hatched], ["whole cell", "nuclear"], loc="upper left",
               bbox_to_anchor=(x, y), ncol=2, fontsize=st.FS_NOTE,
               handlelength=1.3, handleheight=1.0, borderpad=0.0,
               columnspacing=1.1, handletextpad=0.4)


def _paired_receptor_bars(ax, whole, nuclear, ylim):
    """Four receptors measured in one tissue under both preparations."""
    x = np.arange(len(ac.RECEPTORS))
    colors = [st.BAR_BLUE if g == "Oprl1" else st.BAR_GREY for g in ac.RECEPTORS]
    ax.bar(x - 0.20, whole, 0.38, color=colors, edgecolor="black",
           linewidth=0.5, zorder=3)
    ax.bar(x + 0.20, nuclear, 0.38, color=colors, edgecolor="black",
           linewidth=0.5, hatch=NUC_HATCH, zorder=3)
    # Logarithmic: the whole-cell ordering these panels are compared against is
    # a 1.1-fold difference, and a linear axis that fits the 244 CPM nuclear
    # Oprm1 bar renders it thinner than the line weight.
    ax.set_yscale("log")
    ax.set_ylim(*ylim)
    ax.set_xticks(x)
    ax.set_xticklabels(ac.RECEPTORS, rotation=45, ha="right",
                       fontstyle="italic", fontsize=st.FS_TICK)
    ax.set_ylabel("mean expression (CPM)", fontsize=st.FS_LABEL)


def figureS1():
    """What a nuclear preparation does to the receptor ordering."""
    nod = pd.read_csv(ac.RES / "nodose_opioid_levels.csv")
    bias = pd.read_csv(ac.RES / "nuclear_bias_vs_gene_length.csv").set_index("gene")
    gang = pd.read_csv(ac.RES / "receptor_levels_by_ganglion.csv")
    xsp = pd.read_csv(ac.RES / "human_xspecies_drg.csv").set_index("species")

    gw = pd.read_csv(ac.RES / "preparation_bias_genomewide.csv")
    dec = pd.read_csv(ac.RES / "preparation_bias_by_span_decile.csv")
    rec = pd.read_csv(ac.RES / "preparation_bias_receptors.csv").set_index("gene")

    # The four measured panels across the top, the two genome-wide panels below
    # at half the count and so twice the width. Each row gets its own gridspec:
    # the rows differ in column count and only the top row is titled.
    PANEL_H, BELOW = 1.60, 0.38
    TITLED, BARE = 0.36, 0.42
    KEY = 0.20                      # the preparation key, above every panel
    TOP_PAD, BOT_PAD = TITLED + KEY + 0.04, BELOW + 0.06
    height = TOP_PAD + 2 * PANEL_H + BELOW + BARE + BOT_PAD

    st.set_theme()
    plt.rcParams["hatch.linewidth"] = 0.4
    fig = plt.figure(figsize=(st.W_SUPP, height))
    left, right = 0.075, 0.99
    axes = []
    y = height - TOP_PAD
    for ncol, wspace in ((4, 0.62), (2, 0.34)):
        gs = fig.add_gridspec(1, ncol, left=left, right=right, top=y / height,
                              bottom=(y - PANEL_H) / height, wspace=wspace)
        axes += [fig.add_subplot(gs[0, i]) for i in range(ncol)]
        y -= PANEL_H + BELOW + BARE
    for ax, letter in zip(axes, "ABCDEF"):
        st.panel_letter(ax, letter, dx=-0.34 if letter in "ABCD" else -0.17,
                        dy=1.26)

    # (A) One tissue, both preparations, from the same atlas.
    ax = axes[0]
    _paired_receptor_bars(ax, [bias.loc[g, "whole_cell_CPM"] for g in ac.RECEPTORS],
                          [bias.loc[g, "nuclear_CPM"] for g in ac.RECEPTORS],
                          (0.1, 1000))
    wc_n = int(nod[(nod.prep == "whole cell") & (nod.gene == "Oprl1")
                   & (nod.tissue == "nodose+jugular")].n_cells.sum())
    nu_n = int(nod[(nod.prep == "nuclear") & (nod.gene == "Oprl1")].n_cells.iloc[0])
    ax.set_title(f"Vagal ganglia (X)\n{wc_n:,} cells, {nu_n:,} nuclei",
                 fontsize=st.FS_NOTE, pad=3)
    _prep_key(fig, left, 1 - 0.02)

    # (B) A second tissue, a second laboratory, the same inversion.
    ax = axes[1]
    drg = gang[gang.tissue == "dorsal root"].iloc[0]
    _paired_receptor_bars(ax, [getattr(drg, g) for g in ac.RECEPTORS],
                          [xsp.loc["mouse", g.upper()] for g in ac.RECEPTORS],
                          (0.1, 1000))
    ax.set_title(f"Dorsal root ganglion\n{drg.n:,} cells, "
                 f"{int(round(xsp.loc['mouse', 'n'] * xsp.loc['mouse', 'n_samples'])):,}"
                 " nuclei", fontsize=st.FS_NOTE, pad=3)

    # (C and D) Detection rate per gene per deposit, on one scale. Plotted as
    # the two rates rather than their ratio: a ratio hides which of the two
    # moved, and here only one of them does.
    det = (nod[nod.dataset.str.contains(":") & (nod.tissue == "nodose+jugular")]
           .pivot_table(index="dataset", columns="gene", values="pct_detected"))
    prep = (nod[nod.dataset.str.contains(":")].groupby("dataset")["prep"].first())
    det = det.assign(prep=prep)
    det = pd.concat([det[det.prep == "whole cell"].sort_index(),
                     det[det.prep == "nuclear"]])
    x = np.arange(len(det))
    hatch = ["" if p == "whole cell" else NUC_HATCH for p in det.prep]
    # The atlas calls its own deposit "in-house". Named here for its authors, as
    # the other four are, since in-house reads as ours on this page.
    names = [d.split(":")[1].replace("inhouse", "Cheng") for d in det.index]
    top = float(np.ceil(det[["Oprl1", "Oprm1"]].to_numpy().max() / 10) * 10)
    for ax, gene, color in ((axes[2], "Oprl1", st.BAR_BLUE),
                            (axes[3], "Oprm1", st.BAR_GREY)):
        ax.bar(x, det[gene], 0.68, color=color, edgecolor="black",
               linewidth=0.5, zorder=3, hatch=hatch)
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=45, ha="right", fontsize=st.FS_TICK)
        ax.set_ylim(0, top)
        ax.set_ylabel(rf"$\it{{{gene}}}$ detection (%)", fontsize=st.FS_LABEL)

    # (E) Which of the two measurements moves with gene length. Whichever it is
    # is the distorted one: the length of a transcription unit says nothing
    # about how much mature message a neuron carries.
    ax = axes[4]
    for col, label, marker, fill in (("median_whole_cell", "whole cell", "o", "white"),
                                     ("median_nuclear", "nuclear", "s", "black")):
        ax.plot(dec.span_med, dec[col], linestyle="none", marker=marker,
                color="black", markersize=3.0, markerfacecolor=fill,
                markeredgewidth=0.5, label=label, zorder=3)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(5, 200)
    ax.set_xlabel("genomic span (kb), decile median", fontsize=st.FS_LABEL)
    ax.set_ylabel("median expression (CPM)", fontsize=st.FS_LABEL)
    ax.legend(loc="upper left", fontsize=st.FS_NOTE, handlelength=1.8,
              borderpad=0.2, labelspacing=0.3, handletextpad=0.5)

    # (F) Every expressed gene, and where the receptors sit among them.
    ax = axes[5]
    ax.scatter(gw.span_kb, gw.ratio, s=0.6, c="#C8C8C8", linewidths=0,
               rasterized=True, zorder=2)
    ax.plot(dec.span_med, dec.median_ratio, "-", color="black", linewidth=0.8,
            zorder=4)
    ax.axhline(1.0, color="black", lw=0.5, ls=(0, (3, 2)), zorder=3)
    offsets = {"Oprl1": (-4, -11), "Oprm1": (-21, 5), "Oprk1": (4, 3),
               "Oprd1": (-22, 3)}
    for gene in rec.index:
        colour = st.BAR_BLUE if gene == "Oprl1" else "black"
        ax.scatter([rec.loc[gene, "span_kb"]], [rec.loc[gene, "ratio"]], s=14,
                   c=colour, edgecolors="black", linewidths=0.4, zorder=5)
        ax.annotate(gene, (rec.loc[gene, "span_kb"], rec.loc[gene, "ratio"]),
                    fontsize=st.FS_TICK, fontstyle="italic", zorder=5,
                    xytext=offsets.get(gene, (4, 3)), textcoords="offset points")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.15, 3000)
    ax.set_ylim(0.02, 200)
    ax.set_xlabel("genomic span (kb)", fontsize=st.FS_LABEL)
    ax.set_ylabel("nuclear / whole cell", fontsize=st.FS_LABEL)
    return fig


def emit(fig, name, max_w=174.5, max_h=235.0):
    """Write the figure and refuse it if it breaks the journal's limits."""
    OUT.mkdir(parents=True, exist_ok=True)
    # A tight bounding box would crop to the drawn content and hand the journal
    # a figure a millimetre or two off the column width. The canvas is already
    # built to the column width, so it is written as it stands.
    fig.savefig(OUT / f"{name}.pdf", bbox_inches=fig.bbox_inches)
    fig.savefig(OUT / f"{name}.png", bbox_inches=fig.bbox_inches)

    w, h = (v * 25.4 for v in fig.get_size_inches())
    tight = fig.get_tightbbox(fig.canvas.get_renderer())
    plt.close(fig)

    over = []
    if w > max_w:
        over.append(f"width {w:.1f} mm exceeds {max_w:.0f} mm")
    if h > max_h:
        over.append(f"height {h:.1f} mm exceeds {max_h:.0f} mm")
    # Nothing may sit outside the canvas, or writing it untrimmed clips it.
    for edge, past in (("left", -tight.x0), ("bottom", -tight.y0),
                       ("right", tight.x1 - w / 25.4),
                       ("top", tight.y1 - h / 25.4)):
        if past > 0.005:
            over.append(f"clipped {past * 25.4:.2f} mm at the {edge}")
    if over:
        raise ac.SanityCheckError(f"{name}: " + "; ".join(over))
    print(f"  [paper] paper/figures/{name}.pdf|.png  {w:.1f} x {h:.1f} mm")


if __name__ == "__main__":
    emit(figure1(), "Figure1")
    emit(figureS1(), "FigureS1",
         max_w=st.W_SUPP * 25.4 + 0.5, max_h=st.H_SUPP * 25.4)
