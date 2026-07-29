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

import re

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
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


def figure2():
    """What a nuclear preparation does to the receptor ordering."""
    nod = pd.read_csv(ac.RES / "nodose_opioid_levels.csv")
    bias = pd.read_csv(ac.RES / "nuclear_bias_vs_gene_length.csv").set_index("gene")
    gang = pd.read_csv(ac.RES / "receptor_levels_by_ganglion.csv")
    xsp = pd.read_csv(ac.RES / "human_xspecies_drg.csv").set_index("species")

    gw = pd.read_csv(ac.RES / "preparation_bias_genomewide.csv")
    dec = pd.read_csv(ac.RES / "preparation_bias_by_span_decile.csv")
    rec = pd.read_csv(ac.RES / "preparation_bias_receptors.csv").set_index("gene")

    # Three rows of two. As a supplemental figure this was four panels across
    # the top and two below, on a 228 mm landscape canvas that a reviewer reads
    # on one page; a main figure is bound by the 174 mm text column, and four
    # panels across that width leaves each one under 25 mm of drawing area with
    # five rotated deposit names beneath it. Two per row is what the column
    # affords. It also falls out by pairs: the two tissues measured under both
    # preparations, the two detection rates, the two genome-wide panels. Each
    # row gets its own gridspec, since only the top row carries titles.
    PANEL_H, BELOW = 1.60, 0.38
    TITLED, BARE = 0.36, 0.42
    KEY = 0.20                      # the preparation key, above every panel
    TOP_PAD, BOT_PAD = TITLED + KEY + 0.04, BELOW + 0.06
    NROW = 3
    height = (TOP_PAD + NROW * PANEL_H + (NROW - 1) * (BELOW + BARE) + BOT_PAD)

    st.set_theme()
    plt.rcParams["hatch.linewidth"] = 0.4
    fig = plt.figure(figsize=(st.W_2COL, height))
    left, right = 0.085, 0.99
    axes = []
    y = height - TOP_PAD
    for _ in range(NROW):
        gs = fig.add_gridspec(1, 2, left=left, right=right, top=y / height,
                              bottom=(y - PANEL_H) / height, wspace=0.30)
        axes += [fig.add_subplot(gs[0, i]) for i in range(2)]
        y -= PANEL_H + BELOW + BARE
    for ax, letter in zip(axes, "ABCDEF"):
        st.panel_letter(ax, letter, dx=-0.155, dy=1.26)

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
    for col, label, marker, ls, fill in (
            ("median_whole_cell", "whole cell", "o", "-", "white"),
            ("median_nuclear", "nuclear", "s", "--", "black")):
        ax.plot(dec.span_med, dec[col], ls, marker=marker, color="black",
                markersize=2.8, linewidth=0.7, markerfacecolor=fill,
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


# The eight purified dorsal root ganglion subtypes of GSE131230, in the
# deposit's own order: unmyelinated nociceptors and C-LTMRs, then the
# myelinated low-threshold afferents, then the proprioceptor.
SUBTYPE_ORDER = ["Nonpeptidergic Nociceptor", "Peptidergic Nociceptor", "C-LTMR",
                 "Aδ-LTMR", "Aβ RA-LTMR", "Aβ SA1-LTMR", "Aβ Field-LTMR",
                 "Proprioceptor"]


def figureS3():
    """The four receptors across purified dorsal root ganglion subtypes."""
    scr = pd.read_csv(ac.RES / "bulk_composition_screen.csv")
    sub = scr[(scr.gse == "GSE131230") & scr.group.isin(SUBTYPE_ORDER)]
    sub = sub.set_index("group").reindex([s for s in SUBTYPE_ORDER
                                          if s in set(sub.group)])

    PANEL_H, ABOVE, BELOW = 1.70, 0.30, 1.15
    height = ABOVE + PANEL_H + BELOW
    st.set_theme()
    fig = plt.figure(figsize=(st.W_15COL, height))
    gs = fig.add_gridspec(1, 1, left=0.135, right=0.99,
                          top=(height - ABOVE) / height, bottom=BELOW / height)
    ax = fig.add_subplot(gs[0, 0])
    x = np.arange(len(sub))
    for i, gene in enumerate(ac.RECEPTORS):
        ax.bar(x + (i - 1.5) * 0.20, sub[gene], 0.20, color=RECEPTOR_FILL[gene],
               edgecolor="black", linewidth=0.4, zorder=3,
               label=rf"$\it{{{gene}}}$")
    ax.set_yscale("log")
    ax.yaxis.set_minor_formatter(mticker.NullFormatter())
    ax.set_xticks(x)
    ax.set_xticklabels(sub.index, rotation=45, ha="right", fontsize=st.FS_TICK)
    ax.set_xlim(-0.6, len(sub) - 0.4)
    ax.set_ylabel("mean expression (CPM)", fontsize=st.FS_LABEL)
    ax.legend(loc="upper left", fontsize=st.FS_NOTE, ncol=4, frameon=False,
              handlelength=1.1, columnspacing=0.9, handletextpad=0.4)
    return fig


def figure3():
    """Bulk ganglion tissue: composition, nerve injury, and neuronal subtype."""
    scr = pd.read_csv(ac.RES / "bulk_composition_screen.csv")
    axo = pd.read_csv(ac.RES / "bulk_axotomy_pairs.csv")

    # The screen gates on neuronal content, which a brain sample also clears,
    # so tissue identity is applied here: the screen's job is to record every
    # group with its reason, and this figure is about peripheral ganglia.
    NOT_PERIPHERAL = re.compile(r"cortex|hippocamp|brain|spinal|striat|"
                                r"arcuate|\bArc|hypothal", re.I)
    # Series whose Atf3 is raised by an experimental nerve lesion rather than by
    # handling; they carry (A) and (B) and are excluded from (C).
    CHRONIC_LESION_SERIES = {"GSE97090", "GSE149770", "GSE161342", "GSE138769",
                             "GSE67130", "GSE188922"}
    read = scr[(scr.verdict == "read")
               & ~scr.group.str.contains(NOT_PERIPHERAL)].copy()
    # A deposit that publishes FPKM and counts of the same samples appears
    # twice; the length-normalised copy is the one this figure plots.
    read["set"] = read.gse + "|" + read.group.str.replace(
        r"^(FPKM|count|TPM)\.", "", regex=True)
    read = read.sort_values("unit").drop_duplicates("set")

    # The bottom row carries subtype names at 45 degrees, the longest of which
    # is "Nonpeptidergic Nociceptor", so the two rows need different clearance
    # beneath them.
    PANEL_H, ABOVE = 1.50, 0.34
    BELOW_TOP, BELOW_BOTTOM = 0.34, 0.46
    height = ABOVE + PANEL_H + (BELOW_TOP + ABOVE) + PANEL_H + BELOW_BOTTOM

    st.set_theme()
    fig = plt.figure(figsize=(st.W_2COL, height))
    left, right = 0.10, 0.99
    row1 = height - ABOVE
    gs1 = fig.add_gridspec(1, 2, left=left, right=right, top=row1 / height,
                           bottom=(row1 - PANEL_H) / height,
                           wspace=0.34, width_ratios=[1, 1.25])
    axA, axB = fig.add_subplot(gs1[0, 0]), fig.add_subplot(gs1[0, 1])
    row2 = row1 - PANEL_H - BELOW_TOP - ABOVE
    gs2 = fig.add_gridspec(1, 2, left=left, right=right, top=row2 / height,
                           bottom=(row2 - PANEL_H) / height, wspace=0.34)
    axC, axD = fig.add_subplot(gs2[0, 0]), fig.add_subplot(gs2[0, 1])

    # (A) and (B): every quantity is a level, so the lesion and the receptors
    # are read on the same axis rather than as a change in a ratio.
    def paired(ax, genes, colours):
        x = [0, 1]
        for gene, colour in zip(genes, colours):
            g = axo[axo.gene == gene]
            for r in g.itertuples():
                ax.plot(x, [r.control_CPM, r.injury_CPM], "-", color=colour,
                        lw=0.7, alpha=0.85, zorder=3)
            ax.plot([0] * len(g), g.control_CPM, "o", color=colour, ms=3.4,
                    mec="black", mew=0.4, zorder=4, label=rf"$\it{{{gene}}}$")
            ax.plot([1] * len(g), g.injury_CPM, "o", color=colour, ms=3.4,
                    mec="black", mew=0.4, zorder=4)
        ax.set_yscale("log")
        ax.set_xlim(-0.35, 1.35)
        ax.set_xticks(x)
        ax.set_xticklabels(["control", "nerve injury"], fontsize=st.FS_TICK)
        ax.set_ylabel("mean expression (CPM)", fontsize=st.FS_LABEL)

    paired(axA, ["Atf3"], [st.HIGHLIGHT])
    axA.legend(loc="upper left", fontsize=st.FS_NOTE, frameon=False,
               handletextpad=0.3, borderpad=0.1)
    paired(axB, ["Oprl1", "Oprm1"], [st.BAR_BLUE, RECEPTOR_FILL["Oprm1"]])
    axB.legend(loc="upper right", fontsize=st.FS_NOTE, frameon=False, ncol=2,
               handletextpad=0.3, borderpad=0.1, columnspacing=1.0)

    # (C) The acute question, which (A) and (B) do not answer: those are
    # chronic lesions over days, while dissociation is an hour of enzyme. Here
    # Atf3 is read as a continuous measure of how hard a preparation was
    # handled, in groups carrying no experimental lesion.
    lesion = scr.gse.isin(CHRONIC_LESION_SERIES)
    q = scr[(scr.Atf3 > 0) & (scr.Oprl1 > 0) & (scr.Snap25 > 0) & ~lesion
            & ~scr.group.str.contains(NOT_PERIPHERAL)].copy()
    for mask, marker, face, label in (
            (q.Snap25_over_Plp1 < 5, "o", st.BAR_BLUE, "ganglion tissue"),
            (q.Snap25_over_Plp1 >= 5, "^", "white", "purified neurons")):
        axC.plot(q.loc[mask, "Atf3"], q.loc[mask, "Oprl1"], marker, ms=3.6,
                 mfc=face, mec="black", mew=0.4, ls="none", zorder=3, label=label)
    axC.set_xscale("log"); axC.set_yscale("log")
    axC.set_xlabel(r"$\it{Atf3}$ (CPM)", fontsize=st.FS_LABEL)
    axC.set_ylabel(r"$\it{Oprl1}$ (CPM)", fontsize=st.FS_LABEL)
    axC.legend(loc="lower left", fontsize=st.FS_NOTE, frameon=False,
               handletextpad=0.3, borderpad=0.1)

    # (D) The two receptors against each other, so neither axis is a ratio and
    # the diagonal carries the comparison.
    tissue = read.Snap25_over_Plp1 < 5
    for mask, marker, face, label in (
            (tissue, "o", st.BAR_BLUE, "ganglion tissue"),
            (~tissue, "^", "white", "purified neurons")):
        axD.plot(read.loc[mask, "Oprm1"].clip(lower=0.01),
                 read.loc[mask, "Oprl1"].clip(lower=0.01), marker,
                 ms=3.6, mfc=face, mec="black", mew=0.4, ls="none",
                 zorder=3, label=label)
    lim = (0.004, 1000)
    axD.plot(lim, lim, "--", color="black", lw=0.5, zorder=2)
    axD.set_xscale("log"); axD.set_yscale("log")
    axD.set_xlim(*lim); axD.set_ylim(*lim)
    axD.set_xlabel(r"$\it{Oprm1}$ (CPM)", fontsize=st.FS_LABEL)
    axD.set_ylabel(r"$\it{Oprl1}$ (CPM)", fontsize=st.FS_LABEL)
    axD.legend(loc="lower right", fontsize=st.FS_NOTE, frameon=False,
               handletextpad=0.3, borderpad=0.1)

    # On a log axis spanning less than two decades matplotlib labels the minor
    # ticks as well, which fills the axis with 3 x 10^2 and its neighbours.
    for ax in (axA, axB, axC, axD):
        ax.yaxis.set_minor_formatter(mticker.NullFormatter())
        ax.xaxis.set_minor_formatter(mticker.NullFormatter())

    for ax, letter, dx in ((axA, "A", -0.17), (axB, "B", -0.12),
                           (axC, "C", -0.16), (axD, "D", -0.16)):
        st.panel_letter(ax, letter, dx=dx, dy=1.06)
    return fig


# Figure 1 says "Oprl1 blue, the other three grey". The supplemental panels put
# all four side by side, so the other three are separated by value rather than
# by hue, and the sentence still holds.
RECEPTOR_FILL = {"Oprl1": st.BAR_BLUE, "Oprm1": "#595959",
                 "Oprd1": "#8C8C8C", "Oprk1": "#C9C9C9"}
MIN_POSITIVE = 100      # matches src/17_prevalence_decomposition.py


def _population_labels(d):
    """Panel letter and a short tissue name, for an x axis of 15 populations."""
    short = {"enteric submucosal, P24": "enteric P24",
             "enteric submucosal, P7": "enteric P7",
             "superior cervical": "sup. cervical",
             "sphenopalatine (VII)": "sphenopalatine"}
    out = []
    for panel, tissue in d[["panel", "population"]].drop_duplicates().values:
        t = short.get(tissue, tissue.split(" (")[0])
        out.append(f"{panel}  {t}")
    return out


def figureS1():
    """A population mean is a prevalence and a per-cell level, separated here."""
    d = pd.read_csv(ac.RES / "receptor_prevalence_decomposition.csv")
    o = pd.read_csv(ac.RES / "receptor_prevalence_orderings.csv")
    o = o[o.margin_per_cell.notna()]
    panels = list(dict.fromkeys(d.panel))
    labels = _population_labels(d)

    # Laid out in inches from the top down. Both rows carry rotated population
    # names, so each needs BELOW beneath it and ABOVE for its panel letter.
    PANEL_H, BELOW, ABOVE = 1.55, 0.78, 0.34
    height = ABOVE + PANEL_H + (BELOW + ABOVE) + PANEL_H + BELOW

    st.set_theme()
    fig = plt.figure(figsize=(st.W_SUPP, height))
    left, right = 0.055, 0.99
    row1 = height - ABOVE
    gs1 = fig.add_gridspec(1, 1, left=left, right=right, top=row1 / height,
                           bottom=(row1 - PANEL_H) / height)
    axA = fig.add_subplot(gs1[0, 0])
    row2 = row1 - PANEL_H - BELOW - ABOVE
    gs2 = fig.add_gridspec(1, 3, left=left, right=right, top=row2 / height,
                           bottom=(row2 - PANEL_H) / height, wspace=0.38)
    axB = fig.add_subplot(gs2[0, 0:2])
    axC = fig.add_subplot(gs2[0, 2])

    x = np.arange(len(panels))
    w = 0.20
    for i, gene in enumerate(ac.RECEPTORS):
        g = d[d.gene == gene].set_index("panel").loc[panels]
        off = (i - 1.5) * w
        axA.bar(x + off, g.pct_detected, w, color=RECEPTOR_FILL[gene],
                edgecolor="black", linewidth=0.4, zorder=3,
                label=rf"$\it{{{gene}}}$")
        # A conditional mean over few positive cells is floored by the one-count
        # detection limit rather than measured, so it is drawn but faded.
        solid = g.n_positive >= MIN_POSITIVE
        axB.bar(x + off, g.level_if_positive.where(solid), w,
                color=RECEPTOR_FILL[gene], edgecolor="black", linewidth=0.4, zorder=3)
        axB.bar(x + off, g.level_if_positive.where(~solid), w,
                color=RECEPTOR_FILL[gene], edgecolor="black", linewidth=0.4,
                alpha=0.30, zorder=3)

    for ax, ylab in ((axA, "neurons detecting (%)"),
                     (axB, "level in positive cells (CPM)")):
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=st.FS_TICK)
        ax.set_xlim(-0.6, len(panels) - 0.4)
        ax.set_ylabel(ylab, fontsize=st.FS_LABEL)
    axA.set_ylim(0, 100)
    axB.set_yscale("log")
    axA.legend(loc="upper right", fontsize=st.FS_NOTE, ncol=4, frameon=False,
               handlelength=1.1, columnspacing=1.0, handletextpad=0.4)

    lim = (0.9, max(o.margin_by_mean.max(), o.margin_per_cell.max()) * 1.6)
    axC.plot(lim, lim, color="black", lw=0.5, ls="--", zorder=1)
    axC.scatter(o.margin_by_mean, o.margin_per_cell, s=16, color=st.BAR_BLUE,
                edgecolor="black", linewidth=0.4, zorder=3)
    # F and G sit on top of one another at the foot of the identity line, so
    # labels alternate side rather than all trailing right.
    for i, r in enumerate(o.sort_values("margin_by_mean").itertuples()):
        dx, dy = ((4, -1) if i % 2 == 0 else (-4, 3))
        axC.annotate(r.panel, (r.margin_by_mean, r.margin_per_cell),
                     textcoords="offset points", xytext=(dx, dy),
                     ha="left" if dx > 0 else "right", fontsize=st.FS_NOTE)
    axC.set_xscale("log"); axC.set_yscale("log")
    axC.set_xlim(*lim); axC.set_ylim(*lim)
    axC.set_xlabel("margin on the population mean", fontsize=st.FS_LABEL)
    axC.set_ylabel("margin in positive cells", fontsize=st.FS_LABEL)

    for ax, letter in zip((axA, axB, axC), "ABC"):
        st.panel_letter(ax, letter, dx=-0.055 if ax is not axC else -0.20, dy=1.10)
    return fig


def figureS2():
    """Dissociation is a second preparation effect, and it is not this one."""
    g = pd.read_csv(ac.RES / "preparation_bias_genomewide.csv")
    s = pd.read_csv(ac.RES / "dissociation_signature.csv")
    dis, rec = s[s.set == "dissociation"], s[s.set == "receptor"]

    PANEL_H, BELOW, ABOVE = 2.05, 0.42, 0.30
    height = ABOVE + PANEL_H + BELOW + 0.08
    st.set_theme()
    fig = plt.figure(figsize=(st.W_2COL, height))
    # The 99th-percentile marker in (B) sits on the axis edge and its radius
    # falls outside the canvas at a full-bleed right margin.
    gs = fig.add_gridspec(1, 2, left=0.085, right=0.982,
                          top=(height - ABOVE) / height,
                          bottom=BELOW / height, wspace=0.34)
    axA, axB = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])

    axA.scatter(g.span_kb, g.ratio, s=0.6, color=st.FEATURE_LOW, linewidths=0,
                rasterized=True, zorder=1)
    axA.axhline(1.0, color="black", lw=0.5, ls="--", zorder=2)
    axA.scatter(dis.span_kb, dis.ratio, s=14, color=st.HIGHLIGHT,
                edgecolor="black", linewidth=0.4, zorder=4,
                label="dissociation-induced")
    axA.scatter(rec.span_kb, rec.ratio, s=14,
                color=[RECEPTOR_FILL[x] for x in rec.gene],
                edgecolor="black", linewidth=0.4, zorder=4)
    for r in rec.itertuples():
        axA.annotate(rf"$\it{{{r.gene}}}$", (r.span_kb, r.ratio),
                     textcoords="offset points", xytext=(5, -1),
                     fontsize=st.FS_NOTE)
    axA.set_xscale("log"); axA.set_yscale("log")
    axA.set_xlabel("genomic span (kb)", fontsize=st.FS_LABEL)
    axA.set_ylabel("nuclear / whole cell", fontsize=st.FS_LABEL)
    axA.legend(loc="upper left", fontsize=st.FS_NOTE, frameon=False,
               handletextpad=0.3, borderpad=0.1)

    # Percentile within the gene's own length class, which is what separates an
    # induced gene from a long one.
    both = pd.concat([dis.sort_values("pct_of_peers"),
                      rec.sort_values("pct_of_peers")], ignore_index=True)
    ypos = np.arange(len(both))
    colors = [st.HIGHLIGHT if r.set == "dissociation" else RECEPTOR_FILL[r.gene]
              for r in both.itertuples()]
    axB.hlines(ypos, 0, both.pct_of_peers, color=colors, lw=0.9, zorder=3)
    axB.scatter(both.pct_of_peers, ypos, s=13, color=colors, edgecolor="black",
                linewidth=0.4, zorder=4)
    axB.axvline(50, color="black", lw=0.5, ls="--", zorder=2)
    axB.set_yticks(ypos)
    axB.set_yticklabels([rf"$\it{{{n}}}$" for n in both.gene], fontsize=st.FS_TICK)
    axB.set_ylim(-0.8, len(both) - 0.2)
    axB.set_xlim(0, 103)
    axB.set_xlabel("percentile among genes of the same span", fontsize=st.FS_LABEL)

    for ax, letter in zip((axA, axB), "AB"):
        st.panel_letter(ax, letter, dx=-0.16, dy=1.06)
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
    emit(figure2(), "Figure2")
    emit(figure3(), "Figure3")
    # Supplemental figures are supplied as separate files, so the 174 mm text
    # column does not bind them; S1 needs the width for 15 populations.
    emit(figureS1(), "FigureS1",
         max_w=st.W_SUPP * 25.4 + 0.5, max_h=st.H_SUPP * 25.4)
    emit(figureS2(), "FigureS2")
    emit(figureS3(), "FigureS3")
