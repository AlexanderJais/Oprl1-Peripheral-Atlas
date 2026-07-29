"""Plotting helpers, matched to the sibling PNOC-Nodose project.

The figure idiom here is deliberately the same one used in
https://github.com/AlexanderJais/PNOC-Nodose (`src/nodomap_style.py`), which in
turn follows the NodoMap authors' own plotting code: Seurat-style ranked bars,
dot plots sized by percent-expressing and coloured by mean expression, and
FeaturePlot-style UMAP overlays.

Everything in this module describes `Oprl1` where it is expressed. Nothing here
computes a derived ratio.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import numpy as np
import pandas as pd

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
FIG = ROOT / "figures"

# One colour per tissue, used wherever tissues appear side by side.
TISSUE_COLORS = {
    "geniculate": "#1B7837",
    "nodose": "#08306B",
    "nodose+jugular": "#08306B",
    "jugular": "#6BAED6",
    "NTS": "#B2182B",
    "NTS (central)": "#B2182B",
}
FEATURE_LOW = "#DEDEDE"

# Preparation is hatched rather than coloured, so tissue keeps the colour axis.
PREP_HATCH = {"whole cell": "", "nuclear": "///"}


# Cell Press artwork specification. Column widths are fixed by the journal and
# every main figure is built to one of them, so nothing is rescaled at
# submission and no font is resized away from the value set here.
MM = 1 / 25.4
W_1COL, W_15COL, W_2COL = 85 * MM, 114 * MM, 174 * MM
H_MAX = 235 * MM

# Arial or Helvetica, and nothing else. Helvetica is listed first for a
# production system that licenses it; Nimbus Sans is the URW clone with
# identical metrics, and Liberation Sans carries Arial's metrics. DejaVu, the
# matplotlib default, is neither and is last so a missing font is visible rather
# than silently substituted.
SANS = ["Helvetica", "Nimbus Sans", "Arial", "Liberation Sans", "FreeSans",
        "DejaVu Sans"]

# Type sizes in points, as they will appear on the printed page.
FS_PANEL = 8       # panel letter, bold capital
FS_LABEL = 7       # axis labels
FS_TICK = 6.5      # tick labels
FS_NOTE = 6.5      # in-panel identifiers


def set_theme(base: float = FS_TICK) -> None:
    """Cell Press house style: Arial or Helvetica, thin rules, vector text."""
    mpl.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "font.family": "sans-serif",
        "font.sans-serif": SANS,
        "font.size": base,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": "black",
        # Cell Press sets 0.25 pt as the minimum reproducible rule; 0.5 pt
        # survives reduction without thickening the figure.
        "axes.linewidth": 0.5,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "xtick.major.size": 2.0,
        "ytick.major.size": 2.0,
        "axes.titlesize": FS_NOTE,
        "axes.labelsize": FS_LABEL,
        "xtick.labelsize": FS_TICK,
        "ytick.labelsize": FS_TICK,
        "legend.frameon": False,
        "legend.fontsize": FS_TICK,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })
    # Report which family actually resolved, so a missing Helvetica is caught
    # here rather than at proof stage.
    resolved = mpl.font_manager.findfont(
        mpl.font_manager.FontProperties(family=SANS))
    name = mpl.font_manager.FontProperties(fname=resolved).get_name()
    if name not in ("Helvetica", "Nimbus Sans", "Arial", "Liberation Sans"):
        print(f"  [warn] the sans-serif stack resolves to {name}; Cell Press "
              "requires Arial or Helvetica")
    # Point mathtext at the same family, so an italic gene symbol inside an axis
    # label is set in the figure's own font rather than matplotlib's default
    # serif math face.
    mpl.rcParams.update({
        "mathtext.fontset": "custom",
        "mathtext.rm": name,
        "mathtext.it": f"{name}:italic",
        "mathtext.bf": f"{name}:bold",
        "mathtext.cal": name,
        "mathtext.sf": name,
        "mathtext.tt": name,
        "mathtext.default": "regular",
    })
    return name


def save(fig, name: str) -> None:
    """Write both a vector PDF and a PNG, as the sibling project does."""
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / f"{name}.pdf")
    fig.savefig(FIG / f"{name}.png")
    plt.close(fig)
    print(f"  [fig ] figures/{name}.pdf|.png")


# ---------------------------------------------------------------------------
# Per-group summaries
# ---------------------------------------------------------------------------

def percent_expressing(expr: pd.DataFrame, groups) -> pd.DataFrame:
    """Percentage of cells per group with non-zero expression."""
    detected = (expr > 0).astype(float)
    detected["__group__"] = np.asarray(groups)
    pct = detected.groupby("__group__", observed=True).mean() * 100.0
    pct.index.name = getattr(groups, "name", None)
    return pct


def mean_expression(expr: pd.DataFrame, groups) -> pd.DataFrame:
    """Mean expression per group, over all cells including zeros."""
    tmp = expr.copy()
    tmp["__group__"] = np.asarray(groups)
    out = tmp.groupby("__group__", observed=True).mean()
    out.index.name = getattr(groups, "name", None)
    return out


# ---------------------------------------------------------------------------
# Figure primitives
# ---------------------------------------------------------------------------

# Every main-text panel is a bar of mean expression in the unit the dataset was
# measured in, or a violin of the per-cell values behind it. No ratios, no ranks,
# no dot plots: the reader should be able to read the number off the axis.
BAR_BLUE = "#1F7FEF"
BAR_GREY = "#BFBFBF"
HIGHLIGHT = "#D62728"


def panel_letter(ax, letter, dx=-0.16, dy=1.06, fontsize=FS_PANEL):
    """Bold capital, top left of the panel, as Cell Press sets them."""
    ax.text(dx, dy, str(letter).upper(), transform=ax.transAxes,
            fontsize=fontsize, fontweight="bold", va="top", ha="left")


def expression_bars(ax, values, labels, ylabel, colors=None, nd_mask=None,
                    italic=True, rotation=45, annotate=False, fontsize=13,
                    linewidth=1.4):
    """Mean expression per gene, in the dataset's own unit.

    `nd_mask` marks genes that were measured and not detected; they get an
    "n.d." tick rather than a zero-height bar, so an undetected gene cannot be
    mistaken for a missing one.
    """
    v = np.asarray(values, dtype=float)
    x = np.arange(len(v))
    if colors is None:
        colors = [BAR_BLUE] * len(v)
    ax.bar(x, np.nan_to_num(v), color=colors, edgecolor="black",
           linewidth=linewidth, width=0.68, zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=rotation, ha="right", fontsize=fontsize,
                       fontstyle="italic" if italic else "normal")
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.tick_params(axis="y", labelsize=fontsize)
    top = float(np.nanmax(v)) if np.isfinite(v).any() else 1.0
    ax.set_ylim(0, top * 1.15)
    if nd_mask is not None:
        for i, nd in enumerate(nd_mask):
            if nd:
                ax.text(i, top * 0.02, "n.d.", ha="center", va="bottom",
                        fontsize=fontsize - 3)
    if annotate:
        for i, val in enumerate(v):
            if np.isfinite(val):
                ax.text(i, val + top * 0.02, f"{val:.2f}".rstrip("0").rstrip("."),
                        ha="center", va="bottom", fontsize=fontsize - 3)
    return ax


def violin_points(ax, groups, values_by_group, ylabel, color=BAR_BLUE,
                  fontsize=13, seed=0):
    """Violin + every individual cell + a median bar, as in the geniculate figure."""
    data = [np.asarray(values_by_group[g], dtype=float) for g in groups]
    parts = ax.violinplot(data, positions=range(len(groups)), widths=0.8,
                          showextrema=False, showmedians=False)
    for body in parts["bodies"]:
        body.set_facecolor(color)
        body.set_alpha(0.45)
        body.set_edgecolor("#555555")
        body.set_linewidth(1.0)
    rng = np.random.default_rng(seed)
    for i, v in enumerate(data):
        ax.scatter(i + rng.uniform(-0.11, 0.11, v.size), v, s=26, color=color,
                   edgecolors="black", linewidths=0.6, zorder=3, alpha=0.95)
        ax.hlines(np.median(v), i - 0.34, i + 0.34, color="black", lw=3.0, zorder=4)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(groups, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.tick_params(axis="y", labelsize=fontsize)
    ax.set_xlim(-0.6, len(groups) - 0.4)
    return ax


def ranked_bars(ax, values, labels, xlabel, cmap="viridis", threshold=None,
                threshold_label=None, annotate=None, fontsize=8):
    """One gene across many groups, sorted high to low, viridis-coloured.

    The direct analogue of PNOC-Nodose `figure2_oprl1_ranked_bars`: the reader
    should be able to see at a glance which populations carry the transcript.
    """
    v = np.asarray(values, dtype=float)
    labels = list(labels)
    order = np.argsort(-v)
    v, labels = v[order], [labels[i] for i in order]
    if annotate is not None:
        annotate = [list(annotate)[i] for i in order]

    top = float(np.nanmax(v)) if np.isfinite(v).any() and np.nanmax(v) > 0 else 1.0
    colors = plt.get_cmap(cmap)(v / top)
    y = np.arange(len(v))
    ax.barh(y, v, color=colors, edgecolor="black", linewidth=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=fontsize)
    ax.invert_yaxis()
    ax.set_xlabel(xlabel)

    right = top * 1.18 if annotate is not None else top * 1.02
    if threshold is not None:
        right = max(right, threshold * 1.05)
    ax.set_xlim(0, right)

    if threshold is not None:
        ax.axvline(threshold, color="#D62728", lw=1.0, ls="--")
        if threshold_label:
            # Right-align inside the axes: a threshold usually sits near the
            # top of the range, where a left-aligned label would be clipped.
            ax.text(threshold * 0.98, len(v) - 0.4, threshold_label,
                    color="#D62728", fontsize=7, va="bottom", ha="right")
    if annotate is not None:
        pad = top * 0.015
        for i, txt in enumerate(annotate):
            ax.text(v[i] + pad, i, txt, va="center", fontsize=6.5, color="#333333")
    return ax


def strip_by_group(ax, groups, values_by_group, ylabel, color=BAR_BLUE,
                   fontsize=11, seed=0, log=False):
    """One point per cluster, grouped, with a median bar.

    Used where the grouping variable is a property of the cluster rather than of
    the cell: plotting cells there would treat 26,000 correlated observations as
    independent, and the honest n is the number of clusters.
    """
    rng = np.random.default_rng(seed)
    for i, g in enumerate(groups):
        v = np.asarray(values_by_group[g], dtype=float)
        if not v.size:
            continue
        ax.scatter(i + rng.uniform(-0.13, 0.13, v.size), v, s=44, color=color,
                   edgecolors="black", linewidths=0.7, zorder=3, alpha=0.95)
        ax.hlines(np.median(v), i - 0.3, i + 0.3, color="black", lw=2.6, zorder=4)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(groups, fontsize=fontsize - 1)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.set_xlim(-0.6, len(groups) - 0.4)
    if log:
        ax.set_yscale("log")
    return ax


def dim_plot(ax, umap, labels, order, colors, background=None, point_size=1.2):
    """Cells coloured by a categorical label, on the published embedding."""
    if background is not None:
        ax.scatter(umap[background, 0], umap[background, 1], s=0.35, c="#EDEDED",
                   linewidths=0, rasterized=True)
    labels = np.asarray(labels)
    for lab in order:
        m = labels == lab
        if m.any():
            ax.scatter(umap[m, 0], umap[m, 1], s=point_size, c=colors[lab],
                       linewidths=0, rasterized=True, label=f"{lab} ({m.sum():,})")
    ax.set_aspect("equal")
    ax.axis("off")
    return ax


def dot_plot(ax, pct: pd.DataFrame, mean: pd.DataFrame, cmap="viridis",
             max_dot=90.0, vmax=None, xtick_fontsize=7):
    """Seurat DotPlot: rows = genes, cols = groups, size = % expressing."""
    genes = list(pct.columns)
    groups = list(pct.index)
    xs, ys, ss, cs = [], [], [], []
    for gi, g in enumerate(genes):
        for ci, c in enumerate(groups):
            xs.append(ci)
            ys.append(len(genes) - 1 - gi)
            ss.append(float(pct.loc[c, g]) / 100.0 * max_dot)
            cs.append(float(mean.loc[c, g]))
    vmax = vmax if vmax is not None else max(np.max(cs), 1e-6)
    sc = ax.scatter(xs, ys, s=ss, c=cs, cmap=cmap, vmin=0, vmax=vmax,
                    linewidths=0.2, edgecolors="black")
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(groups, rotation=90, fontsize=xtick_fontsize)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(list(reversed(genes)), fontstyle="italic", fontsize=9)
    ax.set_xlim(-0.8, len(groups) - 0.2)
    ax.set_ylim(-0.8, len(genes) - 0.2)
    ax.tick_params(length=2)
    return sc


def dot_size_legend(ax, max_dot=90.0, values=(1, 5, 10, 20, 30),
                    title="% expressing"):
    handles = [
        ax.scatter([], [], s=v / 100.0 * max_dot, c="#4d4d4d", linewidths=0.2,
                   edgecolors="black", label=f"{v}")
        for v in values
    ]
    return ax.legend(handles=handles, title=title, loc="center left",
                     bbox_to_anchor=(1.005, 0.5), labelspacing=1.1,
                     fontsize=7, title_fontsize=7, borderpad=0.8,
                     handletextpad=0.6)


def feature_plot(ax, umap, values, title, high="#08306B", vmax=None,
                 point_size=0.6):
    """Seurat FeaturePlot look: grey background cells, expression in colour."""
    from matplotlib.colors import LinearSegmentedColormap

    cmap = LinearSegmentedColormap.from_list("feat", [FEATURE_LOW, high])
    v = np.asarray(values, dtype=float)
    if vmax is None:
        vmax = np.percentile(v[v > 0], 99) if (v > 0).any() else 1.0
    vmax = max(float(vmax), 1e-6)
    order = np.argsort(v)
    ax.scatter(umap[order, 0], umap[order, 1], c=v[order], cmap=cmap,
               s=point_size, vmin=0, vmax=vmax, linewidths=0, rasterized=True)
    ax.set_title(title, fontstyle="italic")
    ax.set_aspect("equal")
    ax.axis("off")
    return cmap, vmax


def add_colorbar(fig, ax, cmap, vmax, label="mean expression"):
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize

    sm = ScalarMappable(norm=Normalize(0, vmax), cmap=cmap)
    cb = fig.colorbar(sm, ax=ax, fraction=0.030, pad=0.01, shrink=0.55, aspect=18)
    cb.ax.tick_params(labelsize=6, length=2)
    cb.set_label(label, size=6)
    cb.outline.set_linewidth(0.4)
    return cb


def violin_by_group(ax, values, groups, order, colors, width=0.85):
    """Seurat VlnPlot look: one violin per group, filled in the group colour."""
    values = np.asarray(values, dtype=float)
    groups = np.asarray(groups)
    data, present = [], []
    for c in order:
        v = values[groups == c]
        if v.size:
            data.append(v)
            present.append(c)
    parts = ax.violinplot(data, positions=range(len(data)), widths=width,
                          showextrema=False, showmedians=True)
    for body, c in zip(parts["bodies"], present):
        body.set_facecolor(colors.get(c, "#999999"))
        body.set_edgecolor("black")
        body.set_linewidth(0.5)
        body.set_alpha(0.9)
    if "cmedians" in parts:
        parts["cmedians"].set_color("black")
        parts["cmedians"].set_linewidth(1.0)
    ax.set_xticks(range(len(present)))
    ax.set_xticklabels(present, fontsize=8)
    return present
