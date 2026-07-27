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


def set_theme() -> None:
    mpl.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Nimbus Sans", "Helvetica", "Arial"],
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": "black",
        "axes.linewidth": 0.8,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.frameon": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


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
