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

# Panels carry identifiers only: panel letter, tissue, sample size, gene symbol,
# axis value, unit, and the division name. Margins, bootstrap support, dataset
# accessions and preparation live in paper/figure_legends.md.


def figure1():
    """The four opioid receptor genes in 19 peripheral neuronal populations."""
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
    return fig


def emit(fig, name):
    """Write the figure and refuse it if it breaks the journal's limits."""
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{name}.pdf")
    fig.savefig(OUT / f"{name}.png")
    plt.close(fig)

    w, h = (v * 25.4 for v in fig.get_size_inches())
    over = []
    if w > 174.5:
        over.append(f"width {w:.1f} mm exceeds 174 mm")
    if h > 235.0:
        over.append(f"height {h:.1f} mm exceeds 235 mm")
    if over:
        raise ac.SanityCheckError(f"{name}: " + "; ".join(over))
    print(f"  [paper] paper/figures/{name}.pdf|.png  {w:.1f} x {h:.1f} mm")


if __name__ == "__main__":
    emit(figure1(), "Figure1")
