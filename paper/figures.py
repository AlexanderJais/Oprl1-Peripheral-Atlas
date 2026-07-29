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


def emit(fig, name):
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
    if w > 174.5:
        over.append(f"width {w:.1f} mm exceeds 174 mm")
    if h > 235.0:
        over.append(f"height {h:.1f} mm exceeds 235 mm")
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
