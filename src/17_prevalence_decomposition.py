"""Is the receptor ordering a difference in level, or in how many neurons carry it?

A pseudobulk mean over a population is the product of two quantities that a
single number hides: the fraction of neurons in which the gene is detected, and
the mean level in those neurons. Two genes can reach the same population mean
with opposite biology, one present at a low level in most neurons and one at a
high level in a few, and every claim in this atlas so far has rested on the
product rather than on either factor.

This script separates them. For every population where detection is recorded for
all four receptors, it reports the prevalence, the mean among positive cells,
and the ordering each of the two produces.

The conditional mean is biased, and in a direction that matters here. A droplet
library cannot record less than one count, so a gene detected in very few cells
has its conditional mean floored by the quantisation limit rather than measured.
That bias inflates the rare receptors and therefore works against Oprl1, so the
comparison is conservative for the prevalence claim and unreliable wherever the
positive-cell count is small. `n_positive` is reported for every entry so a
reader can see which conditional means rest on a handful of cells, and no
ordering is called from a gene detected in fewer than `MIN_POSITIVE` of them.

Run from the repository root, after 10_peripheral_ganglia.py and 11_quality_panel.py.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import atlas_common as ac

# A conditional mean over fewer positive cells than this is dominated by the
# detection floor and by whatever ambient signal reached those barcodes. It is
# reported, but it is not allowed to decide an ordering.
MIN_POSITIVE = 100


def key(dataset, tissue):
    return f"{dataset}|{str(tissue).split(' (')[0].strip()}"


def gather() -> pd.DataFrame:
    """Level and detection for all four receptors, per population.

    Three tables record detection, because the populations entered the project
    through three different scripts. A population is included only where all
    four receptors carry both a level and a detection rate, since the point of
    the decomposition is the comparison between them.
    """
    fig = pd.read_csv(ac.RES / "receptor_levels_by_ganglion.csv")
    rows = {}

    qual = pd.read_csv(ac.RES / "dataset_quality_panel.csv")
    for r in qual.itertuples():
        vals = {g: (getattr(r, g), getattr(r, f"{g}_pct")) for g in ac.RECEPTORS}
        if all(pd.notna(v) and pd.notna(p) for v, p in vals.values()):
            rows[key(r.dataset, r.tissue)] = (vals, int(r.n))

    for name, ds_col in (("geniculate_opioid_levels.csv", None),
                         ("nodose_opioid_levels.csv", "NodoMap:whole-cell")):
        d = pd.read_csv(ac.RES / name)
        if ds_col:
            d = d[d.dataset == ds_col]
        d = d[d.gene.isin(ac.RECEPTORS)]
        for (dataset, tissue), grp in d.groupby(["dataset", "tissue"]):
            g = grp.set_index("gene")
            if not set(ac.RECEPTORS) <= set(g.index):
                continue
            if g.loc[ac.RECEPTORS, ["mean_level", "pct_detected"]].isna().any().any():
                continue
            label = "NodoMap" if ds_col else dataset
            rows[key(label, tissue)] = (
                {x: (float(g.loc[x, "mean_level"]), float(g.loc[x, "pct_detected"]))
                 for x in ac.RECEPTORS},
                int(g.loc[ac.RECEPTORS[0], "n_cells"]))

    out = []
    for r in fig.itertuples():
        k = key(r.dataset, r.tissue)
        if k not in rows:
            continue
        vals, n = rows[k]
        for gene in ac.RECEPTORS:
            level, pct = vals[gene]
            out.append({
                "panel": str(r.panel).upper(), "division": r.group,
                "population": r.tissue, "dataset": r.dataset, "unit": r.unit,
                "n_neurons": n, "gene": gene,
                "mean_level": round(level, 4),
                "pct_detected": round(pct, 3),
                "n_positive": int(round(pct / 100 * n)),
                "level_if_positive": round(level / (pct / 100), 3) if pct > 0 else np.nan,
            })
    return pd.DataFrame(out)


def orderings(tbl: pd.DataFrame) -> pd.DataFrame:
    """The ordering the population mean gives, against the one per-cell level gives."""
    rows = []
    for panel, g in tbl.groupby("panel", sort=False):
        g = g.set_index("gene")
        by_mean = g.mean_level.sort_values(ascending=False)
        usable = g[(g.n_positive >= MIN_POSITIVE) & g.level_if_positive.notna()]
        by_cond = usable.level_if_positive.sort_values(ascending=False)
        row = {
            "panel": panel, "population": g.population.iloc[0],
            "n_neurons": int(g.n_neurons.iloc[0]),
            "top_by_mean": by_mean.index[0],
            "margin_by_mean": round(by_mean.iloc[0] / by_mean.iloc[1], 3)
                              if by_mean.iloc[1] > 0 else np.inf,
            "n_genes_testable_per_cell": len(by_cond),
        }
        if len(by_cond) >= 2:
            row.update({
                "top_per_cell": by_cond.index[0],
                "margin_per_cell": round(by_cond.iloc[0] / by_cond.iloc[1], 3),
                "ordering_changes": by_cond.index[0] != by_mean.index[0],
            })
        else:
            row.update({"top_per_cell": None, "margin_per_cell": np.nan,
                        "ordering_changes": None})
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> int:
    tbl = gather()
    if tbl.empty:
        raise ac.SanityCheckError("no population records detection for all four receptors")
    ac.save_table(tbl, "receptor_prevalence_decomposition.csv")

    ordr = orderings(tbl)
    ac.save_table(ordr, "receptor_prevalence_orderings.csv")

    pan = sorted(tbl.panel.unique())
    print(f"\n  {len(pan)} of 19 Figure 1 panels record detection for all four "
          f"receptors: {', '.join(pan)}")

    piv = tbl.pivot(index="panel", columns="gene", values="pct_detected")
    print("\n  Prevalence, percent of neurons detecting each receptor:")
    print(piv[ac.RECEPTORS].round(1).to_string())

    o = ordr[ordr.margin_per_cell.notna()]
    print(f"\n  Margin on the population mean, median {o.margin_by_mean.median():.2f}x;"
          f" among positive cells, median {o.margin_per_cell.median():.2f}x.")
    changed = o[o.ordering_changes]
    print(f"  The top receptor differs between the two in {len(changed)} of "
          f"{len(o)} populations: {', '.join(changed.panel)}")

    l1 = tbl[tbl.gene == "Oprl1"].pct_detected
    m1 = tbl[tbl.gene == "Oprm1"].pct_detected
    print(f"\n  Oprl1 is detected in {l1.min():.1f} to {l1.max():.1f}% of neurons, "
          f"Oprm1 in {m1.min():.1f} to {m1.max():.1f}%.")
    most_prevalent = piv[ac.RECEPTORS].idxmax(axis=1)
    print(f"  Oprl1 is the most prevalent of the four in "
          f"{int((most_prevalent == 'Oprl1').sum())} of {len(pan)} populations.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
