"""Where Oprl1 sits among the GPCRs, not among four opioid receptors.

"Highest of the four opioid receptors" is a claim about a family of four, and it
invites the reply that the family is uniformly low. This script replaces the
denominator: Oprl1 is ranked against every G protein-coupled receptor measured
in the same cells, in every population that carries a full transcriptome.

THE DENOMINATOR

The receptor list is the curated non-sensory set from the IUPHAR/BPS Guide to
Pharmacology, 356 mouse symbols. Olfactory, vomeronasal and taste receptors are
excluded, and that exclusion is the point rather than an oversight: the mouse
genome carries over a thousand olfactory receptors that are silent outside the
olfactory epithelium, and counting them would inflate any rank in any tissue
without saying anything. Two denominators are therefore reported for each
population: every GPCR present in that annotation, and every GPCR actually
detected in those neurons. The second is the one a reader should use.

WHAT IS RESAMPLED

The populations, as in section 2. A rank in one ganglion is one measurement; the
statistic reported is the median rank across populations with a 95% interval
from 10,000 draws of the populations with replacement.

Run from the repository root. The receptor list is fetched into data/raw on
first use; the expression matrices come through src/12_oprl1_signature.py, which
applies the same neuron definition as the rest of the atlas.
"""

import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp

import atlas_common as ac
import atlas_style as st

_sig = __import__("12_oprl1_signature")

GTOP = ac.DATA / "gtopdb_targets_and_families.csv"
GTOP_URL = "https://www.guidetopharmacology.org/DATA/targets_and_families.csv"

# The geniculate SMART-seq data is back in: section 2 had to drop it because 88
# of its 96 neurons are Oprl1-positive and there was nothing to compare against,
# but a rank needs no split.
POPULATIONS = _sig.POPULATIONS

N_BOOT = ac.N_BOOT
TOP_SHOWN = 25


def gpcr_symbols():
    """Mouse symbols for the curated non-sensory GPCRs."""
    if not GTOP.exists():
        print(f"  downloading {GTOP.name}")
        GTOP.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(["curl", "-fsSL", "-o", str(GTOP), GTOP_URL], check=True)
    d = pd.read_csv(GTOP, skiprows=1, low_memory=False)
    rows = d[d.Type.str.lower() == "gpcr"]["MGI symbol"].dropna()
    # A few entries list several symbols for one target, pipe-separated.
    symbols = sorted({s.strip() for v in rows for s in str(v).split("|")
                      if s.strip()})
    if "Oprl1" not in symbols:
        raise ac.SanityCheckError(
            "Oprl1 is not in the GPCR list; the rank would be meaningless")
    for g in ac.RECEPTORS:
        if g not in symbols:
            raise ac.SanityCheckError(f"{g} is missing from the GPCR list")
    print(f"  {len(symbols)} mouse GPCR symbols "
          f"(non-sensory, IUPHAR/BPS Guide to Pharmacology)")
    return symbols


def rank_in(acc, tissue, symbols):
    """Oprl1's position among the GPCRs of one population."""
    cpm, counts, genes, _ = _sig.neurons_of(acc, tissue)
    genes = pd.Index(genes)
    mean, det = _sig._column_stats(cpm, genes)

    here = [g for g in symbols if g in genes]
    panel = pd.DataFrame({"mean": mean[here], "detected": det[here]})
    panel = panel.sort_values("mean", ascending=False)
    detected = panel[panel.detected > 0]

    def position(frame, gene):
        return (int(frame.index.get_loc(gene)) + 1 if gene in frame.index
                else np.nan)

    row = {
        "dataset": acc, "tissue": tissue, "n_cells": int(cpm.shape[0]),
        "n_gpcr_in_annotation": len(panel),
        "n_gpcr_detected": len(detected),
        "Oprl1_CPM": round(float(panel.loc["Oprl1", "mean"]), 3),
        "Oprl1_pct_detected": round(float(panel.loc["Oprl1", "detected"]) * 100, 2),
        "rank_of_annotated": position(panel, "Oprl1"),
        "rank_of_detected": position(detected, "Oprl1"),
        "percentile_of_detected": round(
            100 * (1 - (position(detected, "Oprl1") - 1) / len(detected)), 1),
        "transcriptome_percentile": round(
            ac.transcriptome_percentile(mean, "Oprl1"), 1),
    }
    for g in ac.RECEPTORS[1:]:
        row[f"{g}_rank"] = position(detected, g)
    row["higher_than_Oprl1"] = ";".join(
        detected.index[:position(detected, "Oprl1") - 1])
    print(f"    {acc} {tissue}: Oprl1 {row['Oprl1_CPM']:.2f} CPM, rank "
          f"{row['rank_of_detected']} of {len(detected)} detected GPCRs "
          f"({row['percentile_of_detected']:.0f}th percentile), "
          f"{len(panel)} in the annotation")
    return row, panel, detected


def bootstrap_median_rank(values, n_boot=N_BOOT, seed=0):
    """95% interval on the median rank, resampling populations."""
    rng = np.random.default_rng(seed)
    v = np.asarray(values, dtype=float)
    draws = np.median(v[rng.integers(0, len(v), (n_boot, len(v)))], axis=1)
    return float(np.median(v)), float(np.percentile(draws, 2.5)), \
        float(np.percentile(draws, 97.5))


# ---------------------------------------------------------------------- figure
def figure3(table, panels, beaters):
    st.set_theme()
    fig = plt.figure(figsize=(17.5, 8.8))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.0, 0.95], wspace=0.55,
                          left=0.055, right=0.985, top=0.86, bottom=0.30)

    # (a) the ranked GPCR panel of the deepest droplet population, raw units
    name, panel, oprl1_rank = panels
    ax = fig.add_subplot(gs[0])
    # The panel has to reach Oprl1, or the highlight it exists to show is
    # off the end of the axis.
    top = panel[panel.detected > 0].head(max(TOP_SHOWN, oprl1_rank + 3))
    colors = [st.BAR_BLUE if g == "Oprl1" else st.BAR_GREY for g in top.index]
    st.expression_bars(ax, top["mean"].to_numpy(), list(top.index),
                       "Mean expression (CPM)", colors=colors, rotation=90,
                       fontsize=10)
    st.panel_letter(ax, "a", dx=-0.10)
    ax.set_title(f"The {len(top)} highest GPCRs in one population\n{name}",
                 fontsize=11, pad=10)

    # (b) the rank itself, every population
    ax = fig.add_subplot(gs[1])
    t = table.sort_values("rank_of_detected", ascending=False)
    y = np.arange(len(t))
    ax.hlines(y, 0, t.rank_of_detected, color="#999999", lw=1.2, zorder=2)
    ax.scatter(t.rank_of_detected, y, s=90, c=st.BAR_BLUE, edgecolors="black",
               linewidths=0.9, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r.tissue} ({r.dataset})\n{r.n_gpcr_detected} GPCRs "
                        f"detected" for r in t.itertuples()], fontsize=8.5)
    ax.set_xlabel("Rank of Oprl1 among the GPCRs detected", fontsize=11)
    ax.set_xlim(0, float(t.rank_of_detected.max()) * 1.12)
    ax.invert_xaxis()
    st.panel_letter(ax, "b", dx=-0.62)
    ax.set_title("Oprl1's position in each population\n"
                 "lower is higher-expressed", fontsize=11, pad=10)

    # (c) which GPCRs are above it, and how often
    ax = fig.add_subplot(gs[2])
    b = beaters.head(15).iloc[::-1]
    ax.barh(np.arange(len(b)), b.values, color=st.BAR_GREY, edgecolor="black",
            linewidth=1.0, height=0.72, zorder=3)
    ax.set_yticks(np.arange(len(b)))
    ax.set_yticklabels(b.index, fontstyle="italic", fontsize=10)
    ax.set_xlabel(f"Populations of {len(table)} where the GPCR "
                  "exceeds Oprl1", fontsize=10)
    ax.set_xlim(0, len(table))
    st.panel_letter(ax, "c", dx=-0.36)
    ax.set_title("The GPCRs that outrank Oprl1,\nand in how many populations",
                 fontsize=11, pad=10)

    fig.suptitle("Oprl1 among all G protein-coupled receptors, not among four "
                 "opioid receptors", fontsize=15, y=0.955)
    st.save(fig, "figure3_gpcr_rank")


# ------------------------------------------------------------------------ main
def main() -> int:
    ac.RES.mkdir(parents=True, exist_ok=True)
    symbols = gpcr_symbols()

    rows, above, panels = [], [], None
    for acc, tissue, _ in POPULATIONS:
        row, panel, detected = rank_in(acc, tissue, symbols)
        rows.append(row)
        above.append(pd.Series(1, index=detected.index[
            :int(row["rank_of_detected"]) - 1]))
        if acc == "GSE309608":
            panels = (f"{tissue}, {acc}, {row['n_cells']:,} neurons", panel,
                      int(row["rank_of_detected"]))

    table = pd.DataFrame(rows)
    ac.save_table(table, "gpcr_rank_by_population.csv")

    beaters = (pd.concat(above, axis=1).sum(axis=1).sort_values(ascending=False)
               if above else pd.Series(dtype=float))
    ac.save_table(beaters.rename("n_populations_above_Oprl1").reset_index()
                  .rename(columns={"index": "gene"}), "gpcr_above_oprl1.csv")

    med, lo, hi = bootstrap_median_rank(table.rank_of_detected)
    pmed, plo, phi = bootstrap_median_rank(table.percentile_of_detected)
    print(f"\n  Across {len(table)} populations, Oprl1's median rank among "
          f"detected GPCRs is {med:.0f} (95% interval {lo:.0f} to {hi:.0f}), "
          f"the {pmed:.1f}th percentile ({plo:.1f} to {phi:.1f}).")
    print(f"  Detected GPCRs per population: "
          f"{table.n_gpcr_detected.min()} to {table.n_gpcr_detected.max()}, "
          f"of {len(symbols)} in the list.")
    print(f"  Oprl1 is in the top 10 GPCRs in "
          f"{int((table.rank_of_detected <= 10).sum())} of {len(table)} "
          f"populations and the top 20 in "
          f"{int((table.rank_of_detected <= 20).sum())}.")
    for g in ac.RECEPTORS[1:]:
        r = table[f"{g}_rank"].dropna()
        print(f"  {g}: median rank {np.median(r):.0f} where detected "
              f"({len(r)} of {len(table)} populations)")

    print("\n  Per population:")
    print(table[["dataset", "tissue", "n_cells", "Oprl1_CPM",
                 "rank_of_detected", "n_gpcr_detected",
                 "percentile_of_detected", "transcriptome_percentile"]]
          .to_string(index=False))
    print("\n  GPCRs above Oprl1 in the most populations:")
    print(beaters.head(15).to_string())

    ac.save_table(pd.DataFrame([{
        "n_populations": len(table),
        "median_rank": med, "rank_ci_lo": lo, "rank_ci_hi": hi,
        "median_percentile": pmed, "percentile_ci_lo": plo,
        "percentile_ci_hi": phi,
        "gpcr_list_size": len(symbols),
        "in_top_10": int((table.rank_of_detected <= 10).sum()),
        "in_top_20": int((table.rank_of_detected <= 20).sum()),
        "n_boot": N_BOOT}]), "gpcr_rank_summary.csv")

    figure3(table, panels, beaters)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
