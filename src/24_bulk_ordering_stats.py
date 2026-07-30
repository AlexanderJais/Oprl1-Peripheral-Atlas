"""The ordering statistics Figure 3 and its Results paragraph quote.

src/20_bulk_composition_screen.py screens 223 GEO series and records, for every
column group, whether it is neuron-rich enough and whether its unit permits one
gene to be compared with another. This script takes the groups that clear both
gates and asks the question the screen was built for: over independent deposits,
where does each of the four opioid receptors fall?

Three statistics, in increasing conservatism.

  group level     one count per column group. 52 groups, but a deposit that
                  publishes four treatment arms contributes four of them, so
                  this over-counts a deposit's evidence.
  deposit level   one count per series, and a series counts for Oprl1 only if
                  Oprl1 leads every group it contributes. This is the unit the
                  claim is about.
  rank profile    how often each receptor takes each of the four positions.
                  Whether a receptor is ever last matters as much as how often
                  it is first: a receptor can win the most groups by being
                  enormous in a few and absent from the rest, which is what
                  Oprm1 does.

The breadth statistic comes last, on GSE131230, the one deposit that sequenced
purified subtypes rather than whole ganglion. Breadth is the fold range of a
receptor across the eight subtypes: it is what "broadly expressed" means, and it
separates the four receptors more sharply than the mean rank does.

Levels here are as the depositors normalised them, which is a length-normalised
unit by the screen's own gate. They are compared between genes within a group,
never between groups of different deposits.

Run from the repository root, after 20_bulk_composition_screen.py.
"""

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
from scipy import stats

import atlas_common as ac

# The screen gates on neuronal content, which a brain sample also clears, so
# tissue identity is applied here rather than there; paper/figures.py figure3()
# uses the same pattern and the two must not drift apart.
NOT_PERIPHERAL = re.compile(r"cortex|hippocamp|brain|spinal|striat|"
                            r"arcuate|\bArc|hypothal", re.I)
SUBTYPE_SERIES = "GSE131230"
SUBTYPES = ["Nonpeptidergic Nociceptor", "Peptidergic Nociceptor", "C-LTMR",
            "Aδ-LTMR", "Aβ RA-LTMR", "Aβ SA1-LTMR", "Aβ Field-LTMR",
            "Proprioceptor"]
NULL = 0.25  # four receptors, so a coin with four faces


def peripheral_read(scr: pd.DataFrame) -> pd.DataFrame:
    """The groups Figure 3 plots: gate passed, peripheral, one copy each."""
    read = scr[(scr.verdict == "read")
               & ~scr.group.str.contains(NOT_PERIPHERAL)].copy()
    # A deposit that publishes FPKM and counts of the same samples appears
    # twice; the length-normalised copy is the one that clears the unit gate,
    # but the dedup is explicit so a future screen cannot double-count.
    read["set"] = read.gse + "|" + read.group.str.replace(
        r"^(FPKM|count|TPM)\.", "", regex=True)
    read = read.sort_values("unit").drop_duplicates("set")
    if read.empty:
        raise ac.SanityCheckError(
            "no peripheral group clears both gates; the screen may have failed")
    return read


def main() -> int:
    scr = pd.read_csv(ac.RES / "bulk_composition_screen.csv")
    read = peripheral_read(scr)

    rows = []
    k, n = int((read.top_receptor == "Oprl1").sum()), len(read)
    rows.append({"level": "column group", "unit": "group", "n": n,
                 "Oprl1_first": k,
                 "p": f"{stats.binomtest(k, n, NULL, alternative='greater').pvalue:.3g}"})
    dep = read.groupby("gse").top_receptor.agg(
        lambda s: "Oprl1" if (s == "Oprl1").all() else "other")
    kd, nd = int((dep == "Oprl1").sum()), len(dep)
    rows.append({"level": "deposit, leading every group it contributes",
                 "unit": "series", "n": nd, "Oprl1_first": kd,
                 "p": f"{stats.binomtest(kd, nd, NULL, alternative='greater').pvalue:.3g}"})
    counts = pd.DataFrame(rows)
    ac.save_table(counts, "bulk_ordering_counts.csv")

    rk = read[ac.RECEPTORS].rank(axis=1, ascending=False)
    profile = pd.DataFrame({
        "gene": ac.RECEPTORS,
        "first": [int((rk[g] == 1).sum()) for g in ac.RECEPTORS],
        "second": [int((rk[g] == 2).sum()) for g in ac.RECEPTORS],
        "third": [int((rk[g] == 3).sum()) for g in ac.RECEPTORS],
        "last": [int((rk[g] == 4).sum()) for g in ac.RECEPTORS],
        "mean_rank": [round(float(rk[g].mean()), 3) for g in ac.RECEPTORS],
        "median_level": [round(float(read[g].median()), 3) for g in ac.RECEPTORS],
    })
    ac.save_table(profile, "bulk_ordering_rank_profile.csv")

    sub = scr[(scr.gse == SUBTYPE_SERIES) & scr.group.isin(SUBTYPES)]
    sub = sub.set_index("group").reindex([s for s in SUBTYPES
                                          if s in set(sub.group)])
    if len(sub) != len(SUBTYPES):
        raise ac.SanityCheckError(
            f"{SUBTYPE_SERIES} contributes {len(sub)} of {len(SUBTYPES)} "
            "purified subtypes; the breadth statistic would be computed on a "
            "different set than the one it is described on")
    srk = sub[ac.RECEPTORS].rank(axis=1, ascending=False)
    breadth = pd.DataFrame({
        "gene": ac.RECEPTORS,
        "subtypes_first": [int((srk[g] == 1).sum()) for g in ac.RECEPTORS],
        "mean_rank": [round(float(srk[g].mean()), 3) for g in ac.RECEPTORS],
        "min_level": [round(float(sub[g].min()), 3) for g in ac.RECEPTORS],
        "max_level": [round(float(sub[g].max()), 3) for g in ac.RECEPTORS],
        "fold_range": [round(float(sub[g].max() / sub[g].min()), 1)
                       if sub[g].min() > 0 else np.inf for g in ac.RECEPTORS],
    })
    ac.save_table(breadth, "bulk_subtype_breadth.csv")

    axo = pd.read_csv(ac.RES / "bulk_axotomy_pairs.csv")
    inj = []
    for g in ["Atf3"] + ac.RECEPTORS:
        s = axo[axo.gene == g]
        fc = s.injury_CPM / s.control_CPM.replace(0, np.nan)
        inj.append({"gene": g, "pairs": len(s),
                    "rose": int((s.injury_CPM > s.control_CPM).sum()),
                    "median_fold_change": round(float(fc.median()), 3)})
    injury = pd.DataFrame(inj)
    ac.save_table(injury, "bulk_injury_direction.csv")

    print("  Oprl1 first, over independent deposits:")
    print(counts.to_string(index=False))
    print("\n  Where each receptor falls, over the "
          f"{len(read)} groups of {read.gse.nunique()} series:")
    print(profile.to_string(index=False))
    print(f"\n  Median margin over the runner-up where Oprl1 leads: "
          f"{read.loc[read.top_receptor == 'Oprl1', 'margin'].median():.2f}-fold")
    print(f"\n  Breadth across the eight purified subtypes of {SUBTYPE_SERIES}:")
    print(breadth.to_string(index=False))
    print("\n  Matched control against nerve injury:")
    print(injury.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
