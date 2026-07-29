"""A test at the altitude the claim is made.

The cell bootstrap elsewhere in this project resamples the cells of one
population. That is the right question for one panel and the wrong unit for the
claim the atlas makes: cells within a ganglion are not independent replicates of
anything, and four of the nineteen populations come from a single animal. An
interval of 1.15 to 1.28 over 26,047 nodose cells says how much the ordering
depends on which cells were captured, not on which mouse, laboratory or
platform produced them.

The claim is that Oprl1 leads across the peripheral nervous system. Its unit is
therefore the population, and more conservatively the deposit, since the five
autonomic ganglia of one experiment are not independent of each other. This
script counts how many place Oprl1 first and tests that count against two nulls:

  p = 0.25   the four receptors are exchangeable, so any one of them leads a
             given population by chance
  p = 0.5    Oprl1 against the field, which assumes nothing about how the other
             three divide the remainder and is the conservative choice

A deposit counts as leading only if Oprl1 is first in every population it
contributes, so GSE263422 is scored against the atlas on the strength of its P7
enteric sample.

Run from the repository root.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd
from scipy import stats

import atlas_common as ac

NULLS = {0.25: "four receptors exchangeable", 0.50: "Oprl1 against the field"}


def main() -> int:
    fig = pd.read_csv(ac.RES / "receptor_levels_by_ganglion.csv")

    units = {
        "population": fig.assign(first=fig.top == "Oprl1")[["panel", "first"]],
        "deposit": (fig.assign(first=fig.top == "Oprl1")
                    .groupby("dataset")["first"].all().reset_index()),
    }

    rows = []
    for unit, d in units.items():
        k, n = int(d.first.sum()), len(d)
        for null, label in NULLS.items():
            t = stats.binomtest(k, n, null, alternative="greater")
            rows.append({"unit": unit, "k": k, "n": n,
                         "null_p": null, "null": label,
                         "p_value": float(t.pvalue)})
        print(f"  {unit:11} {k} of {n} place Oprl1 first"
              + ("" if unit == "population" else
                 f"; mixed: {', '.join(d.loc[~d.first, 'dataset'])}"))

    out = pd.DataFrame(rows)
    out["p_value"] = out.p_value.map(lambda v: float(f"{v:.3g}"))
    ac.save_table(out, "receptor_ordering_sign_test.csv")
    print()
    print(out.to_string(index=False))
    print("\n  The deposit-level test on the conservative null is the one to "
          "quote: it treats the five ganglia of GSE232789 as a single "
          "observation and assumes nothing about how Oprm1, Oprd1 and Oprk1 "
          "divide the populations Oprl1 does not lead.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
