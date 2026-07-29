"""Dissociation is a second preparation effect, and it is separable from length.

Figure 2 attributes the preparation-dependent receptor reversal to gene length.
The competing account is dissociation: whole-cell libraries are made from live
enzymatic dissociation, which axotomises every neuron and induces an immediate-
early and injury programme, while nuclei here come from tissue that was never
dissociated alive. Axotomy lowers opioid responsiveness and raises nociceptin
responsiveness in dorsal root ganglion neurons (Abdulla and Smith, 1998), which
is the direction that would produce this atlas's result artefactually.

The two accounts are separable on the same table Figure 2 is drawn from. A
dissociation-induced gene should be enriched in whole cell far beyond what its
length predicts; a length artefact should scale with span and have nothing to do
with induction. This script scores both sets of genes the same way
16_preparation_bias.py scores the receptors: the nuclear-to-whole-cell ratio,
against the ratio of every gene within a factor of 1.6 of the same span.

The result constrains the manuscript in both directions and is reported as such.
The dissociation signature is real and large, so no claim here may rest on the
whole-cell libraries being unperturbed; and the receptors sit nowhere near it,
so it does not account for the reversal that Figure 2 attributes to length.

Run from the repository root, after 16_preparation_bias.py.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import atlas_common as ac

# Immediate-early and injury genes. Fos, Fosb, Jun, Junb, Jund, Egr1, Arc and
# Npas4 are the activity-induced set that dissociation switches on; Hspa1a and
# Hspa1b report the heat-shock arm of the same handling; Atf3 and Socs3 are the
# axotomy markers, and Atf3 is the one that speaks directly to the competing
# account, since it is induced by the nerve injury dissociation inflicts.
DISSOCIATION_GENES = ["Fos", "Fosb", "Jun", "Junb", "Jund", "Egr1", "Arc",
                      "Npas4", "Atf3", "Socs3", "Hspa1a", "Hspa1b"]
SPAN_BAND = 1.6           # the peer definition used for the receptors


def peer_scores(d: pd.DataFrame, genes) -> pd.DataFrame:
    rows = []
    for g in genes:
        hit = d[d.gene == g]
        if hit.empty:
            print(f"  [warn] {g} is not in the genome-wide table")
            continue
        r = hit.iloc[0]
        peers = d[(d.span_kb >= r.span_kb / SPAN_BAND)
                  & (d.span_kb <= r.span_kb * SPAN_BAND)]
        rows.append({
            "gene": g, "span_kb": round(float(r.span_kb), 3),
            "whole_cell_CPM": round(float(r.whole_cell_CPM), 3),
            "nuclear_CPM": round(float(r.nuclear_CPM), 3),
            "ratio": round(float(r.ratio), 4),
            "n_peers": int(len(peers)),
            "peer_median_ratio": round(float(peers.ratio.median()), 3),
            "pct_of_peers": round(float((peers.ratio < r.ratio).mean() * 100), 1),
        })
    return pd.DataFrame(rows)


def main() -> int:
    path = ac.RES / "preparation_bias_genomewide.csv"
    if not path.exists():
        raise FileNotFoundError(
            f"{path.name} is missing; run src/16_preparation_bias.py first")
    d = pd.read_csv(path)

    tbl = peer_scores(d, DISSOCIATION_GENES)
    tbl["set"] = "dissociation"
    rec = peer_scores(d, ac.RECEPTORS)
    rec["set"] = "receptor"
    out = pd.concat([tbl, rec], ignore_index=True)[
        ["set", "gene", "span_kb", "whole_cell_CPM", "nuclear_CPM", "ratio",
         "n_peers", "peer_median_ratio", "pct_of_peers"]]
    ac.save_table(out, "dissociation_signature.csv")

    print("\n  Nuclear over whole cell, against genes of the same length:")
    print(out.to_string(index=False))

    dis, recs = tbl.pct_of_peers, rec.pct_of_peers
    print(f"\n  The {len(tbl)} dissociation genes sit at a median "
          f"{dis.median():.1f}th percentile of their own length class, "
          f"{int((dis < 10).sum())} of them below the 10th.")
    print(f"  The four receptors sit at {', '.join(f'{g} {p:.1f}' for g, p in zip(rec.gene, recs))}.")
    below = int((dis < recs.min()).sum())
    print(f"  {below} of {len(tbl)} dissociation genes fall below the lowest of "
          f"the four receptors; the exception is "
          f"{', '.join(tbl.loc[dis >= recs.min(), 'gene'])}.")

    fos = tbl[tbl.gene == "Fos"].iloc[0]
    atf3 = tbl[tbl.gene == "Atf3"].iloc[0]
    print(f"\n  Fos reads {fos.whole_cell_CPM:,.1f} CPM in whole cells against "
          f"{fos.nuclear_CPM:.2f} in nuclei, and Atf3, the axotomy marker, "
          f"{atf3.whole_cell_CPM:,.1f} against {atf3.nuclear_CPM:.2f}. The "
          "whole-cell libraries carry an injury programme, and no claim in this "
          "atlas may assume otherwise.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
