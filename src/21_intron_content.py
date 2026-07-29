"""Is it length, or is it intron content?

src/16_preparation_bias.py shows that nuclear levels rise with genomic span and
whole-cell levels do not, and attributes the difference to pre-mRNA. Genomic
span is only a proxy for that. If the mechanism is what Figure 2 says it is,
the operative variable is the intronic sequence a nucleus holds unspliced, and
exonic sequence should carry the mature message in both preparations without
adding a nuclear gain.

That is testable. This script takes the exon models of every gene in the atlas
comparison from the Ensembl annotation, computes the exonic length as the union
of exons over all transcripts of a gene, and takes intronic length as the
remainder of the span. It then regresses the two preparations, and their ratio,
on each of the three lengths.

The union over isoforms overstates the mature transcript of any single isoform,
and understates introns by the same amount, so the split is an approximation
that runs against the intronic account rather than for it.

Run from the repository root, after 16_preparation_bias.py. The annotation is
fetched into scratch/ on first use and reused.
"""

import gzip
import os
import re
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
from scipy import stats

import atlas_common as ac

SCRATCH = Path(os.environ.get("SCRATCH", "scratch"))
GTF_URL = ("https://ftp.ensembl.org/pub/release-112/gtf/mus_musculus/"
           "Mus_musculus.GRCm39.112.gtf.gz")
GENE_ID = re.compile(r'gene_id "([^"]+)"')
N_DECILE = 10


def gtf_path() -> Path:
    SCRATCH.mkdir(parents=True, exist_ok=True)
    dest = SCRATCH / Path(GTF_URL).name
    if not dest.exists():
        print(f"  fetching {Path(GTF_URL).name}")
        with urllib.request.urlopen(GTF_URL, timeout=600) as r:
            dest.write_bytes(r.read())
    return dest


def union_length(intervals) -> int:
    """Bases covered by a set of exons, counting overlap once."""
    intervals = sorted(intervals)
    total = 0
    start, end = intervals[0]
    for s, e in intervals[1:]:
        if s <= end + 1:
            end = max(end, e)
        else:
            total += end - start + 1
            start, end = s, e
    return total + end - start + 1


def exonic_lengths(wanted) -> pd.DataFrame:
    """Exonic length per gene, as the union of exons over all its transcripts."""
    exons = {}
    with gzip.open(gtf_path(), "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if f[2] != "exon":
                continue
            m = GENE_ID.search(f[8])
            if m and m.group(1) in wanted:
                exons.setdefault(m.group(1), []).append((int(f[3]), int(f[4])))
    if not exons:
        raise ac.SanityCheckError(
            "no exon records matched the atlas gene identifiers; the annotation "
            "release may not use the same gene IDs")
    return pd.DataFrame([{"ensembl_id": g, "exonic_kb": union_length(v) / 1000.0}
                         for g, v in exons.items()])


def main() -> int:
    path = ac.RES / "preparation_bias_genomewide.csv"
    if not path.exists():
        raise FileNotFoundError(
            f"{path.name} is missing; run src/16_preparation_bias.py first")
    gw = pd.read_csv(path)

    ex = exonic_lengths(set(gw.ensembl_id))
    d = gw.merge(ex, on="ensembl_id")
    missing = len(gw) - len(d)
    if missing:
        print(f"  [warn] {missing} of {len(gw)} genes have no exon model")

    # A span shorter than its own exons means the two annotations disagree about
    # the locus; the floor keeps such a gene in the table without a negative
    # intron length, and there are few enough to report rather than drop.
    d["intronic_kb"] = (d.span_kb - d.exonic_kb).clip(lower=0.001)
    d["intron_fraction"] = (d.intronic_kb / d.span_kb).round(4)
    d = d[(d.exonic_kb > 0) & (d.ratio > 0)].copy()
    ac.save_table(d.round(4), "intron_content.csv")

    lr = np.log10(d.ratio)
    logs = {k: np.log10(d[k]) for k in ("span_kb", "intronic_kb", "exonic_kb")}
    rows = []
    for name, x in logs.items():
        for quantity, y in (("nuclear", np.log10(d.nuclear_CPM.clip(lower=1e-3))),
                            ("whole cell", np.log10(d.whole_cell_CPM.clip(lower=1e-3))),
                            ("nuclear / whole cell", lr)):
            r, p = stats.pearsonr(x, y)
            rows.append({"length": name.replace("_kb", ""), "quantity": quantity,
                         "n": len(d), "pearson_r": round(float(r), 4),
                         "r_squared": round(float(r) ** 2, 4), "p": f"{p:.3g}"})

    def partial(x, y, z):
        rx = x - np.polyval(np.polyfit(z, x, 1), z)
        ry = y - np.polyval(np.polyfit(z, y, 1), z)
        return stats.pearsonr(rx, ry)

    for a, b in (("intronic_kb", "exonic_kb"), ("exonic_kb", "intronic_kb")):
        r, p = partial(logs[a], lr, logs[b])
        rows.append({"length": f"{a.replace('_kb','')} | {b.replace('_kb','')} held fixed",
                     "quantity": "nuclear / whole cell", "n": len(d),
                     "pearson_r": round(float(r), 4),
                     "r_squared": round(float(r) ** 2, 4), "p": f"{p:.3g}"})
    tests = pd.DataFrame(rows)
    ac.save_table(tests, "intron_content_tests.csv")

    # Deciles of intronic length, the same summary src/16 makes for span, so the
    # two are read the same way.
    d["decile"] = pd.qcut(d.intronic_kb, N_DECILE, labels=False, duplicates="drop")
    dec = (d.groupby("decile")
             .agg(n=("gene", "size"),
                  intronic_kb=("intronic_kb", "median"),
                  exonic_kb=("exonic_kb", "median"),
                  whole_cell_CPM=("whole_cell_CPM", "median"),
                  nuclear_CPM=("nuclear_CPM", "median"))
             .reset_index().round(4))
    ac.save_table(dec, "intron_content_by_decile.csv")

    print("\n  Correlations on log10 axes:")
    print(tests.to_string(index=False))
    print("\n  Median level by decile of intronic length:")
    print(dec.to_string(index=False))

    rec = d[d.gene.isin(ac.RECEPTORS)].set_index("gene").reindex(ac.RECEPTORS)
    print("\n  The four receptors:")
    print(rec[["span_kb", "exonic_kb", "intronic_kb", "intron_fraction",
               "whole_cell_CPM", "nuclear_CPM", "ratio"]].round(3).to_string())
    print(f"\n  Genome-wide median intron fraction: {d.intron_fraction.median():.3f}")
    by_intron = rec.intronic_kb.sort_values(ascending=False).index.tolist()
    by_exon = rec.exonic_kb.sort_values(ascending=False).index.tolist()
    by_ratio = rec.ratio.sort_values(ascending=False).index.tolist()
    print(f"  ordered by ratio        : {' > '.join(by_ratio)}")
    print(f"  ordered by intronic kb  : {' > '.join(by_intron)}"
          f"   {'(same)' if by_intron == by_ratio else '(differs)'}")
    print(f"  ordered by exonic kb    : {' > '.join(by_exon)}"
          f"   {'(same)' if by_exon == by_ratio else '(differs)'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
