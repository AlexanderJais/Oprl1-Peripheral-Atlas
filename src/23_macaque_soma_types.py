"""The primate question in tissue that was never ischaemic.

The human somata of src/22_human_soma_types.py come from ventilated organ
donors, recovered 1.3 to 6.9 hours after the ventilator. That is a good deal
better than an autopsy, but it is not a live biopsy, and every human peripheral
ganglion dataset located for this project shares the constraint: human ganglia
are not removed from healthy people. Macaque is the species where they are.

Kupari et al. (2021) took lumbar, thoracic and sacral ganglia from rhesus
monkeys under terminal anaesthesia, dissociated them and sequenced full-length
libraries from single cells: 1,040 neurons by SMART-seq2 and 2,518 on a WaferGen
chip, from five animals, typed by the authors on the scheme the mouse atlases
use. Whole cell, so no nuclear pre-mRNA; no ischaemic interval; and an
independent laboratory from the human atlas.

It buys that at the price of dissociation, which is the mouse preparation's own
problem and which Figure S2 shows is separable from gene length. And it is
poorly sampled where the mouse answer lives: enzymatic dissociation loses large
myelinated neurons, and the two deposits together hold 38 A-LTMRs against 3,188
nociceptors, with no proprioceptor class at all. Where the human atlas is 33.6%
myelinated and the mouse atlas 13.4%, these deposits are 2.1% and 0.6%. The
A-fibre contrast is therefore reported with its sample size beside it and should
be read as underpowered rather than as a null.

That composition also governs the pooled ordering, which is why the pooled
ordering is reported here and not leaned on. A ganglion-wide average is a
statement about which neurons a protocol captured as much as about the species.

Counts from a full-length protocol carry a transcript-length term, so levels are
divided by exonic length from the NCBI Mmul_10 annotation, which is the
annotation the deposit was quantified against. Mmul_10 is a sparser annotation
than GRCh38 -- it gives OPRM1 7.3 kb of exon against the human 21.0 kb, mostly
in untranslated sequence -- so TPM values are comparable within this script and
not across species.

Run from the repository root, after 22_human_soma_types.py. Sources are fetched
into scratch/macaque/ on first use; the annotation is ~24 MB.
"""

import gzip
import re
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
from scipy import stats

import atlas_common as ac

SCRATCH = Path("scratch") / "macaque"
SS2 = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165553/suppl"
WG = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165566/suppl"
GTF_URL = ("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/339/765/"
           "GCF_003339765.1_Mmul_10/GCF_003339765.1_Mmul_10_genomic.gtf.gz")

SS2_PLATES = [
    "GSE165553_SS2-16-225_expression_for_R.tab.gz",
    "GSE165553_SS2-16-227_expression_for_R.tab.gz",
    "GSE165553_SS2-16-259_expression_for_R.tab.gz",
    "GSE165553_SS2-16-303_expression_for_R.tab.gz",
    "GSE165553_SS2-16-304_expression_for_R.tab.gz",
    "GSE165553_SS2-16-305_expression_for_R_from_comb_plate.tab.gz",
    "GSE165553_SS2-16-414_expression_for_R_from_comb_plate.tab.gz",
    "GSE165553_SS2-16-415_expression_for_R_from_comb_plate.tab.gz",
    "GSE165553_SS2-16-416_expression_for_R_from_comb_plate.tab.gz",
]
WG_PLATES = [
    "GSE165566_WG17019_expression_for_R.tab.gz",
    "GSE165566_WG17020_expression_for_R.tab.gz",
    "GSE165566_WG18008_expression_for_R.tab.gz",
]

RECEPTORS = ["OPRL1", "OPRM1", "OPRD1", "OPRK1"]
NOT_NEURONS = ["lowq/nonneuron"]
# Kupari et al.'s labels. A-LTMR and TrpM8high are the only non-nociceptive
# types they resolve; there is no proprioceptor class in either deposit.
MYELINATED = ["A-LTMR"]
NOCICEPTIVE = ["PEP1", "PEP2", "PEP3", "NP1", "NP2", "NP3"]
UNGROUPED = ["C-LTMR", "TrpM8high"]

GENE_NAME = re.compile(r'gene_id "([^"]+)"')
PEER_BAND = 1.6
MIN_GROUP = 10  # below this a per-type mean is a handful of cells, not a level


def fetch(url: str, name: str) -> Path:
    SCRATCH.mkdir(parents=True, exist_ok=True)
    dest = SCRATCH / name
    if not dest.exists():
        print(f"  fetching {name}")
        with urllib.request.urlopen(url, timeout=900) as r:
            dest.write_bytes(r.read())
    return dest


def union_length(intervals) -> int:
    intervals = sorted(intervals)
    total, (start, end) = 0, intervals[0]
    for s, e in intervals[1:]:
        if s <= end + 1:
            end = max(end, e)
        else:
            total += end - start + 1
            start, end = s, e
    return total + end - start + 1


def exonic_kb() -> pd.Series:
    """Exonic length per gene symbol, from the annotation the deposit used.

    The matrices carry NCBI symbols including LOC identifiers, so the RefSeq
    annotation joins on more of the transcriptome than Ensembl's would.
    """
    path = SCRATCH / "exonic_kb.csv"
    if path.exists():
        return pd.read_csv(path, index_col=0).exonic_kb
    exons: dict[str, list[tuple[int, int]]] = {}
    with gzip.open(fetch(GTF_URL, "Mmul_10.gtf.gz"), "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if f[2] != "exon":
                continue
            m = GENE_NAME.search(f[8])
            if m:
                exons.setdefault(m.group(1), []).append((int(f[3]), int(f[4])))
    out = pd.Series({g: union_length(v) / 1000.0 for g, v in exons.items()},
                    name="exonic_kb")
    out.to_frame().to_csv(path)
    return out


def read_plates(base: str, names) -> pd.DataFrame:
    """One counts matrix per deposit, genes by cells.

    Each plate is deposited separately and the plates do not carry identical
    gene lists, so they are joined on the gene index rather than concatenated
    blind; a gene absent from a plate is absent, not zero, and is dropped.
    """
    frames = []
    for n in names:
        d = pd.read_csv(fetch(f"{base}/{n}", n), sep="\t", index_col=0)
        d = d[~d.index.duplicated()]
        frames.append(d)
    joined = pd.concat(frames, axis=1, join="inner")
    lost = max(len(f) for f in frames) - len(joined)
    if lost:
        print(f"    {lost} genes are not quantified on every plate and are dropped")
    return joined


def prepare(base: str, plates, meta: pd.DataFrame, label: str):
    """One deposit's neurons, joined to their type labels.

    The two files write the same cell with different punctuation: the matrices
    have SS2-16-225_A01 and WG17019_A01-W01, the metadata SS2.16.225_A01 and
    WG17019_A01.W03. Joining them raw silently keeps only the four plates whose
    names happen to agree, so the separators are normalised and the join is
    checked rather than assumed.
    """
    counts = read_plates(base, plates)
    counts.columns = [c.replace("-", ".") for c in counts.columns]
    shared = [c for c in counts.columns if c in set(meta.index)]
    if len(shared) < 0.9 * len(counts.columns):
        raise ac.SanityCheckError(
            f"{label}: only {len(shared)} of {len(counts.columns)} matrix "
            "columns appear in the metadata; the cell identifiers do not join")
    meta = meta.loc[shared]
    counts = counts[shared]
    meta = meta[~meta["Cell.Identity"].isin(NOT_NEURONS)]
    counts = counts[meta.index]
    for g in RECEPTORS:
        if g not in counts.index:
            raise ac.SanityCheckError(f"{label}: {g} absent from the matrix")
    return counts, meta


def to_tpm(counts: pd.DataFrame, kb: pd.Series) -> pd.DataFrame:
    length = kb.reindex(counts.index)
    modelled = length.notna() & (length > 0)
    print(f"    {int(modelled.sum())} of {len(counts)} rows carry an exon model")
    rate = counts[modelled].div(length[modelled], axis=0)
    return rate.div(rate.sum(axis=0), axis=1) * 1e6


def by_type(tpm, counts, meta, label) -> pd.DataFrame:
    rows = []
    for t, cells in meta.groupby("Cell.Identity").groups.items():
        for g in RECEPTORS:
            rows.append({"deposit": label, "type": t, "n_cells": len(cells),
                         "gene": g,
                         "TPM_mean": tpm.loc[g, cells].mean(),
                         "pct_detected": 100.0 * (counts.loc[g, cells] > 0).mean()})
    return pd.DataFrame(rows)


def class_contrast(tpm, meta, label) -> pd.DataFrame:
    a = meta.index[meta["Cell.Identity"].isin(MYELINATED)]
    b = meta.index[meta["Cell.Identity"].isin(NOCICEPTIVE)]
    rows = []
    for g in RECEPTORS:
        u, p = stats.mannwhitneyu(tpm.loc[g, a], tpm.loc[g, b],
                                  alternative="two-sided")
        rows.append({"deposit": label, "gene": g, "n_myelinated": len(a),
                     "n_nociceptive": len(b),
                     "myelinated_TPM": tpm.loc[g, a].mean(),
                     "nociceptive_TPM": tpm.loc[g, b].mean(),
                     "ratio": (tpm.loc[g, a].mean() / tpm.loc[g, b].mean()
                               if tpm.loc[g, b].mean() else np.nan),
                     "p": f"{p:.3g}"})
    return pd.DataFrame(rows)


def peer_percentiles(tpm, counts, kb, label) -> pd.DataFrame:
    prof = pd.DataFrame({"exonic_kb": kb.reindex(tpm.index),
                         "TPM": tpm.mean(axis=1),
                         "pct": 100.0 * (counts.reindex(tpm.index) > 0).mean(axis=1)})
    prof = prof[prof.TPM > 0]
    rows = []
    for g in RECEPTORS:
        e = prof.exonic_kb[g]
        band = prof[(prof.exonic_kb >= e / PEER_BAND)
                    & (prof.exonic_kb <= e * PEER_BAND)]
        rows.append({"deposit": label, "gene": g, "exonic_kb": e,
                     "n_peers": len(band), "TPM": prof.TPM[g],
                     "TPM_percentile": stats.percentileofscore(band.TPM, prof.TPM[g]),
                     "pct_detected": prof.pct[g],
                     "detection_percentile": stats.percentileofscore(band.pct,
                                                                     prof.pct[g])})
    return pd.DataFrame(rows)


def main() -> int:
    kb = exonic_kb()

    ss2_meta = pd.read_csv(
        fetch(f"{SS2}/GSE165553_SS2_SmartSeq2_all_metadata.tab.gz",
              "SS2_metadata.tab.gz"), sep="\t").set_index("Cell.ID")
    wg_meta = pd.read_csv(
        fetch(f"{WG}/GSE165566_WG_all_metadata.tab.gz", "WG_metadata.tab.gz"),
        sep="\t").rename(columns={"cell": "Cell.ID"}).set_index("Cell.ID")

    deposits = {}
    for label, base, plates, meta in (("GSE165553 SMART-seq2", SS2, SS2_PLATES, ss2_meta),
                                      ("GSE165566 WaferGen", WG, WG_PLATES, wg_meta)):
        print(f"  {label}")
        counts, m = prepare(base, plates, meta, label)
        tpm = to_tpm(counts, kb)
        cpm = counts.div(counts.sum(axis=0), axis=1) * 1e6
        myelin = m["Cell.Identity"].isin(MYELINATED).sum()
        print(f"    {len(m)} neurons, {m['Cell.Identity'].nunique()} types, "
              f"{myelin} myelinated ({100.0 * myelin / len(m):.1f}%)")
        deposits[label] = (counts, cpm, tpm, m)

    types = pd.concat([by_type(t, c, m, k) for k, (c, _, t, m) in deposits.items()])
    ac.save_table(types.round(4), "macaque_receptors_by_type.csv")

    pooled = pd.concat([
        pd.DataFrame({"deposit": k, "gene": RECEPTORS,
                      "n_neurons": len(m),
                      "exonic_kb": kb.reindex(RECEPTORS).values,
                      "CPM_mean": cpm.loc[RECEPTORS].mean(axis=1).values,
                      "TPM_mean": t.loc[RECEPTORS].mean(axis=1).values,
                      "pct_detected": (100.0 * (c.loc[RECEPTORS] > 0)
                                       .mean(axis=1)).values})
        for k, (c, cpm, t, m) in deposits.items()])
    ac.save_table(pooled.round(4), "macaque_receptors_pooled.csv")

    contrast = pd.concat([class_contrast(t, m, k)
                          for k, (_, _, t, m) in deposits.items()])
    ac.save_table(contrast.round(4), "macaque_class_contrast.csv")

    peers = pd.concat([peer_percentiles(t, c, kb, k)
                       for k, (c, _, t, m) in deposits.items()])
    ac.save_table(peers.round(3), "macaque_peer_percentiles.csv")

    for k, (c, _, t, m) in deposits.items():
        print(f"\n  {k}: mean TPM by type "
              f"(types under {MIN_GROUP} cells marked)")
        sub = types[types.deposit == k]
        piv = sub.pivot(index="type", columns="gene", values="TPM_mean")[RECEPTORS]
        piv.insert(0, "n", sub.groupby("type").n_cells.first())
        piv["top"] = piv[RECEPTORS].idxmax(axis=1)
        piv["top"] = np.where(piv.n < MIN_GROUP, piv.top + " (n<10)", piv.top)
        print(piv.sort_values("n", ascending=False).round(3).to_string())

    print("\n  Pooled over neurons:")
    print(pooled.round(3).to_string(index=False))
    print("\n  A-LTMR against nociceptive types:")
    print(contrast.to_string(index=False))
    print("\n  Against genes of the same exonic length:")
    print(peers.round(2).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
