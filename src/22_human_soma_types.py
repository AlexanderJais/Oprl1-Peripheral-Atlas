"""Does the mouse cell-type partition of the opioid receptors hold in human?

src/15_human_ganglia.py could not answer the human question, because every human
peripheral ganglion dataset it found was single-nucleus, and Figure 2 shows what
that preparation does to a 7 kb gene sitting beside a 280 kb one. Yu et al.
(2024) removes that objection: 1,136 human dorsal root ganglion somata, cut out
of the tissue one at a time by laser capture and sequenced full-length. A soma
is a whole cell, so there is no nuclear enrichment of pre-mRNA.

It does not remove every length term. A full-length protocol reported as raw
counts assigns reads in proportion to transcript length, so a long mRNA earns
more counts than a short one at the same molar abundance. That is the reason
TPM exists, and this script divides by the exonic length of each gene before
comparing genes with each other. Pooled over all somata the ordering is the same
either way, so the correction changes the size of the gap and not its sign:

                    raw CPM   per kb of exon
      OPRM1           73.3         8.74
      OPRD1           16.3         4.54
      OPRK1            4.4         1.92
      OPRL1            0.9         0.61

That is the opposite of the mouse ordering, and taken alone it says the mouse
result does not transfer. Taken alone is the problem. In mouse, Oprl1's lead is
not spread evenly over the ganglion: it is carried by the myelinated NF classes,
where Oprl1 sits at 17.1 CPM against 4.0 for Oprm1, while Oprm1 leads in the
peptidergic and SST nociceptors. A pooled average over a sample weighted toward
nociceptors would understate Oprl1 in mouse too.

Yu et al. assign each of their somata to one of 16 neuronal types, and the
assignments are in the GEO sample records, so the same split can be made in
human. This script makes it. It reads the per-cell type labels and soma
diameters from the series matrix, converts the count matrix to TPM against the
Ensembl exon models, and reports each receptor by type.

The answer is that the cell-type partition replicates and the ordering does not.
OPRL1 is the only one of the four enriched in the A-fibre low-threshold
mechanoreceptor and proprioceptive types over the nociceptive types; the other
three are depleted there, OPRM1 by twelve-fold. OPRL1 is also the only one whose
expressing somata are larger than its non-expressing somata. But it reaches
those types at roughly 1 TPM, where OPRM1 reaches nociceptors at 13, so the
pooled ordering still places OPRL1 last.

The tissue is worth stating precisely, because src/23 turns on it. The three
donors were ventilated organ donors, and the ganglia were recovered 1.3, 6.9 and
1.3 hours after the ventilator, then frozen (Yu et al., Supplementary Table 1).
That is far from an autopsy interval, but it is not a live biopsy either, and
the somata carry a stress signature: ATF3 sits at 25 TPM and is detected in half
of them even outside the 23-cell cluster the authors themselves label hATF3.

Run from the repository root. Sources are fetched into scratch/human/ on first
use and reused. The GRCh38 annotation is ~55 MB.
"""

import csv
import gzip
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
from scipy import stats

import atlas_common as ac

SCRATCH = Path("scratch") / "human" / "soma"
GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE249nnn/GSE249746"
COUNTS = f"{GEO}/suppl/GSE249746_Expression_matrix_raw_counts.csv.gz"
SERIES = f"{GEO}/matrix/GSE249746_series_matrix.txt.gz"
GTF_URL = ("https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/"
           "Homo_sapiens.GRCh38.112.gtf.gz")

RECEPTORS = ["OPRL1", "OPRM1", "OPRD1", "OPRK1"]

# Yu et al.'s own labels. hAb.LTMR and hAd.LTMR are the myelinated
# low-threshold mechanoreceptors and hPropr the proprioceptors; these are the
# types that correspond to the mouse NF classes where Oprl1 leads. Glia.cont is
# their label for a soma whose profile is glial contamination.
LTMR_PROPRIOCEPTIVE = ["hAb.LTMR", "hAd.LTMR", "hPropr"]
NOCICEPTIVE = ["hNP1", "hNP2", "hPEP.0", "hPEP.CHRNA7", "hPEP.KIT", "hPEP.NTRK3",
               "hPEP.PIEZOh", "hPEP.SST", "hPEP.TRPV1/A1.1", "hPEP.TRPV1/A1.2",
               "hTRPM8"]
NOT_NEURONS = ["Glia.cont", ""]
# hC.LTMR and hATF3 are in neither group: the C-LTMRs are unmyelinated, so they
# do not belong with the A-fibre classes, and hATF3 is a stress signature rather
# than a modality.

# The mouse counterpart, from results/drg_oprl1_by_subtype.csv, so the two
# species are contrasted on the same grouping rather than by eye.
MOUSE_MYELINATED = ["NF1", "NF2", "NF3"]
MOUSE_NOCICEPTIVE = ["PEP1", "PEP2", "NP", "SST"]

PEER_BAND = 1.6  # a gene's length class, as a multiplicative half-width


def fetch(url: str, name: str) -> Path:
    SCRATCH.mkdir(parents=True, exist_ok=True)
    dest = SCRATCH / name
    if not dest.exists():
        print(f"  fetching {name}")
        with urllib.request.urlopen(url, timeout=900) as r:
            dest.write_bytes(r.read())
    return dest


def cell_metadata() -> pd.DataFrame:
    """Per-soma type, donor, spinal level and diameter, from the GEO records.

    The count matrix carries no annotation at all; the labels live in the
    per-sample characteristics of the series matrix, one row per field.
    """
    fields = {}
    with gzip.open(fetch(SERIES, "series_matrix.txt.gz"), "rt") as fh:
        for line in fh:
            if not line.startswith("!Sample_"):
                continue
            f = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
            if f[0] == "!Sample_title":
                fields["cell"] = f[1:]
            elif f[0] == "!Sample_characteristics_ch1":
                tag = f[1].split(":")[0].strip()
                fields[tag] = [v.split(":", 1)[1].strip() if ":" in v else ""
                               for v in f[1:]]
    for need in ("cell", "cell identity", "dissected soma size"):
        if need not in fields:
            raise ac.SanityCheckError(
                f"the series matrix carries no '{need}' field; GEO may have "
                "re-deposited the records in another layout")
    m = pd.DataFrame({
        "cell": fields["cell"],
        "type": fields["cell identity"],
        "donor": fields["donor id"],
        "region": fields["region"],
        "soma_um": pd.to_numeric(fields["dissected soma size"], errors="coerce"),
    }).set_index("cell")
    unknown = set(m.type) - set(LTMR_PROPRIOCEPTIVE) - set(NOCICEPTIVE) \
        - set(NOT_NEURONS) - {"hC.LTMR", "hATF3"}
    if unknown:
        raise ac.SanityCheckError(
            f"unrecognised cell identities {sorted(unknown)}; the grouping in "
            "this script would silently drop them")
    return m


def exonic_kb() -> pd.Series:
    """Exonic length per gene symbol, as the union of exons over transcripts."""
    path = SCRATCH / "exonic_kb.csv"
    if path.exists():
        return pd.read_csv(path, index_col=0).exonic_kb
    exons: dict[str, list[tuple[int, int]]] = {}
    with gzip.open(fetch(GTF_URL, "GRCh38.112.gtf.gz"), "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if f[2] != "exon":
                continue
            attrs = f[8]
            start = attrs.find('gene_name "')
            if start < 0:
                continue
            name = attrs[start + 11:attrs.index('"', start + 11)]
            exons.setdefault(name, []).append((int(f[3]), int(f[4])))
    out = pd.Series({g: union_length(v) / 1000.0 for g, v in exons.items()},
                    name="exonic_kb")
    out.to_frame().to_csv(path)
    return out


def union_length(intervals) -> int:
    """Bases covered by a set of exons, counting overlap once."""
    intervals = sorted(intervals)
    total, (start, end) = 0, intervals[0]
    for s, e in intervals[1:]:
        if s <= end + 1:
            end = max(end, e)
        else:
            total += end - start + 1
            start, end = s, e
    return total + end - start + 1


def to_tpm(counts: pd.DataFrame, kb: pd.Series):
    """Counts to TPM, and the rows that carry an exon model.

    A full-length library reported as counts carries a transcript-length term.
    Dividing by exonic length and renormalising removes it, which is what makes
    one gene comparable with another rather than with itself across cells.
    """
    length = kb.reindex(counts.index)
    modelled = length.notna() & (length > 0)
    rate = counts[modelled].div(length[modelled], axis=0)
    return rate.div(rate.sum(axis=0), axis=1) * 1e6, modelled


def class_contrast(tpm, counts, meta) -> pd.DataFrame:
    """Each receptor in the A-fibre LTMR and proprioceptive types against the
    nociceptive types, with the mouse value from the same grouping beside it."""
    mouse = pd.read_csv(ac.RES / "drg_oprl1_by_subtype.csv").set_index("subtype")

    def weighted(frame, types):
        cells = meta.index[meta.type.isin(types)]
        return tpm.loc[RECEPTORS, cells].mean(axis=1), cells

    lt, lt_cells = weighted(tpm, LTMR_PROPRIOCEPTIVE)
    no, no_cells = weighted(tpm, NOCICEPTIVE)
    rows = []
    for g in RECEPTORS:
        u, p = stats.mannwhitneyu(tpm.loc[g, lt_cells], tpm.loc[g, no_cells],
                                  alternative="two-sided")
        mg = g.capitalize()
        m_lt = m_no = np.nan
        if mg in mouse.columns:
            n = mouse.n
            m_lt = float((mouse.loc[MOUSE_MYELINATED, mg]
                          * n[MOUSE_MYELINATED]).sum() / n[MOUSE_MYELINATED].sum())
            m_no = float((mouse.loc[MOUSE_NOCICEPTIVE, mg]
                          * n[MOUSE_NOCICEPTIVE]).sum() / n[MOUSE_NOCICEPTIVE].sum())
        rows.append({
            "gene": g,
            "human_LTMR_proprioceptive_TPM": lt[g],
            "human_nociceptive_TPM": no[g],
            "human_ratio": lt[g] / no[g] if no[g] else np.nan,
            "human_p": f"{p:.3g}",
            "mouse_myelinated_CPM": m_lt,
            "mouse_nociceptive_CPM": m_no,
            "mouse_ratio": m_lt / m_no if m_no else np.nan,
        })
    return pd.DataFrame(rows)


def soma_size(tpm, counts, meta) -> pd.DataFrame:
    """Is a soma that expresses the receptor larger than one that does not?

    Diameter is the one property of these cells measured before sequencing, by
    the microdissection itself, so it is independent of the library.
    """
    rows = []
    for g in RECEPTORS:
        v = counts.loc[g, meta.index]
        pos = meta.soma_um[v.values > 0].dropna()
        neg = meta.soma_um[v.values == 0].dropna()
        u, p = stats.mannwhitneyu(pos, neg, alternative="two-sided")
        rows.append({"gene": g, "n_expressing": len(pos),
                     "median_um_expressing": pos.median(),
                     "n_silent": len(neg), "median_um_silent": neg.median(),
                     "p": f"{p:.3g}"})
    return pd.DataFrame(rows)


def peer_percentiles(tpm, counts, kb) -> pd.DataFrame:
    """Where each receptor sits among genes of its own exonic length.

    The length correction is applied to the level, but detection is not
    correctable that way: a longer transcript earns more reads and so clears the
    one-count floor more often. Comparing each receptor with genes of the same
    length holds that constant.
    """
    prof = pd.DataFrame({"exonic_kb": kb.reindex(tpm.index),
                         "TPM": tpm.mean(axis=1),
                         "pct": 100.0 * (counts.reindex(tpm.index) > 0).mean(axis=1)})
    prof = prof[prof.TPM > 0]
    rows = []
    for g in RECEPTORS:
        e = prof.exonic_kb[g]
        band = prof[(prof.exonic_kb >= e / PEER_BAND)
                    & (prof.exonic_kb <= e * PEER_BAND)]
        rows.append({
            "gene": g, "exonic_kb": e, "n_peers": len(band),
            "TPM": prof.TPM[g],
            "TPM_percentile": stats.percentileofscore(band.TPM, prof.TPM[g]),
            "pct_detected": prof.pct[g],
            "detection_percentile": stats.percentileofscore(band.pct, prof.pct[g]),
        })
    return pd.DataFrame(rows)


def main() -> int:
    meta = cell_metadata()
    counts = pd.read_csv(fetch(COUNTS, "raw_counts.csv.gz"), index_col=0)
    counts.index = counts.index.astype(str)
    counts = counts[~counts.index.duplicated()]
    if list(counts.columns) != list(meta.index):
        raise ac.SanityCheckError(
            "the count matrix columns are not the series matrix samples in "
            "order; the per-cell labels cannot be joined by position")

    missing = [g for g in RECEPTORS if g not in counts.index]
    if missing:
        raise ac.SanityCheckError(f"{missing} absent from the matrix")

    meta = meta[~meta.type.isin(NOT_NEURONS)]
    counts = counts[meta.index]
    print(f"  {len(meta)} neuronal somata, {len(set(meta.type))} types, "
          f"{meta.donor.nunique()} donors")

    kb = exonic_kb()
    tpm, modelled = to_tpm(counts, kb)
    print(f"  {int(modelled.sum())} of {len(counts)} rows carry an exon model")
    cpm = counts.div(counts.sum(axis=0), axis=1) * 1e6

    rows = []
    for t, cells in meta.groupby("type").groups.items():
        for g in RECEPTORS:
            rows.append({
                "type": t, "n_cells": len(cells), "gene": g,
                "TPM_mean": tpm.loc[g, cells].mean(),
                "TPM_median": tpm.loc[g, cells].median(),
                "CPM_mean": cpm.loc[g, cells].mean(),
                "pct_detected": 100.0 * (counts.loc[g, cells] > 0).mean(),
                "soma_um_median": meta.loc[cells, "soma_um"].median(),
            })
    by_type = pd.DataFrame(rows)
    ac.save_table(by_type.round(4), "human_soma_receptors_by_type.csv")

    pooled = pd.DataFrame({
        "gene": RECEPTORS,
        "exonic_kb": kb.reindex(RECEPTORS).values,
        "CPM_mean": cpm.loc[RECEPTORS].mean(axis=1).values,
        "TPM_mean": tpm.loc[RECEPTORS].mean(axis=1).values,
        "pct_detected": (100.0 * (counts.loc[RECEPTORS] > 0).mean(axis=1)).values,
    })
    ac.save_table(pooled.round(4), "human_soma_receptors_pooled.csv")

    contrast = class_contrast(tpm, counts, meta)
    ac.save_table(contrast.round(4), "human_soma_class_contrast.csv")
    size = soma_size(tpm, counts, meta)
    ac.save_table(size.round(4), "human_soma_size_by_expression.csv")
    peers = peer_percentiles(tpm, counts, kb)
    ac.save_table(peers.round(3), "human_soma_peer_percentiles.csv")

    if pooled.set_index("gene").TPM_mean.idxmax() != "OPRM1":
        raise ac.SanityCheckError(
            "OPRM1 no longer leads the pooled human somata; the text of this "
            "script and the Results describe a result that has changed")

    piv = by_type.pivot(index="type", columns="gene", values="TPM_mean")[RECEPTORS]
    piv.insert(0, "n", by_type.groupby("type").n_cells.first())
    piv.insert(1, "soma_um", by_type.groupby("type").soma_um_median.first().round(1))
    piv["top"] = piv[RECEPTORS].idxmax(axis=1)
    print("\n  Mean TPM by type:")
    print(piv.sort_values("OPRL1", ascending=False).round(3).to_string())

    pct = by_type.pivot(index="type", columns="gene", values="pct_detected")[RECEPTORS]
    print("\n  Percent of somata detected:")
    print(pct.sort_values("OPRL1", ascending=False).round(1).to_string())

    print("\n  Pooled over all neuronal somata:")
    print(pooled.round(3).to_string(index=False))
    print("\n  A-fibre LTMR and proprioceptive types against nociceptive types:")
    print(contrast.round(4).to_string(index=False))
    print("\n  Soma diameter, expressing against silent:")
    print(size.round(2).to_string(index=False))
    print("\n  Against genes of the same exonic length:")
    print(peers.round(2).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
