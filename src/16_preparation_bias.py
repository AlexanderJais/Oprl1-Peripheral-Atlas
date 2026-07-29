"""Which preparation is distorted, measured on every gene rather than eight.

The receptor ordering in vagal neurons depends on the preparation: whole-cell
libraries put Oprl1 first, nuclear libraries put Oprm1 first. Deciding between
them is the question the rest of the atlas rests on, and it cannot be settled by
asserting a mechanism. It can be settled by asking which of the two measurements
is predicted by a variable that has nothing to do with how much receptor a
neuron carries.

Genomic span is that variable. The length of a transcription unit is a property
of the locus, not of the abundance of its mature message, so a quantification
that faithfully reports abundance should be uncorrelated with it. This script
computes, for every gene in the NodoMap annotation, the mean per-cell CPM in
whole-cell vagal neurons and in nuclear vagal neurons of the same atlas, and
regresses each of them, and their ratio, on genomic span.

Both preparations are taken from the same integration, the same tissue and the
same 10x chemistry, so the comparison is between preparations rather than
between studies.

Run from the repository root, after 00b_fetch_gene_spans.py and 02_nodose.py.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import h5py
import numpy as np
import pandas as pd
from scipy import stats

import atlas_common as ac

NEURON_CLASSES = ["Nodose ganglion neuron", "Jugular ganglion neuron"]
BLOCK = 4096          # cells per read, to keep the CSR stream off the heap
# A floor only to keep ratios off a near-zero denominator. It is set low enough
# to admit all four opioid receptors, which are the subject of this atlas and
# must never be filtered out of an analysis about them. The relationship is
# reported at a range of floors below, and does not depend on the choice.
MIN_CPM = 0.1
SPANS_ALL = "ensembl_gene_spans_all.csv.gz"


def cell_table(path):
    """cell_class, suspension_type and library size, without loading the matrix."""
    import anndata as ad

    adata = ad.read_h5ad(path, backed="r")
    obs = adata.obs
    lib = obs["nCount_RNA"].values.astype(float)
    if not np.isfinite(lib).all():
        raise ac.SanityCheckError("nCount_RNA has non-finite entries")
    return obs, lib, adata.n_obs


def mean_cpm_by_group(path, masks, lib, n_genes):
    """Mean per-cell CPM per gene, for each named mask, in one pass over the CSR.

    The same statistic as `02_nodose.pseudobulk_cpm`, computed for every gene
    instead of a chosen few. Rows are scaled by their own library size before
    summing, so a deep cell does not outweigh a shallow one.
    """
    acc = {k: np.zeros(n_genes) for k in masks}
    scale = np.where(lib > 0, 1e6 / np.where(lib > 0, lib, 1.0), 0.0)
    with h5py.File(path, "r") as fh:
        X = fh["raw/X"]
        indptr = X["indptr"][:]
        n_obs = indptr.size - 1
        for lo in range(0, n_obs, BLOCK):
            hi = min(lo + BLOCK, n_obs)
            if not any(m[lo:hi].any() for m in masks.values()):
                continue
            a, b = int(indptr[lo]), int(indptr[hi])
            data = X["data"][a:b]
            idx = X["indices"][a:b]
            starts = indptr[lo:hi] - a
            ends = indptr[lo + 1:hi + 1] - a
            for key, m in masks.items():
                rows = np.nonzero(m[lo:hi])[0]
                if not rows.size:
                    continue
                for r in rows:
                    s, e = int(starts[r]), int(ends[r])
                    if e > s:
                        np.add.at(acc[key], idx[s:e], data[s:e] * scale[lo + r])
    return {k: acc[k] / max(int(masks[k].sum()), 1) for k in masks}


def main() -> int:
    path = ac.NODOSE_H5AD
    if not path.exists():
        raise SystemExit(f"NodoMap atlas not found at {path}; see 00_download_data.sh")

    obs, lib, n_obs = cell_table(path)
    prep = obs["suspension_type"].values
    is_neuron = obs["cell_class"].isin(NEURON_CLASSES).values
    masks = {"whole_cell": is_neuron & (prep == "cell"),
             "nuclear": is_neuron & (prep == "nucleus")}
    for k, m in masks.items():
        print(f"  {k}: {int(m.sum()):,} vagal neurons")

    with h5py.File(path, "r") as fh:
        var = fh["raw/var"]
        gene_ids = np.array([s.decode() if isinstance(s, bytes) else str(s)
                             for s in var["gene_ids"][:]])
        cats = var["feature_name/categories"][:]
        cats = np.array([s.decode() if isinstance(s, bytes) else str(s)
                         for s in cats])
        symbols = cats[var["feature_name/codes"][:]]
        n_genes = gene_ids.size

    means = mean_cpm_by_group(path, masks, lib, n_genes)
    d = pd.DataFrame({"ensembl_id": gene_ids, "gene": symbols,
                      "whole_cell_CPM": means["whole_cell"],
                      "nuclear_CPM": means["nuclear"]})

    spans = pd.read_csv(ac.DATA / SPANS_ALL)
    d = d.merge(spans[["ensembl_id", "span_kb", "biotype"]], on="ensembl_id",
                how="inner")
    d = d[d.biotype == "protein_coding"]
    print(f"  {len(d):,} protein-coding genes matched to a span")

    # Filter on expression, never on span, so the selection cannot manufacture
    # the relationship the script is testing.
    keep = d[(d.whole_cell_CPM >= MIN_CPM) & (d.nuclear_CPM >= MIN_CPM)].copy()
    keep["ratio"] = keep.nuclear_CPM / keep.whole_cell_CPM
    print(f"  {len(keep):,} of them at >= {MIN_CPM} CPM in both preparations")
    missing = [g for g in ac.RECEPTORS if g not in set(keep.gene)]
    if missing:
        raise ac.SanityCheckError(
            f"the {MIN_CPM} CPM floor excludes {missing}; this atlas is about "
            "those four genes and an analysis of them cannot filter them out")

    # The floor is a judgement call, so its effect is measured rather than
    # asserted. Only the ratios need it: they divide by the whole-cell level.
    sens = []
    for floor in (0.0, 0.05, 0.1, 0.25, 0.5, 1.0):
        # Strictly greater than, so the no-floor row still has a positive
        # denominator for the ratio and a defined logarithm for both levels.
        # Substituting a pseudo-value for zero would put the correlation on
        # numbers that were never measured.
        k = d[(d.whole_cell_CPM > floor) & (d.nuclear_CPM > floor)]
        if len(k) < 50:
            continue
        lx = np.log10(k.span_kb)
        sens.append({"floor_cpm": floor, "n": len(k),
                     "r_ratio": round(float(np.corrcoef(
                         lx, np.log10(k.nuclear_CPM / k.whole_cell_CPM))[0, 1]), 3),
                     "r_whole_cell": round(float(np.corrcoef(
                         lx, np.log10(k.whole_cell_CPM))[0, 1]), 3),
                     "r_nuclear": round(float(np.corrcoef(
                         lx, np.log10(k.nuclear_CPM))[0, 1]), 3),
                     "receptors_included": int(k.gene.isin(ac.RECEPTORS).sum())})
    sens = pd.DataFrame(sens)
    print("\n  Sensitivity to the expression floor:")
    print(sens.to_string(index=False))
    ac.save_table(sens, "preparation_bias_floor_sensitivity.csv")

    x = np.log10(keep.span_kb)
    out = []
    for col, label in (("whole_cell_CPM", "whole-cell level"),
                       ("nuclear_CPM", "nuclear level"),
                       ("ratio", "nuclear / whole cell")):
        y = np.log10(keep[col])
        r = float(np.corrcoef(x, y)[0, 1])
        rho, p = stats.spearmanr(keep.span_kb, keep[col])
        out.append({"quantity": label, "n": len(keep), "log_log_pearson_r": round(r, 4),
                    "spearman_rho": round(float(rho), 4), "p": f"{p:.3g}"})
    res = pd.DataFrame(out)
    print("\n  Against genomic span, over every expressed protein-coding gene:")
    print(res.to_string(index=False))

    slope, intercept = np.polyfit(x, np.log10(keep.ratio), 1)
    cross = 10 ** (-intercept / slope)
    print(f"\n  ratio fit: slope {slope:.3f} per decade of span, crossing 1.0 at "
          f"{cross:.1f} kb")

    # Median ratio per span decile: the same result without a fitted line.
    keep["decile"] = pd.qcut(keep.span_kb, 10, labels=False)
    dec = (keep.groupby("decile")
           .agg(n=("ratio", "size"), span_lo=("span_kb", "min"),
                span_hi=("span_kb", "max"), span_med=("span_kb", "median"),
                median_whole_cell=("whole_cell_CPM", "median"),
                median_nuclear=("nuclear_CPM", "median"),
                median_ratio=("ratio", "median"))
           .round(3).reset_index())
    print("\n  By span decile. The preparation that moves with length is the "
          "one distorted by it:")
    print(dec.to_string(index=False))

    # Where the four receptors sit on that distribution.
    rec = keep[keep.gene.isin(ac.RECEPTORS)].set_index("gene").reindex(ac.RECEPTORS)
    rec["pct_all_genes"] = [
        round(float((keep.ratio < v).mean() * 100), 1) for v in rec.ratio]
    # Also against genes of its own length, which is the comparison that says
    # whether length alone accounts for where a receptor lands.
    band = []
    for sp, v in zip(rec.span_kb, rec.ratio):
        peers = keep[(keep.span_kb >= sp / 1.6) & (keep.span_kb <= sp * 1.6)]
        band.append((len(peers), round(float(peers.ratio.median()), 2),
                     round(float((peers.ratio < v).mean() * 100), 1)))
    rec[["n_peers", "peer_median_ratio", "pct_of_peers"]] = band
    print("\n  The receptors on that distribution, against all genes and "
          "against genes of their own length:")
    print(rec[["span_kb", "whole_cell_CPM", "nuclear_CPM", "ratio",
               "pct_all_genes", "n_peers", "peer_median_ratio",
               "pct_of_peers"]].round(3).to_string())
    ac.save_table(rec.reset_index(), "preparation_bias_receptors.csv")

    ac.save_table(keep[["gene", "ensembl_id", "span_kb", "whole_cell_CPM",
                        "nuclear_CPM", "ratio"]].sort_values("span_kb"),
                  "preparation_bias_genomewide.csv")
    ac.save_table(res, "preparation_bias_tests.csv")
    ac.save_table(dec, "preparation_bias_by_span_decile.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
