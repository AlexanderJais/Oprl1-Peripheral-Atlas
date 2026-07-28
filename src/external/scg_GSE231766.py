"""Superior cervical ganglion, GSE231766: the second sympathetic dataset.

The sympathetic control in section 1 of the README rested on GSE78845, 298
thoracic neurons on a full-length platform. That is the platform class which
produces the largest margins in this project, so the control shared its main
weakness with the datasets it was meant to test. GSE231766 removes the overlap:
a cranial rather than a thoracic sympathetic ganglion, a different laboratory,
and 10x droplet chemistry.

Four samples, two untreated and two from mice with transverse aortic
constriction at 5 and 18 days. The reported row uses the two untreated animals;
the disease samples are reported separately so the condition can be seen not to
carry the result.

Run from the repository root. The GEO tarball is fetched into `scratch/` on
first use and reused afterwards.
"""

import gzip
import glob
import os
import subprocess
import sys
import tarfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

import atlas_common as ac

ACC = "GSE231766"
URL = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE231nnn/GSE231766/suppl/"
       "GSE231766_RAW.tar")
RAW = Path(os.environ.get("SCRATCH", "scratch")) / "GSE231766"

# Disease state read from the GEO sample characteristics, not from the titles,
# which carry only internal animal identifiers.
UNTREATED = ("GSM7300193", "GSM7300194")

# Neurons are called on raw counts so the threshold does not move with library
# size, and glial and immune barcodes are excluded on CPM so the threshold does
# not move with sequencing run.
MIN_UMI, MIN_GENES = 2000, 1000
NEURON_RAW = {"Snap25": 5, "Th": 5}
NOT_NEURON_CPM = {"Sox10": 200, "Plp1": 500, "Ptprc": 50}


def fetch():
    RAW.mkdir(parents=True, exist_ok=True)
    if glob.glob(str(RAW / "*_matrix.mtx.gz")):
        return
    tar = RAW / f"{ACC}_RAW.tar"
    if not tar.exists():
        print(f"  downloading {URL}")
        subprocess.run(["curl", "-fsSL", "-o", str(tar), URL], check=True)
    with tarfile.open(tar) as t:
        t.extractall(RAW)
    tar.unlink()


def load():
    """Every sample stacked into one cells-by-genes matrix of raw counts."""
    samples = sorted({os.path.basename(f).rsplit("_", 1)[0]
                      for f in glob.glob(str(RAW / "*_matrix.mtx.gz"))})
    if not samples:
        raise FileNotFoundError(f"no matrices under {RAW}")

    mats, sample_of, genes = [], [], None
    for s in samples:
        with gzip.open(RAW / f"{s}_features.tsv.gz", "rt") as fh:
            g = [line.split("\t")[1] for line in fh]
        if genes is None:
            genes = g
        elif g != genes:
            raise ValueError(f"{s} has a different feature list; stacking the "
                             "matrices would misalign genes")
        with gzip.open(RAW / f"{s}_matrix.mtx.gz", "rt") as fh:
            m = mmread(fh).tocsc().T.tocsr()          # genes x cells on disk
        with gzip.open(RAW / f"{s}_barcodes.tsv.gz", "rt") as fh:
            n_bc = sum(1 for _ in fh)
        if m.shape[0] != n_bc:
            raise ValueError(f"{s}: {m.shape[0]} rows against {n_bc} barcodes; "
                             "the matrix is not oriented cells-by-genes")
        mats.append(m.astype(np.float32))
        sample_of += [s] * n_bc
        print(f"  {s}: {m.shape[0]:,} barcodes")

    return (sparse.vstack(mats).tocsr(), np.array(sample_of),
            pd.Index(genes), samples)


def report(label, cpm, raw_counts, n_umi):
    print(f"\n  {label}: {len(cpm):,} neurons, "
          f"median library {np.median(n_umi):,.0f} UMI")
    for g in ac.RECEPTORS:
        print(f"    {g:6s} {cpm[g].mean():7.2f} CPM   "
              f"{(raw_counts[g] > 0).mean() * 100:5.1f}% detected")
    b = ac.bootstrap_receptor_support(cpm[ac.RECEPTORS], n_boot=ac.N_BOOT)
    print(f"    top {b['top_gene']} over {b['runner_up']}  "
          f"{b['margin']:.2f}x [{b['margin_lo']:.2f}-{b['margin_hi']:.2f}]  "
          f"support {b['support']:.4f}")
    print(f"    over Oprm1 {cpm['Oprl1'].mean() / cpm['Oprm1'].mean():.1f}x")
    return b


def main() -> int:
    fetch()
    X, sample_of, genes, samples = load()
    print(f"  {X.shape[0]:,} barcodes x {X.shape[1]:,} genes")

    umi = np.asarray(X.sum(axis=1)).ravel()
    keep = (umi >= MIN_UMI) & (np.asarray((X > 0).sum(axis=1)).ravel() >= MIN_GENES)
    X, sample_of, umi = X[keep], sample_of[keep], umi[keep]
    print(f"  QC (>={MIN_UMI} UMI, >={MIN_GENES} genes): {keep.sum():,} cells, "
          f"median {np.median(umi):,.0f} UMI")

    def raw(gene):
        j = np.where(genes == gene)[0]
        return (np.zeros(X.shape[0]) if len(j) == 0
                else np.asarray(X[:, j].sum(axis=1)).ravel())

    counts = {g: raw(g) for g in
              set(NEURON_RAW) | set(NOT_NEURON_CPM) | set(ac.RECEPTORS)
              | set(ac.SANITY_GENES) | {"Dbh", "Prph", "Scn1a", "Scn10a", "Pvalb"}}
    cpm_of = lambda v: v / umi * 1e6

    neuron = np.ones(X.shape[0], bool)
    for g, t in NEURON_RAW.items():
        neuron &= counts[g] >= t
    for g, t in NOT_NEURON_CPM.items():
        neuron &= cpm_of(counts[g]) < t
    print(f"  sympathetic neurons: {neuron.sum():,} "
          f"({neuron.mean() * 100:.1f}% of QC-passing cells)")

    marker_cpm = pd.Series({g: float(cpm_of(v)[neuron].mean())
                            for g, v in counts.items()})
    print("  markers (CPM): " +
          "  ".join(f"{g} {marker_cpm[g]:,.0f}" for g in
                    ("Snap25", "Th", "Dbh", "Prph", "Plp1", "Sox10")))
    ac.save_table(ac.check_markers(marker_cpm, f"{ACC} (SCG neurons)"),
                  "scg_marker_checks.csv")

    cpm = pd.DataFrame({g: cpm_of(v) for g, v in counts.items()})
    untreated = neuron & np.isin([s.split("_")[0] for s in sample_of], UNTREATED)

    rows = []
    for label, mask in (("untreated only", untreated), ("all samples", neuron)):
        b = report(label, cpm[mask], {g: v[mask] for g, v in counts.items()},
                   umi[mask])
        rows.append({"population": label, "n": int(mask.sum()),
                     **{g: round(float(cpm.loc[mask, g].mean()), 3)
                        for g in ac.RECEPTORS},
                     "margin": round(b["margin"], 3),
                     "margin_lo": round(b["margin_lo"], 3),
                     "margin_hi": round(b["margin_hi"], 3),
                     "support": b["support"], "n_boot": b["n_boot"]})

    print(f"\n  Scn1a {cpm.loc[untreated, 'Scn1a'].mean():.2f} CPM   "
          f"Scn10a {cpm.loc[untreated, 'Scn10a'].mean():.2f}   "
          f"Pvalb {cpm.loc[untreated, 'Pvalb'].mean():.2f}")

    print("\n  each sample on its own:")
    for s in samples:
        m = neuron & (sample_of == s)
        if m.sum() < 20:
            print(f"    {s}: {m.sum()} neurons, too few to order")
            continue
        lv = {g: float(cpm.loc[m, g].mean()) for g in ac.RECEPTORS}
        order = sorted(lv, key=lv.get, reverse=True)
        rows.append({"population": s, "n": int(m.sum()),
                     **{g: round(lv[g], 3) for g in ac.RECEPTORS},
                     "margin": round(lv[order[0]] / lv[order[1]], 3)
                     if lv[order[1]] else np.inf,
                     "margin_lo": np.nan, "margin_hi": np.nan,
                     "support": np.nan, "n_boot": np.nan})
        print(f"    {s}: n = {m.sum():5,d}   " +
              "  ".join(f"{g} {lv[g]:.2f}" for g in order))

    ac.save_table(pd.DataFrame(rows), "scg_GSE231766_receptor_levels.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
