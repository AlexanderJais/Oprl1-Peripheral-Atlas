"""Autonomic and enteric ganglia: sympathetic, parasympathetic, enteric.

Section 1 of the README established the receptor ordering across the sensory
series and then, with GSE78845 and GSE231766, across two sympathetic ganglia.
This script extends the comparison to the rest of the peripheral nervous
system, so that "peripheral neurons" in the title is measured rather than
assumed.

  GSE232789  six autonomic ganglia from one laboratory on one platform:
             stellate, coeliac and lumbar chain (sympathetic), sphenopalatine
             (cranial parasympathetic) and two pelvic ganglia (mixed). The
             within-dataset design removes laboratory and platform from the
             comparison between divisions.
  GSE231924  stellate ganglion, cardiac-projecting sympathetic neurons, already
             filtered to neurons by the authors and deposited as Seurat
             log-normalised values at a scale of 10,000. CPM is recovered
             exactly as expm1(x) * 100.
  GSE330884  intrinsic cardiac nervous system, the parasympathetic ganglia of
             the heart.
  GSE263422  enteric neurons of the mouse small intestine at P7 and P24. The
             two ages are reported separately because they do not agree.

Neurons are called on raw counts, so the threshold does not move with library
size, and glial and immune barcodes are excluded on CPM, so that threshold does
not move with sequencing run. Every population passes the same marker gate as
the numbered pipeline before any receptor number is read from it.

Run from the repository root. Source files are fetched into `scratch/` on first
use and reused afterwards.
"""

import glob
import gzip
import os
import subprocess
import sys
import tarfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
from scipy.io import mmread

import atlas_common as ac

SCRATCH = Path(os.environ.get("SCRATCH", "scratch"))
FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"

MIN_UMI, MIN_GENES = 2000, 1000

# Neuronal identity, on raw counts. Snap25 alone admits ambient-dominated
# barcodes in preparations this glia-rich, so each dataset adds the pan-neuronal
# marker its own tissue expresses most strongly.
NOT_NEURON_CPM = {"Sox10": 200, "Plp1": 500, "Ptprc": 50}

# One entry per population that reaches a figure panel or a table row.
#   division: the grouping used in figure 1
#   pos:      raw-count thresholds that define a neuron
#   match:    substring selecting samples out of a multi-ganglion series
DATASETS = [
    dict(acc="GSE232789", kind="mtx",
         url=f"{FTP}/GSE232nnn/GSE232789/suppl/GSE232789_RAW.tar",
         pos={"Snap25": 5, "Tubb3": 5},
         markers=["Th", "Dbh", "Chat", "Slc18a3", "Phox2b", "Prph"],
         populations=[
             ("stellate", "sympathetic", "stellate"),
             ("coeliac", "sympathetic", "coeliac"),
             ("lumbar chain", "sympathetic", "lumbar"),
             ("sphenopalatine", "parasympathetic", "sphenopalatine"),
             ("pelvic", "mixed autonomic", "pelvic"),
         ]),
    dict(acc="GSE330884", kind="mtx",
         url=f"{FTP}/GSE330nnn/GSE330884/suppl/GSE330884_RAW.tar",
         pos={"Snap25": 5, "Phox2b": 2},
         markers=["Chat", "Slc5a7", "Slc18a3", "Prph", "Th"],
         populations=[("intrinsic cardiac", "parasympathetic", "")]),
    dict(acc="GSE263422", kind="mtx",
         url=f"{FTP}/GSE263nnn/GSE263422/suppl/GSE263422_RAW.tar",
         pos={"Snap25": 5, "Elavl4": 3},
         markers=["Chat", "Nos1", "Ret", "Slc18a3", "Phox2b"],
         populations=[("enteric submucosal, P24", "enteric", "P24"),
                      ("enteric submucosal, P7", "enteric", "P7")]),
    dict(acc="GSE231924", kind="seurat",
         url=f"{FTP}/GSE231nnn/GSE231924/suppl/GSE231924_neuron.dge.csv.gz",
         markers=["Th", "Dbh", "Prph"],
         populations=[("stellate", "sympathetic", "")]),
]


# ------------------------------------------------------------------ fetching
def fetch(ds):
    d = SCRATCH / ds["acc"]
    d.mkdir(parents=True, exist_ok=True)
    if ds["kind"] == "mtx":
        if glob.glob(str(d / "*_matrix.mtx.gz")):
            return d
        tar = d / "RAW.tar"
        if not tar.exists():
            print(f"  downloading {ds['acc']}")
            subprocess.run(["curl", "-fsSL", "-o", str(tar), ds["url"]], check=True)
        with tarfile.open(tar) as t:
            t.extractall(d)
        tar.unlink()
        # GSE232789 separates the sample name from the file role with a dot;
        # every other series uses an underscore. Normalise so one glob works.
        for f in d.glob("*.barcodes.tsv.gz"):
            for role in ("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"):
                src = Path(str(f).replace("barcodes.tsv.gz", role))
                src.rename(str(src).replace("." + role, "_" + role))
        return d

    out = d / "dge.csv.gz"
    if not out.exists():
        print(f"  downloading {ds['acc']}")
        subprocess.run(["curl", "-fsSL", "-o", str(out), ds["url"]], check=True)
    return d


# ------------------------------------------------------------------- loading
def load_mtx(d, match, want):
    """Per-sample QC and gene extraction, keeping each sample's own features.

    Samples in a series are not always annotated against the same reference, so
    the matrices are never stacked. Only the columns of interest are pulled out
    and concatenated, and the transcriptome-wide level used for the percentile
    is computed inside each sample.
    """
    samples = sorted({os.path.basename(f).rsplit("_", 1)[0]
                      for f in glob.glob(str(d / "*_matrix.mtx.gz"))})
    samples = [s for s in samples if match in s]
    if not samples:
        raise FileNotFoundError(f"no sample under {d} matches {match!r}")

    parts, umis, sample_of, levels = [], [], [], []
    for s in samples:
        with gzip.open(d / f"{s}_features.tsv.gz", "rt") as fh:
            g = pd.Index([line.split("\t")[1] for line in fh])
        with gzip.open(d / f"{s}_matrix.mtx.gz", "rt") as fh:
            m = mmread(fh).tocsc().T.tocsr().astype(np.float32)
        u = np.asarray(m.sum(axis=1)).ravel()
        n_gene = np.asarray((m > 0).sum(axis=1)).ravel()
        keep = (u >= MIN_UMI) & (n_gene >= MIN_GENES)
        m, u = m[keep], u[keep]

        cols = {}
        for w in want:
            j = np.where(g == w)[0]
            cols[w] = (np.zeros(m.shape[0]) if len(j) == 0
                       else np.asarray(m[:, j].sum(axis=1)).ravel())
        parts.append(pd.DataFrame(cols))
        umis.append(u)
        sample_of += [s] * int(keep.sum())
        tot = np.asarray(m.sum(axis=0)).ravel()
        levels.append(pd.Series(tot / tot.sum() * 1e6, index=g).groupby(level=0).sum())
        print(f"    {s}: {keep.sum():,}/{len(keep):,} cells pass QC, "
              f"median {np.median(u):,.0f} UMI")

    return (pd.concat(parts, ignore_index=True), np.concatenate(umis),
            np.array(sample_of), levels)


def load_seurat(d):
    """Seurat LogNormalize at a scale of 10,000, back to CPM.

    The normalisation is invertible: expm1 recovers counts-per-10,000, and the
    row sums confirm the scale rather than assuming it.
    """
    x = pd.read_csv(d / "dge.csv.gz", index_col=0)
    linear = np.expm1(x.values)
    scale = linear.sum(axis=1)
    if not np.allclose(scale, 10_000, rtol=1e-3):
        raise ValueError(f"expected a scale of 10,000 per cell, found "
                         f"{scale.min():.0f} to {scale.max():.0f}; the "
                         "back-transform to CPM would be wrong")
    return pd.DataFrame(linear * 100, index=x.index, columns=x.columns)


# ------------------------------------------------------------------ analysis
def summarise(acc, tissue, division, cpm, detected, umi, levels, n_samples):
    for g in ac.RECEPTORS:
        print(f"      {g:6s} {cpm[g].mean():8.2f} CPM  "
              f"{detected[g].mean() * 100:5.1f}% detected")
    b = ac.bootstrap_receptor_support(cpm[ac.RECEPTORS], n_boot=ac.N_BOOT)
    pct = [ac.transcriptome_percentile(l, "Oprl1") for l in levels
           if "Oprl1" in l.index]
    print(f"      top {b['top_gene']} over {b['runner_up']}  {b['margin']:.2f}x "
          f"[{b['margin_lo']:.2f}-{b['margin_hi']:.2f}]  "
          f"support {b['support']:.4f}")
    return {
        "dataset": acc, "tissue": tissue, "division": division,
        "n": int(len(cpm)), "n_samples": n_samples,
        "median_umi": (round(float(np.median(umi))) if umi is not None else np.nan),
        **{g: round(float(cpm[g].mean()), 3) for g in ac.RECEPTORS},
        **{f"{g}_pct": round(float(detected[g].mean() * 100), 2)
           for g in ac.RECEPTORS},
        "top_gene": b["top_gene"], "runner_up": b["runner_up"],
        "margin": round(b["margin"], 3),
        "margin_lo": round(b["margin_lo"], 3),
        "margin_hi": round(b["margin_hi"], 3),
        "support": b["support"], "n_boot": b["n_boot"],
        "over_Oprm1": round(float(cpm["Oprl1"].mean()
                                  / max(cpm["Oprm1"].mean(), 1e-9)), 2),
        "Oprl1_percentile": round(float(np.mean(pct)), 1) if pct else np.nan,
    }


def run_mtx(ds, d):
    want = sorted(set(ds["pos"]) | set(NOT_NEURON_CPM) | set(ac.RECEPTORS)
                  | set(ac.SANITY_GENES) | set(ds["markers"]))
    rows, checks = [], []
    for tissue, division, match in ds["populations"]:
        print(f"  {ds['acc']} {tissue}:")
        counts, umi, sample_of, levels = load_mtx(d, match, want)
        cpm = counts.div(umi, axis=0) * 1e6

        neuron = np.ones(len(counts), bool)
        for g, t in ds["pos"].items():
            neuron &= counts[g].values >= t
        for g, t in NOT_NEURON_CPM.items():
            neuron &= cpm[g].values < t
        print(f"    neurons: {neuron.sum():,} ({neuron.mean() * 100:.1f}% of "
              f"QC-passing cells), median {np.median(umi[neuron]):,.0f} UMI")
        print("    markers (CPM): " + "  ".join(
            f"{g} {cpm.loc[neuron, g].mean():,.0f}"
            for g in list(ds["pos"]) + ds["markers"] + list(NOT_NEURON_CPM)))

        label = f"{ds['acc']} ({tissue})"
        checks.append(ac.check_markers(cpm.loc[neuron].mean(), label))
        rows.append(summarise(ds["acc"], tissue, division, cpm.loc[neuron],
                              counts.loc[neuron] > 0, umi[neuron], levels,
                              len(set(sample_of))))

        for s in sorted(set(sample_of)):
            m = neuron & (sample_of == s)
            if m.sum() < 20:
                print(f"      {s}: {m.sum()} neurons, too few to order")
                continue
            lv = {g: float(cpm.loc[m, g].mean()) for g in ac.RECEPTORS}
            order = sorted(lv, key=lv.get, reverse=True)
            print(f"      {s}: n = {m.sum():5,d}   " +
                  "  ".join(f"{g} {lv[g]:.2f}" for g in order))
    return rows, checks


def run_seurat(ds, d):
    tissue, division, _ = ds["populations"][0]
    print(f"  {ds['acc']} {tissue}:")
    cpm = load_seurat(d)
    print(f"    {len(cpm):,} author-filtered neurons")
    print("    markers (CPM): " + "  ".join(
        f"{g} {cpm[g].mean():,.0f}" for g in ds["markers"] if g in cpm.columns))
    checks = [ac.check_markers(cpm.mean(), f"{ds['acc']} ({tissue})")]
    row = summarise(ds["acc"], tissue, division, cpm, cpm > 0, None,
                    [cpm.mean()], np.nan)
    # The deposit carries no library sizes, so median UMI is not recoverable.
    return [row], checks


def main() -> int:
    rows, checks = [], []
    for ds in DATASETS:
        d = fetch(ds)
        r, c = (run_mtx if ds["kind"] == "mtx" else run_seurat)(ds, d)
        rows += r
        checks += c

    table = pd.DataFrame(rows)
    ac.save_table(table, "autonomic_receptor_levels.csv")
    ac.save_table(pd.concat(checks, ignore_index=True),
                  "autonomic_marker_checks.csv")
    print("\n  Summary:")
    print(table[["dataset", "tissue", "division", "n", "Oprl1", "Oprm1",
                 "Oprd1", "Oprk1", "top_gene", "margin", "support"]]
          .to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
