"""Shared definitions for the cross-ganglion Oprl1 comparison.

The central problem this module exists to solve: the three tissues were measured
on assays with very different per-cell sensitivity. Full-length SMART-seq on 96
geniculate neurons detects Oprl1 in 92% of cells; droplet data on nodose neurons
detects the same gene in 9%. Comparing those two numbers would measure the
platform, not the biology.

The subject of this atlas is `Oprl1` itself: where the NOP receptor is expressed,
in which neurons, and how strongly. Levels and detection rates are reported
directly, per dataset and per cluster. The other opioid genes are reported
alongside it for context, and `Pnoc` is carried as a descriptive column only.

Every cross-tissue statistic here is computed WITHIN a sample and compared only
as a rank or percentile:

  receptor_rank()        order of the four opioid receptors measured on the
                         same cells with the same chemistry
  transcriptome_percentile() where a gene sits in its own sample's expression
                         distribution

Raw detection percentages and expression levels are still reported per dataset,
but are never placed on a common axis across assays.

Two guards matter enough to state here, because getting either wrong is how this
comparison would produce a confident wrong answer:

  * A gene absent from a dataset's annotation is NOT a measured zero. Every
    table carries an `in_matrix` column derived from the matrix itself rather
    than assumed, and an unmeasured gene gets a null level, not 0.0.
  * A rank is only reported as determinate when the top receptor is strictly
    above the runner-up. `bootstrap_receptor_support` puts a number on how much
    a given ordering depends on which cells were sampled.
"""

from __future__ import annotations

import gzip
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "raw"
RES = ROOT / "results"
FIG = ROOT / "figures"

# The nodose atlas lives in the sibling PNOC-Nodose project. Override without
# editing this file:  export NODOSE_ROOT=/path/to/PNOC-Nodose
NODOSE_ROOT = Path(os.environ.get("NODOSE_ROOT", "/home/user/PNOC-Nodose"))
NODOSE_H5AD = NODOSE_ROOT / "data" / "nodomap_integrated.h5ad"

# Source filenames, shared so no script re-hardcodes them.
GENICULATE_FPKM = "GSE102443_GEO-ID_Dvoryanchikov_2017_Datatable_FPKM.txt.gz"
GENICULATE_XLSX = "gse135801/GSM4037432_GG_scRNAseq_Phox2b_expressed454.xlsx"

RECEPTORS = ["Oprl1", "Oprm1", "Oprd1", "Oprk1"]
LIGANDS = ["Pnoc", "Penk", "Pdyn", "Pomc"]
OPIOID_GENES = RECEPTORS + LIGANDS

# Genes used to confirm each dataset is what it claims to be before any opioid
# number is read out of it. REQUIRED_MARKERS are the ones whose absence or
# silence means the matrix is not what we think it is, so they are enforced;
# the rest are tissue-specific and are reported without being enforced.
SANITY_GENES = ["Snap25", "Phox2b", "Slc17a6", "Tac1", "Calca", "Actb"]
REQUIRED_MARKERS = ["Snap25", "Actb"]

# Dataset registry: preparation type is the single most important covariate in
# this comparison, so it is declared once, here, rather than per script.
#
# Provenance of the `prep` field. GSE102443 and GSE135801 are whole-cell by the
# publications' own methods (SMART-seq on picked cells; dissociated Phox2b+
# cells). The NodoMap assignments follow the `assay`/`suspension_type` fields of
# the CELLxGENE-hosted atlas: Bai, Buchanan, Kupari and Zhao are cell
# suspensions, the in-house set is nuclei. GSE166648 is snRNA-seq by title.
#
# NOTE (see AUDIT.md 3.2): `prep` is perfectly confounded with laboratory in
# this collection. There is exactly one nuclear nodose dataset and no in-house
# whole-cell dataset, so "nuclear" and "in-house" cannot be separated here.
DATASETS = {
    "GSE102443":       dict(tissue="geniculate",     prep="whole cell"),
    "GSE135801":       dict(tissue="geniculate",     prep="whole cell"),
    "NodoMap:Bai":      dict(tissue="nodose/jugular", prep="whole cell"),
    "NodoMap:Buchanan": dict(tissue="nodose/jugular", prep="whole cell"),
    "NodoMap:Kupari":   dict(tissue="nodose/jugular", prep="whole cell"),
    "NodoMap:Zhao":     dict(tissue="nodose/jugular", prep="whole cell"),
    "NodoMap:inhouse":  dict(tissue="nodose/jugular", prep="nuclear"),
    "GSE166648":       dict(tissue="NTS (central)",  prep="nuclear"),
}
PREP = {k: v["prep"] for k, v in DATASETS.items()}
TISSUE_OF = {k: v["tissue"] for k, v in DATASETS.items()}
WHOLE_CELL_NODOSE = [d for d, v in DATASETS.items()
                     if d.startswith("NodoMap:") and v["prep"] == "whole cell"]

# Platform description only. Cell counts live in each table's `n_cells` column,
# where they refer to the rows actually summarised.
ASSAY_LABEL = {
    "GSE102443": "SMART-seq (full-length)",
    "GSE135801": "scRNA-seq (3')",
    "NodoMap": "10x droplet sc/snRNA-seq",
    "GSE166648": "snRNA-seq",
}


class SanityCheckError(RuntimeError):
    """A dataset failed its marker-gene check, so no opioid number is read."""


# ---------------------------------------------------------------------------
# Within-sample statistics: the only quantities compared across tissues
# ---------------------------------------------------------------------------

def receptor_rank(levels: pd.Series) -> pd.DataFrame:
    """Rank the four opioid receptors within one sample.

    `levels` maps gene symbol to any monotonic measure of abundance (mean FPKM,
    mean CPM, percent detected) computed on the same cells. Rank is invariant to
    the units, which is what makes it comparable between assays; the underlying
    values are not.

    `determinate` is False when the ordering does not actually identify a top
    receptor: either every level is zero, or the top two are tied. Callers must
    not read a `rank == 1` row without checking it.
    """
    present = [g for g in RECEPTORS if g in levels.index]
    if not present:
        raise ValueError("no opioid receptors present in this sample")
    s = levels.loc[present].astype(float)
    if not np.isfinite(s.values).all():
        bad = list(s.index[~np.isfinite(s.values)])
        raise ValueError(f"non-finite receptor levels for {bad}; "
                         "a missing gene must be absent from the index, not NaN")

    top = float(s.max())
    n_at_top = int((s == top).sum())
    determinate = bool(top > 0 and n_at_top == 1)

    out = pd.DataFrame({
        "gene": present,
        "level": s.values,
        "rank": s.rank(ascending=False, method="min").astype(int).values,
        "determinate": determinate,
    })
    out["fraction_of_top"] = (s / top).values if top > 0 else np.nan
    return out.sort_values("rank").reset_index(drop=True)


def bootstrap_receptor_support(per_cell: pd.DataFrame, n_boot: int = 2000,
                               seed: int = 0) -> dict:
    """How much of the receptor ordering survives resampling the cells.

    `per_cell` is cells x receptor genes of per-cell normalised values (CPM or
    FPKM) — the same quantity whose column mean produces the pseudobulk level.
    Cells are resampled with replacement; `support` is the fraction of
    resamples in which the observed top receptor is still top.

    This is the piece the project was missing: a pseudobulk ordering with a 1.8%
    margin and one with a 27% margin are otherwise reported identically.
    """
    cols = [g for g in RECEPTORS if g in per_cell.columns]
    X = per_cell[cols].to_numpy(dtype=float, copy=False)
    n = X.shape[0]
    obs = X.mean(axis=0)
    if n < 2 or not np.isfinite(obs).all() or obs.max() <= 0:
        return {"top_gene": None, "support": np.nan, "margin": np.nan,
                "margin_lo": np.nan, "margin_hi": np.nan, "n_boot": 0}

    order = np.argsort(-obs)
    top_i, second_i = order[0], order[1]
    rng = np.random.default_rng(seed)

    wins = 0
    margins = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, n, n)
        m = X[idx].mean(axis=0)
        wins += int(np.argmax(m) == top_i)
        denom = m[second_i]
        margins[b] = m[top_i] / denom if denom > 0 else np.inf

    finite = margins[np.isfinite(margins)]
    return {
        "top_gene": cols[top_i],
        "runner_up": cols[second_i],
        "support": wins / n_boot,
        "margin": float(obs[top_i] / obs[second_i]) if obs[second_i] > 0 else np.inf,
        "margin_lo": float(np.percentile(finite, 2.5)) if finite.size else np.nan,
        "margin_hi": float(np.percentile(finite, 97.5)) if finite.size else np.nan,
        "n_cells": int(n),
        "n_boot": int(n_boot),
    }


def transcriptome_percentile(levels: pd.Series, gene: str, min_level=0.0):
    """Where a gene sits in its own sample's expression distribution.

    Expressed as a percentile of all genes above `min_level`, so it answers
    "how abundant is this transcript for this tissue" without importing units
    from another assay.
    """
    if gene not in levels.index:
        return np.nan
    pool = levels[levels > min_level].astype(float)
    if len(pool) < 100:
        return np.nan
    return float(100.0 * (pool < float(levels[gene])).mean())


def leave_one_out_pearson(x, y, labels) -> pd.DataFrame:
    """Pearson r on log10 axes, recomputed with each point dropped in turn.

    A correlation over seven genes can be carried entirely by two of them. This
    reports that directly instead of leaving it to the reader to notice.
    """
    x = np.log10(np.asarray(x, dtype=float))
    y = np.log10(np.asarray(y, dtype=float))
    labels = list(labels)
    full = float(np.corrcoef(x, y)[0, 1])
    rows = [{"dropped": "(none)", "n": len(x), "pearson_r": round(full, 4),
             "delta_vs_full": 0.0}]
    for i, lab in enumerate(labels):
        keep = np.arange(len(x)) != i
        r = float(np.corrcoef(x[keep], y[keep])[0, 1])
        rows.append({"dropped": lab, "n": int(keep.sum()),
                     "pearson_r": round(r, 4),
                     "delta_vs_full": round(r - full, 4)})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Marker-gene gate
# ---------------------------------------------------------------------------

def check_markers(levels: pd.Series, label: str, required=None,
                  strict: bool = True) -> pd.DataFrame:
    """Confirm a dataset is what it claims before any opioid number is read.

    Reports every gene in SANITY_GENES, and enforces REQUIRED_MARKERS: those
    must be present in the matrix and non-zero, or `SanityCheckError` is raised.
    The remaining markers are tissue-specific (Phox2b is not expected in every
    preparation), so they are recorded and not enforced.
    """
    required = list(REQUIRED_MARKERS if required is None else required)
    rows = []
    for g in SANITY_GENES:
        here = g in levels.index
        rows.append({
            "dataset": label,
            "gene": g,
            "in_matrix": bool(here),
            "level": round(float(levels[g]), 4) if here else np.nan,
            "enforced": g in required,
        })
    tbl = pd.DataFrame(rows)

    failed = [r["gene"] for r in rows
              if r["enforced"] and not (r["in_matrix"] and r["level"] > 0)]
    shown = ", ".join(f"{r['gene']} {r['level']:.1f}" for r in rows if r["in_matrix"])
    absent = [r["gene"] for r in rows if not r["in_matrix"]]
    print(f"  {label} markers: {shown}" + (f"  [absent: {absent}]" if absent else ""))
    if failed:
        msg = f"{label}: required marker(s) missing or zero: {failed}"
        if strict:
            raise SanityCheckError(msg)
        print(f"  [warn] {msg}")
    return tbl


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def stream_gene_rows(path: Path, wanted, sep=",", totals=False,
                     progress_every: int = 0):
    """Pull gene rows out of a huge genes-by-cells CSV without loading it.

    GSE166648's expression matrix is ~100 MB gzipped and expands to more than
    the session's disk allowance, so it is read one line at a time and only the
    requested rows are kept. With `totals=True` a per-gene sum is accumulated
    for every row as it passes, which costs one pass instead of two.

    Returns (frame, absent, totals_or_None) where `frame` is cells x found
    genes and `absent` lists the requested symbols the file does not contain —
    a reference gap, which callers must not record as a measured zero.

    The first header field is the row-index column and is dropped; a row whose
    value count disagrees with the header is skipped and counted. If every row
    disagrees the header offset is wrong for this file, which is raised as such
    rather than surfacing later as an empty-array error.
    """
    wanted = {w.strip('"') for w in wanted}
    opener = gzip.open if str(path).endswith(".gz") else open
    rows, names, tot = {}, [], []
    n_rows = n_bad = 0

    with opener(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split(sep)
        cells = [c.strip('"') for c in header[1:]]
        n = len(cells)
        for line in fh:
            cut = line.find(sep)
            if cut < 0:
                continue
            n_rows += 1
            name = line[:cut].strip('"')
            vals = np.fromstring(line[cut + 1:], sep=sep)
            if vals.size != n:
                n_bad += 1
                continue
            if totals:
                names.append(name)
                tot.append(vals.sum())
            if name in wanted and name not in rows:
                rows[name] = vals.astype(np.float32)
            if not totals and len(rows) == len(wanted):
                break
            if progress_every and n_rows % progress_every == 0:
                print(f"    {n_rows:,} genes", flush=True)

    if n_rows and n_bad == n_rows:
        raise ValueError(
            f"{path.name}: every one of {n_rows:,} rows had a value count other "
            f"than the {n} implied by the header. The header offset is wrong for "
            f"this file, or the separator is not {sep!r}.")
    if n_bad:
        print(f"  [warn] skipped {n_bad:,}/{n_rows:,} malformed rows in {path.name}")
    if not rows:
        raise ValueError(f"none of {sorted(wanted)} found in {path}")

    found = sorted(rows)
    frame = pd.DataFrame({g: rows[g] for g in found}, index=cells)
    absent = sorted(wanted - set(found))
    return frame, absent, (pd.Series(tot, index=names) if totals else None)


def csr_gene_columns(h5group, columns_per_gene: dict, n_rows: int) -> np.ndarray:
    """Pull a few genes out of an on-disk CSR matrix without loading it all.

    The NodoMap atlas has ~176M non-zeros; materialising it costs well over a
    gigabyte, and we only ever want a dozen genes. The CSR arrays are streamed
    in chunks and only entries falling in the requested columns are kept.

    `columns_per_gene` maps output gene name to the list of source column
    indices carrying that symbol. A symbol appearing on more than one row of the
    annotation has its columns SUMMED, matching how transcript-level rows are
    collapsed elsewhere in this project — taking the first column instead would
    silently undercount the gene.
    """
    genes = list(columns_per_gene)
    out_of = {}
    for j, g in enumerate(genes):
        for c in columns_per_gene[g]:
            out_of.setdefault(int(c), []).append(j)

    wanted = np.zeros(int(h5group.attrs["shape"][1]), dtype=bool)
    wanted[list(out_of)] = True
    indptr = h5group["indptr"][:]
    indices_ds, data_ds = h5group["indices"], h5group["data"]
    nnz = indices_ds.shape[0]

    out = np.zeros((n_rows, len(genes)), dtype=np.float64)
    chunk = 20_000_000
    for start in range(0, nnz, chunk):
        stop = min(start + chunk, nnz)
        idx = indices_ds[start:stop]
        keep = np.nonzero(wanted[idx])[0]
        if not keep.size:
            continue
        vals = data_ds[start:stop][keep]
        cols = idx[keep]
        rows = np.searchsorted(indptr, keep + start, side="right") - 1
        for src, targets in out_of.items():
            m = cols == src
            if not m.any():
                continue
            for j in targets:
                np.add.at(out, (rows[m], j), vals[m])
    return out


def save_table(df: pd.DataFrame, name: str, index=False):
    RES.mkdir(parents=True, exist_ok=True)
    df.to_csv(RES / name, index=index)
    print(f"  [csv ] results/{name}")


# Figures live in src/atlas_style.py, matched to the sibling PNOC-Nodose project.
