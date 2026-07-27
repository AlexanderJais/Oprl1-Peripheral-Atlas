"""Shared definitions for the cross-ganglion Oprl1 comparison.

The central problem this module exists to solve: the three tissues were measured
on assays with very different per-cell sensitivity. Full-length SMART-seq on 96
geniculate neurons detects Oprl1 in 92% of cells; droplet data on nodose neurons
detects the same gene in 9%. Comparing those two numbers would measure the
platform, not the biology.

Every cross-tissue statistic here is therefore computed WITHIN a sample and
compared only as a rank, ratio or percentile:

  receptor_rank()        order of the four opioid receptors measured on the
                         same cells with the same chemistry
  ligand_receptor_ratio() Oprl1 against Pnoc in the same cells
  transcriptome_percentile() where a gene sits in its own sample's expression
                         distribution, among abundance-matched genes

Raw detection percentages and expression levels are still reported per dataset,
but are never placed on a common axis across assays.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "raw"
RES = ROOT / "results"
FIG = ROOT / "figures"

# The nodose atlas lives in the sibling PNOC-Nodose project.
NODOSE_ROOT = Path("/home/user/PNOC-Nodose")
NODOSE_H5AD = NODOSE_ROOT / "data" / "nodomap_integrated.h5ad"

RECEPTORS = ["Oprl1", "Oprm1", "Oprd1", "Oprk1"]
LIGANDS = ["Pnoc", "Penk", "Pdyn", "Pomc"]
OPIOID_GENES = RECEPTORS + LIGANDS

# Genes used to confirm each dataset is what it claims to be before any opioid
# number is read out of it.
SANITY_GENES = ["Snap25", "Phox2b", "Slc17a6", "Tac1", "Calca", "Actb"]

TISSUE_COLORS = {
    "geniculate": "#1B7837",
    "nodose": "#08306B",
    "jugular": "#6BAED6",
    "NTS": "#B2182B",
}

ASSAY_LABEL = {
    "GSE102443": "SMART-seq (full-length, 96 cells)",
    "GSE135801": "scRNA-seq (454 cells)",
    "NodoMap": "10x droplet sc/snRNA-seq (106,436 cells)",
    "GSE166648": "snRNA-seq (49,392 neuronal nuclei)",
}


# ---------------------------------------------------------------------------
# Within-sample statistics: the only quantities compared across tissues
# ---------------------------------------------------------------------------

def receptor_rank(levels: pd.Series) -> pd.DataFrame:
    """Rank the four opioid receptors within one sample.

    `levels` maps gene symbol to any monotonic measure of abundance (mean FPKM,
    mean CPM, percent detected) computed on the same cells. Rank is invariant to
    the units, which is what makes it comparable between assays; the underlying
    values are not.
    """
    present = [g for g in RECEPTORS if g in levels.index]
    s = levels.loc[present].astype(float)
    out = pd.DataFrame({
        "gene": present,
        "level": s.values,
        "rank": s.rank(ascending=False, method="min").astype(int).values,
    })
    top = s.max()
    out["fraction_of_top"] = (s / top).values if top > 0 else np.nan
    return out.sort_values("rank").reset_index(drop=True)


def ligand_receptor_ratio(levels: pd.Series, receptor="Oprl1", ligand="Pnoc"):
    """Receptor against its ligand precursor in the same cells.

    Both genes come from one matrix, so shared depth and shared capture
    efficiency cancel. A ratio far above 1 in every tissue is the
    receptor-without-local-ligand pattern.

    Absence from the quantified annotation and a measured zero are treated the
    same way here, both giving an infinite ratio. They are not the same evidence,
    so callers report `<ligand>_in_matrix` alongside this value: a gene missing
    from the annotation could be a reference gap rather than a biological zero.
    """
    r = float(levels.get(receptor, np.nan))
    l = float(levels.get(ligand, 0.0))
    if not np.isfinite(r):
        return np.nan
    if not np.isfinite(l):
        l = 0.0
    return r / l if l > 0 else np.inf


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


def abundance_matched_percentile(levels: pd.Series, gene: str,
                                 group_a: pd.Series, group_b: pd.Series,
                                 tol=0.25, min_controls=50):
    """Enrichment of `gene` in group A over group B, against matched controls.

    `levels` is overall abundance (used only to pick controls); `group_a` and
    `group_b` are per-gene abundance within each group. Returns the percentile
    of this gene's log2 enrichment among all genes of similar overall
    abundance, plus the one-sided empirical p (Phipson-Smyth).
    """
    if gene not in levels.index:
        return np.nan, np.nan, 0
    target = float(levels[gene])
    if not np.isfinite(target) or target <= 0:
        return np.nan, np.nan, 0
    lo, hi = target * (1 - tol), target * (1 + tol)
    pool = [g for g in levels[(levels >= lo) & (levels <= hi)].index if g != gene]
    if len(pool) < min_controls:
        return np.nan, np.nan, len(pool)

    def enrich(g):
        a = float(group_a.get(g, 0.0))
        b = float(group_b.get(g, 0.0))
        return np.log2((a + 0.01) / (b + 0.01))

    s = enrich(gene)
    ctrl = np.array([enrich(g) for g in pool])
    pct = 100.0 * float((ctrl < s).mean())
    p = float(((ctrl >= s).sum() + 1) / (ctrl.size + 1))
    return pct, p, len(pool)


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def collapse_isoforms(df: pd.DataFrame, symbol_col=None) -> pd.DataFrame:
    """Sum transcript-level rows to gene level, matching the geniculate project."""
    if symbol_col is None:
        return df.groupby(df.index).sum()
    return df.groupby(df[symbol_col]).sum(numeric_only=True)


def stream_gene_rows(path: Path, wanted, sep=","):
    """Pull a few gene rows out of a huge genes-by-cells CSV without loading it.

    GSE166648's expression matrix is ~100 MB gzipped and expands to more than
    the session's disk allowance, so it is read one line at a time and only the
    requested rows are kept.
    """
    wanted = {w.strip('"') for w in wanted}
    opener = gzip.open if str(path).endswith(".gz") else open
    rows, header = {}, None
    with opener(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split(sep)
        # A leading empty field means the first column is the row index.
        cells = [c.strip('"') for c in header[1:]]
        for line in fh:
            cut = line.find(sep)
            if cut < 0:
                continue
            name = line[:cut].strip('"')
            if name in wanted:
                vals = np.fromstring(line[cut + 1:], sep=sep)
                if vals.size == len(cells):
                    rows[name] = vals
                if len(rows) == len(wanted):
                    break
    if not rows:
        raise ValueError(f"none of {sorted(wanted)} found in {path}")
    return pd.DataFrame(rows, index=cells)


def save_table(df: pd.DataFrame, name: str, index=False):
    RES.mkdir(parents=True, exist_ok=True)
    df.to_csv(RES / name, index=index)
    print(f"  [csv ] results/{name}")


def save_fig(fig, name: str):
    FIG.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"{name}.{ext}", bbox_inches="tight")
    print(f"  [fig ] figures/{name}.pdf|.png")
    import matplotlib.pyplot as plt
    plt.close(fig)


def set_theme():
    import matplotlib as mpl
    mpl.use("Agg")
    mpl.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "font.family": "sans-serif",
        "font.sans-serif": ["Nimbus Sans", "Helvetica", "Arial", "DejaVu Sans"],
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.frameon": False,
        "pdf.fonttype": 42,
    })
