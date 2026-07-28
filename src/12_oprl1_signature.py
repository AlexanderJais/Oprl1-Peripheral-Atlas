"""The Oprl1 signature: genes that travel with Oprl1 across the whole atlas.

Section 5 finds that Oprl1 tracks Nav1.1 in the nodose ganglion and the DRG.
That was a candidate-gene result read off cluster annotations. This script asks
the unbiased version of the same question in every dataset that carries a full
transcriptome, and keeps only what reproduces across them.

WHAT IS CORRELATED, AND WHY NOT SPEARMAN

Per-cell Spearman between two sparsely detected transcripts is dominated by
capture depth: a cell with more UMI detects more of everything, so almost every
gene correlates positively with almost every other gene. Section 5 of the README
records what that produced the first time, a correlate list whose top entry was
a Schwann-cell transcript riding in on ambient RNA.

The statistic here is instead the log2 fold change of each gene between
Oprl1-positive and Oprl1-negative neurons, computed WITHIN depth strata and
pooled across them. Depth is held fixed by construction, the comparison is
between cells of the same dataset sequenced together, and the quantity is the
one a reader already knows how to interpret.

WHAT IS RESAMPLED, AND WHY

The question is not whether a gene reaches significance in one ganglion; it is
whether it reproduces across the peripheral nervous system. So the unit of
resampling is the dataset, not the cell. The 13 populations are drawn with
replacement 10,000 times, and each draw produces a new median fold change for
every gene. That gives a 95% interval on the median and, for each gene, the
fraction of draws in which it stays in the top 30. A gene carried by one deep
dataset falls apart under this; a gene present in most of them does not.

Two things are excluded from the signature rather than ranked in it. Oprl1
itself is positive by construction and is reported as a positive control. Genes
more abundant in the glia of the same ganglion than in its neurons are flagged
from `results/peripheral_ambient_checks.csv`, since this pipeline applies no
ambient correction.

Run from the repository root. Sources are fetched into `scratch/` on first use.
"""

import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp

import atlas_common as ac
import atlas_style as st

_peri = __import__("10_peripheral_ganglia")
_gen = __import__("01_geniculate")

# A gene has to be present often enough for a fold change to mean anything.
MIN_DETECTION = 0.05          # fraction of neurons detecting the gene
MIN_CPM = 1.0                 # mean level in neurons
MIN_PER_GROUP = 30            # Oprl1-positive and negative cells per stratum
MAX_STRATA = 5
CELLS_PER_STRATUM = 200
PSEUDO_CPM = 0.1              # keeps a zero group mean from giving -inf

# A median over few datasets is unstable, and ranking on it puts the genes with
# the fewest observations on top. Requiring 10 of the 13 makes breadth a
# condition of entry rather than something the reader has to check afterwards.
MIN_DATASETS = 10             # of 13, for a gene to enter the signature
TOP_N = 30                    # the signature itself
N_BOOT = ac.N_BOOT

# Populations contributing a full transcriptome. The GSE231766 heart-disease
# pair and the two pelvic samples are already pooled into their rows above; the
# disease pair is left out here because it shares its ganglion and deposit with
# the untreated pair and would count that tissue twice.
POPULATIONS = [
    ("GSE102443", "geniculate (VII)", "sensory"),
    ("GSE135801", "geniculate (VII)", "sensory"),
    ("GSE309608", "vestibular (VIII)", "sensory"),
    ("GSE232789", "stellate", "sympathetic"),
    ("GSE232789", "coeliac", "sympathetic"),
    ("GSE232789", "lumbar chain", "sympathetic"),
    ("GSE231766", "superior cervical", "sympathetic"),
    ("GSE231924", "stellate", "sympathetic"),
    ("GSE232789", "sphenopalatine", "parasympathetic"),
    ("GSE330884", "intrinsic cardiac", "parasympathetic"),
    ("GSE232789", "pelvic", "mixed autonomic"),
    ("GSE263422", "enteric submucosal, P24", "enteric"),
    ("GSE263422", "enteric submucosal, P7", "enteric"),
]


# --------------------------------------------------------------------- loading
def neurons_of(acc, tissue):
    """Full transcriptome of one population's neurons, as (CPM, counts, genes).

    Returns cells x genes CPM, the raw counts needed for detection, and the gene
    index. The neuron definition is the one src/10_peripheral_ganglia.py already
    applies, so the cells here are exactly the cells behind the section 1 rows.
    """
    if acc == "GSE102443":
        m = _gen.load_dvoryanchikov().T              # cells x genes, FPKM
        m = m.loc[:, ~m.columns.duplicated()]
        return m, m, m.columns, None
    if acc == "GSE135801":
        m = _gen.load_zuker().T                      # cells x genes, CPM
        m = m.loc[:, ~m.columns.duplicated()]
        return m, m, m.columns, None

    ds = next(d for d in _peri.DATASETS if d["acc"] == acc)
    d = _peri.fetch(ds)
    if ds["kind"] == "seurat":
        m = _peri.load_seurat(d)
        m = m.loc[:, ~m.columns.duplicated()]
        return m, m, m.columns, None

    match = next(m for t, _, m in ds["populations"] if t == tissue)
    frames, glia = [], []
    for s in _peri._samples(d, ds["kind"], match):
        x, genes = _peri._sample_matrix(d, s, ds["kind"])
        umi = np.asarray(x.sum(axis=1)).ravel()
        n_gene = np.asarray((x > 0).sum(axis=1)).ravel()
        keep = (umi >= _peri.MIN_UMI) & (n_gene >= _peri.MIN_GENES)
        x, umi = x[keep], umi[keep]

        # The neuron call needs a handful of genes at raw-count and CPM scale.
        want = set(ds["pos"]) | set(_peri.NOT_NEURON_CPM)
        col = {}
        for w in want:
            j = np.where(genes == w)[0]
            col[w] = (np.zeros(x.shape[0]) if len(j) == 0
                      else np.asarray(x[:, j].sum(axis=1)).ravel())
        neuron = np.ones(x.shape[0], bool)
        for g, t in ds["pos"].items():
            neuron &= col[g] >= t
        for g, t in _peri.NOT_NEURON_CPM.items():
            neuron &= (col[g] / umi * 1e6) < t

        # The barcodes this population excludes are the ambient reference: what
        # a transcript reads in the glia and immune cells of the same
        # dissociation. Kept per sample so the signature can be scored for
        # contamination gene by gene rather than only for Oprl1.
        if (~neuron).sum():
            other = sp.diags(1e6 / umi[~neuron]) @ x[~neuron]
            glia.append((np.asarray(other.mean(axis=0)).ravel(),
                         pd.Index(genes), int((~neuron).sum())))

        x, umi = x[neuron], umi[neuron]
        if x.shape[0] == 0:
            continue
        # Collapse duplicate symbols before anything is read out of the matrix:
        # a symbol on two annotation rows is one gene, and leaving both in place
        # would put the same transcript into the signature twice.
        codes, uniq = pd.factorize(pd.Index(genes))
        collapse = sp.csr_matrix(
            (np.ones(len(codes), dtype=np.float32),
             (np.arange(len(codes)), codes)), shape=(len(codes), len(uniq)))
        frames.append((sp.csr_matrix(x) @ collapse, pd.Index(uniq), umi))
        print(f"    {s}: {x.shape[0]:,} neurons")

    if not frames:
        raise ac.SanityCheckError(f"{acc} {tissue}: no neurons survived QC")

    shared = frames[0][1]
    for _, g, _ in frames[1:]:
        if not g.equals(shared):
            shared = shared.intersection(g)
    counts, umis = [], []
    for x, g, u in frames:
        idx = g.get_indexer(shared)
        counts.append(x[:, idx])
        umis.append(u)
    X = sp.vstack(counts).tocsr()
    umi = np.concatenate(umis)
    cpm = sp.diags(1e6 / umi) @ X

    other_mean = None
    if glia:
        total = sum(n for _, _, n in glia)
        acc_v = np.zeros(len(shared))
        for v, g, n in glia:
            acc_v += pd.Series(v, index=g).groupby(level=0).sum().reindex(
                shared).fillna(0.0).to_numpy() * n
        other_mean = pd.Series(acc_v / total, index=shared)
    return cpm.tocsr(), X.tocsr(), shared, other_mean


# ------------------------------------------------------------------- statistic
def _column_stats(mat, genes):
    """Mean level and detection fraction per gene, for sparse or dense input."""
    if sp.issparse(mat):
        n = mat.shape[0]
        mean = np.asarray(mat.sum(axis=0)).ravel() / n
        det = np.asarray((mat > 0).sum(axis=0)).ravel() / n
    else:
        mean = mat.to_numpy().mean(axis=0)
        det = (mat.to_numpy() > 0).mean(axis=0)
    return pd.Series(mean, index=genes), pd.Series(det, index=genes)


def stratified_log2fc(cpm, counts, genes, label):
    """log2 fold change between Oprl1-positive and negative neurons.

    Cells are split into depth strata by the number of genes they detect, the
    comparison is made inside each stratum, and the strata are pooled weighted
    by how much each one can say. Depth therefore cannot produce the difference:
    every comparison is between cells that detected a similar number of genes.
    """
    if "Oprl1" not in genes:
        raise ac.SanityCheckError(f"{label}: Oprl1 absent from the annotation")

    dense = not sp.issparse(cpm)
    n_cells = cpm.shape[0]
    j = list(genes).index("Oprl1")
    oprl1 = (counts.iloc[:, j].to_numpy() if dense
             else np.asarray(counts[:, j].todense()).ravel())
    level = (cpm.iloc[:, j].to_numpy() if dense
             else np.asarray(cpm[:, j].todense()).ravel())

    depth = (np.asarray((counts > 0).sum(axis=1)).ravel() if not dense
             else (counts.to_numpy() > 0).sum(axis=1))
    n_strata = int(np.clip(n_cells // CELLS_PER_STRATUM, 1, MAX_STRATA))
    stratum = ac.depth_strata(depth, n_bins=n_strata)

    # Which cells count as Oprl1-positive depends on what the platform can
    # resolve. Where Oprl1 is detected in 10 to 90 percent of neurons, detection
    # is the split. Full-length data detects it in 92 percent, where a
    # detected-against-not split has almost nothing on one side, so the contrast
    # there is high against low level within the same depth stratum.
    det = float((oprl1 > 0).mean())
    if 0.10 <= det <= 0.90:
        pos, rule = oprl1 > 0, "detected"
    else:
        pos = np.zeros(n_cells, bool)
        for s in np.unique(stratum):
            in_s = stratum == s
            pos[in_s] = level[in_s] > np.median(level[in_s])
        rule = "above the within-stratum median"
    print(f"    {label}: Oprl1 detected in {det * 100:.1f}% of neurons, "
          f"split on {rule}")
    if pos.sum() < MIN_PER_GROUP or (~pos).sum() < MIN_PER_GROUP:
        print(f"    {label}: {pos.sum()} positive of {n_cells}, "
              "too few on one side to compare")
        return None

    mean, det = _column_stats(cpm, genes)
    keep = (det >= MIN_DETECTION) & (mean >= MIN_CPM)
    keep_idx = np.where(keep.to_numpy())[0]
    print(f"    {label}: {n_cells:,} neurons, {pos.sum():,} Oprl1-positive, "
          f"{n_strata} depth strata, {len(keep_idx):,} genes tested")

    num = np.zeros(len(keep_idx))
    den = 0.0
    used = 0
    for s in np.unique(stratum):
        in_s = stratum == s
        a, b = in_s & pos, in_s & ~pos
        if a.sum() < MIN_PER_GROUP or b.sum() < MIN_PER_GROUP:
            continue
        if dense:
            m1 = cpm.to_numpy()[a][:, keep_idx].mean(axis=0)
            m0 = cpm.to_numpy()[b][:, keep_idx].mean(axis=0)
        else:
            m1 = np.asarray(cpm[a][:, keep_idx].mean(axis=0)).ravel()
            m0 = np.asarray(cpm[b][:, keep_idx].mean(axis=0)).ravel()
        w = float(min(a.sum(), b.sum()))
        num += w * np.log2((m1 + PSEUDO_CPM) / (m0 + PSEUDO_CPM))
        den += w
        used += 1
    if den == 0:
        print(f"    {label}: no depth stratum has {MIN_PER_GROUP} cells "
              "on both sides")
        return None
    print(f"      {used} of {n_strata} strata usable")
    return pd.Series(num / den, index=np.asarray(genes)[keep_idx], name=label)


# ------------------------------------------------------------------ aggregation
def bootstrap_signature(wide, n_boot=N_BOOT, seed=0):
    """Resample the datasets, not the cells.

    `wide` is genes x populations of log2 fold changes with gaps where a gene
    was not testable. Each draw takes a population sample with replacement and
    recomputes every gene's median over whatever it has in that draw, so a gene
    carried by one dataset loses its median as soon as that dataset is missed.
    """
    rng = np.random.default_rng(seed)
    cols = wide.columns.to_numpy()
    values = wide.to_numpy()
    medians = np.empty((n_boot, values.shape[0]), dtype=np.float32)
    top = np.zeros(values.shape[0], dtype=np.int64)
    for b in range(n_boot):
        pick = rng.integers(0, len(cols), len(cols))
        draw = values[:, pick]
        with np.errstate(invalid="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            med = np.nanmedian(draw, axis=1)
        medians[b] = med
        order = np.argsort(-np.nan_to_num(med, nan=-np.inf))[:TOP_N]
        top[order] += 1
    lo, hi = np.nanpercentile(medians, [2.5, 97.5], axis=0)
    return pd.DataFrame({
        "median_log2fc": np.nanmedian(values, axis=1).round(3),
        "ci_lo": np.round(lo, 3), "ci_hi": np.round(hi, 3),
        "top30_support": np.round(top / n_boot, 4),
    }, index=wide.index)


# ---------------------------------------------------------------------- figure
def figure2(sig, wide, meta, raw_panel):
    """Four views of one result, in the order the argument is made.

    Panel a is the signature. Panel b shows it holding population by population.
    Panel c is the control that decides whether to believe it: the negative tail
    is the non-neuronal compartment, and the positive tail is not. Panel d gives
    the raw levels the fold changes were computed from.
    """
    st.set_theme()
    up = sig.head(20)
    down = sig.tail(10)
    fig = plt.figure(figsize=(18.5, 12.4))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.0, 1.15, 1.05],
                          height_ratios=[1.0, 0.62],
                          wspace=0.62, hspace=0.42,
                          left=0.055, right=0.985, top=0.90, bottom=0.115)

    def bars(ax, frame, color, letter, title):
        f = frame.iloc[::-1]
        y = np.arange(len(f))
        # A gene more abundant in the glia of the same ganglion is drawn open:
        # it is a contamination candidate rather than a partner.
        glial = (f.neuron_over_glia_log2 < 0).to_numpy()
        ax.barh(y, f.median_log2fc,
                color=[("white" if g else color) for g in glial],
                edgecolor="black", linewidth=1.0, height=0.72, zorder=3,
                hatch=None)
        ax.errorbar(f.median_log2fc, y,
                    xerr=[f.median_log2fc - f.ci_lo, f.ci_hi - f.median_log2fc],
                    fmt="none", ecolor="black", elinewidth=1.0, capsize=2.5,
                    zorder=4)
        ax.set_yticks(y)
        ax.set_yticklabels(f.index, fontstyle="italic", fontsize=10)
        ax.axvline(0, color="black", lw=0.8)
        ax.set_xlabel("log2 fold change, Oprl1+ over Oprl1- neurons", fontsize=10)
        st.panel_letter(ax, letter, dx=-0.40)
        ax.set_title(title, fontsize=11, pad=10)

    bars(fig.add_subplot(gs[0, 0]), up, st.BAR_BLUE, "a",
         "The 20 genes highest in Oprl1+ neurons\nmedian of 13 populations, 95% "
         "interval from\n10,000 resamples of the populations")
    bars(fig.add_subplot(gs[1, 0]), down, st.BAR_GREY, "c",
         "The 10 lowest\nopen bars are more abundant in glia")

    ax = fig.add_subplot(gs[:, 1])
    keys = list(up.index) + list(down.index)
    block = wide.loc[keys, meta.index]
    v = float(np.nanpercentile(np.abs(block.to_numpy()), 98))
    im = ax.imshow(block.to_numpy(), cmap="RdBu_r", vmin=-v, vmax=v,
                   aspect="auto")
    ax.set_xticks(range(len(meta)))
    ax.set_xticklabels(meta.label, rotation=90, fontsize=8.5)
    ax.set_yticks(range(len(block)))
    ax.set_yticklabels(block.index, fontstyle="italic", fontsize=9.5)
    ax.axhline(len(up) - 0.5, color="black", lw=2.0)
    for x, div in enumerate(meta.division):
        if x and div != meta.division.iloc[x - 1]:
            ax.axvline(x - 0.5, color="black", lw=1.6)
    fig.colorbar(im, ax=ax, fraction=0.030, pad=0.02, label="log2 fold change")
    st.panel_letter(ax, "b", dx=-0.26)
    ax.set_title("Every population on its own, divisions ruled apart\n"
                 "white where the gene was not testable there", fontsize=11,
                 pad=10)

    ax = fig.add_subplot(gs[0, 2])
    x = sig.median_log2fc.to_numpy()
    y = sig.neuron_over_glia_log2.to_numpy()
    ax.scatter(x, y, s=5, c="#CCCCCC", linewidths=0, rasterized=True)
    ax.scatter(up.median_log2fc, up.neuron_over_glia_log2, s=34,
               c=st.BAR_BLUE, edgecolors="black", linewidths=0.6, zorder=3)
    ax.scatter(down.median_log2fc, down.neuron_over_glia_log2, s=34,
               c=st.BAR_GREY, edgecolors="black", linewidths=0.6, zorder=3)
    for g in ["Dcn", "Col3a1", "Igfbp7", "Mpz"]:
        if g in sig.index:
            ax.annotate(g, (sig.median_log2fc[g], sig.neuron_over_glia_log2[g]),
                        fontsize=8.5, fontstyle="italic", xytext=(4, -9),
                        textcoords="offset points")
    for g in ["Grm7", "Epha10", "Caln1"]:
        if g in sig.index:
            ax.annotate(g, (sig.median_log2fc[g], sig.neuron_over_glia_log2[g]),
                        fontsize=8.5, fontstyle="italic", xytext=(4, 4),
                        textcoords="offset points")
    ax.axhline(0, color="black", lw=0.8, ls="--")
    ax.axvline(0, color="black", lw=0.8, ls="--")
    ax.set_xlabel("log2 fold change, Oprl1+ over Oprl1- neurons", fontsize=10)
    ax.set_ylabel("log2 neuron over non-neuron", fontsize=10)
    st.panel_letter(ax, "d", dx=-0.24)
    ax.set_title(f"All {len(sig):,} genes tested\n"
                 "the negative tail is the non-neuronal compartment",
                 fontsize=11, pad=10)

    ax = fig.add_subplot(gs[1, 2])
    r = raw_panel.iloc[::-1]
    y = np.arange(len(r))
    floor = 0.05
    ax.barh(y - 0.19, np.maximum(r.neg_cpm, floor), height=0.36,
            color=st.BAR_GREY, edgecolor="black", linewidth=0.9,
            label="Oprl1- neurons", zorder=3)
    ax.barh(y + 0.19, np.maximum(r.pos_cpm, floor), height=0.36,
            color=st.BAR_BLUE, edgecolor="black", linewidth=0.9,
            label="Oprl1+ neurons", zorder=3)
    ax.set_xscale("log")
    ax.set_xlim(floor, float(max(r.pos_cpm.max(), r.neg_cpm.max())) * 2.2)
    ax.set_yticks(y)
    ax.set_yticklabels(r.index, fontstyle="italic", fontsize=10)
    ax.set_xlabel("Mean expression (CPM, log scale)", fontsize=10)
    ax.legend(fontsize=9, loc="lower right")
    st.panel_letter(ax, "e", dx=-0.30)
    ax.set_title("Raw levels behind the fold changes\n"
                 "vestibular ganglion, 6,596 neurons\n"
                 "Snap25 and Prph are the flat controls", fontsize=11, pad=10)

    fig.suptitle("The Oprl1 signature across the peripheral nervous system",
                 fontsize=16, y=0.965)
    st.save(fig, "figure2_oprl1_signature")


# ------------------------------------------------------------------------ main
def main() -> int:
    ac.RES.mkdir(parents=True, exist_ok=True)
    series, meta, ambient, raw_store = [], [], [], {}
    for acc, tissue, division in POPULATIONS:
        label = f"{acc} {tissue}"
        print(f"  {label}")
        cpm, counts, genes, other = neurons_of(acc, tissue)
        fc = stratified_log2fc(cpm, counts, genes, label)
        if fc is None:
            continue
        series.append(fc)
        if other is not None:
            neuron_mean, _ = _column_stats(cpm, genes)
            ambient.append(np.log2((neuron_mean + PSEUDO_CPM)
                                   / (other.reindex(genes) + PSEUDO_CPM))
                           .rename(label))
        meta.append({"label": f"{tissue} ({acc})", "acc": acc,
                     "tissue": tissue, "division": division,
                     "n": int(cpm.shape[0])})
        if acc == "GSE309608":
            raw_store["cpm"], raw_store["counts"] = cpm, counts
            raw_store["genes"] = genes

    wide = pd.concat(series, axis=1)
    meta = pd.DataFrame(meta)
    meta.index = wide.columns
    order = {"sensory": 0, "sympathetic": 1, "parasympathetic": 2,
             "mixed autonomic": 3, "enteric": 4}
    meta = meta.sort_values("division", key=lambda s: s.map(order))
    wide = wide[meta.index]
    print(f"\n  {wide.shape[0]:,} genes tested in at least one of "
          f"{wide.shape[1]} populations")

    n_seen = wide.notna().sum(axis=1)
    wide = wide[n_seen >= MIN_DATASETS]
    print(f"  {wide.shape[0]:,} genes testable in at least {MIN_DATASETS}")

    print(f"  bootstrapping {N_BOOT:,} dataset resamples ...")
    sig = bootstrap_signature(wide)
    sig["n_datasets"] = wide.notna().sum(axis=1)
    sig["n_positive"] = (wide > 0).sum(axis=1)
    sig["n_deposits"] = [
        meta.loc[wide.columns[~wide.loc[g].isna().to_numpy()], "acc"].nunique()
        for g in wide.index]
    # Neuron over non-neuron for every gene in the signature, median across the
    # datasets that have a non-neuronal compartment. A signature gene that is
    # more abundant in glia is a contamination candidate, not a partner.
    amb = pd.concat(ambient, axis=1).median(axis=1).round(3)
    ac.save_table(amb.rename("log2_neuron_over_glia").reset_index()
                  .rename(columns={"index": "gene"}),
                  "oprl1_signature_ambient.csv")
    sig["neuron_over_glia_log2"] = amb.reindex(sig.index).to_numpy()
    sig = sig.sort_values("median_log2fc", ascending=False)

    control = sig.loc[["Oprl1"]] if "Oprl1" in sig.index else sig.iloc[:0]
    sig = sig.drop(index="Oprl1", errors="ignore")
    ac.save_table(pd.concat([control, sig]), "oprl1_signature.csv", index=True)
    ac.save_table(wide.round(3), "oprl1_signature_by_dataset.csv", index=True)

    if len(control):
        print(f"\n  positive control, Oprl1 itself: "
              f"{float(control.median_log2fc.iloc[0]):+.2f} log2, "
              f"support {float(control.top30_support.iloc[0]):.3f}")
    print("\n  The Oprl1 signature, top 25:")
    print(sig.head(25).round(3).to_string())
    print("\n  The other end, 10 genes lowest in Oprl1-positive neurons:")
    print(sig.tail(10).round(3).to_string())

    # Genes this project already has a position on, wherever they land.
    landmarks = ["Scn1a", "Scn10a", "Pvalb", "Trpv1", "Piezo2", "Glp1r",
                 "Cckar", "Calca", "Prph", "Snap25", "Plp1", "Sox10", "Pnoc",
                 "Oprm1", "Oprk1", "Oprd1"]
    here = sig.reindex([g for g in landmarks if g in sig.index])
    here.insert(0, "rank", [int(sig.index.get_loc(g)) + 1 for g in here.index])
    print(f"\n  Landmarks, out of {len(sig):,} genes in the signature:")
    print(here.round(3).to_string())

    n_glial = int((sig.head(TOP_N).neuron_over_glia_log2 < 0).sum())
    print(f"\n  Ambient: {TOP_N - n_glial} of the top {TOP_N} are more "
          f"abundant in neurons than in the non-neuronal cells of the same "
          f"ganglion; {n_glial} are not.")

    top8 = list(sig.head(8).index) + ["Snap25", "Prph"]
    raw = []
    if raw_store:
        cpm, counts, genes = (raw_store["cpm"], raw_store["counts"],
                              raw_store["genes"])
        j = list(genes).index("Oprl1")
        pos = np.asarray(counts[:, j].todense()).ravel() > 0
        for g in top8:
            k = list(genes).index(g)
            col = np.asarray(cpm[:, k].todense()).ravel()
            raw.append({"gene": g, "pos_cpm": round(float(col[pos].mean()), 2),
                        "neg_cpm": round(float(col[~pos].mean()), 2)})
    raw_panel = pd.DataFrame(raw).set_index("gene")
    ac.save_table(raw_panel.reset_index(), "oprl1_signature_raw_levels.csv")

    figure2(sig, wide, meta, raw_panel)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
