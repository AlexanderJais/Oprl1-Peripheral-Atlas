"""Does the mouse result hold in human peripheral ganglia? Not yet decidable.

Every claim in this atlas so far is a mouse claim. This script asks the same
question of human tissue and finds that the available data cannot answer it, for
a reason that is itself worth reporting.

THE PREPARATION PROBLEM, MEASURED A SECOND TIME

Human peripheral ganglia come from post-mortem or surgical tissue that is
snap-frozen, so essentially all human single-cell data from them is
single-nucleus. Figure S1 measured what that does in mouse nodose tissue, the
one ganglion for which both preparations exist: nuclear libraries retain
unspliced pre-mRNA, so Oprm1 at 250 kb gains 29-fold while Oprl1 at 6 kb loses
half its signal, which moves Oprm1 from second to first.

GSE201654 tests that in a second tissue and a second laboratory. It is a
cross-species atlas of dorsal root ganglion nuclei: human, macaque, mouse and
guinea pig, dissected and sequenced by one group on one protocol. The mouse arm
is the control, because mouse DRG has whole-cell data too. Whole-cell mouse DRG
places Oprl1 first at 9.31 CPM against Oprm1 8.23. Nuclear mouse DRG, same
tissue, places Oprm1 first at 94.92 against Oprl1 7.60. The inversion is
reproduced.

That is what the human rows have to be read against. A human nuclear dataset
placing OPRM1 first has not established a species difference, because the mouse
does the same thing under the same preparation.

WHAT THE UNBIASED HUMAN DATA SAYS

Bulk RNA-seq of whole ganglion tissue escapes the problem, since it sequences
cytoplasmic RNA along with nuclear. GSE231763 is bulk RNA-seq of six human
superior cervical ganglia and is the only non-nuclear human peripheral ganglion
source located for this project. It places OPRL1 last of the four, at 0.45 TPM
against OPRM1 3.01, and at exactly zero in three of six donors. It cannot
separate neurons from glia, which is the trade, and in mouse superior cervical
neurons Oprl1 sits only 1.7-fold above the non-neuronal compartment, so dilution
is a live concern. It is one ganglion from one laboratory and it is the only
evidence of its kind.

The honest summary is that the mouse result is not established in human, the
single-nucleus data cannot settle it either way, and the one unbiased human
dataset points against it.

  GSE201654   dorsal root ganglion nuclei, four species, one protocol
  GSE241386   human stellate ganglion, single-nucleus 10x
  GSE231763   human superior cervical ganglion, six samples, bulk RNA-seq

Run from the repository root. Sources are fetched into `scratch/` on first use.
"""

import glob
import gzip
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

import atlas_common as ac
import atlas_style as st

SCRATCH = Path("scratch") / "human"
RECEPTORS = ["OPRL1", "OPRM1", "OPRD1", "OPRK1"]

# Nuclei carry less RNA than whole cells, so the droplet thresholds used
# elsewhere in this project would reject most of a good nuclear library.
MIN_UMI, MIN_GENES = 1000, 500
# Sensory neuron identity across four species, on raw counts. Any one of the
# three is enough: SNAP25 alone recovered 0.2% of human DRG nuclei, which is
# below the neuronal fraction the source publication reports.
NEURON_ANY = {"SNAP25": 1, "PRPH": 1, "SCN9A": 1}
NOT_NEURON_CPM = {"PLP1": 500}

# Mouse DRG whole cell, from the iPain rows of section 1, for the one comparison
# that decides how the human rows can be read.
MOUSE_WHOLE_CELL = {"OPRL1": 9.31, "OPRM1": 8.23, "OPRD1": 0.20, "OPRK1": 2.78}

# Explicit accession ranges, not string prefixes: the macaque samples are
# GSM6068077 to GSM6068092 and the human samples GSM6068093 to GSM6068110, so a
# prefix rule on "GSM60680" silently labels human nuclei as macaque.
SPECIES_OF = {
    **{f"GSM{i}": "mouse" for i in range(6069069, 6069079)},
    **{f"GSM{i}": "macaque" for i in range(6068077, 6068093)},
    **{f"GSM{i}": "guinea pig" for i in range(6067401, 6067409)},
    **{f"GSM{i}": "human" for i in range(6068093, 6068111)},
}


def species_of(path):
    stem = os.path.basename(path).split("_")[0]
    if stem not in SPECIES_OF:
        raise ac.SanityCheckError(f"{stem} is not in the species map; a sample "
                                  "of unknown species cannot enter the "
                                  "comparison")
    return SPECIES_OF[stem]


# ------------------------------------------------------------------- loading
def load_h5(path):
    with h5py.File(path, "r") as f:
        g = f["matrix"]
        genes = pd.Index([s.decode().upper() for s in g["features"]["name"][:]])
        x = sparse.csc_matrix((g["data"][:], g["indices"][:], g["indptr"][:]),
                              shape=g["shape"][:]).T.tocsr().astype(np.float32)
    return x, genes


def load_mtx(d, sample):
    with gzip.open(d / f"{sample}_features.tsv.gz", "rt") as fh:
        genes = pd.Index([line.split("\t")[1].upper() for line in fh])
    with gzip.open(d / f"{sample}_matrix.mtx.gz", "rt") as fh:
        return mmread(fh).tocsc().T.tocsr().astype(np.float32), genes


def nuclei(x, genes):
    """QC, neuron call and per-cell CPM for one nuclear sample."""
    umi = np.asarray(x.sum(axis=1)).ravel()
    n_gene = np.asarray((x > 0).sum(axis=1)).ravel()
    keep = (umi >= MIN_UMI) & (n_gene >= MIN_GENES)
    x, umi = x[keep], umi[keep]

    def col(w):
        j = np.where(genes == w)[0]
        return (np.zeros(x.shape[0]) if len(j) == 0
                else np.asarray(x[:, j].sum(axis=1)).ravel())

    neuron = np.zeros(x.shape[0], bool)
    for g, t in NEURON_ANY.items():
        neuron |= col(g) >= t
    for g, t in NOT_NEURON_CPM.items():
        neuron &= (col(g) / umi * 1e6) < t

    want = sorted(set(RECEPTORS) | set(NEURON_ANY) | set(NOT_NEURON_CPM)
                  | {"TUBB3", "ACTB"})
    cpm = pd.DataFrame({w: col(w) / umi * 1e6 for w in want})
    counts = pd.DataFrame({w: col(w) for w in want})
    return cpm, counts, neuron, umi, int(keep.sum())


def population(files, label, loader=load_h5):
    """One row per dataset, plus the per-sample ordering."""
    cpms, countss, neurons, umis, sample_of, n_qc = [], [], [], [], [], 0
    for name, path in files:
        x, genes = loader(path) if loader is load_h5 else loader(*path)
        cpm, counts, neuron, umi, kept = nuclei(x, genes)
        cpms.append(cpm)
        countss.append(counts)
        neurons.append(neuron)
        umis.append(umi)
        sample_of += [name] * len(cpm)
        n_qc += kept
        print(f"    {name}: {kept:,} nuclei pass QC, {neuron.sum():,} neuronal "
              f"({neuron.mean() * 100:.1f}%), median {np.median(umi):,.0f} UMI")

    cpm = pd.concat(cpms, ignore_index=True)
    counts = pd.concat(countss, ignore_index=True)
    neuron = np.concatenate(neurons)
    umi = np.concatenate(umis)
    sample_of = np.array(sample_of)

    print(f"  {label}: {neuron.sum():,} neuronal nuclei of {n_qc:,}, markers "
          "(CPM) " + "  ".join(f"{g} {cpm.loc[neuron, g].mean():,.0f}"
                               for g in ("SNAP25", "PRPH", "SCN9A", "PLP1")))
    b = ac.bootstrap_receptor_support(
        cpm.loc[neuron, RECEPTORS].rename(columns=dict(zip(RECEPTORS, ac.RECEPTORS))),
        n_boot=ac.N_BOOT)
    back = dict(zip(ac.RECEPTORS, RECEPTORS))

    first = {}
    for s in sorted(set(sample_of)):
        m = neuron & (sample_of == s)
        if m.sum() < 30:
            continue
        lv = {g: float(cpm.loc[m, g].mean()) for g in RECEPTORS}
        first[s] = max(lv, key=lv.get)
    tops = pd.Series(first).value_counts()

    row = {"dataset": label, "prep": "single nucleus", "n": int(neuron.sum()),
           "n_samples": len(set(sample_of)),
           "median_umi": int(np.median(umi[neuron])),
           **{g: round(float(cpm.loc[neuron, g].mean()), 3) for g in RECEPTORS},
           **{f"{g}_pct": round(float((counts.loc[neuron, g] > 0).mean() * 100), 2)
              for g in RECEPTORS},
           "top_gene": back[b["top_gene"]], "runner_up": back[b["runner_up"]],
           "margin": round(b["margin"], 3),
           "margin_lo": round(b["margin_lo"], 3),
           "margin_hi": round(b["margin_hi"], 3),
           "support": b["support"], "n_boot": b["n_boot"],
           "samples_by_top": "; ".join(f"{k} {v}" for k, v in tops.items())}
    print(f"    top {row['top_gene']} over {row['runner_up']} "
          f"{row['margin']:.2f}x, support {row['support']:.4f}, "
          f"per sample: {row['samples_by_top']}")
    return row


def load_bulk():
    """STAR gene abundance tables from bulk human superior cervical ganglia."""
    files = sorted(glob.glob(str(SCRATCH / "bulk" / "*.gtf.gz")))
    if not files:
        raise FileNotFoundError(f"no bulk tables under {SCRATCH / 'bulk'}")
    frames = {}
    for f in files:
        d = pd.read_csv(f, sep="\t").groupby("Gene Name")["TPM"].sum()
        frames[os.path.basename(f).split("_")[0]] = d
    return pd.DataFrame(frames)


# ---------------------------------------------------------------------- figure
def figure(human, xspecies, bulk):
    st.set_theme()
    fig = plt.figure(figsize=(17.0, 9.6))
    gs = fig.add_gridspec(2, 3, wspace=0.46, hspace=0.62,
                          left=0.06, right=0.985, top=0.88, bottom=0.09)

    def bars(ax, values, title, unit="Mean expression (CPM)"):
        order = np.argsort(-np.array([values[g] for g in RECEPTORS]))
        genes = [RECEPTORS[i] for i in order]
        st.expression_bars(
            ax, [values[g] for g in genes], genes, unit,
            colors=[st.BAR_BLUE if g == "OPRL1" else st.BAR_GREY for g in genes],
            fontsize=11)
        ax.set_title(title, fontsize=10.5, pad=8)

    for i, r in enumerate(human.itertuples()):
        ax = fig.add_subplot(gs[0, i])
        bars(ax, {g: getattr(r, g) for g in RECEPTORS},
             f"{r.dataset}\n{r.n:,} nuclei, single nucleus")
        st.panel_letter(ax, "abc"[i], dx=-0.30)

    ax = fig.add_subplot(gs[0, 2])
    b = bulk.loc[[g for g in RECEPTORS if g in bulk.index]]
    x = np.arange(len(b))
    ax.bar(x, b.mean(axis=1),
           color=[st.BAR_BLUE if g == "OPRL1" else st.BAR_GREY for g in b.index],
           edgecolor="black", linewidth=1.4, width=0.66, zorder=2)
    for j, c in enumerate(b.columns):
        ax.scatter(x + (j - len(b.columns) / 2) * 0.07, b[c], s=26, c="black",
                   zorder=4)
    ax.set_xticks(x)
    ax.set_xticklabels(b.index, rotation=45, ha="right", fontsize=11,
                       fontstyle="italic")
    ax.set_ylabel("TPM", fontsize=11)
    st.panel_letter(ax, "c", dx=-0.28)
    ax.set_title("Human superior cervical ganglion\n"
                 f"{b.shape[1]} donors, bulk RNA-seq, not nuclear",
                 fontsize=10.5, pad=8)

    # The control: one tissue, two preparations, same species.
    ax = fig.add_subplot(gs[1, 0])
    bars(ax, MOUSE_WHOLE_CELL, "Mouse DRG, whole cell\niPain atlas, 31,802 cells")
    st.panel_letter(ax, "d", dx=-0.30)
    ax = fig.add_subplot(gs[1, 1])
    mouse = xspecies[xspecies.species == "mouse"].iloc[0]
    bars(ax, {g: mouse[g] for g in RECEPTORS},
         f"Mouse DRG, single nucleus\nGSE201654, {int(mouse['n']):,} nuclei")
    st.panel_letter(ax, "e", dx=-0.30)

    ax = fig.add_subplot(gs[1, 2])
    order = ["guinea pig", "macaque", "mouse", "human"]
    present = [s for s in order if s in set(xspecies.species)]
    width = 0.2
    for k, g in enumerate(RECEPTORS):
        vals = [float(xspecies.loc[xspecies.species == s, g].iloc[0])
                for s in present]
        ax.bar(np.arange(len(present)) + (k - 1.5) * width, vals, width=width,
               color=st.BAR_BLUE if g == "OPRL1" else ["#B2182B", "#4D4D4D",
                                                       "#B8B8B8"][k - 1],
               edgecolor="black", linewidth=0.8, label=g, zorder=3)
    ax.set_xticks(np.arange(len(present)))
    ax.set_xticklabels(present, fontsize=10)
    ax.set_ylabel("Mean expression (CPM)", fontsize=11)
    ax.set_yscale("log")
    ax.legend(fontsize=8.5, ncol=2)
    st.panel_letter(ax, "f", dx=-0.26)
    ax.set_title("Four species, one laboratory, one nuclear protocol\n"
                 "GSE201654 dorsal root ganglion", fontsize=10.5, pad=8)

    fig.suptitle("Human peripheral ganglia, and why single-nucleus data cannot "
                 "settle the ordering", fontsize=14, y=0.955)
    st.save(fig, "figure8_human_ganglia")


# ------------------------------------------------------------------------ main
def main() -> int:
    ac.RES.mkdir(parents=True, exist_ok=True)

    print("\n  Human dorsal root ganglion (GSE201654, single nucleus):")
    drg = [(os.path.basename(f).split("_")[0], f)
           for f in sorted(glob.glob(str(SCRATCH / "drg" / "*.h5")))]
    human_drg = population(drg, "GSE201654 human DRG")

    print("\n  Human stellate ganglion (GSE241386, single nucleus):")
    human_sg = population([("GSM7728058", (SCRATCH / "sg", "GSM7728058"))],
                          "GSE241386 human stellate", loader=load_mtx)

    print("\n  The same protocol in three other species (GSE201654):")
    rows, per_sample = [], None
    for f in sorted(glob.glob(str(SCRATCH / "xspecies" / "*.h5"))):
        x, genes = load_h5(f)
        cpm, counts, neuron, umi, kept = nuclei(x, genes)
        if neuron.sum() < 30:
            continue
        rows.append({"species": species_of(f),
                     "sample": os.path.basename(f).split("_")[0],
                     "n": int(neuron.sum()),
                     **{g: round(float(cpm.loc[neuron, g].mean()), 3)
                        for g in RECEPTORS},
                     **{g: round(float(cpm.loc[neuron, g].mean()), 1)
                        for g in ("SNAP25", "PRPH", "SCN9A", "PLP1")},
                     # Neuron over non-neuron inside each sample, the check that
                     # says whether the called population is what it claims.
                     **{f"{g}_enr": round(float(np.log2(
                         (cpm.loc[neuron, g].mean() + 0.01)
                         / (cpm.loc[~neuron, g].mean() + 0.01))), 2)
                        for g in RECEPTORS}})
    for f in sorted(glob.glob(str(SCRATCH / "drg" / "*.h5"))):
        x, genes = load_h5(f)
        cpm, counts, neuron, umi, kept = nuclei(x, genes)
        if neuron.sum() < 30:
            continue
        rows.append({"species": species_of(f),
                     "sample": os.path.basename(f).split("_")[0],
                     "n": int(neuron.sum()),
                     **{g: round(float(cpm.loc[neuron, g].mean()), 3)
                        for g in RECEPTORS},
                     **{g: round(float(cpm.loc[neuron, g].mean()), 1)
                        for g in ("SNAP25", "PRPH", "SCN9A", "PLP1")},
                     **{f"{g}_enr": round(float(np.log2(
                         (cpm.loc[neuron, g].mean() + 0.01)
                         / (cpm.loc[~neuron, g].mean() + 0.01))), 2)
                        for g in RECEPTORS}})
    per_sample = pd.DataFrame(rows)
    if per_sample.empty:
        raise ac.SanityCheckError("no cross-species sample passed QC")
    per_sample["top"] = per_sample[RECEPTORS].idxmax(axis=1)
    ac.save_table(per_sample, "human_xspecies_per_sample.csv")

    marker = per_sample.groupby("species")[
        ["SNAP25", "PRPH", "SCN9A", "PLP1"] +
        [f"{g}_enr" for g in RECEPTORS]].mean().round(2)
    ac.save_table(marker.reset_index(), "human_xspecies_markers.csv")
    print("\n  What the called neurons look like in each species:")
    print(marker.to_string())
    snap = marker["SNAP25"]
    if snap.get("human", np.inf) < 0.25 * snap.drop("human", errors="ignore").median():
        print("\n  [warn] SNAP25 in the human population is far below the other "
              "species under the\n  same protocol, so the human neuron call is "
              "not trustworthy and no human\n  receptor ordering is reported "
              "from it.")

    xspecies = per_sample.groupby("species")[RECEPTORS + ["n"]].mean().round(3)
    xspecies["n_samples"] = per_sample.groupby("species").size()
    xspecies["top_in_samples"] = per_sample.groupby("species")["top"].apply(
        lambda s: "; ".join(f"{k} {v}" for k, v in s.value_counts().items()))
    order = [s for s in ("guinea pig", "macaque", "mouse", "human")
             if s in xspecies.index]
    xspecies = xspecies.loc[order]
    xspecies.index.name = "species"
    ac.save_table(xspecies.reset_index(), "human_xspecies_drg.csv")
    print("\n  Dorsal root ganglion nuclei, four species, one protocol:")
    print(xspecies.to_string())
    print(f"\n  Mouse DRG whole cell, for comparison: "
          + "  ".join(f"{g} {v}" for g, v in MOUSE_WHOLE_CELL.items()))
    print("  The preparation inverts the mouse ordering in the DRG as it does "
          "in the nodose\n  ganglion (figure S1), so no human nuclear row can "
          "establish a species difference.")

    print("\n  Human superior cervical ganglion (GSE231763, bulk RNA-seq):")
    bulk = load_bulk()
    present = [g for g in RECEPTORS if g in bulk.index]
    b = bulk.loc[present]
    ac.save_table(b.reset_index().rename(columns={"Gene Name": "gene"}),
                  "human_scg_bulk_tpm.csv")
    print(b.round(3).to_string())
    top = b.mean(axis=1).idxmax()
    runner = b.mean(axis=1).drop(top).idxmax()
    n_first = int((b.idxmax(axis=0) == top).sum())
    print(f"  identity of the tissue: " + "  ".join(
        f"{g} {bulk.loc[g].mean():,.0f} TPM"
        for g in ("SNAP25", "TH", "DBH", "PRPH", "NPY") if g in bulk.index))
    print(f"  top {top} over {runner} "
          f"{b.loc[top].mean() / max(b.loc[runner].mean(), 1e-9):.2f}x, first in "
          f"{n_first} of {b.shape[1]} donors; OPRL1 mean {b.loc['OPRL1'].mean():.3f} "
          f"TPM and zero in {int((b.loc['OPRL1'] == 0).sum())} of {b.shape[1]}")

    bulk_row = {"dataset": "GSE231763 human SCG", "prep": "bulk RNA-seq",
                "n": np.nan, "n_samples": b.shape[1], "median_umi": np.nan,
                **{g: round(float(b.loc[g].mean()), 3) for g in present},
                "top_gene": top, "runner_up": runner,
                "margin": round(float(b.loc[top].mean()
                                      / max(b.loc[runner].mean(), 1e-9)), 3),
                "samples_by_top": f"{top} {n_first} of {b.shape[1]}"}
    table = pd.DataFrame([human_drg, human_sg, bulk_row])
    ac.save_table(table, "human_receptor_levels.csv")
    print("\n  Summary:")
    print(table[["dataset", "prep", "n", "n_samples"] + RECEPTORS +
                ["top_gene", "margin", "samples_by_top"]].to_string(index=False))

    figure(pd.DataFrame([human_drg, human_sg]), xspecies.reset_index(), bulk)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
