"""When the Oprl1 ordering appears, and which spiral neuron carries it.

Sections 1 to 3 describe adult tissue. Two questions they do not answer:

  WHEN. Is Oprl1 dominance present as soon as the neurons differentiate, or is
  it acquired later? The otic lineage is the one place this can be followed end
  to end in a single tissue. GSE178931 covers otic neurogenesis at E9.5, E11.5
  and E13.5 on droplet; GSE165502 covers the cochlea from E14.5 to P3 on
  SMART-seq2; GSE114997 is the same ganglion at P25 on SMART-seq. The enteric
  result in section 1 already shows that age can reverse the ordering, so this
  is not a formality.

  WHICH NEURON. The spiral ganglion splits into type Ia, Ib and Ic, which differ
  in spontaneous rate and threshold. GSE114997 has 186 wild-type neurons at a
  median 2.7 million reads, enough to assign subtypes and ask whether Oprl1
  tracks Scn1a inside one ganglion rather than across ganglia.

Subtypes are assigned in two independent ways and the assignment is only used
where they agree: unsupervised clustering of the 186 neurons, and the published
markers Calb2 (Ia), Calb1 (Ib) and Lypd1 with Pou4f1 (Ic). The GEO deposit
carries no subtype labels, so an assignment made here has to show its working.

Run from the repository root. Sources are fetched into `scratch/` on first use.
"""

import gzip
import io
import subprocess
import sys
import tarfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import atlas_common as ac
import atlas_style as st

SCRATCH = Path("scratch") / "otic"
FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"

SOURCES = {
    "GSE165502_counts.csv.gz":
        f"{FTP}/GSE165nnn/GSE165502/suppl/GSE165502_counts.csv.gz",
    "GSE114997_Shrestha_raw_counts.csv.gz":
        f"{FTP}/GSE114nnn/GSE114997/suppl/GSE114997_Shrestha_raw_counts.csv.gz",
    "GSE178931_RAW.tar":
        f"{FTP}/GSE178nnn/GSE178931/suppl/GSE178931_RAW.tar",
}

# Plate to age, read from the GEO sample titles of GSE165502.
PLATE_AGE = {
    "SS2_15_0135": "P3", "SS2_15_0144": "P3",
    "SS2_18_349": "E16.5", "SS2_18_350": "E16.5", "SS2_20_139": "E16.5",
    "SS2_19_166": "E18.5", "SS2_19_168": "E18.5",
    "SS2_20_133": "E14.5", "SS2_20_135": "E14.5",
    "SS2_20_151": "E15.5", "SS2_20_153": "E15.5",
    "SS2_20_143": "E17.5",
}
AGE_ORDER = ["E9.5", "E11.5", "E13.5", "E14.5", "E15.5", "E16.5", "E17.5",
             "E18.5", "P3", "P25"]

# Neuronal identity, on raw counts, per dataset. Snap25 is still low in the
# embryonic material, where Neurod1 marks the otic neuroblast, and by P3 the
# sorted SGNs carry Snap25 but have turned Neurod1 down. One rule across all
# three ages would either admit the E9.5 otocyst or exclude the P3 ganglion.
OTIC_NEURON = ["Neurod1", "Isl1", "Snap25", "Tubb3"]
NEURON_RULE = {
    "GSE178931": {"Neurod1": 3, "Tubb3": 3},
    "GSE165502": {"Tubb3": 5, "Snap25": 1},
    "GSE114997": {"Tubb3": 5},
}
MIN_GENES = 1000
MIN_CELLS_PER_AGE = 30

SUBTYPE_MARKERS = {"Ia": ["Calb2"], "Ib": ["Calb1"], "Ic": ["Lypd1", "Pou4f1"]}


def fetch(name):
    SCRATCH.mkdir(parents=True, exist_ok=True)
    path = SCRATCH / name
    if not path.exists():
        print(f"  downloading {name}")
        subprocess.run(["curl", "-fsSL", "-o", str(path), SOURCES[name]],
                       check=True)
    return path


# --------------------------------------------------- GSE165502, E14.5 to P3
def load_gse165502():
    """SMART-seq2 cochlea, cells x genes counts, with age from the plate id."""
    x = pd.read_csv(fetch("GSE165502_counts.csv.gz"), index_col=0)
    x = x.loc[:, ~x.columns.duplicated()]
    plate = pd.Index(x.index).str.split(":").str[0]
    age = pd.Series([PLATE_AGE.get(p) for p in plate], index=x.index)
    if age.isna().any():
        # The counts file carries one plate, SS2_20_141, that the series sample
        # list does not. Without an age it cannot enter a developmental series,
        # so it is dropped rather than guessed at.
        unknown = sorted(set(plate[age.isna()]))
        print(f"  [note] dropping {int(age.isna().sum()):,} cells from "
              f"{unknown}: no age in the GEO sample list")
        x, age = x[age.notna().to_numpy()], age.dropna()
    print(f"  GSE165502: {x.shape[0]:,} cells x {x.shape[1]:,} genes, "
          f"ages {sorted(set(age), key=AGE_ORDER.index)}")
    return x, age


# ------------------------------------------------- GSE178931, E9.5 to E13.5
def stream_targets(handle, targets, chunk=2000):
    """Target gene rows and per-cell totals from a genes-by-cells CSV.

    The E9.5 sample alone is 24,344 cells by 30,000 genes, which does not fit in
    memory as a dense frame and does not need to: this analysis reads four
    receptors and four identity markers. One pass keeps those rows and
    accumulates the library size and detected-gene count every cell needs for
    QC, so nothing is held that is not used.
    """
    kept, lib, n_gene, cells = [], None, None, None
    for block in pd.read_csv(handle, index_col=0, chunksize=chunk):
        if cells is None:
            cells = block.columns
            lib = np.zeros(len(cells))
            n_gene = np.zeros(len(cells), dtype=np.int64)
        values = block.to_numpy(dtype=np.float32)
        lib += values.sum(axis=0)
        n_gene += (values > 0).sum(axis=0)
        here = block.index.intersection(targets)
        if len(here):
            kept.append(block.loc[here])
    counts = (pd.concat(kept).groupby(level=0).sum().T if kept
              else pd.DataFrame(index=cells))
    counts = counts.reindex(columns=[t for t in targets if t in counts.columns])
    return counts, pd.Series(lib, index=cells), pd.Series(n_gene, index=cells)


def load_gse178931(targets):
    """Droplet otic tissue at three embryonic ages, one CSV per sample."""
    tar = fetch("GSE178931_RAW.tar")
    counts, libs, ngenes, ages = [], [], [], []
    with tarfile.open(tar) as t:
        for m in sorted(t.getmembers(), key=lambda x: x.name):
            if not m.name.endswith("_raw_count.csv.gz"):
                continue
            # GSM5401073_E11_5_1_raw_count.csv.gz -> E11.5
            parts = m.name.split("_")
            age = f"{parts[1]}.{parts[2]}"
            with gzip.open(io.BytesIO(t.extractfile(m).read()), "rt") as fh:
                c, lib, ng = stream_targets(fh, targets)
            c.index = [f"{m.name}:{i}" for i in c.index]
            lib.index, ng.index = c.index, c.index
            counts.append(c)
            libs.append(lib)
            ngenes.append(ng)
            ages += [age] * len(c)
            print(f"    {m.name}: {len(c):,} cells, {age}")
    counts = pd.concat(counts)
    return (counts, pd.concat(libs), pd.concat(ngenes),
            pd.Series(ages, index=counts.index))


# ---------------------------------------------------------- shared analysis
def neurons_and_cpm(x, label, lib=None, n_gene=None):
    """QC, neuron call and CPM for a cells-by-genes count matrix.

    `lib` and `n_gene` are passed in when the matrix was streamed and only
    carries the target columns, since a library size computed over those
    columns alone would not be a library size.
    """
    lib = (x.sum(axis=1).to_numpy(dtype=float) if lib is None
           else lib.reindex(x.index).to_numpy(dtype=float))
    n_gene = ((x > 0).sum(axis=1).to_numpy() if n_gene is None
              else n_gene.reindex(x.index).to_numpy())
    keep = (n_gene >= MIN_GENES) & (lib > 0)
    x, lib = x[keep], lib[keep]
    cpm = x.div(lib, axis=0) * 1e6

    rule = NEURON_RULE[label]
    neuron = np.ones(len(x), bool)
    for g, t in rule.items():
        if g not in x.columns:
            raise ac.SanityCheckError(f"{label}: {g} absent, neurons cannot "
                                      "be called")
        neuron &= x[g].to_numpy() >= t
    print(f"  {label}: {keep.sum():,} cells pass QC, {neuron.sum():,} neuronal "
          f"({neuron.mean() * 100:.0f}%) by "
          + " and ".join(f"{g} >= {t}" for g, t in rule.items())
          + f", median {np.median(lib):,.0f} counts")
    print("    markers in the neurons (CPM): " + "  ".join(
        f"{g} {cpm.loc[neuron, g].mean():,.0f}"
        for g in OTIC_NEURON if g in cpm.columns))
    return cpm[neuron], x[neuron], keep


def by_age(cpm, counts, age, label, rows):
    for a in sorted(set(age), key=AGE_ORDER.index):
        m = (age.reindex(cpm.index) == a).to_numpy()
        if m.sum() < MIN_CELLS_PER_AGE:
            print(f"    {a}: {m.sum()} neurons, below {MIN_CELLS_PER_AGE}")
            continue
        sub = cpm.loc[m, [g for g in ac.RECEPTORS if g in cpm.columns]]
        b = ac.bootstrap_receptor_support(sub, n_boot=ac.N_BOOT)
        row = {"dataset": label, "age": a, "n": int(m.sum()),
               **{g: round(float(cpm.loc[m, g].mean()), 3) if g in cpm else np.nan
                  for g in ac.RECEPTORS},
               **{f"{g}_pct": round(float((counts.loc[m, g] > 0).mean() * 100), 1)
                  if g in counts else np.nan for g in ac.RECEPTORS},
               "top_gene": b["top_gene"], "runner_up": b["runner_up"],
               "margin": round(b["margin"], 3),
               "margin_lo": round(b["margin_lo"], 3),
               "margin_hi": round(b["margin_hi"], 3),
               "support": b["support"], "n_boot": b["n_boot"]}
        rows.append(row)
        print(f"    {a}: n = {m.sum():4d}  top {b['top_gene']} over "
              f"{b['runner_up']} {b['margin']:.2f}x  support {b['support']:.3f}  "
              + "  ".join(f"{g} {row[g]:.2f}" for g in ac.RECEPTORS))
    return rows


# ------------------------------------------- GSE114997, the adult subtypes
def load_gse114997():
    """186 wild-type P25 spiral ganglion neurons, genes x cells raw counts."""
    x = pd.read_csv(fetch("GSE114997_Shrestha_raw_counts.csv.gz"), index_col=0)
    x = x.groupby(level=0).sum().T                     # cells x genes
    titles = geo_titles("GSE114997")
    genotype = pd.Series({c: titles.get(str(c), "") for c in x.index})
    wt = genotype.str.contains("Wildtype").to_numpy()
    print(f"  GSE114997: {x.shape[0]} cells, {wt.sum()} wild type "
          f"({(~wt).sum()} Vglut3-/- excluded)")
    return x[wt]


def geo_titles(acc):
    """Cell number to sample title, from the GEO sample records."""
    cache = SCRATCH / f"{acc}_titles.txt"
    if not cache.exists():
        subprocess.run(
            ["curl", "-fsSL", "-o", str(cache),
             f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}"
             "&targ=gsm&form=text&view=brief"], check=True)
    out = {}
    for line in cache.read_text(errors="ignore").splitlines():
        if line.startswith("!Sample_title"):
            t = line.split("= ", 1)[1].strip()
            out[t.split("_")[0]] = t
    return out


def spiral_subtypes(cpm, counts):
    """Assign Ia, Ib and Ic two ways and keep only the cells that agree."""
    genes = [g for m in SUBTYPE_MARKERS.values() for g in m if g in cpm.columns]
    missing = [g for m in SUBTYPE_MARKERS.values() for g in m
               if g not in cpm.columns]
    if missing:
        raise ac.SanityCheckError(f"subtype markers absent: {missing}")

    # 1. Marker call: the subtype whose markers are highest after scaling each
    # marker to its own spread, so one abundant marker cannot win by scale.
    z = np.log1p(cpm[genes])
    z = (z - z.mean()) / z.std(ddof=0)
    score = pd.DataFrame({k: z[[g for g in m]].mean(axis=1)
                          for k, m in SUBTYPE_MARKERS.items()})
    marker_call = score.idxmax(axis=1)

    # 2. Unsupervised: three clusters on the top principal components of the
    # 2,000 most variable genes, with no knowledge of the markers.
    lx = np.log1p(cpm)
    hv = lx.var(axis=0).sort_values(ascending=False).head(2000).index
    m = lx[hv].to_numpy()
    m = (m - m.mean(0)) / (m.std(0) + 1e-9)
    u, s, _ = np.linalg.svd(m - m.mean(0), full_matrices=False)
    pcs = u[:, :10] * s[:10]
    cluster = kmeans(pcs, k=3, seed=0)

    # Name each cluster by the marker score it maximises, then keep agreement.
    names = {}
    for c in range(3):
        names[c] = score[cluster == c].mean().idxmax()
    if len(set(names.values())) < 3:
        print("  [warn] two clusters map to the same subtype; the unsupervised "
              "split does not reproduce the marker split")
    cluster_call = pd.Series([names[c] for c in cluster], index=cpm.index)

    agree = marker_call == cluster_call
    print(f"  subtype assignment: {agree.sum()} of {len(agree)} neurons agree "
          f"between markers and clustering ({agree.mean() * 100:.0f}%)")
    print(pd.crosstab(marker_call, cluster_call).to_string())
    return marker_call, cluster_call, agree


def kmeans(x, k=3, seed=0, n_iter=100):
    """Plain Lloyd's algorithm; the project has no sklearn dependency."""
    rng = np.random.default_rng(seed)
    centres = x[rng.choice(len(x), k, replace=False)]
    for _ in range(n_iter):
        d = ((x[:, None, :] - centres[None]) ** 2).sum(-1)
        lab = d.argmin(1)
        new = np.array([x[lab == c].mean(0) if (lab == c).any() else centres[c]
                        for c in range(k)])
        if np.allclose(new, centres):
            break
        centres = new
    return lab


# ---------------------------------------------------------------------- figure
def figure(dev, sub_table, per_cell, kruskal_p):
    st.set_theme()
    fig = plt.figure(figsize=(16.5, 5.2))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.3, 1.0, 1.0], wspace=0.42,
                          left=0.06, right=0.985, top=0.80, bottom=0.20)

    ax = fig.add_subplot(gs[0])
    d = dev.sort_values("age", key=lambda s: s.map(AGE_ORDER.index))
    x = np.arange(len(d))
    for g, colour in zip(ac.RECEPTORS,
                         [st.BAR_BLUE, "#B2182B", "#4D4D4D", "#B8B8B8"]):
        ax.plot(x, d[g], marker="o", ms=6, lw=1.8, color=colour, label=g,
                zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels(d.age, fontsize=10, rotation=45, ha="right")
    ax.set_ylabel("Mean expression (CPM)", fontsize=11)
    ax.set_yscale("symlog", linthresh=0.1)
    ax.set_ylim(0, float(d[ac.RECEPTORS].to_numpy().max()) * 4)
    # The deposit changes twice along this axis, and absolute level is not
    # comparable across the changes. Only the ordering is.
    for i in range(1, len(d)):
        if d.dataset.iloc[i] != d.dataset.iloc[i - 1]:
            ax.axvline(i - 0.5, color="black", lw=1.0, ls="--")
    for lab, pos in [("droplet", 1.0), ("SMART-seq2", 5.5), ("SMART-seq", 9.0)]:
        ax.text(pos, float(d[ac.RECEPTORS].to_numpy().max()) * 2.2, lab,
                ha="center", fontsize=8.5, color="#444444")
    ax.legend(fontsize=9, ncol=4, loc="lower left")
    st.panel_letter(ax, "a", dx=-0.10)
    ax.set_title("The otic lineage from neurogenesis to the adult ganglion\n"
                 "level is not comparable across the dashed lines, ordering is",
                 fontsize=11, pad=10)

    ax = fig.add_subplot(gs[1])
    s = sub_table.set_index("subtype")
    order = [t for t in ["Ia", "Ib", "Ic"] if t in s.index]
    st.expression_bars(ax, s.loc[order, "Oprl1"].to_numpy(), order,
                       "Oprl1 (mean CPM)", italic=False, rotation=0,
                       fontsize=12)
    for i, t in enumerate(order):
        ax.text(i, 0, f"n = {int(s.loc[t, 'n'])}", ha="center", va="bottom",
                fontsize=9, color="white", fontweight="bold")
    st.panel_letter(ax, "b", dx=-0.26)
    ax.set_title(f"Spiral ganglion subtypes at P25\nGSE114997, "
                 f"Kruskal-Wallis p = {kruskal_p:.3f}", fontsize=11, pad=10)

    ax = fig.add_subplot(gs[2])
    colours = {"Ia": st.BAR_BLUE, "Ib": "#7F7F7F", "Ic": st.HIGHLIGHT}
    for t in order:
        m = per_cell.subtype == t
        ax.scatter(per_cell.loc[m, "Scn1a"], per_cell.loc[m, "Oprl1"], s=34,
                   c=colours[t], edgecolors="black", linewidths=0.5, label=t,
                   zorder=3)
    rho, p = stats.spearmanr(per_cell.Scn1a, per_cell.Oprl1)
    ax.set_xlabel("Scn1a (CPM)", fontsize=11)
    ax.set_ylabel("Oprl1 (CPM)", fontsize=11)
    ax.legend(fontsize=9, title="subtype", title_fontsize=9)
    st.panel_letter(ax, "c", dx=-0.26)
    ax.set_title(f"Within one ganglion, per neuron\nrho = {rho:.2f}, "
                 f"p = {p:.1e}", fontsize=11, pad=10)

    fig.suptitle("When the Oprl1 ordering appears, and which spiral neuron "
                 "carries it", fontsize=14, y=0.965)
    st.save(fig, "figure1b_otic_lineage")


# ------------------------------------------------------------------------ main
def main() -> int:
    ac.RES.mkdir(parents=True, exist_ok=True)
    rows = []

    targets = ac.RECEPTORS + OTIC_NEURON
    print("\n  Otic neurogenesis, E9.5 to E13.5 (GSE178931, droplet):")
    x, lib, ngene, age = load_gse178931(targets)
    cpm, counts, keep = neurons_and_cpm(x, "GSE178931", lib, ngene)
    rows = by_age(cpm, counts, age[keep], "GSE178931", rows)
    del x, cpm, counts

    print("\n  Cochlea, E14.5 to P3 (GSE165502, SMART-seq2):")
    x, age = load_gse165502()
    cpm, counts, keep = neurons_and_cpm(x, "GSE165502")
    rows = by_age(cpm, counts, age[keep], "GSE165502", rows)
    del x

    print("\n  Spiral ganglion at P25 (GSE114997, SMART-seq):")
    adult = load_gse114997()
    a_cpm, a_counts, a_keep = neurons_and_cpm(adult, "GSE114997")
    rows = by_age(a_cpm, a_counts, pd.Series("P25", index=a_cpm.index),
                  "GSE114997", rows)

    dev = pd.DataFrame(rows)
    ac.save_table(dev, "otic_lineage_by_age.csv")
    print("\n  Receptor ordering along the otic lineage:")
    print(dev[["dataset", "age", "n", "Oprl1", "Oprm1", "Oprd1", "Oprk1",
               "top_gene", "margin", "support"]].to_string(index=False))

    ac.save_table(ac.check_markers(a_cpm.mean(), "GSE114997 (P25 wild type)"),
                  "spiral_marker_checks.csv")
    marker_call, cluster_call, agree = spiral_subtypes(a_cpm, a_counts)

    per_cell = pd.DataFrame({
        "subtype": marker_call[agree],
        **{g: a_cpm.loc[agree, g] for g in
           ["Oprl1", "Scn1a", "Calb2", "Calb1", "Lypd1", "Pou4f1", "Pvalb"]
           if g in a_cpm.columns}})
    ac.save_table(per_cell.reset_index().rename(columns={"index": "cell"}),
                  "spiral_subtype_per_cell.csv")

    sub = per_cell.groupby("subtype").agg(
        n=("Oprl1", "size"),
        **{g: (g, "mean") for g in per_cell.columns if g != "subtype"}
    ).round(2).reset_index()
    ac.save_table(sub, "spiral_subtype_summary.csv")
    print("\n  Spiral ganglion subtypes, agreeing cells only:")
    print(sub.to_string(index=False))

    groups = [per_cell.loc[per_cell.subtype == t, "Oprl1"].to_numpy()
              for t in sub.subtype]
    h, p = stats.kruskal(*groups)
    rho, prho = stats.spearmanr(per_cell.Scn1a, per_cell.Oprl1)
    print(f"\n  Oprl1 across subtypes: Kruskal-Wallis H = {h:.2f}, p = {p:.4f}")
    print(f"  Oprl1 against Scn1a per neuron: rho = {rho:.3f}, p = {prho:.2e} "
          f"(n = {len(per_cell)})")
    ac.save_table(pd.DataFrame([{
        "n_neurons": len(per_cell), "kruskal_H": round(h, 3),
        "kruskal_p": p, "spearman_rho_Oprl1_Scn1a": round(rho, 3),
        "spearman_p": prho,
        "agreement_marker_vs_cluster": round(float(agree.mean()), 3)}]),
        "spiral_subtype_tests.csv")

    figure(dev, sub, per_cell, p)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
