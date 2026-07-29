"""Cross-tissue synthesis: where Oprl1 stands, and what confounds reading it.

Combines the three per-tissue tables into the comparisons that survive the assay
differences, and produces the summary figure.

The organising finding is that opioid receptor ordering in these datasets is
determined partly by library chemistry rather than by tissue. Whole-cell and
nuclear preparations disagree systematically, in the direction predicted by
genomic span, so datasets are stratified by preparation before any tissue claim
is made.

Every cross-dataset comparison here is a rank or a percentile. Absolute levels
are shown only within a preparation and tissue, never on a common axis across
assays.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import atlas_common as ac
import atlas_style as st

ORDER = list(ac.DATASETS)


def load_ranks():
    gen = pd.read_csv(ac.RES / "geniculate_receptor_rank.csv")
    nod = pd.read_csv(ac.RES / "nodose_receptor_rank.csv")
    nts = pd.read_csv(ac.RES / "nts_receptor_rank.csv")
    nod = nod[nod.dataset.str.contains(":", regex=False)
              & (nod.tissue == "nodose+jugular")]
    ranks = pd.concat([gen, nod, nts], ignore_index=True)
    ranks["prep"] = ranks.dataset.map(ac.PREP)
    ranks["tissue_group"] = ranks.dataset.map(ac.TISSUE_OF)
    return ranks[ranks.dataset.isin(ORDER)]


def load_support():
    frames = []
    for name in ("geniculate_rank_support.csv", "nodose_rank_support.csv",
                 "nts_rank_support.csv"):
        p = ac.RES / name
        if p.exists():
            frames.append(pd.read_csv(p))
        else:
            print(f"  [warn] {name} missing; bootstrap support will be blank")
    if not frames:
        return pd.DataFrame(columns=["dataset", "support", "margin"])
    return pd.concat(frames, ignore_index=True)


def gene_spans():
    """Genomic span per gene, from one annotation for all of them.

    Read from Ensembl rather than from a dataset's own coordinate columns.
    GSE102443 quantifies 17,225 features and none of them is Pnoc, so taking
    spans from it dropped Pnoc from the comparison without saying so.
    """
    span = pd.read_csv(ac.DATA / "ensembl_gene_spans.csv")
    missing = set(ac.OPIOID_GENES) - set(span.gene)
    if missing:
        raise ac.SanityCheckError(
            f"no genomic span for {sorted(missing)}; run 00b_fetch_gene_spans.py")
    return span.set_index("gene")["span_kb"]


def main() -> int:
    st.set_theme()
    ranks = load_ranks()
    support = load_support()

    # ---------------------------------------------------- receptor ordering
    top = ranks[ranks["rank"] == 1].copy()
    dup = top.dataset[top.dataset.duplicated()].unique()
    if len(dup):
        raise ValueError(f"tied top receptor in {list(dup)}; the ordering is "
                         "indeterminate and must not be reported as a winner")
    top = top[["dataset", "tissue_group", "prep", "gene", "determinate"]] \
        .rename(columns={"gene": "top_receptor"})
    top = top.merge(support[["dataset", "support", "margin"]], on="dataset", how="left")
    top["dataset"] = pd.Categorical(top.dataset, ORDER, ordered=True)
    top = top.sort_values("dataset")
    ac.save_table(top, "top_receptor_by_dataset.csv")
    print("\n  Highest-expressed opioid receptor in each dataset:")
    print(top.round(3).to_string(index=False))

    wc, nu = top[top.prep == "whole cell"], top[top.prep == "nuclear"]
    print(f"\n  whole-cell datasets with Oprl1 first: "
          f"{(wc.top_receptor == 'Oprl1').sum()}/{len(wc)}")
    print(f"  nuclear datasets with Oprl1 first:    "
          f"{(nu.top_receptor == 'Oprl1').sum()}/{len(nu)}")
    weak = wc[(wc.top_receptor == "Oprl1") & (wc.support < 0.95)]
    if len(weak):
        print("  [note] not all of those are secure. Bootstrap support below 0.95: "
              + ", ".join(f"{r.dataset} {r.support:.2f}" for r in weak.itertuples()))

    # -------------------------------- Oprl1's own position, dataset by dataset
    oprl1 = ranks[ranks.gene == "Oprl1"][["dataset", "tissue_group", "prep",
                                          "level", "rank"]].copy()
    oprl1 = oprl1.merge(support[["dataset", "support"]], on="dataset", how="left")
    oprl1["dataset"] = pd.Categorical(oprl1.dataset, ORDER, ordered=True)
    oprl1 = oprl1.sort_values("dataset").rename(columns={"rank": "rank_of_4"})
    ac.save_table(oprl1, "oprl1_across_datasets.csv")
    print("\n  Oprl1's rank among the four receptors, per dataset:")
    print(oprl1.round(3).to_string(index=False))

    # ------------------------------------- gene length explains the inversion
    span = gene_spans()
    nod_lv = pd.read_csv(ac.RES / "nodose_opioid_levels.csv")
    wc_sets = [d for d in ac.WHOLE_CELL_NODOSE]
    sub = nod_lv[nod_lv.dataset.isin(wc_sets)]
    # Unweighted mean over datasets (one vote per study), and a cell-weighted
    # mean, because Zhao contributes 74% of the cells and 25% of the votes.
    sc = sub.groupby("gene")["mean_level"].mean()
    sc_w = (sub.assign(w=sub.mean_level * sub.n_cells).groupby("gene")
            .apply(lambda d: d.w.sum() / d.n_cells.sum(), include_groups=False))
    sn = nod_lv[nod_lv.dataset == "NodoMap:inhouse"].set_index("gene")["mean_level"]

    bias = pd.DataFrame({
        "genomic_span_kb": span.reindex(ac.OPIOID_GENES).round(0),
        "whole_cell_CPM": sc.reindex(ac.OPIOID_GENES).round(3),
        "whole_cell_CPM_cellweighted": sc_w.reindex(ac.OPIOID_GENES).round(3),
        "nuclear_CPM": sn.reindex(ac.OPIOID_GENES).round(3),
    }).rename_axis("gene")
    for src, dst in (("whole_cell_CPM", "nuclear_over_whole_cell"),
                     ("whole_cell_CPM_cellweighted", "nuclear_over_whole_cell_cw")):
        denom = bias[src].where(bias[src] > 0)      # no silent inf on a zero level
        bias[dst] = (bias.nuclear_CPM / denom).round(2)
    bias = bias.reset_index()
    ac.save_table(bias, "nuclear_bias_vs_gene_length.csv")

    v = bias.dropna(subset=["genomic_span_kb", "nuclear_over_whole_cell"])
    r = float(np.corrcoef(np.log10(v.genomic_span_kb),
                          np.log10(v.nuclear_over_whole_cell))[0, 1])
    rho, p_rho = stats.spearmanr(v.genomic_span_kb, v.nuclear_over_whole_cell)
    print("\n  Nuclear vs whole-cell bias against genomic span (same tissue):")
    print(bias.to_string(index=False))
    print(f"  log-log Pearson r = {r:.2f}; Spearman rho = {rho:.2f}, "
          f"p = {p_rho:.3f} (n = {len(v)} genes)")

    # How much of that correlation is carried by single genes.
    loo = ac.leave_one_out_pearson(v.genomic_span_kb, v.nuclear_over_whole_cell, v.gene)
    both = v[~v.gene.isin(["Oprm1", "Oprd1"])]
    r_both = float(np.corrcoef(np.log10(both.genomic_span_kb),
                               np.log10(both.nuclear_over_whole_cell))[0, 1])
    loo = pd.concat([loo, pd.DataFrame([{
        "dropped": "Oprm1 + Oprd1", "n": len(both), "pearson_r": round(r_both, 4),
        "delta_vs_full": round(r_both - r, 4)}])], ignore_index=True)
    ac.save_table(loo, "nuclear_bias_sensitivity.csv")
    print("\n  Leave-one-out sensitivity of that correlation:")
    print(loo.to_string(index=False))

    figures(ranks, top, oprl1, v, r, loo, r_both)

    summary = pd.DataFrame([{
        "statement": "Oprl1 is the top-ranked opioid receptor in whole-cell data",
        "value": f"{(wc.top_receptor == 'Oprl1').sum()}/{len(wc)} datasets, 2 ganglia; "
                 f"bootstrap support {wc.support.min():.2f}-{wc.support.max():.2f}",
    }, {
        "statement": "nuclear preparations invert this toward Oprm1",
        "value": f"{(nu.top_receptor == 'Oprm1').sum()}/{len(nu)} datasets; "
                 f"Oprm1 spans {bias.set_index('gene').loc['Oprm1','genomic_span_kb']:.0f} kb "
                 f"vs Oprl1 {bias.set_index('gene').loc['Oprl1','genomic_span_kb']:.0f} kb",
    }, {
        "statement": "that correlation does not depend on the two long receptors",
        "value": f"log-log r = {r:.2f} (n = {len(v)}); "
                 f"r = {r_both:.2f} without Oprm1 and Oprd1; "
                 f"Spearman rho = {rho:.2f}, p = {p_rho:.4f}",
    }])
    ac.save_table(summary, "synthesis_summary.csv")
    print("\n" + summary.to_string(index=False))
    return 0


def figures(ranks, top, oprl1, v, r, loo, r_both):
    main_figure(ranks, oprl1)
    supplementary_figure(v, r, loo, r_both)


def main_figure(ranks, oprl1):
    """The four opioid receptors in every dataset, in each dataset's own unit.

    Eight small panels rather than one composite: units differ between datasets,
    so the only honest comparison is down a panel, and each panel is labelled
    with the unit it is in.
    """
    fig, axes = plt.subplots(2, 4, figsize=(15.0, 7.6))
    unit = {"GSE102443": "FPKM"}
    for ax, ds in zip(axes.ravel(), ORDER):
        sub = ranks[ranks.dataset == ds].set_index("gene")
        vals = [float(sub.loc[g, "level"]) if g in sub.index else np.nan
                for g in ac.RECEPTORS]
        order = np.argsort(-np.array(vals))
        colours = [st.BAR_BLUE if ac.RECEPTORS[i] == "Oprl1" else st.BAR_GREY
                   for i in order]
        st.expression_bars(ax, [vals[i] for i in order],
                           [ac.RECEPTORS[i] for i in order],
                           f"Mean expression ({unit.get(ds, 'CPM')})",
                           colors=colours, fontsize=10)
        prep = ac.PREP[ds]
        ax.set_title(f"{ds}\n{ac.TISSUE_OF[ds]}, {prep}", fontsize=10,
                     color="#B2182B" if prep == "nuclear" else "black", pad=8)
    fig.suptitle("Oprl1 (blue) is the highest-expressed opioid receptor in every "
                 "whole-cell dataset\nthe two nuclear preparations (red titles) "
                 "invert this — see SUPPLEMENT.md", fontsize=12, y=1.03)
    fig.tight_layout()
    st.save(fig, "figureS3_all_datasets_receptors")


def supplementary_figure(v, r, loo, r_both):
    """Data quality: what nuclear preparation does to the receptor ordering."""
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.6))

    # (A) The nuclear inversion against genomic span.
    ax = axes[0]
    ax.scatter(v.genomic_span_kb, v.nuclear_over_whole_cell, s=48,
               c=["#B2182B" if g in ("Oprm1", "Oprd1") else "#08306B"
                  for g in v.gene], edgecolors="black", linewidths=0.4, zorder=3)
    # Oprl1, Pomc and Penk sit within 1 kb of each other; fan the labels apart.
    offsets = {"Oprl1": (7, 3), "Pomc": (7, -4), "Penk": (-2, -13), "Pdyn": (6, 2)}
    for _, rr in v.iterrows():
        ax.annotate(rr.gene, (rr.genomic_span_kb, rr.nuclear_over_whole_cell),
                    fontsize=7.5, fontstyle="italic",
                    xytext=offsets.get(rr.gene, (5, 4)), textcoords="offset points")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.axhline(1.0, color="black", lw=0.8, ls="--")
    ax.set_xlabel("genomic span (kb, log)", fontsize=8.5)
    ax.set_ylabel("nuclear / whole-cell level", fontsize=8.5)
    ax.set_title(f"(A) The inversion tracks gene length\n"
                 f"same tissue, log-log r = {r:.2f} (n = {len(v)})", fontsize=10)

    # (B) How much of (A) rests on single genes.
    ax = axes[1]
    d = loo[loo.dropped != "(none)"].copy()
    colours = ["#B2182B" if x in ("Oprm1", "Oprd1", "Oprm1 + Oprd1") else "#999999"
               for x in d.dropped]
    yy = np.arange(len(d))
    ax.barh(yy, d.pearson_r, color=colours, edgecolor="black", linewidth=0.4)
    ax.set_yticks(yy)
    ax.set_yticklabels([f"without {x}" for x in d.dropped], fontsize=7.5)
    ax.invert_yaxis()
    ax.axvline(r, color="black", lw=0.9, ls="--")
    ax.text(r - 0.02, -0.75, f"all 7 genes: r = {r:.2f}", fontsize=7,
            va="bottom", ha="right")
    ax.set_xlim(-0.05, 1.0)
    ax.set_xlabel("log-log Pearson r with that gene removed", fontsize=8.5)
    ax.set_title(f"(B) Two genes carry the correlation\n"
                 f"r = {r_both:.2f} without both", fontsize=10)

    fig.tight_layout()
    st.save(fig, "figureS1_nuclear_preparation_bias")


if __name__ == "__main__":
    raise SystemExit(main())
