"""Cross-tissue synthesis: what the Oprl1 signature is, and what it is not.

Combines the three per-tissue tables into the comparisons that survive the
assay differences, and produces the figures.

The organising finding is that opioid receptor ordering in these datasets is
determined more by library chemistry than by tissue. Whole-cell and nuclear
preparations disagree systematically, in the direction predicted by genomic
span, so datasets are stratified by preparation before any tissue claim is made.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import atlas_common as ac

# Which datasets are whole-cell and which are nuclear. This is the single most
# important covariate in the comparison.
PREP = {
    "GSE102443": "whole cell", "GSE135801": "whole cell",
    "NodoMap:Bai": "whole cell", "NodoMap:Buchanan": "whole cell",
    "NodoMap:Kupari": "whole cell", "NodoMap:Zhao": "whole cell",
    "NodoMap:inhouse": "nuclear", "GSE166648": "nuclear",
}
TISSUE_OF = {
    "GSE102443": "geniculate", "GSE135801": "geniculate",
    "NodoMap:Bai": "nodose/jugular", "NodoMap:Buchanan": "nodose/jugular",
    "NodoMap:Kupari": "nodose/jugular", "NodoMap:Zhao": "nodose/jugular",
    "NodoMap:inhouse": "nodose/jugular", "GSE166648": "NTS (central)",
}
ORDER = list(PREP)


def load_all():
    gen = pd.read_csv(ac.RES / "geniculate_receptor_rank.csv")
    nod = pd.read_csv(ac.RES / "nodose_receptor_rank.csv")
    nts = pd.read_csv(ac.RES / "nts_receptor_rank.csv")
    nod = nod[nod.dataset.str.contains(":") & (nod.tissue == "nodose+jugular")]
    ranks = pd.concat([gen, nod, nts], ignore_index=True)
    ranks["prep"] = ranks.dataset.map(PREP)
    ranks["tissue_group"] = ranks.dataset.map(TISSUE_OF)
    return ranks


def main() -> int:
    ac.set_theme()
    ranks = load_all()
    ranks = ranks[ranks.dataset.isin(ORDER)]

    # ---------------------------------------------------- receptor ordering
    top = (ranks[ranks["rank"] == 1][["dataset", "tissue_group", "prep", "gene"]]
           .rename(columns={"gene": "top_receptor"}))
    top["dataset"] = pd.Categorical(top.dataset, ORDER, ordered=True)
    top = top.sort_values("dataset")
    ac.save_table(top, "top_receptor_by_dataset.csv")
    print("\n  Highest-expressed opioid receptor in each dataset:")
    print(top.to_string(index=False))

    wc = top[top.prep == "whole cell"]
    nu = top[top.prep == "nuclear"]
    print(f"\n  whole-cell datasets with Oprl1 first: "
          f"{(wc.top_receptor == 'Oprl1').sum()}/{len(wc)}")
    print(f"  nuclear datasets with Oprl1 first:    "
          f"{(nu.top_receptor == 'Oprl1').sum()}/{len(nu)}")

    # ------------------------------------- gene length explains the inversion
    # Genomic coordinates come from the GSE102443 transcript table, the only
    # source here that carries them.
    span = pd.read_csv(
        ac.DATA / "GSE102443_GEO-ID_Dvoryanchikov_2017_Datatable_FPKM.txt.gz",
        sep="\t", usecols=["Gene", "Begin", "End"], low_memory=False)
    span = (span.assign(kb=(span.End - span.Begin).abs() / 1000)
            .groupby("Gene")["kb"].max())
    nod_lv = pd.read_csv(ac.RES / "nodose_opioid_levels.csv")
    sc = (nod_lv[nod_lv.dataset.isin([f"NodoMap:{d}" for d in
                                      ["Bai", "Buchanan", "Kupari", "Zhao"]])]
          .groupby("gene")["mean_level"].mean())
    sn = nod_lv[nod_lv.dataset == "NodoMap:inhouse"].set_index("gene")["mean_level"]
    bias = pd.DataFrame({
        "genomic_span_kb": span.reindex(ac.OPIOID_GENES).round(0),
        "whole_cell_CPM": sc.reindex(ac.OPIOID_GENES).round(3),
        "nuclear_CPM": sn.reindex(ac.OPIOID_GENES).round(3),
    })
    bias["nuclear_over_whole_cell"] = (bias.nuclear_CPM / bias.whole_cell_CPM).round(2)
    bias = bias.reset_index().rename(columns={"index": "gene"})
    ac.save_table(bias, "nuclear_bias_vs_gene_length.csv")
    v = bias.dropna(subset=["genomic_span_kb", "nuclear_over_whole_cell"])
    r = float(np.corrcoef(np.log10(v.genomic_span_kb),
                          np.log10(v.nuclear_over_whole_cell))[0, 1])
    rho, p_rho = stats.spearmanr(v.genomic_span_kb, v.nuclear_over_whole_cell)
    print("\n  Nuclear vs whole-cell bias against genomic span (same tissue):")
    print(bias.to_string(index=False))
    print(f"  log-log Pearson r = {r:.2f}; Spearman rho = {rho:.2f}, "
          f"p = {p_rho:.3f} (n = {len(v)} genes)")

    # ------------------------------------------- receptor without local ligand
    lr = pd.concat([
        pd.read_csv(ac.RES / "geniculate_ligand_receptor.csv"),
        pd.read_csv(ac.RES / "nodose_ligand_receptor.csv"),
        pd.read_csv(ac.RES / "nts_ligand_receptor.csv"),
    ], ignore_index=True)
    lr["prep"] = lr.dataset.map(PREP)
    keep = lr.dataset.isin(ORDER) | lr.tissue.isin(["nodose", "jugular"])
    lr = lr[keep].copy()
    lr.loc[lr.dataset == "NodoMap", "prep"] = "mixed (pooled)"
    ac.save_table(lr, "ligand_receptor_by_tissue.csv")
    print("\n  Oprl1 against Pnoc, within the same cells:")
    print(lr[["tissue", "dataset", "prep", "Oprl1", "Pnoc",
              "Oprl1_over_Pnoc", "Pnoc_in_matrix"]].round(3).to_string(index=False))

    # --------------------------------------------------------------- figures
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.4))

    ax = axes[0]
    y = np.arange(len(ORDER))
    colour = {"Oprl1": "#08306B", "Oprm1": "#B2182B",
              "Oprd1": "#999999", "Oprk1": "#7F7F7F"}
    tt = top.set_index("dataset")
    for i, ds in enumerate(ORDER):
        g = tt.loc[ds, "top_receptor"]
        frac = ranks[(ranks.dataset == ds) & (ranks["rank"] == 2)]["fraction_of_top"]
        ax.barh(i, 1.0, color=colour.get(g, "#CCCCCC"), edgecolor="black", linewidth=0.4)
        if len(frac):
            ax.barh(i, float(frac.iloc[0]), color="white", alpha=0.0,
                    edgecolor="black", linewidth=0.4, hatch="///")
        ax.text(1.03, i, g, va="center", fontsize=8, fontstyle="italic")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{d}  ({PREP[d]})" for d in ORDER], fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.35)
    ax.set_xlabel("top receptor (bar) and runner-up as a fraction of it (hatched)",
                  fontsize=8)
    ax.set_title("(A) Highest-expressed opioid receptor\nby dataset and preparation",
                 fontsize=9.5)

    ax = axes[1]
    v2 = v[v.gene != "Pnoc"]
    ax.scatter(v2.genomic_span_kb, v2.nuclear_over_whole_cell, s=45,
               c=["#B2182B" if g in ("Oprm1", "Oprd1") else "#08306B"
                  for g in v2.gene], edgecolors="black", linewidths=0.4, zorder=3)
    # Penk, Pomc and Oprl1 sit almost on top of each other; fan the labels out.
    offsets = {"Penk": (5, -9), "Pomc": (5, 3), "Oprl1": (5, 6)}
    for _, rr in v2.iterrows():
        ax.annotate(rr.gene, (rr.genomic_span_kb, rr.nuclear_over_whole_cell),
                    fontsize=7.5, fontstyle="italic",
                    xytext=offsets.get(rr.gene, (4, 4)), textcoords="offset points")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.axhline(1.0, color="black", lw=0.8, ls="--")
    ax.set_xlabel("genomic span (kb, log)", fontsize=8.5)
    ax.set_ylabel("nuclear / whole-cell level", fontsize=8.5)
    ax.set_title(f"(B) The inversion is a gene-length effect\n"
                 f"same tissue, log-log r = {r:.2f} (n = {len(v2)})", fontsize=9.5)

    ax = axes[2]
    sel = lr[lr.dataset.isin(["GSE135801", "NodoMap:Bai", "NodoMap:Buchanan",
                              "NodoMap:Kupari", "NodoMap:Zhao",
                              "NodoMap:inhouse", "GSE166648"])
             | lr.tissue.isin(["jugular"])].copy()
    sel = sel[np.isfinite(sel.Oprl1_over_Pnoc)]
    lbl = sel.apply(lambda r_: f"{r_.tissue} / {r_.dataset}", axis=1)
    cols = ["#1B7837" if t == "geniculate" else "#6BAED6" if t == "jugular"
            else "#B2182B" if t == "NTS" else "#08306B" for t in sel.tissue]
    yy = np.arange(len(sel))
    ax.barh(yy, sel.Oprl1_over_Pnoc, color=cols, edgecolor="black", linewidth=0.4)
    ax.set_yticks(yy); ax.set_yticklabels(lbl, fontsize=7.5)
    ax.invert_yaxis()
    ax.axvline(1.0, color="black", lw=0.8, ls="--")
    ax.set_xscale("log")
    ax.set_xlabel("Oprl1 : Pnoc, same cells (log)", fontsize=8.5)
    ax.set_title("(C) Receptor exceeds its own ligand\nin every peripheral dataset",
                 fontsize=9.5)

    fig.tight_layout()
    ac.save_fig(fig, "figure1_cross_tissue_signature")

    # Summary object for the README.
    summary = pd.DataFrame([{
        "statement": "Oprl1 is the top-ranked opioid receptor in whole-cell data",
        "value": f"{(wc.top_receptor == 'Oprl1').sum()}/{len(wc)} datasets, "
                 f"2 ganglia",
    }, {
        "statement": "nuclear preparations invert this toward Oprm1",
        "value": f"{(nu.top_receptor == 'Oprm1').sum()}/{len(nu)} datasets; "
                 f"Oprm1 spans {bias.set_index('gene').loc['Oprm1','genomic_span_kb']:.0f} kb "
                 f"vs Oprl1 {bias.set_index('gene').loc['Oprl1','genomic_span_kb']:.0f} kb",
    }, {
        "statement": "Oprl1 exceeds Pnoc in every peripheral dataset",
        "value": f"ratios {sel[sel.tissue!='NTS'].Oprl1_over_Pnoc.min():.1f} to "
                 f"{sel[sel.tissue!='NTS'].Oprl1_over_Pnoc.max():.1f}",
    }])
    ac.save_table(summary, "synthesis_summary.csv")
    print("\n" + summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
