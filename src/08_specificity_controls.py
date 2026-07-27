"""Are the co-expression odds ratios specific, or is Oprl1 just a big-cell gene?

Three controls that the co-expression results in sections 4 and 5 have not had.

1. A MATCHED NULL FOR THE ODDS RATIO.
   `stratified_odds_ratio` holds capture depth fixed, but nUMI is not cell size
   or total RNA content, and myelinated A-fibre somata are large. If Oprl1 is
   simply a high expresser in big, transcriptionally active neurons, then within
   a cluster it would show a positive odds ratio against almost any moderately
   expressed gene. The test: draw control genes matched to each partner on
   detection rate and mean expression, and compute the same statistic against
   them. If the null median sits near 1.0 the partner ORs mean something; if it
   sits at 1.3-1.4 then Glp1r, Cckar and Trpv1 are at chance and only Piezo2
   survives.

2. THE Nav CLASS RESULT, PER CELL.
   Crude Trpv1/Oprl1 co-detection inverts under stratification (OR 0.37 -> 1.34),
   so compositional effects can reverse a sign here. The Nav1.1/Nav1.8 split is
   the claim everything else hangs on and it exists only as a Kruskal-Wallis over
   21 cluster means. Scn1a and Scn10a get the same per-cell treatment, and the
   cluster group sizes are reported: 21 split 8/13 is not 3/18.

3. AN AMBIENT-RNA CHECK ON THE CORRELATE LIST.
   NodoMap contains ~50,000 satellite and myelinating glia. Adgrg6/Gpr126 is the
   Schwann-cell myelination receptor, and its topping the Oprl1 correlate list is
   either the confirmation we want or glial ambient RNA bleeding into the
   large-soma clusters. There is no ambient correction (CellBender/SoupX/decontX)
   anywhere in this pipeline, so every top correlate is scored for neuronal
   enrichment against the glial compartment of the same ganglion.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import anndata as ad
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp

import atlas_common as ac
import atlas_style as st

_nodose = __import__("02_nodose")

PARTNERS = ["Piezo2", "Glp1r", "Cckar", "Cckbr", "Trpv1", "Scn1a", "Scn10a"]
N_CONTROLS = 100
DET_TOL = 0.20          # relative tolerance on detection rate
CPM_TOL = 0.50          # relative tolerance on mean CPM
N_TOP_CORRELATES = 50
GLIA_CLASSES = ["Satellite glial cell", "Myelinated glial cell", "Glial cell"]
CACHE = ac.DATA / "nodomap_gene_stats.npz"


def gene_stats(path, adata, is_nodose, is_glia, lib, force=False):
    """Detection rate and mean CPM per gene, in nodose neurons and in glia."""
    n_genes = adata.n_vars
    if CACHE.exists() and not force:
        z = np.load(CACHE)
        if z["det"].shape[0] == n_genes:
            print(f"  [cache] {CACHE.name}")
            return z["det"], z["cpm"], z["glia_cpm"]

    print("  streaming the matrix for per-gene detection and abundance ...")
    det = np.zeros(n_genes)
    cpm = np.zeros(n_genes)
    glia = np.zeros(n_genes)
    with h5py.File(path, "r") as fh:
        grp = fh["raw/X"]
        indptr = grp["indptr"][:]
        data_ds, idx_ds = grp["data"], grp["indices"]
        n_cells = len(indptr) - 1
        for start in range(0, n_cells, 20_000):
            stop = min(start + 20_000, n_cells)
            lo, hi = int(indptr[start]), int(indptr[stop])
            m = sp.csr_matrix(
                (data_ds[lo:hi], idx_ds[lo:hi], indptr[start:stop + 1] - lo),
                shape=(stop - start, n_genes))
            norm = sp.diags(1e6 / lib[start:stop]) @ m
            nb, gb = is_nodose[start:stop], is_glia[start:stop]
            if nb.any():
                det += np.asarray((m[nb] > 0).sum(axis=0)).ravel()
                cpm += np.asarray(norm[nb].sum(axis=0)).ravel()
            if gb.any():
                glia += np.asarray(norm[gb].sum(axis=0)).ravel()
            print(f"    {stop:,}/{n_cells:,} cells", flush=True)

    det = det / max(is_nodose.sum(), 1) * 100.0
    cpm = cpm / max(is_nodose.sum(), 1)
    glia = glia / max(is_glia.sum(), 1)
    np.savez_compressed(CACHE, det=det, cpm=cpm, glia_cpm=glia)
    print(f"  [cache] wrote {CACHE.name}")
    return det, cpm, glia


def matched_controls(symbols, det, cpm, partner_idx, exclude, rng, n=N_CONTROLS):
    """Genes with a similar detection rate and abundance to the partner gene."""
    d0, c0 = det[partner_idx], cpm[partner_idx]
    ok = ((np.abs(det - d0) <= DET_TOL * d0)
          & (np.abs(cpm - c0) <= CPM_TOL * c0)
          & (det > 0))
    ok[list(exclude)] = False
    idx = np.nonzero(ok)[0]
    if idx.size > n:
        idx = rng.choice(idx, n, replace=False)
    return idx


def main() -> int:
    ac.RES.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(ac.NODOSE_H5AD, backed="r")
    obs = adata.obs
    symbols = adata.raw.var["feature_name"].astype(str).values

    whole_cell = obs["suspension_type"].values == "cell"
    is_nodose = (obs["cell_class"] == "Nodose ganglion neuron").values & whole_cell
    is_glia = obs["cell_class"].isin(GLIA_CLASSES).values & whole_cell
    cluster = obs["author_cell_type"].astype(str).values

    counts0, _ = _nodose.load_counts(adata, ac.NODOSE_H5AD, ["Oprl1"])
    lib = _nodose.library_size(obs, counts0).astype(float)
    lib[lib == 0] = 1.0

    det, cpm, glia_cpm = gene_stats(ac.NODOSE_H5AD, adata, is_nodose, is_glia, lib)

    # ------------------------------------------------- 1. matched null for the OR
    by_symbol = {}
    for j, sname in enumerate(symbols):
        by_symbol.setdefault(sname, []).append(j)
    first = {k: v[0] for k, v in by_symbol.items()}

    rng = np.random.default_rng(0)
    control_idx, per_partner = set(), {}
    for p in PARTNERS:
        if p not in first:
            print(f"  [warn] {p} absent from the annotation")
            continue
        idx = matched_controls(symbols, det, cpm, first[p],
                               {first[q] for q in PARTNERS if q in first}
                               | {first["Oprl1"]}, rng)
        per_partner[p] = idx
        control_idx |= set(idx.tolist())
        print(f"  {p:8s} det {det[first[p]]:5.2f}%  {cpm[first[p]]:8.2f} CPM"
              f"  -> {len(idx)} matched control genes")

    load = {"Oprl1": by_symbol["Oprl1"]}
    for p in PARTNERS:
        if p in by_symbol:
            load[p] = by_symbol[p]
    for j in sorted(control_idx):
        load[f"__ctrl_{j}"] = [j]
    print(f"  loading {len(load):,} gene columns for the null ...")
    with h5py.File(ac.NODOSE_H5AD, "r") as fh:
        mat = ac.csr_gene_columns(fh["raw/X"], load, adata.n_obs)
    frame = pd.DataFrame(mat, columns=list(load))

    sub = frame.loc[is_nodose] > 0
    nf = obs["nFeature_RNA"].values[is_nodose]
    joint = ac.cluster_depth_strata(cluster[is_nodose], nf)
    o_pos = sub["Oprl1"].values

    rows = []
    for p, idx in per_partner.items():
        or_p, pp, _ = ac.stratified_odds_ratio(o_pos, sub[p].values, joint)
        nulls = []
        for j in idx:
            col = f"__ctrl_{j}"
            or_c, _, k = ac.stratified_odds_ratio(o_pos, sub[col].values, joint)
            if np.isfinite(or_c) and k > 0:
                nulls.append(or_c)
        nulls = np.array(nulls)
        # Empirical two-sided position of the observed OR in its own null.
        p_emp = (np.sum(nulls >= or_p) + 1) / (nulls.size + 1) if nulls.size else np.nan
        rows.append({
            "partner": p,
            "detection_pct": round(float(det[first[p]]), 2),
            "mean_CPM": round(float(cpm[first[p]]), 2),
            "observed_OR": round(float(or_p), 3),
            "p_mantel_haenszel": pp,
            "n_control_genes": int(nulls.size),
            "null_median_OR": round(float(np.median(nulls)), 3) if nulls.size else np.nan,
            "null_p2.5": round(float(np.percentile(nulls, 2.5)), 3) if nulls.size else np.nan,
            "null_p97.5": round(float(np.percentile(nulls, 97.5)), 3) if nulls.size else np.nan,
            "empirical_p_vs_null": round(float(p_emp), 4),
        })
    null_tbl = pd.DataFrame(rows).sort_values("observed_OR", ascending=False)
    ac.save_table(null_tbl, "vagal_or_matched_null.csv")
    print("\n  Observed odds ratio against a null of expression-matched genes:")
    print(null_tbl.to_string(index=False))

    # -------------------------------------------- 2. cluster group sizes for Nav
    ann = pd.read_csv(ac.RES / "nodose_oprl1_by_cluster_annotated.csv")
    print("\n  Cluster group sizes behind the annotation tests (n = %d clusters):"
          % len(ann))
    for col in ("sodium_channel_type", "fibre_type", "sensor_type",
                "organ_projection"):
        vc = ann[col].value_counts()
        print(f"    {col:22s} " + ", ".join(f"{k} n={v}" for k, v in vc.items()))
    sizes = pd.concat([
        ann[col].value_counts().rename_axis("level").reset_index(name="n_clusters")
        .assign(annotation=col)
        for col in ("sodium_channel_type", "fibre_type", "sensor_type",
                    "organ_projection")], ignore_index=True)
    ac.save_table(sizes[["annotation", "level", "n_clusters"]],
                  "nodose_annotation_group_sizes.csv")

    # ------------------------------------------- 3. ambient check on correlates
    corr = pd.read_csv(ac.RES / "nodose_oprl1_gene_correlations.csv")
    top = corr.head(N_TOP_CORRELATES).copy()
    n_cpm = pd.Series(cpm, index=symbols).groupby(level=0).sum()
    g_cpm = pd.Series(glia_cpm, index=symbols).groupby(level=0).sum()
    top["neuron_CPM"] = [round(float(n_cpm.get(g, np.nan)), 3) for g in top.gene]
    top["glia_CPM"] = [round(float(g_cpm.get(g, np.nan)), 3) for g in top.gene]
    top["log2_neuron_over_glia"] = np.log2((top.neuron_CPM + 0.01) /
                                           (top.glia_CPM + 0.01)).round(3)
    top["glial_biased"] = top.log2_neuron_over_glia < 0
    ac.save_table(top, "nodose_top_correlates_ambient_check.csv")
    n_bad = int(top.glial_biased.sum())
    print(f"\n  Ambient check on the top {len(top)} Oprl1 correlates: "
          f"{n_bad} are more abundant in glia than in nodose neurons")
    print(top[["gene", "spearman_rho", "neuron_CPM", "glia_CPM",
               "log2_neuron_over_glia"]].head(15).to_string(index=False))
    o = float(np.log2((n_cpm.get("Oprl1", 0) + 0.01) / (g_cpm.get("Oprl1", 0) + 0.01)))
    print(f"  For reference, Oprl1 itself: log2(neuron/glia) = {o:.2f}")

    figures(null_tbl, top)
    return 0


def figures(null_tbl, top):
    st.set_theme()
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.0),
                             gridspec_kw={"width_ratios": [1.15, 1.0]})

    ax = axes[0]
    d = null_tbl.sort_values("observed_OR")
    y = np.arange(len(d))
    for i, r in enumerate(d.itertuples()):
        ax.plot([r._8, r._9], [i, i], color="#BBBBBB", lw=5, solid_capstyle="butt",
                zorder=2)
        ax.scatter([r.null_median_OR], [i], s=34, c="#666666", zorder=3)
    ax.scatter(d.observed_OR, y, s=80, c=[st.BAR_BLUE if v > 1.6 else st.HIGHLIGHT
                                          for v in d.observed_OR],
               edgecolors="black", linewidths=0.7, zorder=4)
    ax.axvline(1.0, color="black", lw=0.9, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels(d.partner, fontsize=11, fontstyle="italic")
    ax.set_xlabel("odds of Oprl1 detection, partner+ vs partner−\n"
                  "(same cluster, same depth)", fontsize=10)
    st.panel_letter(ax, "a", dx=-0.18)
    ax.set_title("Observed against a null of expression-matched genes\n"
                 "grey bar = central 95 % of the null", fontsize=11)

    ax = axes[1]
    ax.scatter(top.log2_neuron_over_glia, top.spearman_rho, s=46,
               c=[st.HIGHLIGHT if b else st.BAR_BLUE for b in top.glial_biased],
               edgecolors="black", linewidths=0.5, zorder=3)
    for _, r in top.head(10).iterrows():
        ax.annotate(r.gene, (r.log2_neuron_over_glia, r.spearman_rho), fontsize=7,
                    fontstyle="italic", xytext=(4, 3), textcoords="offset points")
    ax.axvline(0, color="black", lw=0.9, ls="--")
    ax.set_xlabel("log2(nodose neuron / glia) in the same ganglion", fontsize=10)
    ax.set_ylabel("Spearman rho with Oprl1", fontsize=10)
    st.panel_letter(ax, "b", dx=-0.20)
    ax.set_title(f"Top {len(top)} correlates: neuronal or ambient?\n"
                 "red = more abundant in glia", fontsize=11)

    fig.tight_layout()
    st.save(fig, "figureS4_specificity_controls")


if __name__ == "__main__":
    raise SystemExit(main())
