"""Nucleus of the solitary tract: Oprl1 in the first central relay.

Both peripheral ganglia in this atlas project to the NTS — nodose afferents
terminate there, and geniculate gustatory afferents in its rostral pole — so it
is the natural place to ask whether the Oprl1 signature continues centrally.

GSE166648 is snRNA-seq of the dorsal vagal complex (NTS, area postrema, DMV),
72,128 nuclei in 25 annotated neuronal subtypes. The expression matrix is a
dense genes-by-cells CSV that expands well beyond the session's disk allowance,
so it is streamed once: target gene rows are kept and a per-gene total is
accumulated in the same pass. The result is cached so later runs are cheap.

This is a nuclear preparation. Long-intron genes gain signal in nuclei (see
src/04_synthesis.py), so the receptor ordering here is reported but should not
be read as a tissue difference.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import atlas_common as ac
import atlas_style as st

MATRIX = ac.DATA / "GSE166648_snRNA_unnormdata.csv.gz"
META = ac.DATA / "GSE166648_snRNA_metadata.csv.gz"
CACHE = ac.DATA / "gse166648_targets.npz"

TARGETS = ac.OPIOID_GENES + ac.SANITY_GENES + ["Gad1", "Slc17a7", "Glp1r", "Calcr"]
MIN_SUBTYPE_CELLS = 30
N_BOOT = 2000


def stream_matrix(force=False):
    """One pass over the dense matrix: target rows plus every gene's total."""
    if CACHE.exists() and not force:
        z = np.load(CACHE)
        cached = [str(s) for s in z["target_names"]]
        absent = [str(s) for s in z["absent"]]
        # `absent` is an array, so it needs converting to a set before union:
        # `set(...) | list` is a TypeError, not an empty union.
        if set(TARGETS).issubset(set(cached) | set(absent)):
            print(f"  [cache] {CACHE.name}")
            return (pd.DataFrame(z["target_mat"], columns=cached,
                                 index=[str(c) for c in z["cells"]]),
                    pd.Series(z["gene_total"],
                              index=[str(g) for g in z["gene_names"]]),
                    absent)
        print(f"  [cache] {CACHE.name} predates the current target list; re-streaming")

    print(f"  streaming {MATRIX.name} (one pass, then cached) ...")
    frame, absent, totals = ac.stream_gene_rows(MATRIX, TARGETS, totals=True,
                                                progress_every=5000)
    np.savez_compressed(
        CACHE,
        target_mat=frame.to_numpy(dtype=np.float32),
        target_names=np.array(frame.columns, dtype=object).astype(str),
        cells=np.array(frame.index, dtype=object).astype(str),
        gene_names=np.array(totals.index, dtype=object).astype(str),
        gene_total=totals.to_numpy(dtype=np.float64),
        absent=np.array(absent, dtype=object).astype(str))
    print(f"  [cache] wrote {CACHE.name}: {len(totals):,} genes x "
          f"{frame.shape[0]:,} nuclei")
    return frame, totals, absent


def main() -> int:
    ac.RES.mkdir(parents=True, exist_ok=True)
    X, gene_total, absent = stream_matrix()
    if absent:
        print(f"  [warn] not in the GSE166648 annotation: {absent}")

    meta = pd.read_csv(META, index_col=0)
    if meta.index.has_duplicates:
        raise ValueError("GSE166648 metadata has duplicate cell barcodes; "
                         "aligning on it would silently duplicate nuclei")
    common = X.index.intersection(meta.index)
    if len(common) != len(X):
        print(f"  [warn] {len(X) - len(common):,} nuclei lack metadata")
    X, meta = X.loc[common], meta.loc[common]
    lib = meta["nCount_RNA"].values.astype(float)
    lib[lib == 0] = 1.0
    cpm = X.div(lib, axis=0) * 1e6
    print(f"  GSE166648: {len(X):,} nuclei, median library {np.median(lib):,.0f}")

    neuron = (meta["cell.type"] == "Neurons").values
    print(f"  neurons: {neuron.sum():,}")

    ac.save_table(ac.check_markers(cpm.loc[neuron].mean(), "GSE166648 (neurons)"),
                  "nts_marker_checks.csv")

    # ------------------------------------------- levels, neurons vs everything
    rows = []
    for g in ac.OPIOID_GENES:
        here = g in cpm.columns
        rows.append({
            "gene": g, "in_matrix": here,
            "neuron_CPM": round(float(cpm.loc[neuron, g].mean()), 3) if here else np.nan,
            "neuron_pct_detected":
                round(float((X.loc[neuron, g] > 0).mean() * 100), 2) if here else np.nan,
            "non_neuron_CPM":
                round(float(cpm.loc[~neuron, g].mean()), 3) if here else np.nan,
            "non_neuron_pct_detected":
                round(float((X.loc[~neuron, g] > 0).mean() * 100), 2) if here else np.nan,
            "n_neurons": int(neuron.sum()),
            "n_non_neurons": int((~neuron).sum()),
        })
    overall = pd.DataFrame(rows)
    ac.save_table(overall, "nts_opioid_levels.csv")
    print("\n  NTS/DVC opioid genes:")
    print(overall.drop(columns=["n_neurons", "n_non_neurons"]).to_string(index=False))

    # ------------------------------------------------------- receptor ordering
    neu_mean = cpm.loc[neuron].mean()
    rank = ac.receptor_rank(neu_mean)
    rank.insert(0, "dataset", "GSE166648")
    rank.insert(0, "tissue", "NTS")
    ac.save_table(rank, "nts_receptor_rank.csv")
    print("\n  Receptor rank in NTS neurons (nuclear prep — see 04_synthesis):")
    print(rank.to_string(index=False))

    b = ac.bootstrap_receptor_support(cpm.loc[neuron], n_boot=N_BOOT)
    b.update({"dataset": "GSE166648", "tissue": "NTS"})
    sup = pd.DataFrame([b])[["tissue", "dataset", "top_gene", "runner_up",
                             "margin", "margin_lo", "margin_hi", "support",
                             "n_cells", "n_boot"]]
    ac.save_table(sup, "nts_rank_support.csv")
    print(sup.round(3).to_string(index=False))

    # --------------------------------------------------- Oprl1 by neural subtype
    if "Oprl1" not in cpm.columns:
        raise ac.SanityCheckError("Oprl1 is absent from the GSE166648 annotation; "
                                  "there is no Oprl1 result to report for the NTS")
    sub = meta["neuronal.subtype"].astype(str)
    srows = []
    for s in sub[neuron].value_counts().index:
        m = neuron & (sub == s).values
        if m.sum() < MIN_SUBTYPE_CELLS:
            continue
        row = {"subtype": s, "n": int(m.sum()),
               "Oprl1_CPM": round(float(cpm.loc[m, "Oprl1"].mean()), 2),
               "Oprl1_pct": round(float((X.loc[m, "Oprl1"] > 0).mean() * 100), 2)}
        for g in ("Slc17a6", "Gad1", "Pnoc"):
            row[f"{g}_CPM"] = (round(float(cpm.loc[m, g].mean()), 2)
                               if g in cpm.columns else np.nan)
        srows.append(row)
    subtbl = pd.DataFrame(srows).sort_values("Oprl1_CPM", ascending=False)
    n_dropped = int(neuron.sum() - subtbl.n.sum())
    if n_dropped:
        print(f"  [note] {n_dropped:,} neurons sit in subtypes below "
              f"{MIN_SUBTYPE_CELLS} cells and are not shown per subtype")
    ac.save_table(subtbl, "nts_by_subtype.csv")
    print("\n  NTS neuronal subtypes by Oprl1:")
    print(subtbl.head(10).to_string(index=False))

    figures(X, cpm, neuron, sub, subtbl, overall)
    return 0


def figures(X, cpm, neuron, sub, subtbl, overall):
    st.set_theme()
    fig, axes = plt.subplots(1, 2, figsize=(15.0, 4.8),
                             gridspec_kw={"width_ratios": [1.0, 2.4]})

    ax = axes[0]
    neu = cpm.loc[neuron]
    vals = [float(neu[g].mean()) if g in neu.columns else np.nan
            for g in ac.RECEPTORS]
    order = np.argsort(-np.array(vals))
    st.expression_bars(ax, [vals[i] for i in order],
                       [ac.RECEPTORS[i] for i in order],
                       "Mean expression (CPM)")
    st.panel_letter(ax, "a", dx=-0.28)
    ax.set_title(f"NTS neurons\\n{int(neuron.sum()):,} nuclei (nuclear prep)",
                 fontsize=11, pad=10)

    ax = axes[1]
    s2 = subtbl.sort_values("Oprl1_CPM", ascending=False)
    st.expression_bars(ax, s2.Oprl1_CPM.values, s2.subtype.values,
                       "Oprl1 (mean CPM)", italic=False, rotation=90, fontsize=11)
    ax.axhline(float(neu["Oprl1"].mean()), color="black", lw=1.2, ls="--")
    ax.text(len(s2) - 0.4, float(neu["Oprl1"].mean()), " NTS mean", fontsize=9,
            va="bottom", ha="right")
    st.panel_letter(ax, "b", dx=-0.06)
    ax.set_title("Oprl1 across all 25 NTS neuronal subtypes", fontsize=11, pad=10)

    fig.tight_layout()
    st.save(fig, "figure5_nts_oprl1")


if __name__ == "__main__":
    raise SystemExit(main())
