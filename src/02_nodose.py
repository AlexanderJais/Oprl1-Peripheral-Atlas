"""Nodose and jugular ganglion: the same statistics, from the NodoMap atlas.

Reads the integrated atlas maintained in the sibling PNOC-Nodose project. Two
differences from that project's own analysis, both required for comparability:

  * levels are pseudobulk mean CPM computed from raw counts, so the geniculate
    FPKM/CPM means and these are the same kind of quantity (they are still not
    interchangeable in magnitude, which is why only ranks and ratios cross
    tissues);
  * the five constituent datasets are summarised separately as well as pooled,
    because the in-house snRNA-seq inflates long-intron genes such as Oprm1 and
    would otherwise distort the receptor ordering.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

import atlas_common as ac

NEURON_CLASSES = ["Nodose ganglion neuron", "Jugular ganglion neuron"]


def load_counts(adata, genes):
    """Raw integer counts for a few genes, as a cells x genes frame."""
    sym = adata.raw.var["feature_name"].astype(str)
    lookup = pd.Series(np.arange(adata.raw.n_vars), index=sym.values)
    out = {}
    for g in genes:
        if g not in lookup.index:
            print(f"  [warn] {g} absent from the atlas")
            continue
        hit = lookup.loc[g]
        j = int(hit if np.isscalar(hit) else hit.iloc[0])
        col = adata.raw.X[:, j]
        out[g] = np.asarray(col.todense()).ravel() if sp.issparse(col) else np.asarray(col).ravel()
    return pd.DataFrame(out, index=adata.obs_names)


def pseudobulk_cpm(counts, total_counts, mask):
    """Mean per-cell CPM over `mask`, matching how the geniculate CPM is built."""
    sub = counts.loc[mask]
    lib = total_counts[mask].astype(float)
    lib[lib == 0] = 1.0
    return (sub.div(lib, axis=0) * 1e6).mean()


def main() -> int:
    ac.RES.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(ac.NODOSE_H5AD, backed="r")
    obs = adata.obs
    genes = ac.OPIOID_GENES + ac.SANITY_GENES
    counts = load_counts(adata, genes)

    # nCount_RNA is the library size the authors recorded; fall back to the sum
    # over the genes present only if it is missing.
    total = obs["nCount_RNA"].values.astype(float)
    print(f"  NodoMap: {adata.n_obs:,} cells, median library {np.median(total):,.0f}")

    is_neuron = obs["cell_class"].isin(NEURON_CLASSES).values
    is_nodose = (obs["cell_class"] == "Nodose ganglion neuron").values
    is_jugular = (obs["cell_class"] == "Jugular ganglion neuron").values

    rows, means = [], {}
    groups = [
        ("nodose", "NodoMap", is_nodose),
        ("jugular", "NodoMap", is_jugular),
        ("nodose+jugular", "NodoMap", is_neuron),
    ]
    # Per-dataset, neurons only: the ordering must not rest on one study.
    for ds in obs["dataset"].cat.categories:
        groups.append((f"nodose+jugular", f"NodoMap:{ds}",
                       is_neuron & (obs["dataset"] == ds).values))

    for tissue, label, mask in groups:
        if mask.sum() == 0:
            continue
        m = pseudobulk_cpm(counts, total, mask)
        det = (counts.loc[mask] > 0).mean() * 100.0
        means[(tissue, label)] = m
        for g in ac.OPIOID_GENES:
            rows.append({
                "dataset": label, "tissue": tissue,
                "assay": ac.ASSAY_LABEL["NodoMap"], "gene": g, "unit": "CPM",
                "mean_level": round(float(m.get(g, 0.0)), 4),
                "pct_detected": round(float(det.get(g, 0.0)), 2),
                "n_cells": int(mask.sum()),
                "in_matrix": True,
            })
    tbl = pd.DataFrame(rows)
    ac.save_table(tbl, "nodose_opioid_levels.csv")

    print("\n  Opioid gene levels, mean CPM:")
    piv = tbl[tbl.dataset == "NodoMap"].pivot(index="gene", columns="tissue",
                                              values="mean_level")
    print(piv.reindex(ac.OPIOID_GENES).to_string())

    ranks = []
    for (tissue, label), m in means.items():
        r = ac.receptor_rank(m)
        r.insert(0, "dataset", label)
        r.insert(0, "tissue", tissue)
        ranks.append(r)
    rank_tbl = pd.concat(ranks, ignore_index=True)
    ac.save_table(rank_tbl, "nodose_receptor_rank.csv")
    print("\n  Receptor rank (neurons), pooled and per dataset:")
    for label in rank_tbl.dataset.unique():
        b = rank_tbl[(rank_tbl.dataset == label) & (rank_tbl.tissue.str.contains("nodose"))]
        if b.empty:
            continue
        b = b[b.tissue == b.tissue.iloc[0]]
        print(f"    {label:20s} " + " > ".join(b.sort_values('rank').gene))

    lr = []
    for (tissue, label), m in means.items():
        lr.append({
            "tissue": tissue, "dataset": label,
            "Oprl1": round(float(m.get("Oprl1", 0.0)), 4),
            "Pnoc": round(float(m.get("Pnoc", 0.0)), 4),
            "Oprl1_over_Pnoc": ac.ligand_receptor_ratio(m),
            "Pnoc_in_matrix": True,
        })
    lr_tbl = pd.DataFrame(lr)
    ac.save_table(lr_tbl, "nodose_ligand_receptor.csv")
    print("\n  Receptor vs ligand in the same cells:")
    print(lr_tbl[lr_tbl.dataset == "NodoMap"].to_string(index=False))

    # Neuron-vs-non-neuron enrichment, for the same abundance-matched statistic
    # the geniculate data cannot supply (it contains only neurons).
    m_neu = pseudobulk_cpm(counts, total, is_neuron)
    m_non = pseudobulk_cpm(counts, total, ~is_neuron)
    enr = pd.DataFrame({
        "gene": ac.OPIOID_GENES,
        "neuron_CPM": [round(float(m_neu.get(g, 0)), 4) for g in ac.OPIOID_GENES],
        "non_neuron_CPM": [round(float(m_non.get(g, 0)), 4) for g in ac.OPIOID_GENES],
    })
    enr["log2_enrichment"] = np.log2((enr.neuron_CPM + 0.01) /
                                     (enr.non_neuron_CPM + 0.01)).round(3)
    ac.save_table(enr, "nodose_neuron_enrichment.csv")
    print("\n  Neuronal enrichment (mean CPM):")
    print(enr.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
