"""Nucleus of the solitary tract: where the nociceptin ligand is made.

Both peripheral ganglia carry the NOP receptor without the ligand, so the
obvious question is where nociceptin acting on those afferents comes from. The
NTS is the first central relay for both: nodose afferents terminate there, and
geniculate gustatory afferents terminate in its rostral pole.

GSE166648 is snRNA-seq of the dorsal vagal complex (NTS, area postrema, DMV),
72,128 nuclei. The expression matrix is a dense genes-by-cells CSV that expands
well beyond the session's disk allowance, so it is streamed once: target gene
rows are kept and a per-gene total is accumulated for abundance matching. The
result is cached so later runs are cheap.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import gzip

import numpy as np
import pandas as pd

import atlas_common as ac

MATRIX = ac.DATA / "GSE166648_snRNA_unnormdata.csv.gz"
META = ac.DATA / "GSE166648_snRNA_metadata.csv.gz"
CACHE = ac.DATA / "gse166648_targets.npz"

TARGETS = ac.OPIOID_GENES + ac.SANITY_GENES + ["Gad1", "Slc17a7", "Glp1r", "Calcr"]


def stream_matrix(force=False):
    """One pass over the dense matrix: target rows plus every gene's total."""
    if CACHE.exists() and not force:
        z = np.load(CACHE, allow_pickle=True)
        cached = list(z["target_names"])
        if set(TARGETS).issubset(set(cached) | z["absent"].tolist()):
            print(f"  [cache] {CACHE.name}")
            return (pd.DataFrame(z["target_mat"], columns=cached,
                                 index=z["cells"]),
                    pd.Series(z["gene_total"], index=z["gene_names"]),
                    list(z["absent"]))

    print(f"  streaming {MATRIX.name} (one pass, then cached) ...")
    want = set(TARGETS)
    rows, gene_names, gene_total = {}, [], []
    with gzip.open(MATRIX, "rt") as fh:
        cells = [c.strip('"') for c in fh.readline().rstrip("\n").split(",")[1:]]
        n = len(cells)
        for i, line in enumerate(fh):
            cut = line.find(",")
            if cut < 0:
                continue
            name = line[:cut].strip('"')
            vals = np.fromstring(line[cut + 1:], sep=",")
            if vals.size != n:
                continue
            gene_names.append(name)
            gene_total.append(vals.sum())
            if name in want:
                rows[name] = vals.astype(np.float32)
            if i % 5000 == 0 and i:
                print(f"    {i:,} genes", flush=True)

    absent = sorted(want - set(rows))
    names = sorted(rows)
    mat = np.vstack([rows[g] for g in names]).T.astype(np.float32)
    np.savez_compressed(CACHE, target_mat=mat, target_names=np.array(names),
                        cells=np.array(cells), gene_names=np.array(gene_names),
                        gene_total=np.array(gene_total), absent=np.array(absent))
    print(f"  [cache] wrote {CACHE.name}: {len(gene_names):,} genes x {n:,} nuclei")
    return (pd.DataFrame(mat, columns=names, index=cells),
            pd.Series(gene_total, index=gene_names), absent)


def main() -> int:
    ac.RES.mkdir(parents=True, exist_ok=True)
    X, gene_total, absent = stream_matrix()
    if absent:
        print(f"  [warn] not in the GSE166648 annotation: {absent}")

    meta = pd.read_csv(META, index_col=0)
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

    # Overall: is the ligand made here, and is it neuronal?
    rows = []
    for g in ac.OPIOID_GENES:
        if g not in cpm.columns:
            rows.append({"gene": g, "in_matrix": False})
            continue
        rows.append({
            "gene": g, "in_matrix": True,
            "neuron_CPM": round(float(cpm.loc[neuron, g].mean()), 3),
            "neuron_pct_detected": round(float((X.loc[neuron, g] > 0).mean() * 100), 2),
            "non_neuron_CPM": round(float(cpm.loc[~neuron, g].mean()), 3),
            "non_neuron_pct_detected": round(float((X.loc[~neuron, g] > 0).mean() * 100), 2),
        })
    overall = pd.DataFrame(rows)
    ac.save_table(overall, "nts_opioid_levels.csv")
    print("\n  NTS/DVC opioid genes:")
    print(overall.to_string(index=False))

    # The comparison that matters: in the NTS the ligand is present alongside
    # the receptor, which is the opposite of both peripheral ganglia.
    neu_mean = cpm.loc[neuron].mean()
    lr = pd.DataFrame([{
        "tissue": "NTS", "dataset": "GSE166648",
        "Oprl1": round(float(neu_mean.get("Oprl1", 0)), 4),
        "Pnoc": round(float(neu_mean.get("Pnoc", 0)), 4),
        "Oprl1_over_Pnoc": ac.ligand_receptor_ratio(neu_mean),
        "Pnoc_in_matrix": bool("Pnoc" in cpm.columns),
    }])
    ac.save_table(lr, "nts_ligand_receptor.csv")

    rank = ac.receptor_rank(neu_mean)
    rank.insert(0, "dataset", "GSE166648")
    rank.insert(0, "tissue", "NTS")
    ac.save_table(rank, "nts_receptor_rank.csv")
    print("\n  Receptor rank in NTS neurons:")
    print(rank.to_string(index=False))

    # Per neuronal subtype, to identify which NTS populations carry the ligand.
    sub = meta["neuronal.subtype"].astype(str)
    srows = []
    for s in sub[neuron].value_counts().index:
        m = neuron & (sub == s).values
        if m.sum() < 30:
            continue
        srows.append({
            "subtype": s, "n": int(m.sum()),
            "Pnoc_CPM": round(float(cpm.loc[m, "Pnoc"].mean()), 2)
                if "Pnoc" in cpm else np.nan,
            "Pnoc_pct": round(float((X.loc[m, "Pnoc"] > 0).mean() * 100), 2)
                if "Pnoc" in X else np.nan,
            "Oprl1_CPM": round(float(cpm.loc[m, "Oprl1"].mean()), 2),
            "Oprl1_pct": round(float((X.loc[m, "Oprl1"] > 0).mean() * 100), 2),
            "Slc17a6_CPM": round(float(cpm.loc[m, "Slc17a6"].mean()), 2)
                if "Slc17a6" in cpm else np.nan,
            "Gad1_CPM": round(float(cpm.loc[m, "Gad1"].mean()), 2)
                if "Gad1" in cpm else np.nan,
        })
    subtbl = pd.DataFrame(srows).sort_values("Pnoc_CPM", ascending=False)
    ac.save_table(subtbl, "nts_by_subtype.csv")
    print("\n  Top NTS neuronal subtypes by Pnoc:")
    print(subtbl.head(10).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
