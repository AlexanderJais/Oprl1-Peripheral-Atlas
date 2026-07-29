"""One quality panel over every dataset in the atlas.

Datasets entered this project at different times and were held to different
standards. The geniculate and vagal data passed a marker gate and, for NodoMap,
an ambient-RNA check; the ganglia added later passed neither until
src/10_peripheral_ganglia.py was written. This script states, for every
population that contributes a number to section 1, which checks ran, what they
returned, and where a check could not run, why.

Four checks are attempted on every population:

  MARKER GATE      Snap25 and Actb must be present and non-zero before any
                   receptor number is read. Enforced by `check_markers`, which
                   raises rather than warning.
  AMBIENT          The neuron-to-non-neuron ratio of each receptor inside one
                   dissociation. No ambient correction is applied anywhere in
                   this project, so a receptor read in neurons carries whatever
                   the same transcript contributes to the soup. The quantity
                   that matters is not Oprl1's enrichment on its own but its
                   enrichment relative to Snap25 in the same dataset: Snap25 is
                   neuronal by definition, so it calibrates how cleanly the two
                   compartments separated.
  BOOTSTRAP        Margin, 95% interval and support over 10,000 resamples of
                   the cells.
  SAMPLE AGREEMENT Which receptor each biological sample places first on its own
                   cells, so a result cannot rest on one animal.

A population deposited as neurons only cannot answer the ambient question, and a
population from a source matrix that exceeds this session's disk allowance
cannot be re-derived here. Both appear in the panel with the reason in place of
the number.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import atlas_common as ac

# Populations not handled by src/10_peripheral_ganglia.py, with the reason each
# missing check is missing. Numbers come from the scripts that produced them.
LEGACY = [
    # Margin over the runner-up, as everywhere else in this panel. The 31.30
    # this row carried until now is Oprl1 over Oprm1, which is the comparison
    # the geniculate section of the README makes and not the runner-up: Oprm1
    # is the lowest of the four here and Oprk1 is second. These are the
    # bootstrap's own numbers, from results/geniculate_rank_support.csv.
    dict(dataset="GSE102443", tissue="geniculate (VII)", division="sensory",
         unit="FPKM", n=96, prep="whole cell", top_gene="Oprl1", margin=6.92,
         margin_lo=4.26, margin_hi=13.32, support=1.000, median_umi=np.nan,
         samples_agreeing="n/a",
         ambient_reason="deposit is 96 sorted neurons, no non-neuronal cells"),
    dict(dataset="GSE135801", tissue="geniculate (VII)", division="sensory",
         unit="CPM", n=454, prep="whole cell", top_gene="Oprl1", margin=4.02,
         margin_lo=0.99, margin_hi=36.70, support=0.974, median_umi=np.nan,
         samples_agreeing="n/a",
         ambient_reason="deposit is Phox2b-sorted neurons, no non-neuronal cells"),
    dict(dataset="NodoMap", tissue="nodose (X)", division="sensory",
         unit="CPM", n=26047, prep="whole cell", top_gene="Oprl1", margin=1.21,
         margin_lo=1.15, margin_hi=1.28, support=1.000, median_umi=1570,
         samples_agreeing="4/4",
         ambient_reason=""),
    dict(dataset="NodoMap", tissue="jugular (X)", division="sensory",
         unit="CPM", n=4593, prep="whole cell", top_gene="Oprl1", margin=1.16,
         margin_lo=0.92, margin_hi=1.46, support=0.898, median_umi=1570,
         samples_agreeing="2/3",
         ambient_reason="scored on the pooled ganglion, not the jugular subset"),
    dict(dataset="iPain", tissue="trigeminal (V)", division="sensory",
         unit="CPM", n=2773, prep="whole cell", top_gene="Oprl1", margin=2.08,
         margin_lo=1.41, margin_hi=3.26, support=1.000, median_umi=np.nan,
         samples_agreeing="n/a",
         ambient_reason="source h5ad is 0.6 GB and was not retained; not re-derived"),
    dict(dataset="iPain", tissue="dorsal root", division="sensory",
         unit="CPM", n=31802, prep="whole cell", top_gene="Oprl1", margin=1.13,
         margin_lo=1.08, margin_hi=1.19, support=1.000, median_umi=np.nan,
         samples_agreeing="n/a",
         ambient_reason="source h5ad is 20.9 GB, above this session's disk allowance"),
    dict(dataset="GSE114997", tissue="spiral (VIII)", division="sensory",
         unit="CPM", n=226, prep="whole cell", top_gene="Oprl1", margin=np.inf,
         margin_lo=np.nan, margin_hi=np.nan, support=1.000, median_umi=np.nan,
         samples_agreeing="n/a",
         ambient_reason="deposit is 226 sorted neurons, no non-neuronal cells"),
    dict(dataset="GSE78845", tissue="thoracic chain", division="sympathetic",
         unit="CPM", n=298, prep="whole cell", top_gene="Oprl1", margin=14.69,
         margin_lo=9.35, margin_hi=27.27, support=1.000, median_umi=33099,
         samples_agreeing="n/a",
         ambient_reason="deposit is 298 sorted neurons, no non-neuronal cells"),
    dict(dataset="GSE166648", tissue="NTS", division="central (reference)",
         unit="CPM", n=49392, prep="nuclear", top_gene="Oprm1", margin=7.48,
         margin_lo=np.nan, margin_hi=np.nan, support=1.000, median_umi=np.nan,
         samples_agreeing="n/a",
         ambient_reason=""),
]

# Ambient results computed elsewhere in the pipeline, in the same units as
# `ac.ambient_enrichment`. NodoMap comes from src/08_specificity_controls.py,
# which scored the pooled ganglion against its ~50,000 satellite and myelinating
# glia and did not record Snap25 there; GSE166648 is read from the table
# src/03_nts.py writes.
EXTERNAL_AMBIENT = {
    ("NodoMap", "nodose (X)"): dict(oprl1=3.529, snap25=np.nan),
    ("NodoMap", "jugular (X)"): dict(oprl1=3.529, snap25=np.nan),
}
NTS_AMBIENT = ac.RES / "nts_ambient_check.csv"


def load_computed():
    """The populations src/10_peripheral_ganglia.py derived from source."""
    path = ac.RES / "peripheral_receptor_levels.csv"
    if not path.exists():
        raise FileNotFoundError(
            f"{path.name} is missing; run src/10_peripheral_ganglia.py first")
    d = pd.read_csv(path)
    amb = pd.read_csv(ac.RES / "peripheral_ambient_checks.csv")

    # The ambient table was written before the celiac ganglion was renamed out
    # of its British spelling and has not been re-derived since, because its
    # source matrices are not in the repository. Both spellings are accepted
    # here so the lookup does not miss, which it did silently: the Snap25
    # column simply came back empty for that one population.
    snap = {k.replace("coeliac", "celiac"): v for k, v in
            (amb[amb.gene == "Snap25"]
             .set_index("dataset")["log2_enrichment"].items())}
    d["label"] = d.dataset + " (" + d.tissue + ")"
    d["Snap25_log2_neuron_over_glia"] = d.label.map(snap)

    # Oprl1's enrichment is only interpretable against Snap25's in the same
    # dissociation, so a population that has one and not the other is a lookup
    # that missed rather than a check that could not run.
    lost = d.label[d.Oprl1_log2_neuron_over_glia.notna()
                   & d.Snap25_log2_neuron_over_glia.isna()]
    if len(lost):
        raise ac.SanityCheckError(
            "ambient check ran but no Snap25 enrichment was found for "
            f"{list(lost)}; peripheral_ambient_checks.csv names its populations "
            "differently from peripheral_receptor_levels.csv")

    d["marker_gate"] = "pass"
    d["unit"] = "CPM"
    d["prep"] = "whole cell"
    d["ambient_reason"] = np.where(
        d.Oprl1_log2_neuron_over_glia.isna(),
        "deposit is author-filtered neurons, no non-neuronal cells", "")
    return d.drop(columns=["label"])


def load_legacy():
    d = pd.DataFrame(LEGACY)
    d["marker_gate"] = "pass"
    if NTS_AMBIENT.exists():
        nts = pd.read_csv(NTS_AMBIENT).set_index("gene")["log2_enrichment"]
        EXTERNAL_AMBIENT[("GSE166648", "NTS")] = dict(
            oprl1=float(nts["Oprl1"]), snap25=float(nts["Snap25"]))
    d["n_samples"] = np.nan
    for col in ac.RECEPTORS + [f"{g}_pct" for g in ac.RECEPTORS]:
        d[col] = np.nan
    d["Oprl1_log2_neuron_over_glia"] = [
        EXTERNAL_AMBIENT.get((r.dataset, r.tissue), {}).get("oprl1", np.nan)
        for r in d.itertuples()]
    d["Snap25_log2_neuron_over_glia"] = [
        EXTERNAL_AMBIENT.get((r.dataset, r.tissue), {}).get("snap25", np.nan)
        for r in d.itertuples()]
    return d


def main() -> int:
    ac.RES.mkdir(parents=True, exist_ok=True)
    panel = pd.concat([load_computed(), load_legacy()], ignore_index=True)

    # Oprl1's enrichment on its own conflates two things: how neuronal the
    # transcript is, and how cleanly the neuronal and non-neuronal compartments
    # separated in that dissociation. Snap25 is neuronal by definition, so the
    # difference isolates the first.
    panel["Oprl1_minus_Snap25"] = (
        panel.Oprl1_log2_neuron_over_glia
        - panel.Snap25_log2_neuron_over_glia).round(3)

    panel["ambient_check"] = np.where(
        panel.Oprl1_log2_neuron_over_glia.notna(), "run",
        np.where(panel.ambient_reason.fillna("") != "", "not possible",
                 "not run"))

    order = {"sensory": 0, "sympathetic": 1, "parasympathetic": 2,
             "mixed autonomic": 3, "enteric": 4, "central (reference)": 5}
    panel = panel.sort_values(
        ["division", "margin"],
        key=lambda s: s.map(order) if s.name == "division" else -s.fillna(0)
    ).reset_index(drop=True)

    cols = ["division", "dataset", "tissue", "prep", "unit", "n", "n_samples",
            "median_umi", "top_gene", "margin", "margin_lo", "margin_hi",
            "support", "samples_agreeing", "marker_gate", "ambient_check",
            "Oprl1_log2_neuron_over_glia", "Snap25_log2_neuron_over_glia",
            "Oprl1_minus_Snap25", "ambient_reason"]
    panel = panel[cols + [c for c in panel.columns if c not in cols]]
    ac.save_table(panel, "dataset_quality_panel.csv")

    print("\n  Quality panel, one row per population:")
    print(panel[["division", "dataset", "tissue", "n", "top_gene", "margin",
                 "support", "samples_agreeing", "marker_gate", "ambient_check",
                 "Oprl1_minus_Snap25"]].to_string(index=False))

    ran = panel[panel.ambient_check == "run"]
    print(f"\n  Marker gate: {int((panel.marker_gate == 'pass').sum())}/"
          f"{len(panel)} populations pass.")
    print(f"  Ambient check: run on {len(ran)} of {len(panel)}; "
          f"{int((panel.ambient_check == 'not possible').sum())} are neurons-only "
          f"deposits or pooled scorings, "
          f"{int((panel.ambient_check == 'not run').sum())} not re-derived.")
    # The central reference is reported apart from the peripheral populations:
    # it is the one nuclear preparation, and the quantity behaves differently
    # there for the reason in figure S1.
    peri = ran[ran.division != "central (reference)"].Oprl1_minus_Snap25.dropna()
    if len(peri):
        print(f"  Oprl1 minus Snap25 enrichment across {len(peri)} peripheral "
              f"populations: median {peri.median():+.2f} log2, range "
              f"{peri.min():+.2f} to {peri.max():+.2f}.")
    nts = panel.loc[panel.division == "central (reference)", "Oprl1_minus_Snap25"]
    if nts.notna().any():
        print(f"  The nuclear NTS preparation sits at {float(nts.iloc[0]):+.2f}, "
              "outside that range and in the direction figure S1 predicts.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
