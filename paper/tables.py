"""Manuscript tables, built to Cell Press specification.

Everything under paper/tables/ is a table that appears in the paper, named as
the journal will number it. Table S1 is the inventory behind Figure 1: one row
per panel, giving what was measured, in which preparation and at what depth,
and what the receptor ordering came out at.

Every number is read from results/. What this module holds itself is manuscript
metadata that no analysis produces -- which study a deposit comes from and which
chemistry it was run on -- and two library-depth and two detection values that
the manuscript states from a deposit's own libraries but no results table
records. Each of those is commented with where it came from. A number typed
into a script is how the geniculate margin in src/11_quality_panel.py came to
disagree with the bootstrap that produced it, so `check_margins` recomputes
every margin from the levels Figure 1 is drawn from and refuses to write the
table if a recorded margin disagrees with it.

Run from the repository root: python3 paper/tables.py
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import numpy as np
import pandas as pd

import atlas_common as ac

OUT = Path(__file__).resolve().parent / "tables"

# Deposit to study and chemistry. Provenance for the manuscript rather than
# analysis output, taken from the Provenance section of README.md. A deposit
# whose publication this survey did not establish carries an em dash rather
# than a guess.
PROVENANCE = {
    "GSE114997": ("Shrestha et al., 2018", "SMART-seq, full-length"),
    "GSE102443": ("Dvoryanchikov et al., 2017", "SMART-seq, full-length"),
    "GSE309608": ("--", "10x droplet, 3'"),
    "GSE135801": ("Zhang et al., 2019", "droplet, 3'"),
    "iPain": ("Bhuiyan et al., 2024 (iPain Atlas)", "10x droplet, 3'"),
    "NodoMap": ("Cheng et al., 2026 (NodoMap)", "10x droplet, 3'"),
    "GSE232789": ("--", "10x droplet, 3'"),
    "GSE231766": ("Ziegler et al., 2023", "10x droplet, 3'"),
    "GSE78845": ("Furlan et al., 2016", "full-length"),
    "GSE231924": ("--", "10x droplet, 3'"),
    "GSE330884": ("--", "10x droplet, 3'"),
    "GSE263422": ("--", "10x droplet, 3'"),
}

# Values the manuscript states that results/dataset_quality_panel.csv leaves
# empty, because the populations they belong to predate
# src/10_peripheral_ganglia.py and were never re-derived through it. Both are
# read from results/EXTERNAL_GANGLIA_NOTES.md, which records them from the
# deposits' own libraries. Anything not here stays empty and is footnoted as
# not recorded, rather than being filled in from a plausible neighbour.
DEPTH_FROM_NOTES = {
    "GSE114997": 2_678_701,     # median reads per neuron, full-length
}
DETECTION_FROM_NOTES = {
    "GSE114997|spiral": 74.8,
    "GSE78845|thoracic chain": 82.2,
}
SAMPLES_AGREEING_FROM_NOTES = {
    "GSE231924|stellate": "8/8",     # all eight mice place Oprl1 first
}

# A margin computed from the two-decimal levels in the figure table cannot
# reproduce a bootstrap margin computed on full precision exactly. Anything
# further apart than this is a disagreement about which pair of genes the
# margin is over, which is what the check exists to catch.
MARGIN_TOL = 0.02

COLUMNS = [
    "panel", "division", "population", "dataset", "study", "platform",
    "preparation", "unit", "n_neurons", "n_samples", "median_library",
    "Oprl1", "Oprm1", "Oprd1", "Oprk1", "Oprl1_pct_detected",
    "top_receptor", "runner_up", "margin", "margin_lo", "margin_hi",
    "bootstrap_support", "samples_agreeing", "marker_gate", "ambient_check",
    "notes",
]


def key(dataset, tissue):
    """Join key shared by the figure table and the quality panel.

    The two disagree on whether a cranial nerve numeral belongs in the tissue
    name -- "sphenopalatine (VII)" against "sphenopalatine" -- and on nothing
    else. Dataset and stripped tissue are unique in both, which `build` asserts
    rather than assumes.
    """
    return f"{dataset}|{str(tissue).split(' (')[0].strip()}"


def ordering(levels: pd.Series):
    """Top receptor, runner-up and margin, from the levels Figure 1 draws.

    A population whose second receptor reads exactly zero has no runner-up to
    take a ratio against: the spiral ganglion quantifies all four genes and
    detects one, so its margin is infinite rather than large.
    """
    s = levels.sort_values(ascending=False)
    top, second = s.index[0], s.index[1]
    if s.iloc[1] <= 0:
        return top, None, np.inf
    return top, second, float(s.iloc[0] / s.iloc[1])


def detection_rates():
    """Per-population Oprl1 detection, from whichever table recorded it.

    The quality panel carries it for the populations src/10_peripheral_ganglia.py
    derived; the geniculate and vagal populations predate that script and keep
    theirs in the per-tissue tables their own scripts wrote.
    """
    out = {}

    qual = pd.read_csv(ac.RES / "dataset_quality_panel.csv")
    for r in qual[qual.Oprl1_pct.notna()].itertuples():
        out[key(r.dataset, r.tissue)] = float(r.Oprl1_pct)

    gen = pd.read_csv(ac.RES / "geniculate_opioid_levels.csv")
    for r in gen[gen.gene == "Oprl1"].itertuples():
        out[key(r.dataset, r.tissue)] = float(r.pct_detected)

    nod = pd.read_csv(ac.RES / "nodose_opioid_levels.csv")
    nod = nod[(nod.gene == "Oprl1") & (nod.dataset == "NodoMap:whole-cell")]
    for r in nod.itertuples():
        out[key("NodoMap", r.tissue)] = float(r.pct_detected)

    out.update(DETECTION_FROM_NOTES)
    return out


def check_margins(tbl: pd.DataFrame, recorded: pd.Series):
    """Refuse to write a table whose margins disagree with its own levels.

    `recorded` is what results/ says the margin is; `tbl.margin` is the ratio of
    the top receptor to the runner-up in the same levels Figure 1 is drawn from.
    A recorded margin taken over a gene that is not the runner-up reads as a
    large relative error here and stops the build.
    """
    bad = []
    for r in tbl.itertuples():
        rec, got = float(recorded[r.Index]), float(r.margin)
        if np.isinf(rec) and np.isinf(got):
            continue
        if np.isnan(rec):
            continue
        if not np.isfinite(rec) or not np.isfinite(got) or \
                abs(rec - got) / got > MARGIN_TOL:
            bad.append(f"panel {r.panel.upper()} {r.population} "
                       f"({r.dataset}): results/ records {rec:.2f}, the levels "
                       f"give {got:.2f} over {r.runner_up}")
    if bad:
        raise ac.SanityCheckError(
            "recorded margins disagree with the levels Figure 1 is drawn from:\n  "
            + "\n  ".join(bad))


def build() -> pd.DataFrame:
    fig = pd.read_csv(ac.RES / "receptor_levels_by_ganglion.csv")
    # `samples_agreeing` carries "n/a" for a deposit that does not resolve into
    # biological samples and an empty field for one whose agreement was never
    # recorded. Read as a converter so the two stay distinguishable: pandas
    # reads "n/a" as a null by default and the distinction is lost.
    qual = pd.read_csv(ac.RES / "dataset_quality_panel.csv",
                       converters={"samples_agreeing": str})
    det = detection_rates()

    qual["key"] = [key(r.dataset, r.tissue) for r in qual.itertuples()]
    if qual.key.duplicated().any():
        raise ac.SanityCheckError(
            f"dataset_quality_panel.csv: duplicate keys "
            f"{sorted(qual.key[qual.key.duplicated()])}")
    qual = qual.set_index("key")

    rows = []
    for r in fig.itertuples():
        k = key(r.dataset, r.tissue)
        if k not in qual.index:
            raise ac.SanityCheckError(
                f"panel {r.panel.upper()} ({r.dataset}, {r.tissue}) has no row "
                "in dataset_quality_panel.csv; run src/11_quality_panel.py")
        q = qual.loc[k]
        levels = pd.Series({g: getattr(r, g) for g in ac.RECEPTORS})
        top, second, margin = ordering(levels)
        study, platform = PROVENANCE[r.dataset]

        rows.append({
            "panel": str(r.panel).upper(),
            "division": r.group,
            "population": r.tissue,
            "dataset": r.dataset,
            "study": study,
            "platform": platform,
            "preparation": q.prep,
            "unit": r.unit,
            "n_neurons": int(r.n),
            "n_samples": q.n_samples,
            "median_library": (q.median_umi if pd.notna(q.median_umi)
                               else DEPTH_FROM_NOTES.get(r.dataset, np.nan)),
            **{g: levels[g] for g in ac.RECEPTORS},
            "Oprl1_pct_detected": det.get(k, np.nan),
            "top_receptor": top,
            "runner_up": second,
            "margin": margin,
            "margin_lo": q.margin_lo,
            "margin_hi": q.margin_hi,
            "bootstrap_support": q.support,
            "samples_agreeing": (q.samples_agreeing
                                 or SAMPLES_AGREEING_FROM_NOTES.get(k, "")),
            "marker_gate": q.marker_gate,
            "ambient_check": q.ambient_check,
            "notes": q.ambient_reason if isinstance(q.ambient_reason, str) else "",
            "_recorded_margin": q.margin,
        })

    tbl = pd.DataFrame(rows)
    check_margins(tbl, tbl.pop("_recorded_margin"))
    # The recorded margin is the bootstrap's, computed on full precision; the
    # check above is what licenses reporting it against these levels.
    tbl["margin"] = [q for q in
                     (qual.loc[key(d, t), "margin"]
                      for d, t in zip(tbl.dataset, tbl.population))]

    tbl["runner_up"] = tbl.runner_up.fillna("none detected")
    for col in ("n_samples", "median_library"):
        tbl[col] = tbl[col].round().astype("Int64")
    return tbl[COLUMNS]


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

def _num(v, dp=2, thousands=False, dash="--"):
    """A number, or the dash that says the source tables do not record one."""
    if v is None or pd.isna(v):
        return dash
    return f"{v:,.{dp}f}" if thousands else f"{v:.{dp}f}"


def _markdown(df: pd.DataFrame, headers) -> str:
    rows = [[str(v) for v in row] for row in df.itertuples(index=False)]
    widths = [max(len(h), *(len(r[i]) for r in rows))
              for i, h in enumerate(headers)]
    out = ["| " + " | ".join(h.ljust(w) for h, w in zip(headers, widths)) + " |",
           "|" + "|".join("-" * (w + 2) for w in widths) + "|"]
    for r in rows:
        out.append("| " + " | ".join(v.ljust(w) for v, w in zip(r, widths)) + " |")
    return "\n".join(out)


def render(tbl: pd.DataFrame) -> str:
    """Table S1 as two blocks over one panel key.

    The table is 26 columns wide and belongs in one sheet, which is what the
    CSV is. Split here only because 26 columns do not read on a page: (a) is
    where each population came from, (b) is what it measured.
    """
    t = tbl.copy()
    panel = t.panel.str.upper()

    prov = pd.DataFrame({
        "Panel": panel,
        "Division": t.division,
        "Population": t.population,
        "Deposit": t.dataset,
        "Study": t.study,
        "Platform": t.platform,
        "Preparation": t.preparation,
        "Neurons": t.n_neurons.map("{:,}".format),
        "Samples": [_num(v, 0) for v in t.n_samples],
        "Median library": [_num(v, 0, thousands=True) for v in t.median_library],
    })

    ci = [f"{_num(lo)} to {_num(hi)}" if pd.notna(lo) and pd.notna(hi) else "--"
          for lo, hi in zip(t.margin_lo, t.margin_hi)]
    meas = pd.DataFrame({
        "Panel": panel,
        "Population": t.population,
        "Unit": t.unit,
        "Oprl1": [_num(v) for v in t.Oprl1],
        "Oprm1": [_num(v) for v in t.Oprm1],
        "Oprd1": [_num(v) for v in t.Oprd1],
        "Oprk1": [_num(v) for v in t.Oprk1],
        "Oprl1 detected (%)": [_num(v, 1) for v in t.Oprl1_pct_detected],
        "Top": t.top_receptor,
        # A population with no second receptor to divide by has no margin, and
        # writing one in would turn three measured zeros into a large number.
        "Runner-up": t.runner_up.fillna("none detected"),
        "Margin": [("undefined" if not np.isfinite(v) else f"{v:.2f}x")
                   for v in t.margin],
        "95% CI": ci,
        "Support": [_num(v, 3) for v in t.bootstrap_support],
        "Samples agreeing": t.samples_agreeing.replace("", "--"),
        "Marker gate": t.marker_gate,
        "Ambient check": t.ambient_check,
    })

    notes = [f"({r.panel.upper()}) {r.notes}" for r in t.itertuples()
             if isinstance(r.notes, str) and r.notes]
    return "\n\n".join([
        "# Table S1",
        "*Related to Figure 1. The legend is in "
        "[`paper/table_legends.md`](../table_legends.md); the table the journal "
        "receives is [`TableS1.csv`](TableS1.csv), of which this is a rendering.*",
        "## (a) Provenance and sample", _markdown(prov, list(prov.columns)),
        "## (b) Levels, ordering and quality control",
        _markdown(meas, list(meas.columns)),
        "## Why an ambient check could not run", "\n".join(f"- {n}" for n in notes),
    ]) + "\n"


def main() -> int:
    tbl = build()
    OUT.mkdir(parents=True, exist_ok=True)
    tbl.to_csv(OUT / "TableS1.csv", index=False)
    (OUT / "TableS1.md").write_text(render(tbl))
    print(f"  [table] paper/tables/TableS1.csv|.md  "
          f"{len(tbl)} populations x {len(tbl.columns)} columns")

    missing = tbl.loc[tbl.Oprl1_pct_detected.isna(), "panel"].str.upper().tolist()
    if missing:
        print(f"  detection rate not recorded for panels {', '.join(missing)}")
    missing = tbl.loc[tbl.median_library.isna(), "panel"].str.upper().tolist()
    if missing:
        print(f"  median library not recorded for panels {', '.join(missing)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
