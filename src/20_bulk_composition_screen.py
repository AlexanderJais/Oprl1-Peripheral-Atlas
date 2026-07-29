"""Can bulk RNA-seq of intact ganglia test the receptor ordering? Screen and find out.

The single-cell atlas rests on tissue that was enzymatically dissociated at 37
degrees, which axotomises every neuron. Bulk RNA-seq of intact ganglia avoids
that, and avoids the pre-mRNA length effect of Figure 2 as well, so it looks
like the orthogonal measurement the ordering needs. It is not, for two reasons
that this screen exists to measure rather than assert.

  COMPOSITION  A ganglion is mostly not neurons. Whole-tissue RNA is dominated
               by satellite glia, Schwann cells and connective tissue, and how
               much of that a dissection includes varies enormously between
               deposits: Snap25 runs from single-digit CPM to over a thousand
               across series that all describe themselves as whole ganglion.
               The four receptors differ in how neuron-restricted they are, so
               composition alone moves the apparent ordering.

  GENE LENGTH  Counts scale with mature transcript length, so an ordering read
               from counts or CPM compares transcript sizes as much as
               abundances. Only a length-normalised unit answers the question,
               and most deposits supply counts.

Every candidate is therefore reported with Snap25, Plp1 and Ptprc beside its
receptor values, so composition is visible rather than assumed, and with Atf3
and Fos, so the dissociation state of the tissue is visible too. A series is
passed for reading only if it clears the neuronal gate AND supplies a
length-normalised unit; everything else is kept in the table with the reason it
failed, in the same way src/11_quality_panel.py records why a check did not run.

Multi-tissue deposits are split by column-name prefix before anything is
computed. GSE248462 deposits arcuate nucleus and nodose ganglion in one matrix,
and read whole it reports the brain.

Run from the repository root: python3 src/20_bulk_composition_screen.py
Source files are cached under scratch/bulk/ and reused.
"""

import gzip
import io
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import atlas_common as ac

CACHE = Path(os.environ.get("SCRATCH", "scratch")) / "bulk"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"

MARKERS = ["Snap25", "Tubb3", "Plp1", "Sox10", "Ptprc"]
INJURY = ["Atf3", "Fos", "Jun"]

# Neuronal gate. GSE248462's nodose samples read Snap25 at 6 CPM against Plp1
# at 893, and every receptor there sits at the noise floor; the arcuate samples
# in the same matrix read Snap25 954. The gate is set to separate those two and
# is reported alongside the values so a reader can move it.
MIN_SNAP25_CPM = 100.0
MIN_SNAP25_OVER_PLP1 = 0.05

# Counts scale with transcript length, so an ordering read from them compares
# gene sizes. Only these units answer the question the screen is asking.
LENGTH_NORMALISED = re.compile(r"tpm|fpkm|rpkm", re.I)

# Columns that are statistics about genes rather than measurements of samples.
NOT_A_SAMPLE = re.compile(
    r"fold|log2|logfc|lfc|pval|p[_.]val|padj|fdr|qval|adj|stat|basemean|"
    r"length|width|chr|start|end|strand|biotype|description|symbol|id$|gene",
    re.I)


def _get(url, timeout=120, tries=3):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code in (403, 429) and i < tries - 1:
                time.sleep(2 ** i * 3)
                continue
            raise
        except Exception:
            if i == tries - 1:
                raise
            time.sleep(2 ** i * 2)
    raise RuntimeError(url)


def candidates(limit=None):
    """Mouse RNA-seq series naming a peripheral ganglion, without isolation language.

    Searched on the tissue rather than on neurons: a whole-tissue study need
    never use the word, and restricting to neuronal terms misses exactly the
    deposits this screen wants.
    """
    cache = CACHE / "candidates.json"
    if cache.exists():
        rows = json.loads(cache.read_text())
    else:
        term = ('(ganglion[All Fields] OR ganglia[All Fields] OR plexus[All Fields])'
                ' AND "Mus musculus"[Organism]'
                ' AND "expression profiling by high throughput sequencing"[DataSet Type]'
                ' AND gse[Entry Type]')
        ids = json.loads(_get(f"{EUTILS}/esearch.fcgi?" + urllib.parse.urlencode(
            {"db": "gds", "term": term, "retmax": 3000, "retmode": "json"}))
        )["esearchresult"]["idlist"]
        recs = []
        for i in range(0, len(ids), 200):
            d = json.loads(_get(f"{EUTILS}/esummary.fcgi?" + urllib.parse.urlencode(
                {"db": "gds", "id": ",".join(ids[i:i + 200]), "retmode": "json"})))["result"]
            recs += [d[k] for k in d.get("uids", [])]
            time.sleep(0.35)

        # Ganglia that are not peripheral ganglia, and preparations that isolate
        # cells before sequencing, which is what this screen is trying to avoid.
        wrong = re.compile(r"retinal? gangli|ganglion cell layer|\bRGC\b|basal gangli|"
                           r"gangliosid|ganglioneur|neuroblastoma", re.I)
        periph = re.compile(
            r"dorsal root|trigeminal|nodose|vagal|jugular|sympathetic|superior cervical|"
            r"stellate|celiac|coeliac|pelvic|sphenopalatine|petrosal|geniculate|"
            r"vestibular gangli|spiral gangli|myenteric|submucosal|enteric|"
            r"peripheral gangli|autonomic gangli|sensory gangli|paravertebral|"
            r"prevertebral|cardiac gangli|\bDRG\b", re.I)
        isolated = re.compile(
            r"single[- ]cell|single[- ]nucle|scRNA|snRNA|10x|chromium|smart-seq|drop-?seq|"
            r"FACS|flow[- ]sort|sorted|retrograde|laser[- ]capture|\bLCM\b|microdissect|"
            r"picked|purified|isolated neuron|cultur|dissociat|ribotag|TRAP|nuclei", re.I)

        rows = []
        for r in recs:
            blob = f"{r.get('title','')} {r.get('summary','')}"
            if int(r.get("n_samples") or 0) < 3:
                continue
            if wrong.search(blob) or not periph.search(blob) or isolated.search(blob):
                continue
            rows.append({"gse": r["accession"], "n": int(r["n_samples"]),
                         "title": r.get("title", "").strip()})
        CACHE.mkdir(parents=True, exist_ok=True)
        cache.write_text(json.dumps(rows, indent=1))
    return rows[:limit] if limit else rows


def matrix_file(gse):
    """The one supplementary file worth reading, or None with the reason."""
    stem = f"{gse[:-3]}nnn"
    try:
        html = _get(f"{FTP}/{stem}/{gse}/suppl/", timeout=60).decode("utf-8", "replace")
    except Exception as e:
        return None, f"supplementary listing unavailable ({type(e).__name__})"
    names = [f for f in re.findall(r'href="([^"?/][^"]*)"', html)
             if not f.startswith(("..", "http"))]
    flat = [f for f in names
            if re.search(r"\.(txt|csv|tsv)(\.gz)?$", f, re.I)
            and "filelist" not in f.lower()]
    if not flat:
        return None, "no flat matrix; per-sample RAW archive only"
    # A length-normalised file is preferred over counts from the same series.
    flat.sort(key=lambda f: (not LENGTH_NORMALISED.search(f), len(f)))
    return flat[0], ""


def read_matrix(gse, fname):
    cached = CACHE / f"{gse}__{fname}"
    if not cached.exists():
        CACHE.mkdir(parents=True, exist_ok=True)
        cached.write_bytes(_get(f"{FTP}/{gse[:-3]}nnn/{gse}/suppl/{fname}"))
    raw = cached.read_bytes()
    if fname.endswith(".gz"):
        raw = gzip.decompress(raw)
    sep = "," if ".csv" in fname.lower() else "\t"
    return pd.read_csv(io.BytesIO(raw), sep=sep, low_memory=False)


ENSEMBL = re.compile(r"ENSMUSG\d{11}")


def _symbol_map():
    """Ensembl gene ID to symbol, from the table src/16 already built."""
    path = ac.RES / "preparation_bias_genomewide.csv"
    if not path.exists():
        return {}
    d = pd.read_csv(path, usecols=["gene", "ensembl_id"])
    return dict(zip(d.ensembl_id, d.gene))


def gene_column(d, symbols=None):
    """The column carrying gene identity, as symbols or as Ensembl IDs.

    A deposit keyed on ENSMUSG identifiers is not a deposit without gene names,
    and dropping those loses several of the few series that reach this far. The
    map comes from the atlas's own annotation, so a symbol resolved here is the
    same symbol used everywhere else in this project.
    """
    for c in d.columns:
        if d[c].astype(str).str.fullmatch("Snap25|Actb|Gapdh", case=False).any():
            return c
    if not symbols:
        return None
    for c in d.columns:
        col = d[c].astype(str)
        if col.str.contains(ENSEMBL, regex=True).mean() > 0.5:
            d["_symbol"] = (col.str.extract(f"({ENSEMBL.pattern})", expand=False)
                            .map(symbols))
            if d["_symbol"].notna().sum() > 1000:
                return "_symbol"
            d.drop(columns=["_symbol"], inplace=True)
    return None


def sample_columns(d):
    return [c for c in d.select_dtypes("number").columns
            if not NOT_A_SAMPLE.search(str(c))]


def column_groups(cols):
    """Split a matrix by column-name prefix, so a two-tissue deposit separates.

    GSE248462 holds arcuate nucleus and nodose ganglion in one file; read whole
    it reports the brain. Trailing replicate and timepoint digits are stripped
    and what remains is the group.
    """
    groups = {}
    for c in cols:
        groups.setdefault(re.sub(r"[\W_]*\d+[\W_]*\d*$", "", str(c)) or "all", []).append(c)
    # A prefix per sample carries no grouping information; treat as one group.
    if len(groups) > max(2, len(cols) // 2):
        return {"all": list(cols)}
    return groups


def screen_one(gse, title, n, symbols=None):
    fname, why = matrix_file(gse)
    base = {"gse": gse, "n_samples": n, "title": title, "file": fname or ""}
    if fname is None:
        return [{**base, "group": "", "verdict": "no matrix", "reason": why}]
    try:
        d = read_matrix(gse, fname)
    except Exception as e:
        return [{**base, "group": "", "verdict": "unreadable",
                 "reason": f"{type(e).__name__}"}]

    gc = gene_column(d, symbols)
    if gc is None:
        return [{**base, "group": "", "verdict": "no gene symbols",
                 "reason": "no gene symbol or Ensembl ID column"}]
    cols = sample_columns(d)
    if len(cols) < 3:
        return [{**base, "group": "", "verdict": "no samples",
                 "reason": f"{len(cols)} numeric non-statistic columns"}]

    unit = "length-normalised" if LENGTH_NORMALISED.search(fname) else "counts"
    d = d.copy()
    d[gc] = d[gc].astype(str)
    out = []
    for group, cc in column_groups(cols).items():
        vals = d[cc].apply(pd.to_numeric, errors="coerce")
        total = vals.sum()
        if not (total > 0).all():
            continue
        cpm = vals.div(total, axis=1) * 1e6

        def level(g):
            m = d[gc].str.fullmatch(g, case=False).values
            return float(cpm[m].values.mean()) if m.any() else np.nan

        row = {**base, "group": group, "n_columns": len(cc), "unit": unit}
        row.update({g: round(level(g), 3) for g in MARKERS + INJURY + ac.RECEPTORS})

        snap, plp = row["Snap25"], row["Plp1"]
        row["Snap25_over_Plp1"] = (round(snap / plp, 4)
                                   if pd.notna(snap) and pd.notna(plp) and plp > 0 else np.nan)
        rec = {g: row[g] for g in ac.RECEPTORS if pd.notna(row.get(g))}
        if len(rec) >= 2:
            s = sorted(rec.items(), key=lambda x: -x[1])
            row["top_receptor"], row["runner_up"] = s[0][0], s[1][0]
            row["margin"] = round(s[0][1] / s[1][1], 3) if s[1][1] > 0 else np.inf
        else:
            row["top_receptor"] = row["runner_up"] = None
            row["margin"] = np.nan

        fails = []
        if pd.isna(snap) or snap < MIN_SNAP25_CPM:
            fails.append(f"Snap25 {snap if pd.notna(snap) else 'absent'} below "
                         f"{MIN_SNAP25_CPM:.0f} CPM")
        if pd.notna(row["Snap25_over_Plp1"]) and row["Snap25_over_Plp1"] < MIN_SNAP25_OVER_PLP1:
            fails.append(f"Snap25/Plp1 {row['Snap25_over_Plp1']:.3f} below "
                         f"{MIN_SNAP25_OVER_PLP1}")
        if unit != "length-normalised":
            fails.append("counts, so an ordering would compare transcript lengths")
        if row["top_receptor"] is None:
            fails.append("fewer than two receptors quantified")
        row["verdict"] = "read" if not fails else "hold"
        row["reason"] = "; ".join(fails)
        out.append(row)
    return out or [{**base, "group": "", "verdict": "no samples",
                    "reason": "no column group had a positive library size"}]


def main() -> int:
    limit = int(os.environ.get("BULK_LIMIT", "0")) or None
    cands = candidates(limit)
    symbols = _symbol_map()
    print(f"  {len(cands)} candidate series to screen\n")

    rows = []
    for i, c in enumerate(cands, 1):
        try:
            got = screen_one(c["gse"], c["title"], c["n"], symbols)
        except Exception as e:
            got = [{"gse": c["gse"], "n_samples": c["n"], "title": c["title"],
                    "group": "", "verdict": "error", "reason": f"{type(e).__name__}: {e}"}]
        rows += got
        v = ", ".join(sorted({g["verdict"] for g in got}))
        print(f"  [{i:>3}/{len(cands)}] {c['gse']:<11} {v}")
        time.sleep(0.2)

    tbl = pd.DataFrame(rows)
    cols = (["gse", "group", "verdict", "reason", "unit", "n_samples", "n_columns"]
            + ac.RECEPTORS + ["top_receptor", "runner_up", "margin"]
            + MARKERS + ["Snap25_over_Plp1"] + INJURY + ["title", "file"])
    tbl = tbl.reindex(columns=[c for c in cols if c in tbl.columns])
    ac.save_table(tbl, "bulk_composition_screen.csv")

    print("\n  Verdicts over every column group:")
    print(tbl.verdict.value_counts().to_string())

    read = tbl[tbl.verdict == "read"]
    print(f"\n  {len(read)} group(s) clear the neuronal gate on a length-normalised unit.")
    if len(read):
        show = ["gse", "group", "Snap25", "Plp1", "Snap25_over_Plp1", "Atf3"] \
               + ac.RECEPTORS + ["top_receptor", "margin"]
        print(read[show].to_string(index=False))

    held = tbl[tbl.verdict == "hold"]
    if len(held):
        why = held.reason.str.split("; ").explode().value_counts()
        print("\n  Why the rest are held (a group can fail on more than one):")
        for reason, k in why.head(8).items():
            print(f"    {k:>4}  {reason}")
        near = held[(held.Snap25 >= MIN_SNAP25_CPM)
                    & (held.Snap25_over_Plp1 >= MIN_SNAP25_OVER_PLP1)]
        print(f"\n  {len(near)} group(s) are neuron-rich enough but supply counts only; "
              "their orderings are not length-corrected and are reported unread:")
        if len(near):
            show = ["gse", "group", "Snap25", "Snap25_over_Plp1", "Atf3",
                    "top_receptor", "runner_up", "margin"]
            print(near.sort_values("Snap25", ascending=False)[show]
                  .head(20).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
