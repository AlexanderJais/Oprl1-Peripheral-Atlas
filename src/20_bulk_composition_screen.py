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

Matrices are streamed once for the sixteen genes this screen reads and then
deleted: 223 supplementary files do not fit in a session's disk allowance, and
reading isoform-level tables whole exhausted memory and killed the first run.

Run from the repository root: python3 src/20_bulk_composition_screen.py
"""

import csv
import gzip
import itertools
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


# Only these genes are ever read out of a matrix, so a matrix is streamed and
# discarded rather than loaded. Some deposits publish isoform-level tables of
# several hundred thousand rows, and reading those whole exhausted memory and
# killed the first full run at series 83 of 223.
WANTED = set(MARKERS + INJURY + ac.RECEPTORS)
MAX_DOWNLOAD_MB = 400


def stream_matrix(path, sep, symbols):
    """One pass over a matrix: per-column totals, and only the wanted rows.

    Returns (rows, totals, columns). `totals` is the column sum over every gene
    in the file, which is the library size CPM needs; `rows` holds just the
    sixteen genes this screen reads. Memory is bounded by the header width
    rather than by the file.
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.reader(fh, delimiter=sep)
        try:
            header = next(reader)
        except StopIteration:
            return {}, None, []
        head = []
        for row in reader:
            head.append(row)
            if len(head) >= 200:
                break
        if not head:
            return {}, None, []

        width = len(header)
        # A leading index column leaves the header one field short of the rows.
        if len(head[0]) == width + 1:
            header = [""] + header
            width += 1

        def looks_numeric(j):
            ok = 0
            for row in head:
                if len(row) != width:
                    continue
                try:
                    float(row[j]); ok += 1
                except ValueError:
                    pass
            return ok >= 0.8 * len(head)

        cols = [j for j in range(width)
                if not NOT_A_SAMPLE.search(str(header[j])) and looks_numeric(j)]
        if len(cols) < 3:
            return {}, None, []

        # The gene column is whichever non-numeric column carries symbols or
        # Ensembl identifiers in the rows already buffered.
        gene_j = None
        for j in range(width):
            if j in cols:
                continue
            vals = [row[j] for row in head if len(row) == width]
            if any(re.fullmatch("Snap25|Actb|Gapdh", v, re.I) for v in vals):
                gene_j = j; mapper = None; break
            if symbols and sum(bool(ENSEMBL.search(v)) for v in vals) > 0.5 * len(vals):
                gene_j = j; mapper = symbols; break
        if gene_j is None:
            return {}, None, []

        def symbol(v):
            if mapper is None:
                return v
            m = ENSEMBL.search(v)
            return mapper.get(m.group(0)) if m else None

        totals = np.zeros(len(cols))
        rows = {}
        for row in itertools.chain(head, reader):
            if len(row) != width:
                continue
            try:
                vals = np.array([float(row[j]) for j in cols])
            except ValueError:
                continue
            totals += vals
            g = symbol(row[gene_j])
            if g in WANTED and g not in rows:
                rows[g] = vals
        return rows, totals, [str(header[j]) for j in cols]


def fetch(gse, fname):
    """Download to a scratch file, to be streamed and then deleted."""
    CACHE.mkdir(parents=True, exist_ok=True)
    dest = CACHE / f"{gse}__{fname}"
    if not dest.exists():
        dest.write_bytes(_get(f"{FTP}/{gse[:-3]}nnn/{gse}/suppl/{fname}"))
    return dest


ENSEMBL = re.compile(r"ENSMUSG\d{11}")


def _symbol_map():
    """Ensembl gene ID to symbol, from the table src/16 already built."""
    path = ac.RES / "preparation_bias_genomewide.csv"
    if not path.exists():
        return {}
    d = pd.read_csv(path, usecols=["gene", "ensembl_id"])
    return dict(zip(d.ensembl_id, d.gene))


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
    path = None
    try:
        path = fetch(gse, fname)
        size_mb = path.stat().st_size / 1e6
        if size_mb > MAX_DOWNLOAD_MB:
            return [{**base, "group": "", "verdict": "too large",
                     "reason": f"{size_mb:.0f} MB supplementary file"}]
        sep = "," if ".csv" in fname.lower() else "\t"
        rows, totals, colnames = stream_matrix(path, sep, symbols)
    except Exception as e:
        return [{**base, "group": "", "verdict": "unreadable",
                 "reason": f"{type(e).__name__}"}]
    finally:
        # Streamed once and discarded: 223 cached matrices do not fit in the
        # session's disk allowance.
        if path is not None and path.exists():
            path.unlink()

    if totals is None or len(colnames) < 3:
        return [{**base, "group": "", "verdict": "unusable layout",
                 "reason": "no gene column, or fewer than 3 sample columns"}]
    if not rows:
        return [{**base, "group": "", "verdict": "no markers",
                 "reason": "none of the screened genes appear in the matrix"}]

    unit = "length-normalised" if LENGTH_NORMALISED.search(fname) else "counts"
    index = {c: i for i, c in enumerate(colnames)}
    out = []
    for group, cc in column_groups(colnames).items():
        idx = [index[c] for c in cc]
        tot = totals[idx]
        if not (tot > 0).all():
            continue

        def level(g):
            v = rows.get(g)
            return float(np.mean(v[idx] / tot * 1e6)) if v is not None else np.nan

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
