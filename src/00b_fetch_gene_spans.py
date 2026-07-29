"""Fetch genomic spans from Ensembl, once, into data/raw.

Two tables, both keyed on Ensembl gene ID:

  ensembl_gene_spans.csv        the eight opioid genes, for 04_synthesis.py
  ensembl_gene_spans_all.csv.gz every mouse gene, for 16_preparation_bias.py

The per-dataset expression tables carry coordinates for the features they
quantify, but not for every gene and not on one assembly. GSE102443 has no Pnoc
row at all, which silently dropped Pnoc from the gene-length comparison, and it
puts Pdyn at 2 kb against a 13.5 kb locus. One annotation for every gene removes
both problems and makes the spans comparable to each other.

Run from the repository root, before 04_synthesis.py and 16_preparation_bias.py.
"""

import io
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd

import atlas_common as ac

LOOKUP = "https://rest.ensembl.org/lookup/symbol/mus_musculus/{}?expand=0"
BIOMART = "https://www.ensembl.org/biomart/martservice?query="
# BioMart rejects a reflowed query: the Query attributes must stay on one line.
QUERY = """<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" count="" datasetConfigVersion="0.6">
<Dataset name="mmusculus_gene_ensembl" interface="default">
<Attribute name="ensembl_gene_id"/><Attribute name="external_gene_name"/>
<Attribute name="start_position"/><Attribute name="end_position"/>
<Attribute name="gene_biotype"/>
</Dataset></Query>"""


def get_json(url, tries=5):
    """One REST call, retried: Ensembl answers 429 and 500 under load."""
    req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(req, timeout=60) as fh:
                return json.load(fh)
        except urllib.error.HTTPError as err:
            if err.code not in (429, 500, 502, 503) or attempt == tries - 1:
                raise
            wait = 2 ** attempt
            print(f"  [note] Ensembl returned {err.code}; retrying in {wait}s")
            time.sleep(wait)


def opioid_spans():
    """The eight opioid genes, one REST call each, with the symbol checked."""
    rows = []
    for g in ac.OPIOID_GENES:
        d = get_json(LOOKUP.format(g))
        if d["display_name"] != g:
            raise ac.SanityCheckError(
                f"Ensembl returned {d['display_name']} for {g}")
        rows.append({"gene": g, "ensembl_id": d["id"],
                     "chromosome": d["seq_region_name"], "start": d["start"],
                     "end": d["end"], "strand": d["strand"],
                     "span_kb": round((d["end"] - d["start"]) / 1000, 1),
                     "assembly": d["assembly_name"]})
    return pd.DataFrame(rows)


HEADER = ["Gene stable ID", "Gene name", "Gene start (bp)", "Gene end (bp)",
          "Gene type"]


def all_spans(tries=5):
    """Every mouse gene, in one BioMart query.

    BioMart serves its downtime page with HTTP 200, so the body is checked
    against the header we asked for. Parsed as TSV that page becomes a
    one-column frame, and a silent one-column frame here would poison every
    span downstream.
    """
    for attempt in range(tries):
        with urllib.request.urlopen(BIOMART + urllib.parse.quote(QUERY),
                                    timeout=600) as fh:
            txt = fh.read().decode()
        first = txt.split("\n", 1)[0].strip()
        if first.split("\t") == HEADER:
            break
        if attempt == tries - 1:
            raise ac.SanityCheckError(
                "BioMart did not return the requested table; the body began "
                f"{first[:80]!r}")
        wait = 2 ** attempt
        print(f"  [note] BioMart returned {first[:40]!r}; retrying in {wait}s")
        time.sleep(wait)

    d = pd.read_csv(io.StringIO(txt), sep="\t")
    d.columns = ["ensembl_id", "gene", "start", "end", "biotype"]
    d = d.dropna(subset=["gene"])
    d["span_kb"] = (d.end - d.start) / 1000
    return d


def main() -> int:
    ac.DATA.mkdir(parents=True, exist_ok=True)

    few = opioid_spans()
    few.to_csv(ac.DATA / "ensembl_gene_spans.csv", index=False)
    print(few.to_string(index=False))

    every = all_spans()
    missing = set(few.ensembl_id) - set(every.ensembl_id)
    if missing:
        raise ac.SanityCheckError(
            f"BioMart is missing {sorted(missing)}, which the REST lookup found; "
            "the two tables disagree and must not be used together")
    every.to_csv(ac.DATA / "ensembl_gene_spans_all.csv.gz", index=False)
    print(f"\n  {len(every):,} genes, "
          f"{int((every.biotype == 'protein_coding').sum()):,} protein-coding "
          f"-> data/raw/ensembl_gene_spans_all.csv.gz")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
