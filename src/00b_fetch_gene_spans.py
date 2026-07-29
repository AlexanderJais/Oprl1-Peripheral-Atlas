"""Fetch the genomic span of each opioid gene from Ensembl, once, into data/raw.

The per-dataset expression tables carry coordinates for the features they
quantify, but not for every gene: GSE102443 has no Pnoc row at all, which
silently dropped Pnoc from the gene-length comparison. One annotation for all
eight genes removes that gap and makes the spans comparable to each other.
"""
import json
import urllib.request

import pandas as pd

GENES = ["Oprl1", "Oprm1", "Oprd1", "Oprk1", "Pnoc", "Penk", "Pdyn", "Pomc"]
URL = "https://rest.ensembl.org/lookup/symbol/mus_musculus/{}?expand=0"

rows = []
for g in GENES:
    req = urllib.request.Request(URL.format(g),
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=30) as fh:
        d = json.load(fh)
    if d["display_name"] != g:
        raise SystemExit(f"Ensembl returned {d['display_name']} for {g}")
    rows.append({"gene": g, "ensembl_id": d["id"],
                 "chromosome": d["seq_region_name"], "start": d["start"],
                 "end": d["end"], "strand": d["strand"],
                 "span_kb": round((d["end"] - d["start"]) / 1000, 1),
                 "assembly": d["assembly_name"]})

out = pd.DataFrame(rows)
out.to_csv("data/raw/ensembl_gene_spans.csv", index=False)
print(out.to_string(index=False))
