#!/usr/bin/env bash
# Fetch every source matrix this project reads.
#
#   data/raw/   GEO matrices for the geniculate and NTS analyses  (~150 MB)
#   $NODOSE_ROOT/data/  the NodoMap integrated atlas from CZ CELLxGENE (~830 MB)
#
# Both directories are git-ignored. Re-running skips anything already present
# and verifies checksums, so a truncated or error-page download is caught here
# rather than surfacing later as a confusing parse error.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DIR="$ROOT/data/raw"
NODOSE_ROOT="${NODOSE_ROOT:-/home/user/PNOC-Nodose}"
mkdir -p "$DIR" "$NODOSE_ROOT/data"

# `-f` makes curl fail on an HTTP error instead of writing the error body to the
# output file and exiting 0; `-L` follows redirects. Without -f, `set -e` cannot
# see a 404.
fetch() {  # fetch <url> <destination>
  local url="$1" dest="$2"
  if [ -s "$dest" ]; then
    echo "[skip] $(basename "$dest") already present"
    return
  fi
  echo "[get ] $(basename "$dest")"
  curl -fsSL --retry 3 --retry-delay 2 -o "$dest.partial" "$url"
  mv "$dest.partial" "$dest"
}

GEO=https://ftp.ncbi.nlm.nih.gov/geo/series

# GSE102443 - Dvoryanchikov et al. 2017, 96 geniculate neurons, SMART-seq
fetch "$GEO/GSE102nnn/GSE102443/suppl/GSE102443_GEO-ID_Dvoryanchikov_2017_Datatable_FPKM.txt.gz" \
      "$DIR/GSE102443_GEO-ID_Dvoryanchikov_2017_Datatable_FPKM.txt.gz"

# GSE135801 - Zhang et al. 2019 (Zuker lab), 454 Phox2b+ geniculate neurons
if [ ! -s "$DIR/gse135801/GSM4037432_GG_scRNAseq_Phox2b_expressed454.xlsx" ]; then
  fetch "$GEO/GSE135nnn/GSE135801/suppl/GSE135801_RAW.tar" "$DIR/GSE135801_RAW.tar"
  mkdir -p "$DIR/gse135801"
  tar -xf "$DIR/GSE135801_RAW.tar" -C "$DIR/gse135801"
  rm -f "$DIR/GSE135801_RAW.tar"
fi

# GSE166648 - dorsal vagal complex snRNA-seq (NTS), 72,128 nuclei
fetch "$GEO/GSE166nnn/GSE166648/suppl/GSE166648_snRNA_metadata.csv.gz" \
      "$DIR/GSE166648_snRNA_metadata.csv.gz"
fetch "$GEO/GSE166nnn/GSE166648/suppl/GSE166648_snRNA_unnormdata.csv.gz" \
      "$DIR/GSE166648_snRNA_unnormdata.csv.gz"

# NodoMap integrated atlas - Cheng et al., via CZ CELLxGENE Discover.
# Collection https://cellxgene.cziscience.com/collections/982f9f44-031c-4c8c-91ee-dcaa53b10151
fetch "https://datasets.cellxgene.cziscience.com/a503329b-7bac-4a9f-be64-03cefc978452.h5ad" \
      "$NODOSE_ROOT/data/nodomap_integrated.h5ad"

# Gene coordinates for the preparation-bias analysis come from Ensembl rather
# than from any deposit's own columns; see 00b_fetch_gene_spans.py, which writes
# ensembl_gene_spans.csv, and the BioMart pull that writes the genome-wide table.
echo
echo "Verifying checksums ..."
sha256sum -c --quiet - <<EOF
a9d0fbe064dec4ba0840d07f6fbfc3c2113e6c939271791281cd2a20edd08dfe  $DIR/GSE102443_GEO-ID_Dvoryanchikov_2017_Datatable_FPKM.txt.gz
944d53647f5dac05c1c21ba1abcf00b7dd5580a006b5f9f80cd77172dbe9f143  $DIR/gse135801/GSM4037432_GG_scRNAseq_Phox2b_expressed454.xlsx
d3d4175e43d56587fac0c2f06ab636657bc32425501405b55d9416a423ceb676  $DIR/GSE166648_snRNA_metadata.csv.gz
8490d9ca14820be3c9a70fe1c46df359e8011aa8202a8ae43ee2721d517e8db2  $DIR/GSE166648_snRNA_unnormdata.csv.gz
7a4bfa05b79c076a824c02a03510db427873101eba2101cd58d46e27bd369ba0  $NODOSE_ROOT/data/nodomap_integrated.h5ad
EOF

echo "Done. GEO files in $DIR; NodoMap atlas in $NODOSE_ROOT/data"
