#!/usr/bin/env bash
# Fetch the GEO matrices into data/raw/ (git-ignored). ~110 MB total.
set -euo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)/data/raw"
mkdir -p "$DIR"; cd "$DIR"

# GSE102443 - Dvoryanchikov et al. 2017, 96 geniculate neurons, SMART-seq
b=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE102nnn/GSE102443/suppl
curl -sS -O "$b/GSE102443_GEO-ID_Dvoryanchikov_2017_All_Transcripts_Read_Counts.txt.gz"
curl -sS -O "$b/GSE102443_GEO-ID_Dvoryanchikov_2017_Datatable_FPKM.txt.gz"

# GSE135801 - Zhang et al. 2019 (Zuker lab), 454 Phox2b+ geniculate neurons
curl -sS -O "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135801/suppl/GSE135801_RAW.tar"
mkdir -p gse135801 && tar -xf GSE135801_RAW.tar -C gse135801 && rm -f GSE135801_RAW.tar

# GSE166648 - dorsal vagal complex snRNA-seq (NTS), 72,128 nuclei
b=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166648/suppl
curl -sS -O "$b/GSE166648_snRNA_metadata.csv.gz"
curl -sS -O "$b/GSE166648_snRNA_unnormdata.csv.gz"

echo "Done. Files in $DIR"
