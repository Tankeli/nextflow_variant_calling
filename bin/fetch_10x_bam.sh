#!/usr/bin/env bash
# Reconstruct cellranger-ready FASTQ from the ORIGINAL 10X Genomics BAM.
#
# Some GEO/SRA submissions store only the cDNA read in the normalized .sra
# (Statistics nreads="1") -- the cell-barcode read was stripped, so fasterq-dump
# yields no R1 and bin/fetch_sra_10x.sh cannot work. The barcodes survive only in
# the "Original" 10X BAM that the submitter uploaded, kept in the SRA source bucket
# (sra-pub-src-1). This script downloads that BAM and runs 10x `bamtofastq` to
# regenerate R1 (barcode), R2 (cDNA) and I1, concatenating all emitted lanes/chunks
# into a single <sample>_S1_L001_R{1,2,I1}_001.fastq.gz (one pair per sample, as the
# pipeline stages). Affected here: GSE154109 healthy-donor BM (GSM4664009-12).
#
# Usage: fetch_10x_bam.sh <outdir> <SAMPLE:SRR:BAMNAME> [...]
#   BAMNAME is the original file name from the SRA record (e.g. 0064.bam); the
#   script appends the ".1" version suffix used by the source bucket.
# Run on the LOGIN / data-transfer node (compute nodes usually have no internet), e.g.
#   THREADS=8 nohup bash bin/fetch_10x_bam.sh data/controls \
#       HD_BM_1:SRR12185508:0064.bam HD_BM_2:SRR12185509:3958.bam ... &
set -euo pipefail

# 10x standalone bamtofastq (CellRanger 9 dropped the `cellranger bamtofastq` subcommand).
BAMTOFASTQ="${BAMTOFASTQ:-/opt/apps/eb/software/CellRanger/9.0.0/lib/bin/bamtofastq}"
[[ -x "$BAMTOFASTQ" ]] || { echo "ERROR: bamtofastq not found at $BAMTOFASTQ (set BAMTOFASTQ=...)"; exit 1; }

OUTDIR="$1"; shift
THREADS="${THREADS:-8}"
mkdir -p "$OUTDIR"

for spec in "$@"; do
    IFS=':' read -r sample srr bamname <<< "$spec"
    dest="$OUTDIR/$sample"
    final_r2="$dest/${sample}_S1_L001_R2_001.fastq.gz"
    echo "==== $sample ($srr / $bamname) ===="
    if [[ -s "$final_r2" ]]; then echo "  already done; skipping"; continue; fi
    mkdir -p "$dest"; tmp="$dest/tmp"; mkdir -p "$tmp"

    bam="$tmp/${srr}.bam"
    url="https://sra-pub-src-1.s3.amazonaws.com/${srr}/${bamname}.1"
    echo "  downloading original 10X BAM: $url"
    curl -fL --retry 5 --retry-delay 10 -o "$bam" "$url"

    echo "  bamtofastq..."
    rm -rf "$tmp/out"
    "$BAMTOFASTQ" --nthreads="$THREADS" "$bam" "$tmp/out"

    # bamtofastq writes <out>/<flowcell-dir>/bamtofastq_S1_L00N_R{1,2,I1}_001.fastq.gz.
    # Bash sorts globs, so R1/R2/I1 arrays are in matching lane order -> pairs stay aligned.
    shopt -s nullglob
    r1s=( "$tmp"/out/*/*_R1_001.fastq.gz )
    r2s=( "$tmp"/out/*/*_R2_001.fastq.gz )
    i1s=( "$tmp"/out/*/*_I1_001.fastq.gz )
    [[ ${#r1s[@]} -gt 0 && ${#r2s[@]} -gt 0 ]] || { echo "  ERROR: bamtofastq produced no R1/R2 for $srr"; exit 1; }

    echo "  concatenating ${#r1s[@]} R1 / ${#r2s[@]} R2 chunks -> cellranger naming..."
    cat "${r1s[@]}" > "$dest/${sample}_S1_L001_R1_001.fastq.gz"
    cat "${r2s[@]}" > "$dest/${sample}_S1_L001_R2_001.fastq.gz"
    [[ ${#i1s[@]} -gt 0 ]] && cat "${i1s[@]}" > "$dest/${sample}_S1_L001_I1_001.fastq.gz" || true

    rm -f "$bam"; rm -rf "$tmp"
    echo "  done: $dest/${sample}_S1_L001_R{1,2}_001.fastq.gz"
done

echo "All BAM->FASTQ conversions complete."
