#!/usr/bin/env bash
# Download 10x Genomics SRA runs and convert to cellranger-ready FASTQ.
# Recovers the technical reads and renames by detected read length to the
# cellranger convention <sample>_S1_L001_R{1,2}_001.fastq.gz (+ I1 if present).
#
# Usage: fetch_sra_10x.sh <outdir> <SAMPLE:SRR> [<SAMPLE:SRR> ...]
# Run on the LOGIN / data-transfer node (compute nodes usually have no internet), e.g.
#   THREADS=8 nohup bash bin/fetch_sra_10x.sh data/controls PBMMC_1:SRR9264351 ... &
set -euo pipefail

module load SRA-Toolkit/3.2.0 2>/dev/null || true
command -v pigz >/dev/null 2>&1 && GZIP="pigz -p ${THREADS:-8}" || GZIP="gzip"

OUTDIR="$1"; shift
THREADS="${THREADS:-8}"
mkdir -p "$OUTDIR"

for pair in "$@"; do
    sample="${pair%%:*}"; srr="${pair##*:}"
    dest="$OUTDIR/$sample"
    final_r2="$dest/${sample}_S1_L001_R2_001.fastq.gz"
    echo "==== $sample ($srr) ===="
    if [[ -s "$final_r2" ]]; then echo "  already done; skipping"; continue; fi
    mkdir -p "$dest"
    tmp="$dest/tmp"; mkdir -p "$tmp"

    echo "  prefetch..."
    prefetch -O "$tmp" --max-size 200G "$srr"

    echo "  fasterq-dump (split + technical)..."
    fasterq-dump --split-files --include-technical -e "$THREADS" \
        -t "$tmp" -O "$tmp" "$tmp/$srr/$srr.sra"

    # Classify each split file by modal read length: I1<=12, R1<=30 (v2 26bp), else R2 (cDNA).
    R1=""; R2=""; I1=""
    for f in "$tmp/${srr}"_*.fastq; do
        [[ -f "$f" ]] || continue
        len=$(awk 'NR==2{print length($0); exit}' "$f")
        echo "    $(basename "$f"): read length ${len}"
        if   [[ "$len" -le 12 ]]; then I1="$f"
        elif [[ "$len" -le 30 ]]; then R1="$f"
        else                          R2="$f"; fi
    done
    [[ -n "$R1" && -n "$R2" ]] || { echo "  ERROR: could not identify R1/R2 for $srr"; exit 1; }

    echo "  compressing -> cellranger naming..."
    $GZIP -c "$R1" > "$dest/${sample}_S1_L001_R1_001.fastq.gz"
    $GZIP -c "$R2" > "$dest/${sample}_S1_L001_R2_001.fastq.gz"
    [[ -n "$I1" ]] && $GZIP -c "$I1" > "$dest/${sample}_S1_L001_I1_001.fastq.gz" || true

    rm -rf "$tmp"
    echo "  done: $dest/${sample}_S1_L001_R{1,2}_001.fastq.gz"
done

echo "All SRA downloads complete."
