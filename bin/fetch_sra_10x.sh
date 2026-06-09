#!/usr/bin/env bash
# Download 10x Genomics SRA runs and convert to cellranger-ready FASTQ.
# Recovers the technical reads and renames by detected read length to the
# cellranger convention <sample>_S1_L001_R{1,2}_001.fastq.gz (+ I1 if present).
#
# A sample may map to ONE run or to SEVERAL runs (comma-separated) when the
# library was sequenced across multiple SRA runs/lanes — those are concatenated
# into a single L001 R1/R2 pair (cellranger dedups on CB+UMI, so lane info is
# not needed; this matches the one-pair-per-sample staging the pipeline expects).
#
# Usage: fetch_sra_10x.sh <outdir> <SAMPLE:SRR[,SRR2,...]> [...]
# Run on the LOGIN / data-transfer node (compute nodes usually have no internet), e.g.
#   THREADS=8 nohup bash bin/fetch_sra_10x.sh data/controls PBMMC_1:SRR9264351 ... &
#   THREADS=8 nohup bash bin/fetch_sra_10x.sh data/controls \
#       PBM_1:SRR12338699,SRR12338700,SRR12338701,SRR12338702,SRR12338703,SRR12338704,SRR12338705,SRR12338706 ... &
set -euo pipefail

module load SRA-Toolkit/3.2.0-gompi-2024a 2>/dev/null || module load SRA-Toolkit 2>/dev/null || true
command -v pigz >/dev/null 2>&1 && GZIP="pigz -p ${THREADS:-8}" || GZIP="gzip"

OUTDIR="$1"; shift
THREADS="${THREADS:-8}"
mkdir -p "$OUTDIR"

for pair in "$@"; do
    sample="${pair%%:*}"; srrs="${pair##*:}"
    dest="$OUTDIR/$sample"
    final_r2="$dest/${sample}_S1_L001_R2_001.fastq.gz"
    echo "==== $sample (${srrs}) ===="
    if [[ -s "$final_r2" ]]; then echo "  already done; skipping"; continue; fi
    mkdir -p "$dest"
    tmp="$dest/tmp"; mkdir -p "$tmp"

    # Per-sample accumulators (one per read type) concatenated across all runs.
    : > "$tmp/R1.fastq"; : > "$tmp/R2.fastq"; : > "$tmp/I1.fastq"; have_i1=0
    IFS=',' read -ra srr_arr <<< "$srrs"
    for srr in "${srr_arr[@]}"; do
        echo "  -- $srr: prefetch..."
        prefetch -O "$tmp" --max-size 200G "$srr"
        echo "  -- $srr: fasterq-dump (split + technical)..."
        fasterq-dump --split-files --include-technical -e "$THREADS" \
            -t "$tmp" -O "$tmp" "$tmp/$srr/$srr.sra"

        # Classify each split file by modal read length: I1<=12, R1<=30 (v2 26 / v3 28), else R2 (cDNA).
        R1=""; R2=""; I1=""
        for f in "$tmp/${srr}"_*.fastq; do
            [[ -f "$f" ]] || continue
            len=$(awk 'NR==2{print length($0); exit}' "$f")
            echo "       $(basename "$f"): read length ${len}"
            if   [[ "$len" -le 12 ]]; then I1="$f"
            elif [[ "$len" -le 30 ]]; then R1="$f"
            else                          R2="$f"; fi
        done
        [[ -n "$R1" && -n "$R2" ]] || { echo "  ERROR: could not identify R1/R2 for $srr"; exit 1; }

        # Append this run to the per-sample accumulators (R1/R2 in lockstep keeps pairs aligned).
        cat "$R1" >> "$tmp/R1.fastq"; cat "$R2" >> "$tmp/R2.fastq"
        [[ -n "$I1" ]] && { cat "$I1" >> "$tmp/I1.fastq"; have_i1=1; }
        rm -f "$tmp/${srr}"_*.fastq; rm -rf "$tmp/$srr"
    done

    echo "  compressing -> cellranger naming..."
    $GZIP -c "$tmp/R1.fastq" > "$dest/${sample}_S1_L001_R1_001.fastq.gz"
    $GZIP -c "$tmp/R2.fastq" > "$dest/${sample}_S1_L001_R2_001.fastq.gz"
    [[ "$have_i1" -eq 1 ]] && $GZIP -c "$tmp/I1.fastq" > "$dest/${sample}_S1_L001_I1_001.fastq.gz" || true

    rm -rf "$tmp"
    echo "  done: $dest/${sample}_S1_L001_R{1,2}_001.fastq.gz"
done

echo "All SRA downloads complete."
