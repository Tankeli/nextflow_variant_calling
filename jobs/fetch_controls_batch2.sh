#!/usr/bin/env bash
# Fetch DDE_33 healthy-control batch 2 FASTQ on the login/data-transfer node.
#   - HD_BM_1..4 (GSE154109 healthy-donor BM): SRA dropped the barcode read, so we
#     pull the original 10X BAM and run bamtofastq (bin/fetch_10x_bam.sh).
#   - PBM_1,PBM_2 (paediatric BM, 10x v3, 8 runs each): multi-lane fasterq-dump
#     concatenated to one pair (bin/fetch_sra_10x.sh).
# Outputs land in data/controls/<sample>/<sample>_S1_L001_R{1,2}_001.fastq.gz and
# are consumed by assets/controls_samplesheet_batch2.csv.
set -euo pipefail
cd "$(dirname "$0")/.."
export THREADS="${THREADS:-8}"
OUT=data/controls

echo "### [$(date)] healthy-donor BM via original-BAM + bamtofastq ###"
bash bin/fetch_10x_bam.sh "$OUT" \
    HD_BM_1:SRR12185508:0064.bam \
    HD_BM_2:SRR12185509:3958.bam \
    HD_BM_3:SRR12185510:5286.bam \
    HD_BM_4:SRR12185511:7903.bam

echo "### [$(date)] paediatric BM via multi-lane fasterq-dump ###"
bash bin/fetch_sra_10x.sh "$OUT" \
    PBM_1:SRR12338699,SRR12338700,SRR12338701,SRR12338702,SRR12338703,SRR12338704,SRR12338705,SRR12338706 \
    PBM_2:SRR12338707,SRR12338708,SRR12338709,SRR12338710,SRR12338711,SRR12338712,SRR12338713,SRR12338714

echo "### [$(date)] all batch-2 controls fetched ###"
ls -la "$OUT"/HD_BM_* "$OUT"/PBM_* 2>/dev/null
