#!/usr/bin/env bash
# Prefix CB/CR cell-barcode tags in a BAM with "<sample>__" so cells stay unique when
# Dx + Rel BAMs are merged for a joint per-patient souporcell run.
# Ported from DDE_24 clean_run_for_grant/scripts/08g_souporcell_demux_paired_k2to20_noNK.sh.
#
# Usage: souporcell_retag.sh <sample_prefix> <in_bam_or_cram> <out_bam> [ref_fasta_for_cram]
set -euo pipefail

prefix="$1"
in_aln="$2"
out_bam="$3"
ref="${4:-}"

if [[ "${in_aln}" == *.cram ]]; then
    [[ -n "${ref}" ]] || { echo "CRAM input needs a reference FASTA" >&2; exit 1; }
    view_cmd=(samtools view -h -T "${ref}" "${in_aln}")
else
    view_cmd=(samtools view -h "${in_aln}")
fi

"${view_cmd[@]}" | awk -v pfx="${prefix}__" 'BEGIN{OFS="\t"}
    /^@/ { print; next }
    {
        for (i = 12; i <= NF; i++) {
            if ($i ~ /^CB:Z:/) { sub(/^CB:Z:/, "", $i); $i = "CB:Z:" pfx $i }
            if ($i ~ /^CR:Z:/) { sub(/^CR:Z:/, "", $i); $i = "CR:Z:" pfx $i }
        }
        print
    }' | samtools view -b -o "${out_bam}" -
