#!/usr/bin/env bash
# One-time Numbat setup: pull the container and (optionally) download the 1000G reference data.
# Ported from DDE_32 scripts/setup_numbat.sh.
#
# The pkharchenkolab/numbat-rbase container normally ships the SNP VCF / phasing panel / genetic
# map at /data and /Eagle_v2.4.1 (the defaults in nextflow.config). Only download host-side copies
# if your container build lacks them — then point --numbat_snpvcf/--numbat_paneldir/--numbat_gmap at
# the downloaded paths and bind-mount them.
#
# Run inside a SLURM allocation, NOT on the login node (large download + container pull).
#   sbatch --account=biol-stem-2022 --mem=16G --time=06:00:00 --wrap "bash bin/setup_numbat.sh"
set -euo pipefail

cd "$(dirname "$0")/.."

SIF="containers/numbat.sif"
REFDIR="references/numbat"
DOWNLOAD_REFS="${DOWNLOAD_REFS:-false}"   # set true to also fetch host-side 1000G refs

echo "=== Numbat setup ==="

# 1. Pull container as apptainer .sif
if [[ ! -f "${SIF}" ]]; then
    echo "Pulling numbat container -> ${SIF}"
    mkdir -p containers
    apptainer pull "${SIF}" docker://pkharchenkolab/numbat-rbase:latest
else
    echo "Container already present: ${SIF}"
fi

# 2. (Optional) 1000G reference data (hg38)
if [[ "${DOWNLOAD_REFS}" == "true" ]]; then
    mkdir -p "${REFDIR}"
    if [[ ! -f "${REFDIR}/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz" ]]; then
        echo "Downloading 1000G SNP VCF (hg38)..."
        wget -q -P "${REFDIR}" \
            https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz
    fi
    if [[ ! -d "${REFDIR}/1000G_hg38" ]]; then
        echo "Downloading 1000G phasing panel (hg38)..."
        wget -q -P "${REFDIR}" http://pklab.med.harvard.edu/teng/data/1000G_hg38.zip
        ( cd "${REFDIR}" && unzip -q 1000G_hg38.zip && rm -f 1000G_hg38.zip )
    fi
    echo "References in ${REFDIR}/"
fi

echo "=== Numbat setup complete ==="
