#!/usr/bin/env bash
# Build all DDE_27 per-stage container images as Apptainer .sif files.
#
# RUN ON THE VIKING LOGIN / DATA-TRANSFER NODE: compute nodes have no internet, so all image
# building and registry pulls must happen here (same rule as DDE_33). The resulting .sif files are
# what conf/viking.config points the *_container params at (containers/sif/*.sif).
#
# Why Apptainer .def (not the Dockerfiles): podman's overlay/vfs storage hits xattr / user-namespace
# limits on Lustre, and there is no docker daemon. Apptainer builds straight from the docker:// base
# in a root-mapped namespace (no fakeroot needed) and installs everything via micromamba — including
# procps (for `ps`, used by Nextflow's task metrics) — so apt-get is avoided. The per-stage
# Dockerfiles are kept for portability/CI where a real Docker/Podman builder is available.
#
# Usage (after `module load Apptainer/latest`):
#   bash containers/build_all.sh            # build every stage (light -> heavy)
#   bash containers/build_all.sh protein    # build a single stage
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SIF_DIR="$REPO/containers/sif"
DEF_DIR="$REPO/containers/def"
mkdir -p "$SIF_DIR" "$DEF_DIR"

# Keep Apptainer's cache + build tmp on scratch (home is small; build temp is large).
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-/mnt/scratch/users/$USER/.apptainer/cache}"
export APPTAINER_TMPDIR="${APPTAINER_TMPDIR:-/mnt/scratch/users/$USER/.apptainer/tmp}"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"

command -v apptainer >/dev/null 2>&1 || { echo "apptainer not on PATH — run: module load Apptainer/latest"; exit 1; }

# stage_dir:source_yml:image_tag   (ordered light -> heavy; R/Bioconductor envs last)
STAGES=(
  "rna_pseudotime:pseudotime.yml:dde27-rna-pseudotime-0.1"
  "rna_velocity:vlocity.yml:dde27-rna-velocity-0.1"
  "rna_clustering:clustering.yml:dde27-rna-clustering-0.1"
  "composition:composition.yml:dde27-composition-0.1"
  "protein:surfaceprot.yml:dde27-protein-0.1"
  "rna_annotation:annotation.yml:dde27-rna-annotation-0.1"
  "rna_integration:integration.yml:dde27-rna-integration-0.1"
  "rna_de:DE.yml:dde27-rna-de-0.1"
  "rna_preprocessing:preprocessing.yml:dde27-rna-preprocessing-0.1"
)

want="${1:-all}"
built=0; failed=0
for entry in "${STAGES[@]}"; do
  dir="${entry%%:*}"; rest="${entry#*:}"; yml="${rest%%:*}"; tag="${rest##*:}"
  if [ "$want" != "all" ] && [ "$want" != "$dir" ]; then continue; fi

  def="$DEF_DIR/${dir}.def"
  cat > "$def" <<EOF
Bootstrap: docker
From: mambaorg/micromamba:1.5.8

%files
    $REPO/containers/dde27_envs/$yml /opt/env.yml

%post
    # NB: not /tmp — Apptainer bind-mounts the host /tmp over the container's during %post,
    # which would hide a yml copied there by %files.
    micromamba install -y -n base -f /opt/env.yml
    # procps -> \`ps\`, required by Nextflow's task metric collector; via conda to avoid apt/fakeroot.
    micromamba install -y -n base -c conda-forge procps-ng
    micromamba clean -a -y

%environment
    export PATH=/opt/conda/bin:\$PATH
EOF

  echo "==> [$(date +%H:%M:%S)] building $dir ($yml) -> $SIF_DIR/$tag.sif"
  if apptainer build --force "$SIF_DIR/$tag.sif" "$def"; then
    echo "    OK $tag.sif ($(du -h "$SIF_DIR/$tag.sif" | cut -f1))"
    built=$((built+1))
  else
    echo "    FAILED: $dir"
    failed=$((failed+1))
  fi
done

echo "Done. built=$built failed=$failed; .sif in $SIF_DIR"
[ "$failed" -eq 0 ]
