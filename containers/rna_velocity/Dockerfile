# dde27-rna-velocity-0.1 — built from ../../vlocity.yml (the stage's pinned conda env).
# Single source of truth: the env file at the repo root. Build context = repo root.
#   docker build -f containers/rna_velocity/Dockerfile -t dde27-rna-velocity-0.1 .
# See containers/README.md for the Apptainer (.sif) build used on Viking.
FROM mambaorg/micromamba:1.5.8

USER root
# procps -> Nextflow's task metric launcher needs `ps`; ca-certificates for conda channels.
RUN apt-get update && apt-get install -y --no-install-recommends procps ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY containers/dde27_envs/vlocity.yml /tmp/env.yml
RUN micromamba install -y -n base -f /tmp/env.yml && micromamba clean -a -y

# Put the base env on PATH so the bin/*.py shebangs (#!/usr/bin/env python3) resolve to it.
ENV PATH=/opt/conda/bin:$PATH
ENV MAMBA_DOCKERFILE_ACTIVATE=1
