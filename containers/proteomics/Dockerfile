# Python container for the bulk-proteomics stages (PROT_MS_* except PROT_MS_DESP).
# Build locally / push to a registry, then point params.proteomics_container at it. On Viking the
# viking profile instead uses the `aml_proteomics` conda env (envs/proteomics.yml) and ignores this.
FROM mambaorg/micromamba:1.5.8

COPY --chown=$MAMBA_USER:$MAMBA_USER envs/proteomics.yml /tmp/proteomics.yml
RUN micromamba install -y -n base -f /tmp/proteomics.yml && micromamba clean --all --yes
ENV PATH=/opt/conda/bin:$PATH
