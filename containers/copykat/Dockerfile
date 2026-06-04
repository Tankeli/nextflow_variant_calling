# CopyKAT has no public biocontainer, so build this once and push to a registry your
# cluster can pull (or build an apptainer .sif). Mirrors DDE_32 scripts/setup_copykat.sh.
#
# Build:    docker build -t <your-registry>/copykat:1.1.0 containers/copykat
# Apptainer: apptainer build copykat.sif docker-daemon://<your-registry>/copykat:1.1.0
# Then set --copykat_container accordingly.

FROM mambaorg/micromamba:1.5.8

USER root
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        r-base=4.3 \
        r-matrix \
        r-rcpp \
        r-data.table \
        r-kernlab \
        r-mcmcpack \
        r-mixtools \
        r-paralleldist \
        r-igraph \
        r-remotes && \
    micromamba clean -a -y

ENV PATH=/opt/conda/bin:$PATH
RUN Rscript -e 'remotes::install_github("navinlabcode/copykat", dependencies = FALSE)' && \
    Rscript -e 'library(copykat); cat("CopyKAT OK\n")'
