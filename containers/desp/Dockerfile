# R container carrying the DESP package for the single demix step (PROT_MS_DESP).
# DESP is not packaged on conda/CRAN — install from GitHub at build time. Point
# params.proteomics_r_container at the built image; the viking profile uses the `desp_r` conda env.
FROM rocker/r-ver:4.3.2

RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev libssl-dev libxml2-dev && rm -rf /var/lib/apt/lists/*
RUN R -e "install.packages(c('optparse','remotes','quadprog'), repos='https://cloud.r-project.org')" \
 && R -e "remotes::install_github('AhmedYoussef95/DESP')"   # Cell Reports Methods 2024 DESP repo
