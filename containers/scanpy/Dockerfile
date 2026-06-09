# scanpy + scrublet image for the QC / reference-mapping (annotation) branch.
# Mirrors the DDE_23 scanpy stack. No single biocontainer bundles scanpy + scrublet, so build once:
#   docker build -t <registry>/scanpy-scrublet:1.11 containers/scanpy && push
#   apptainer build scanpy.sif docker-daemon://<registry>/scanpy-scrublet:1.11
# Then set --scanpy_container accordingly.

FROM python:3.11-slim

RUN apt-get update && apt-get install -y --no-install-recommends build-essential procps && \
    rm -rf /var/lib/apt/lists/* && \
    pip install --no-cache-dir \
        "scanpy==1.11.*" \
        scrublet \
        leidenalg \
        igraph \
        numpy scipy pandas h5py anndata \
        matplotlib \
        plotly kaleido
