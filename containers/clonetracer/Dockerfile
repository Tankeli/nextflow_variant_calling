# CloneTracer (veltenlab) Bayesian clonal-inference image — runs bin/run_clonetracer.py.
# Built from the upstream conda spec (envs/clonetracer.yml: pyro-ppl 1.8.4 + pytorch 1.13).
# This is the CPU image (the default path); set clonetracer_gpu=true + a GPU partition for the
# CUDA path. Build once and point --clonetracer_container at it:
#   docker build -t <registry>/clonetracer:1.0 containers/clonetracer && push
#   apptainer build clonetracer.sif docker-daemon://<registry>/clonetracer:1.0
# On Viking the `-profile viking` path instead uses a conda env named `clonetracer`
# (conda env create -f envs/clonetracer.yml -n clonetracer, on the LOGIN node).

FROM python:3.11-slim

RUN apt-get update && apt-get install -y --no-install-recommends build-essential procps && \
    rm -rf /var/lib/apt/lists/* && \
    pip install --no-cache-dir \
        "pyro-ppl==1.8.4" \
        "torch==1.13.1" \
        numpy scipy pandas \
        seaborn matplotlib \
        dill scanpy anndata
