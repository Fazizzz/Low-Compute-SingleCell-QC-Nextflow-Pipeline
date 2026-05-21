FROM mambaorg/micromamba:1.5

# procps provides `ps` which Nextflow uses for task metrics; gzip is needed for some
# downstream tools (kb-python tries to gunzip ondisk fastqs).
USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
        procps gzip ca-certificates \
    && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

COPY envs/scq.yml /tmp/scq.yml
RUN micromamba install -y -n base -f /tmp/scq.yml && micromamba clean --all --yes

ENV PATH=/opt/conda/bin:$PATH
ENV MAMBA_DOCKERFILE_ACTIVATE=1

WORKDIR /work
