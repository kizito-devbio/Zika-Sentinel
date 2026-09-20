# Zika-Sentinel v1.1 — reproducible viral tool stack
FROM ubuntu:22.04

LABEL maintainer="Kizito Sylvester-Ali Ibeojo"
LABEL description="Zika-Sentinel ZIKV genomic surveillance"
LABEL version="1.1.0"

ENV DEBIAN_FRONTEND=noninteractive
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

WORKDIR /workspace

# System dependencies and viral genomics tools.
# MAFFT is installed separately because Ubuntu 22.04 does not
# currently provide the required MAFFT 7.505 package version.
RUN apt-get update && apt-get install -y --no-install-recommends \
        wget \
        curl \
        git \
        unzip \
        zip \
        build-essential \
        autoconf \
        automake \
        libtool \
        ca-certificates \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libncurses5-dev \
        libhts-dev \
        libcurl4-gnutls-dev \
        libssl-dev \
        python3 \
        python3-pip \
        python3-venv \
        iqtree=2.0.7+dfsg-1 \
        fastqc \
        fastp \
        samtools \
        bedtools \
        bwa \
        minimap2 \
        bcftools \
        graphviz \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install MAFFT 7.505 from the official MAFFT source release.
ARG MAFFT_VERSION=7.505

RUN set -eux; \
    wget -q \
        "https://mafft.cbrc.jp/alignment/software/mafft-${MAFFT_VERSION}-with-extensions-src.tgz" \
        -O /tmp/mafft.tgz; \
    mkdir -p /tmp/mafft-src; \
    tar -xzf /tmp/mafft.tgz \
        -C /tmp/mafft-src \
        --strip-components=1; \
    cd /tmp/mafft-src/core; \
    make clean || true; \
    make -j"$(nproc)"; \
    make install PREFIX=/usr/local; \
    cd /; \
    rm -rf /tmp/mafft.tgz /tmp/mafft-src; \
    mafft --version

# Python environment for reporting, visualization and sequence analysis.
RUN python3 -m venv /opt/venv

ENV PATH="/opt/venv/bin:$PATH"

RUN pip install --no-cache-dir \
        "pandas==2.2.2" \
        "numpy==1.26.4" \
        "matplotlib==3.8.4" \
        "seaborn==0.13.2" \
        "biopython==1.83" \
        "requests==2.32.3" \
        "pysam==0.22.1"

# Install iVar 1.4.2 from the tagged upstream release.
ARG IVAR_VERSION=1.4.2

RUN set -eux; \
    git clone \
        --depth 1 \
        --branch "v${IVAR_VERSION}" \
        https://github.com/andersen-lab/ivar.git \
        /tmp/ivar; \
    cd /tmp/ivar; \
    ./autogen.sh; \
    ./configure --prefix=/usr/local; \
    make -j"$(nproc)"; \
    make install; \
    cd /; \
    rm -rf /tmp/ivar; \
    ivar version

# Record installed tool versions inside the image.
RUN set -eux; \
    { \
      echo "Zika-Sentinel: 1.1.0"; \
      echo "fastp: $(fastp --version 2>&1 | head -1)"; \
      echo "FastQC: $(fastqc --version 2>&1 | head -1)"; \
      echo "samtools: $(samtools --version 2>&1 | head -1)"; \
      echo "bcftools: $(bcftools --version 2>&1 | head -1)"; \
      echo "bwa: $(bwa 2>&1 | head -3 | tr '\n' ' ')"; \
      echo "minimap2: $(minimap2 --version 2>&1 | head -1)"; \
      echo "bedtools: $(bedtools --version 2>&1 | head -1)"; \
      echo "mafft: $(mafft --version 2>&1 | head -1)"; \
      echo "iqtree: $(iqtree2 --version 2>&1 | head -1 || iqtree --version 2>&1 | head -1)"; \
      echo "ivar: $(ivar version 2>&1 | head -1)"; \
      echo "python: $(python3 --version 2>&1)"; \
      echo "biopython: $(python3 -c 'import Bio; print(Bio.__version__)')"; \
      echo "pandas: $(python3 -c 'import pandas; print(pandas.__version__)')"; \
      echo "numpy: $(python3 -c 'import numpy; print(numpy.__version__)')"; \
      echo "matplotlib: $(python3 -c 'import matplotlib; print(matplotlib.__version__)')"; \
      echo "pysam: $(python3 -c 'import pysam; print(pysam.__version__)')"; \
    } > /opt/zika-sentinel-tool-versions.txt; \
    cat /opt/zika-sentinel-tool-versions.txt

# Zika-Sentinel runtime scripts.
# Used by the public FASTA pathway to retrieve normalized NCBI metadata.
COPY bin/fetch_ncbi_metadata.py /opt/zika-sentinel/bin/fetch_ncbi_metadata.py

RUN chmod +x /opt/zika-sentinel/bin/fetch_ncbi_metadata.py

WORKDIR /workspace



