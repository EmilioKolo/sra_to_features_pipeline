# Use Ubuntu as base
FROM ubuntu:22.04

# Set non-interactive mode for apt
ENV DEBIAN_FRONTEND=noninteractive

# Install base tools
RUN apt-get update && \
    apt-get install -y \
        build-essential \
        ca-certificates \
        curl \
        default-jre \
        git \
        libx11-dev \
        libgl1 \
        libhdf5-dev \
        python3 \
        python3-pip \
        wget \
        zlib1g-dev \
        && apt-get clean

# Install specific tools
RUN apt-get update && \
    apt-get install -y \
        bwa \
        bcftools \
        bedtools \
        kraken2 \
        samtools \
        sra-toolkit \
        tabix \
        unzip \
        && apt-get clean

# Install pip packages
RUN pip3 install --no-cache-dir \
    numpy==1.26.4 \
    pandas==2.1.4 \
    requests \
    cnvpytor==1.3.1

# Create working directory for pipeline
WORKDIR /pipeline

# Copy files needed for setup
COPY amplicons.bed .
COPY config.ini .
COPY gene_regions.bed .
COPY regions.bed .
COPY setup.sh .
COPY scripts/get_config_value.py .

# Make setup script executable
RUN chmod +x setup.sh
# Run setup script at build time
RUN ./setup.sh

# Copy install.py file
COPY install.py .
# Run python install script
RUN python3 ./install.py

# Copy basic scripts
COPY scripts/ ./scripts/
COPY __init__.py .
COPY main.py .

# Run the Python script at container start
ENTRYPOINT ["python3", "pipeline.py"]
