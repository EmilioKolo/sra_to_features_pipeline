#!/bin/bash

# Exit on error and undefined variables
set -euo pipefail

# Function to log messages with timestamp
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Get the environment variables
source config.env

# Define BASE_DIR with realpath
REAL_BASE=$(realpath $BASE_DIR)

# Make sure the base directory exists
mkdir -p "$REAL_BASE"
# Create data, install, logs, bin and tmp directories
DATA_DIR="$REAL_BASE/data"
INSTALL_DIR="$REAL_BASE/install"
LOGS_DIR="$DATA_DIR/logs"
BIN_DIR="$DATA_DIR/bin"
TMP_DIR="$DATA_DIR/tmp"
mkdir -p "$DATA_DIR"
mkdir -p "$INSTALL_DIR"
mkdir -p "$LOGS_DIR"
mkdir -p "$BIN_DIR"
mkdir -p "$TMP_DIR"
# Add the bin directory to PATH
export PATH="$BIN_DIR:$PATH"

# Try downloading with either curl or wget (whichever exists)
download() {
    if command -v curl >/dev/null; then
        curl -L "$1" -o "$2"
    elif command -v wget >/dev/null; then
        wget "$1" -O "$2"
    else
        echo "ERROR: Need either curl or wget to download files."
        exit 1
    fi
}

# Check that python3 is installed
if ! command -v python3 &> /dev/null; then
    log "Python3 is not installed. Attempting to install Python3."
    download "$PYTHON_URL" "$TMP_DIR/python3.tar.gz"
    tar -xzf "$TMP_DIR/python3.tar.gz" -C "$INSTALL_DIR"
    ln -sf "$INSTALL_DIR/python/bin/python3" "$BIN_DIR/python3"
    if ! command -v python3 &> /dev/null; then
        log "Python3 installation failed. Please install Python3 manually."
        exit 1
    else
        log "Python3 installed successfully."
    fi
fi
# Check that python3 version is at least 3.9
PY_VER_MIN=$(python3 -c "import sys; print(sys.version_info.minor)")
if [[ "$PY_VER_MIN" -lt 9 ]]; then
    log "Python3 version is too old. Minimum required version is 3.9."
    download "$PYTHON_URL" "$TMP_DIR/python3.tar.gz"
    tar -xzf "$TMP_DIR/python3.tar.gz" -C "$INSTALL_DIR"
    ln -sf "$INSTALL_DIR/python/bin/python3" "$BIN_DIR/python3"
    python3 --version
fi
# Update PY_VER_MIN after installation
PY_VER_MIN=$(python3 -c "import sys; print(sys.version_info.minor)")
# Check that pip is installed
if ! command -v pip &> /dev/null; then
    log "pip is not installed. Attempting to install pip..."
    PIP_URL="https://bootstrap-pypa-io.ingress.us-east-2.psfhosted.computer/pip/zipapp/pip-25.1.1.pyz"
    download "$PIP_URL" "$TMP_DIR/pip.pyz"
    chmod +x "$TMP_DIR/pip.pyz"
    mv "$TMP_DIR/pip.pyz" "$INSTALL_DIR/pip"
    ln -sf "$INSTALL_DIR/pip" "$BIN_DIR/pip"
    python3 -m pip --version || {
        log "pip installation failed. Please install pip manually with the following command:\nsudo apt-get install python3-pip"
        exit 1
    }
    log "pip installed successfully. Attempting to upgrade pip..."
    python3 -m pip install --upgrade pip
fi
# Check that unzip is installed
if ! command -v unzip &> /dev/null; then
    log "unzip is not installed. Attempting to install unzip..."
    download "$BUSYBOX_URL" "$TMP_DIR/unzip"
    # Validate binary
    chmod +x "$TMP_DIR/unzip"
    file "$TMP_DIR/unzip" | grep -q 'statically linked' || {
    echo "Downloaded busybox is not static–check architecture or URL"; exit 1
    }
    # Move to install dir
    mv "$TMP_DIR/unzip" "$INSTALL_DIR/unzip"
    ln -sf "$INSTALL_DIR/unzip" "$BIN_DIR/unzip"
    if ! command -v unzip &> /dev/null; then
        log "unzip installation failed. Please install unzip manually with the following command:\nsudo apt-get install unzip"
        exit 1
    else
        log "unzip installed successfully."
    fi
fi
# Check that java is installed
if ! command -v java &> /dev/null; then
    log "Java is not installed. Attemptting to install Java..."
    download "$JDE_URL" "$TMP_DIR/java.tar.gz"
    mkdir -p "$INSTALL_DIR/java"
    tar -xzf "$TMP_DIR/java.tar.gz" -C "$INSTALL_DIR/java" --strip-components=1
    ln -s "$INSTALL_DIR/java/bin/java" "$BIN_DIR/java"
    if ! command -v java &> /dev/null; then
        log "Java installation failed. Please install Java manually with the following command:\nsudo apt-get install default-jre"
        exit 1
    else
        log "Java installed successfully."
    fi
fi

# Install required Python packages
log "Checking and installing required Python packages..."
python3 -m pip install numpy==1.26.4
python3 -m pip install pandas==2.1.4
python3 -m pip install requests
python3 -m pip install cnvpytor==1.3.1

# Install required packages manually

# Download and install samtools if not installed
if ! command -v samtools &> /dev/null; then
    log "samtools is not installed. Installing samtools..."
    download "$SAMTOOLS_URL" "$TMP_DIR/samtools.tar.gz"
    tar -xzf "$TMP_DIR/samtools.tar.gz" -C "$TMP_DIR"
    (cd "$TMP_DIR/samtools-$SAMTOOLS_VER" && ./configure --prefix="$INSTALL_DIR" && make && make install)
    ln -sf "$INSTALL_DIR/bin/samtools" "$BIN_DIR/samtools"
else
    log "samtools is already installed, skipping."
fi
# Download and install bcftools if not installed
if ! command -v bcftools &> /dev/null; then
    log "bcftools is not installed. Installing bcftools..."
    download "$BCFTOOLS_URL" "$TMP_DIR/bcftools.tar.gz"
    tar -xzf "$TMP_DIR/bcftools.tar.gz" -C "$TMP_DIR"
    (cd "$TMP_DIR/bcftools-$BCFTOOLS_VER" && ./configure --prefix="$INSTALL_DIR" && make && make install)
    ln -sf "$INSTALL_DIR/bin/bcftools" "$BIN_DIR/bcftools"
else
    log "bcftools is already installed, skipping."
fi
# Download and install bwa if not installed
if ! command -v bwa &> /dev/null; then
    log "bwa is not installed. Installing bwa..."
    download "$BWA_URL" "$TMP_DIR/bwa.tar.gz"
    tar -xzf "$TMP_DIR/bwa.tar.gz" -C "$TMP_DIR"
    (cd "$TMP_DIR/bwa-$BWA_VER" && make)
    cp "$TMP_DIR/bwa-$BWA_VER/bwa" "$BIN_DIR/bwa"
else
    log "bwa is already installed, skipping."
fi
# Download and install sra-toolkit if not installed
if ! command -v fastq-dump &> /dev/null; then
    log "sra-toolkit is not installed. Installing sra-toolkit..."
    download "$SRA_URL" "$TMP_DIR/sratoolkit.tar.gz"
    tar -xzf "$TMP_DIR/sratoolkit.tar.gz" -C "$TMP_DIR"
    cp -r "$TMP_DIR/sratoolkit.$SRA_VER-$SRA_ARCH/bin/"* "$BIN_DIR"
else
    log "sra-toolkit is already installed, skipping."
fi
# Download and install tabix if not installed
if ! command -v tabix &> /dev/null; then
    log "tabix is not installed. Installing tabix..."
    download "$HTSLIB_URL" "$TMP_DIR/htslib.tar.bz2"
    tar -xjf "$TMP_DIR/htslib.tar.bz2" -C "$TMP_DIR"
    (cd "$TMP_DIR/htslib-$HTSLIB_VER" && ./configure --prefix="$INSTALL_DIR" && make && make install)
    ln -sf "$INSTALL_DIR/bin/tabix" "$BIN_DIR/tabix"
else
    log "tabix is already installed, skipping."
fi
# Download and install bedtools if not installed
if ! command -v bedtools &> /dev/null; then
    log "bedtools is not installed. Installing bedtools..."
    download "$BEDTOOLS_URL" "$TMP_DIR/bedtools.tar.gz"
    tar -xzf "$TMP_DIR/bedtools.tar.gz" -C "$TMP_DIR"
    (cd "$TMP_DIR/bedtools2" && make)
    cp "$TMP_DIR/bedtools2/bin/"* "$BIN_DIR"
else
    log "bedtools is already installed, skipping."
fi
# Download and install Kraken2 if not installed
if ! command -v kraken2 &> /dev/null; then
    log "Kraken2 is not installed. Installing Kraken2..."
    download "$KRAKEN2_URL" "$TMP_DIR/kraken2.tar.gz"
    tar -xzf "$TMP_DIR/kraken2.tar.gz" -C "$TMP_DIR"
    (cd "$TMP_DIR/kraken2-$KRAKEN2_VER" && ./install_kraken2.sh "$INSTALL_DIR/kraken2")
    ln -sf "$INSTALL_DIR/kraken2/kraken2" "$BIN_DIR/kraken2"
    ln -sf "$INSTALL_DIR/kraken2/kraken2-build" "$BIN_DIR/kraken2-build"
else
    log "Kraken2 is already installed, skipping."
fi

# Attempt to install pybedtools
python3 -m pip install pybedtools || {
    log "pybedtools installation failed. Please install pybedtools manually with the following command:\npip install pybedtools"
    exit 1
}

# Cleanup
rm -rf "$TMP_DIR/*"
echo "Installation complete!"

# Add to PATH
echo "Adding binaries to PATH temporarily. To make this permanent, add the following lines to your ~/.bashrc or ~/.bash_profile:"
echo "export PATH=\"$BIN_DIR:\$PATH\""
echo "export KRAKEN2_DB_PATH=\"$INSTALL_DIR/kraken2_db\" # Recommended for Kraken2"
# Add kraken2_db to kraken2 paths
export KRAKEN2_DB_PATH="$INSTALL_DIR/kraken2_db"

# Download reference genome if it does not exist
if [[ ! -f "$DATA_DIR/reference.fasta" ]]; then
    log "Downloading reference genome..."
    download "$FASTA_URL" "$DATA_DIR/reference.fasta.gz"
    gunzip -f "$DATA_DIR/reference.fasta.gz"
else
    log "Reference genome already exists, skipping download."
fi
if [[ ! -f "$DATA_DIR/reference.gff" ]]; then
    log "Downloading reference gff..."
    download "$GFF_URL" "$DATA_DIR/reference.gff.gz"
    gunzip -f "$DATA_DIR/reference.gff.gz"
else
    log "Reference gff already exists, skipping download."
fi

# Install snpEff if it does not exist
snpeff_dir="$BIN_DIR"
genome_name="$GENOME_NAME"
snpeff_jar="$snpeff_dir/snpEff/snpEff.jar"

if [[ ! -f "$snpeff_jar" ]]; then
    log "Downloading and installing snpEff..."
    mkdir -p "$snpeff_dir"
    download "$SNPEFF_URL" "$snpeff_dir/snpEff.zip"
    unzip -o "$snpeff_dir/snpEff.zip" -d "$snpeff_dir"
    rm -f "$snpeff_dir/snpEff.zip"
else
    log "snpEff already installed, skipping download."
fi

# Create custom genome if it does not exist
if [[ ! -d "$snpeff_dir/snpEff/data/$genome_name" ]]; then
    log "Creating snpEff custom genome..."
    mkdir -p "$snpeff_dir/snpEff/data/$genome_name"
    cp "$DATA_DIR/reference.fasta" "$snpeff_dir/snpEff/data/$genome_name/sequences.fa"
    cp "$DATA_DIR/reference.gff" "$snpeff_dir/snpEff/data/$genome_name/genes.gff"
    echo "${genome_name}.genome : Custom genome" >> "$snpeff_dir/snpEff/snpEff.config"
    java -Xmx4g -jar "$snpeff_jar" build -gff3 -v "$genome_name" > "$LOGS_DIR/snpeff.log" 2>&1
else
    log "snpEff custom genome already exists, skipping creation."
fi

# Index reference genome if not indexed
if [[ ! -f "$DATA_DIR/reference.fasta.bwt" ]]; then
    log "Indexing reference genome with bwa..."
    bwa index "$DATA_DIR/reference.fasta" > "$LOGS_DIR/bwa_index.log" 2>&1
else
    log "Reference genome already indexed, skipping."
fi

# Run install.py if output does not exist
if [[ ! -f "$DATA_DIR/genome.sizes" ]]; then
    log "Running install.py..."
    python3 install.py "$DATA_DIR" > "$LOGS_DIR/install_py.log" 2>&1
else
    log "install.py was already run, skipping."
fi

# Download cnvpytor data if it does not exist
if [[ ! -d "$BIN_DIR/.cnvpytor/data" ]]; then
    log "Downloading cnvpytor data..."
    cnvpytor -download > "$LOGS_DIR/cnvpytor.log" 2>&1
else
    log "CNVpytor data already exists, skipping download."
fi

# Download Kraken2 database if it does not exist
kraken_db="$INSTALL_DIR/kraken2-db"
if [[ ! -f "$kraken_db/hash.k2d" ]]; then
    log "Downloading Kraken2 database..."
    mkdir -p "$kraken_db"
    download "$KRAKEN2_DB_URL" "$kraken_db/k2.tar.gz"
    tar -xvzf "$kraken_db/k2.tar.gz" -C "$kraken_db"
    rm -f "$kraken_db/k2.tar.gz"
else
    log "Kraken2 database already exists, skipping download."
fi

# Run the main script
log "Starting main script execution..."
python3 2_run_multiple.py
