# Installation Guide

This guide covers all the steps needed to install and set up the SRA to Features Pipeline.

## 📋 Prerequisites

### System Requirements

- **Operating System**: Linux, macOS, or Windows (with WSL)
- **Python**: 3.8 or higher
- **Memory**: 8GB RAM minimum (16GB+ recommended)
- **Storage**: 50GB+ disk space for reference genomes and intermediate files
- **Network**: Internet connection for downloading SRA data and tools

### External Bioinformatics Tools

The pipeline requires the following bioinformatics tools to be installed and available in your system PATH:

| Tool | Version | Purpose |
|------|---------|---------|
| **BWA** | v0.7.17+ | Sequence alignment |
| **SAMtools** | v1.10+ | SAM/BAM file manipulation |
| **BCFtools** | v1.10+ | Variant calling |
| **BEDtools** | v2.29+ | Genome arithmetic |
| **FastQC** | v0.11+ | Quality control |
| **fastq-dump** (SRA Toolkit) | v2.10+ | SRA data download |
| **tabix/bgzip** (HTSlib) | v1.10+ | VCF compression and indexing |
| **Java** | v8+ | For snpEff |
| **CNVpytor** | v1.3+ | Copy number variation analysis |

## 🚀 Installation Methods

### Method 1: Automated Installation (Recommended)

We provide an automated installation script that handles most of the setup:

```bash
# Clone the repository
git clone https://github.com/your-org/sra-to-features-pipeline.git
cd sra-to-features-pipeline

# Run the automated installation script
./scripts/install/install_tools.sh

# Install the Python package
pip install -e .
```

### Method 2: Manual Installation

If you prefer to install tools manually or need custom configurations:

#### Step 1: Install Bioinformatics Tools

**Ubuntu/Debian:**
```bash
sudo apt update
sudo apt install -y bwa samtools bcftools bedtools fastqc default-jre
```

**CentOS/RHEL/Fedora:**
```bash
sudo yum install -y bwa samtools bcftools bedtools fastqc java-1.8.0-openjdk
# or for newer versions:
sudo dnf install -y bwa samtools bcftools bedtools fastqc java-1.8.0-openjdk
```

**macOS:**
```bash
# Install Homebrew if not already installed
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install tools
brew install bwa samtools bcftools bedtools fastqc openjdk@8
```

#### Step 2: Install SRA Toolkit

```bash
# Download and install SRA Toolkit
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -xzf sratoolkit.current-ubuntu64.tar.gz
sudo cp sratoolkit.*/bin/* /usr/local/bin/
```

#### Step 3: Install CNVpytor

```bash
pip install cnvpytor
```

#### Step 4: Install the Pipeline

```bash
# Clone the repository
git clone https://github.com/your-org/sra-to-features-pipeline.git
cd sra-to-features-pipeline

# Install in development mode
pip install -e ".[dev]"
```

### Method 3: Using Conda

If you prefer using Conda for environment management:

```bash
# Create a new conda environment
conda create -n sra-pipeline python=3.9
conda activate sra-pipeline

# Install bioinformatics tools
conda install -c bioconda bwa samtools bcftools bedtools fastqc sra-tools htslib

# Install the pipeline
pip install sra-to-features-pipeline
```

### Method 4: Docker (Legacy)

For users who prefer Docker (not recommended for production):

```bash
# Build the Docker image
docker build -t sra-pipeline .

# Run the pipeline
docker run -v /path/to/output:/output sra-pipeline --sra-id SRR123456
```

## ✅ Verification

After installation, verify that everything is working correctly:

```bash
# Check if all tools are available
which bwa samtools bcftools bedtools fastqc fastq-dump

# Validate the pipeline setup
sra-pipeline validate

# Run a quick test
sra-pipeline run --sra-id SRR123456 --dry-run
```

## 🔧 Configuration

### Environment Variables

Set up environment variables for the pipeline:

```bash
# Add to your ~/.bashrc or ~/.zshrc
export SRA_PIPELINE_BASE_DIR="/path/to/data"
export SRA_PIPELINE_OUTPUT_DIR="/path/to/output"
export SRA_PIPELINE_THREADS=8
export SRA_PIPELINE_LOG_LEVEL="INFO"
```

### Configuration File

Create a `.env` file in your working directory:

```ini
# Base paths
SRA_PIPELINE_BASE_DIR=/path/to/data
SRA_PIPELINE_OUTPUT_DIR=/path/to/output

# Reference files
SRA_PIPELINE_REFERENCE_FASTA=/path/to/reference.fasta
SRA_PIPELINE_REFERENCE_GFF=/path/to/reference.gff
SRA_PIPELINE_BED_GENES=/path/to/genes.bed
SRA_PIPELINE_GENOME_SIZES=/path/to/genome.sizes

# External tools
SRA_PIPELINE_KRAKEN_DB=/path/to/kraken/db
SRA_PIPELINE_SNPEFF_DIR=/path/to/snpeff

# Parameters
SRA_PIPELINE_THREADS=8
SRA_PIPELINE_BIN_SIZE_GVS=100000
SRA_PIPELINE_BIN_SIZE_CNV=100000

# Quality control
SRA_PIPELINE_MIN_QUALITY_SCORE=20
SRA_PIPELINE_MIN_COVERAGE=10

# Performance
SRA_PIPELINE_MAX_MEMORY_GB=16
SRA_PIPELINE_TIMEOUT_SECONDS=3600

# Logging
SRA_PIPELINE_LOG_LEVEL=INFO
```

## 🐛 Troubleshooting

### Common Installation Issues

**Permission Errors:**
```bash
# Install with user permissions
pip install --user sra-to-features-pipeline

# Or use a virtual environment
python -m venv venv
source venv/bin/activate
pip install sra-to-features-pipeline
```

**Missing Tools:**
```bash
# Check if tools are in PATH
echo $PATH
which bwa

# Add tools to PATH if needed
export PATH="/path/to/tools:$PATH"
```

**Python Version Issues:**
```bash
# Check Python version
python --version

# Use Python 3 explicitly
python3 -m pip install sra-to-features-pipeline
```

### Getting Help

If you encounter issues during installation:

1. Check the [Troubleshooting Guide](troubleshooting.md)
2. Search existing [GitHub Issues](https://github.com/your-org/sra-to-features-pipeline/issues)
3. Create a new issue with detailed error information

## 📚 Next Steps

After successful installation:

1. Read the **[Quick Start Guide](quick_start.md)** to run your first analysis
2. Explore the **[Configuration Guide](configuration.md)** for advanced setup
3. Check the **[User Guide](user_guide/)** for detailed usage instructions

---

**Need help?** Check our [Support Guide](support.md) or create an issue on GitHub. 