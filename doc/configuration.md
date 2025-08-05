# Configuration Guide

This guide explains how to configure the SRA to Features Pipeline for optimal performance and customization.

## 🔧 Configuration Methods

The pipeline supports multiple configuration methods, listed in order of precedence:

1. **Command-line arguments** (highest priority)
2. **Environment variables**
3. **Configuration file** (`.env`)
4. **Default values** (lowest priority)

## 📝 Configuration File

Create a `.env` file in your working directory to set configuration options:

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

## 🌍 Environment Variables

Set environment variables in your shell profile (`~/.bashrc`, `~/.zshrc`, etc.):

```bash
# Base configuration
export SRA_PIPELINE_BASE_DIR="/path/to/data"
export SRA_PIPELINE_OUTPUT_DIR="/path/to/output"
export SRA_PIPELINE_THREADS=8

# Reference files
export SRA_PIPELINE_REFERENCE_FASTA="/path/to/reference.fasta"
export SRA_PIPELINE_REFERENCE_GFF="/path/to/reference.gff"

# Quality parameters
export SRA_PIPELINE_MIN_QUALITY_SCORE=20
export SRA_PIPELINE_MIN_COVERAGE=10

# Logging
export SRA_PIPELINE_LOG_LEVEL="INFO"
```

## ⚙️ Configuration Options

### Base Paths

| Option | Default | Description |
|--------|---------|-------------|
| `BASE_DIR` | `./data` | Base directory for all data files |
| `OUTPUT_DIR` | `./output` | Output directory for results |
| `TMP_DIR` | `./tmp` | Temporary directory for intermediate files |

### Reference Files

| Option | Required | Description |
|--------|----------|-------------|
| `REFERENCE_FASTA` | Yes | Reference genome FASTA file |
| `REFERENCE_GFF` | Yes | Reference genome GFF annotation file |
| `BED_GENES` | Yes | BED file with gene coordinates |
| `GENOME_SIZES` | Yes | File with chromosome sizes |

### External Tools

| Option | Required | Description |
|--------|----------|-------------|
| `KRAKEN_DB` | No | Kraken2 database directory |
| `SNPEFF_DIR` | No | snpEff installation directory |

### Analysis Parameters

| Option | Default | Description |
|--------|---------|-------------|
| `THREADS` | 1 | Number of threads for parallel processing |
| `BIN_SIZE_GVS` | 100000 | Bin size for genomic variant statistics |
| `BIN_SIZE_CNV` | 100000 | Bin size for CNV analysis |

### Quality Control

| Option | Default | Description |
|--------|---------|-------------|
| `MIN_QUALITY_SCORE` | 20 | Minimum quality score for variants |
| `MIN_COVERAGE` | 10 | Minimum coverage for variant calling |

### Performance

| Option | Default | Description |
|--------|---------|-------------|
| `MAX_MEMORY_GB` | 16 | Maximum memory usage in GB |
| `TIMEOUT_SECONDS` | 3600 | Timeout for individual operations |

### Logging

| Option | Default | Description |
|--------|---------|-------------|
| `LOG_LEVEL` | INFO | Logging level (DEBUG, INFO, WARNING, ERROR) |
| `LOG_FILE` | None | Log file path (optional) |

## 🎯 Configuration Examples

### Basic Configuration

```ini
# Basic setup for small datasets
SRA_PIPELINE_THREADS=4
SRA_PIPELINE_MIN_QUALITY_SCORE=20
SRA_PIPELINE_MIN_COVERAGE=10
SRA_PIPELINE_LOG_LEVEL=INFO
```

### High-Performance Configuration

```ini
# Optimized for large datasets
SRA_PIPELINE_THREADS=16
SRA_PIPELINE_MAX_MEMORY_GB=32
SRA_PIPELINE_TIMEOUT_SECONDS=7200
SRA_PIPELINE_MIN_QUALITY_SCORE=30
SRA_PIPELINE_MIN_COVERAGE=20
```

### Development Configuration

```ini
# For development and debugging
SRA_PIPELINE_LOG_LEVEL=DEBUG
SRA_PIPELINE_THREADS=2
SRA_PIPELINE_TIMEOUT_SECONDS=1800
```

### Production Configuration

```ini
# For production environments
SRA_PIPELINE_LOG_LEVEL=WARNING
SRA_PIPELINE_THREADS=8
SRA_PIPELINE_MAX_MEMORY_GB=16
SRA_PIPELINE_MIN_QUALITY_SCORE=25
SRA_PIPELINE_MIN_COVERAGE=15
```

## 🔍 Configuration Validation

Validate your configuration:

```bash
# Validate configuration
sra-pipeline validate

# Validate with custom config file
sra-pipeline validate --config my_config.env
```

## 📊 Performance Tuning

### Memory Optimization

```ini
# For memory-constrained systems
SRA_PIPELINE_MAX_MEMORY_GB=8
SRA_PIPELINE_THREADS=4
```

### Speed Optimization

```ini
# For speed-optimized runs
SRA_PIPELINE_THREADS=16
SRA_PIPELINE_MAX_MEMORY_GB=32
SRA_PIPELINE_TIMEOUT_SECONDS=7200
```

### Quality Optimization

```ini
# For high-quality results
SRA_PIPELINE_MIN_QUALITY_SCORE=30
SRA_PIPELINE_MIN_COVERAGE=20
```

## 🐛 Troubleshooting Configuration

### Common Issues

**Configuration not found:**
```bash
# Check if .env file exists
ls -la .env

# Create default configuration
sra-pipeline validate --create-config
```

**Invalid paths:**
```bash
# Check if paths exist
ls -la /path/to/reference.fasta

# Use absolute paths
SRA_PIPELINE_REFERENCE_FASTA=/absolute/path/to/reference.fasta
```

**Permission errors:**
```bash
# Check file permissions
ls -la /path/to/output

# Fix permissions
chmod 755 /path/to/output
```

### Configuration Debugging

```bash
# Show current configuration
sra-pipeline validate --verbose

# Test configuration without running
sra-pipeline run --sra-id SRR123456 --dry-run
```

## 📚 Advanced Configuration

### Custom Reference Genomes

```ini
# Human genome (GRCh38)
SRA_PIPELINE_REFERENCE_FASTA=/path/to/GRCh38.fasta
SRA_PIPELINE_REFERENCE_GFF=/path/to/GRCh38.gff
SRA_PIPELINE_BED_GENES=/path/to/GRCh38_genes.bed
SRA_PIPELINE_GENOME_SIZES=/path/to/GRCh38.sizes

# Mouse genome (GRCm39)
SRA_PIPELINE_REFERENCE_FASTA=/path/to/GRCm39.fasta
SRA_PIPELINE_REFERENCE_GFF=/path/to/GRCm39.gff
SRA_PIPELINE_BED_GENES=/path/to/GRCm39_genes.bed
SRA_PIPELINE_GENOME_SIZES=/path/to/GRCm39.sizes
```

### Batch Processing Configuration

```ini
# Optimized for batch processing
SRA_PIPELINE_THREADS=8
SRA_PIPELINE_MAX_MEMORY_GB=16
SRA_PIPELINE_TIMEOUT_SECONDS=3600
SRA_PIPELINE_LOG_LEVEL=INFO
```

## 🔄 Configuration Migration

### From Old Version (v0.9.0)

Old `config.ini` format:
```ini
[Paths]
BASE_DIR = /content

[Parameters]
THREADS = 8
```

New `.env` format:
```ini
SRA_PIPELINE_BASE_DIR=/content
SRA_PIPELINE_THREADS=8
```

## 📖 Next Steps

After configuring the pipeline:

1. **Validate your configuration**: `sra-pipeline validate`
2. **Run a test analysis**: See [Quick Start Guide](quick_start.md)
3. **Explore advanced options**: See [User Guide](user_guide/)
4. **Optimize performance**: Monitor resource usage and adjust settings

---

**Need help?** Check our [Support Guide](support.md) or create an issue on GitHub. 