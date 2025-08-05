# Command Line Interface Reference

This guide provides a complete reference for the SRA to Features Pipeline command-line interface, including all commands, options, and examples.

## 🚀 Overview

The pipeline provides a modern command-line interface with the following main commands:

```bash
sra-pipeline [OPTIONS] COMMAND [ARGS]...
```

## 📋 Main Commands

### `run` - Execute Pipeline

Run the complete SRA to Features Pipeline analysis.

```bash
sra-pipeline run [OPTIONS]
```

#### Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--sra-id` | string | None | SRA Run Accession ID (e.g., SRR123456) |
| `--fastq` | path | None | FASTQ file paths (can specify multiple) |
| `--output-dir` | path | `./output` | Output directory for results |
| `--config` | path | None | Configuration file path |
| `--threads` | integer | 1 | Number of threads to use |
| `--log-level` | choice | `INFO` | Logging level (DEBUG, INFO, WARNING, ERROR) |
| `--log-file` | path | None | Log file path |
| `--dry-run` | flag | False | Perform a dry run without executing |

#### Examples

**Run with SRA ID:**
```bash
sra-pipeline run --sra-id SRR123456 --output-dir ./results --threads 8
```

**Run with FASTQ files:**
```bash
sra-pipeline run \
  --fastq sample_R1.fastq.gz sample_R2.fastq.gz \
  --output-dir ./results \
  --threads 8
```

**Run with custom configuration:**
```bash
sra-pipeline run \
  --sra-id SRR123456 \
  --config my_config.env \
  --output-dir ./results
```

**Dry run (test configuration):**
```bash
sra-pipeline run --sra-id SRR123456 --dry-run
```

### `validate` - Validate Setup

Validate that your system is properly configured for running the pipeline.

```bash
sra-pipeline validate [OPTIONS]
```

#### Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--config` | path | None | Configuration file path |
| `--verbose` | flag | False | Verbose output |

#### Examples

**Basic validation:**
```bash
sra-pipeline validate
```

**Validate with custom config:**
```bash
sra-pipeline validate --config my_config.env
```

**Verbose validation:**
```bash
sra-pipeline validate --verbose
```

### `setup-docs` - Generate Documentation

Generate API documentation for the pipeline.

```bash
sra-pipeline setup-docs [OPTIONS]
```

#### Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--output-dir` | path | `./doc` | Output directory for documentation |

#### Examples

**Generate documentation:**
```bash
sra-pipeline setup-docs
```

**Generate to custom directory:**
```bash
sra-pipeline setup-docs --output-dir ./my_docs
```

## 🔧 Global Options

These options can be used with any command:

| Option | Description |
|--------|-------------|
| `--help` | Show help message and exit |
| `--version` | Show version and exit |
| `--verbose` | Enable verbose output |
| `--quiet` | Suppress output |

## 📊 Configuration Options

### Environment Variables

The pipeline can be configured using environment variables with the `SRA_PIPELINE_` prefix:

```bash
# Set environment variables
export SRA_PIPELINE_THREADS=8
export SRA_PIPELINE_OUTPUT_DIR="/path/to/output"
export SRA_PIPELINE_LOG_LEVEL="INFO"
export SRA_PIPELINE_MIN_QUALITY_SCORE=20
export SRA_PIPELINE_MIN_COVERAGE=10
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

## 🎯 Common Usage Patterns

### Basic Analysis

```bash
# Simple SRA analysis
sra-pipeline run --sra-id SRR123456

# With custom output directory
sra-pipeline run --sra-id SRR123456 --output-dir ./my_results

# With multiple threads
sra-pipeline run --sra-id SRR123456 --threads 8
```

### Advanced Analysis

```bash
# With custom configuration
sra-pipeline run \
  --sra-id SRR123456 \
  --config production.env \
  --output-dir ./results \
  --threads 16 \
  --log-level DEBUG

# With FASTQ files
sra-pipeline run \
  --fastq sample_R1.fastq.gz sample_R2.fastq.gz \
  --output-dir ./results \
  --threads 8
```

### Batch Processing

```bash
# Process multiple SRA IDs
for sra_id in SRR123456 SRR123457 SRR123458; do
    sra-pipeline run --sra-id $sra_id --output-dir ./batch_results
done

# Process with error handling
for sra_id in SRR123456 SRR123457 SRR123458; do
    echo "Processing $sra_id..."
    if sra-pipeline run --sra-id $sra_id --output-dir ./batch_results; then
        echo "Success: $sra_id"
    else
        echo "Failed: $sra_id"
    fi
done
```

### Development and Testing

```bash
# Validate setup
sra-pipeline validate

# Dry run to test configuration
sra-pipeline run --sra-id SRR123456 --dry-run

# Debug mode
sra-pipeline run --sra-id SRR123456 --log-level DEBUG --log-file debug.log
```

## 📝 Output and Logging

### Log Levels

| Level | Description | Use Case |
|-------|-------------|----------|
| `DEBUG` | Detailed debug information | Development, troubleshooting |
| `INFO` | General information | Normal operation |
| `WARNING` | Warning messages | Potential issues |
| `ERROR` | Error messages only | Production monitoring |

### Log Files

```bash
# Save logs to file
sra-pipeline run --sra-id SRR123456 --log-file pipeline.log

# Combine with log level
sra-pipeline run \
  --sra-id SRR123456 \
  --log-level DEBUG \
  --log-file debug.log
```

### Progress Monitoring

The pipeline provides real-time progress updates:

```bash
sra-pipeline run --sra-id SRR123456
# Shows progress bars and status updates
```

## 🔍 Error Handling

### Common Error Scenarios

**Invalid SRA ID:**
```bash
sra-pipeline run --sra-id INVALID_ID
# Error: SRA ID not found or invalid
```

**Missing FASTQ files:**
```bash
sra-pipeline run --fastq nonexistent.fastq.gz
# Error: FASTQ file not found
```

**Insufficient permissions:**
```bash
sra-pipeline run --sra-id SRR123456 --output-dir /protected/dir
# Error: Permission denied
```

### Error Recovery

```bash
# Check for errors in logs
cat results/SRR123456/logs/pipeline.log | grep ERROR

# Retry with different parameters
sra-pipeline run --sra-id SRR123456 --threads 4  # Reduce threads

# Check system resources
free -h  # Check memory
df -h    # Check disk space
```

## 🐛 Troubleshooting

### Getting Help

```bash
# Show general help
sra-pipeline --help

# Show command-specific help
sra-pipeline run --help
sra-pipeline validate --help

# Show version
sra-pipeline --version
```

### Debugging

```bash
# Enable debug logging
sra-pipeline run --sra-id SRR123456 --log-level DEBUG

# Validate setup with verbose output
sra-pipeline validate --verbose

# Dry run to test configuration
sra-pipeline run --sra-id SRR123456 --dry-run
```

### Performance Issues

```bash
# Monitor resource usage
sra-pipeline run --sra-id SRR123456 --threads 4  # Reduce threads

# Check for memory issues
export SRA_PIPELINE_MAX_MEMORY_GB=8
sra-pipeline run --sra-id SRR123456
```

## 📚 Examples by Use Case

### Research Analysis

```bash
# High-quality analysis
sra-pipeline run \
  --sra-id SRR123456 \
  --output-dir ./research_results \
  --threads 16 \
  --log-level INFO
```

### Production Pipeline

```bash
# Production-ready analysis
sra-pipeline run \
  --sra-id SRR123456 \
  --config production.env \
  --output-dir /data/results \
  --threads 8 \
  --log-level WARNING \
  --log-file /var/log/pipeline.log
```

### Development Testing

```bash
# Development testing
sra-pipeline run \
  --sra-id SRR123456 \
  --output-dir ./test_results \
  --threads 2 \
  --log-level DEBUG \
  --log-file test.log
```

### Batch Processing

```bash
# Batch processing script
#!/bin/bash
for sra_id in $(cat sra_list.txt); do
    echo "Processing $sra_id..."
    sra-pipeline run \
      --sra-id $sra_id \
      --output-dir ./batch_results \
      --threads 8 \
      --log-file batch_${sra_id}.log
done
```

## 🔗 Related Topics

- **[Basic Usage](basic_usage.md)** - Getting started with the CLI
- **[Configuration Guide](../configuration.md)** - Advanced configuration options
- **[Troubleshooting](troubleshooting.md)** - Solving CLI-related issues
- **[Output Format](output_format.md)** - Understanding pipeline outputs

---

**Need help with the CLI?** Check the [Troubleshooting](troubleshooting.md) section or create an issue on GitHub. 