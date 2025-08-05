# Basic Usage Guide

This guide covers the essential commands and workflows you need to get started with the SRA to Features Pipeline.

## 🚀 Quick Start

### 1. Validate Your Setup

Before running any analysis, validate that your system is properly configured:

```bash
# Check pipeline version
sra-pipeline --version

# Validate setup
sra-pipeline validate
```

You should see output indicating that all required tools are available.

### 2. Run Your First Analysis

#### Using an SRA ID

```bash
# Basic run with SRA ID
sra-pipeline run --sra-id SRR123456 --output-dir ./results
```

#### Using FASTQ Files

```bash
# Single-end sequencing
sra-pipeline run --fastq sample.fastq.gz --output-dir ./results

# Paired-end sequencing
sra-pipeline run --fastq sample_R1.fastq.gz sample_R2.fastq.gz --output-dir ./results
```

### 3. Check Your Results

```bash
# Navigate to results
cd results

# List output files
ls -la

# View summary
cat SRR123456/summary.txt

# Examine features
cat SRR123456/features.json
```

## 📋 Essential Commands

### Basic Pipeline Commands

| Command | Description | Example |
|---------|-------------|---------|
| `run` | Run the pipeline | `sra-pipeline run --sra-id SRR123456` |
| `validate` | Validate setup | `sra-pipeline validate` |
| `--help` | Show help | `sra-pipeline --help` |
| `--version` | Show version | `sra-pipeline --version` |

### Common Options

| Option | Description | Default | Example |
|--------|-------------|---------|---------|
| `--sra-id` | SRA accession ID | None | `--sra-id SRR123456` |
| `--fastq` | FASTQ file paths | None | `--fastq sample.fastq.gz` |
| `--output-dir` | Output directory | `./output` | `--output-dir ./results` |
| `--threads` | Number of threads | 1 | `--threads 8` |
| `--log-level` | Logging level | INFO | `--log-level DEBUG` |
| `--dry-run` | Test without running | False | `--dry-run` |

## 🔧 Basic Configuration

### Environment Variables

Set basic configuration using environment variables:

```bash
# Set in your shell profile (~/.bashrc, ~/.zshrc)
export SRA_PIPELINE_THREADS=8
export SRA_PIPELINE_OUTPUT_DIR="/path/to/output"
export SRA_PIPELINE_LOG_LEVEL="INFO"
```

### Configuration File

Create a `.env` file in your working directory:

```ini
# Basic configuration
SRA_PIPELINE_THREADS=8
SRA_PIPELINE_MIN_QUALITY_SCORE=20
SRA_PIPELINE_MIN_COVERAGE=10
SRA_PIPELINE_LOG_LEVEL=INFO
```

## 📊 Understanding Output

### Output Structure

```
results/
├── SRR123456/
│   ├── features.json      # Main results (JSON format)
│   ├── summary.txt        # Human-readable summary
│   └── logs/              # Processing logs
└── logs/                  # Pipeline logs
```

### Key Output Files

#### `features.json`
Contains all extracted features in structured JSON format:

```json
{
  "sra_id": "SRR123456",
  "fragment_stats": {
    "mean": 150.0,
    "median": 150.0,
    "std": 25.0
  },
  "genomic_bins": [...],
  "gene_stats": [...],
  "cnv_regions": [...],
  "quality_metrics": {...},
  "processing_time": 3600.5,
  "pipeline_version": "1.0.0"
}
```

#### `summary.txt`
Human-readable summary of the analysis:

```
SRA to Features Pipeline Summary
========================================

Sample ID: SRR123456
Total Variants: 15,000
Genes with Variants: 500
CNV Regions: 25
Mapping Rate: 95.0%
Mean Coverage: 30.0x
Processing Time: 1h 0m 5s
```

## 🎯 Common Workflows

### Workflow 1: Single SRA Analysis

```bash
# 1. Validate setup
sra-pipeline validate

# 2. Run analysis
sra-pipeline run --sra-id SRR123456 --output-dir ./results --threads 8

# 3. Check results
ls -la results/SRR123456/
cat results/SRR123456/summary.txt
```

### Workflow 2: FASTQ File Analysis

```bash
# 1. Prepare FASTQ files
ls *.fastq.gz

# 2. Run analysis
sra-pipeline run \
  --fastq sample_R1.fastq.gz sample_R2.fastq.gz \
  --output-dir ./results \
  --threads 8

# 3. Examine results
cat results/sample/summary.txt
```

### Workflow 3: Multiple Samples (Batch Processing)

```bash
# Process multiple samples with VCF merging (recommended)
sra-pipeline batch --sra-ids SRR123456,SRR123457,SRR123458 --output-dir ./batch_results

# Process multiple samples without VCF merging
sra-pipeline batch --sra-ids SRR123456,SRR123457,SRR123458 --output-dir ./batch_results --no-merge-vcfs

# Process multiple samples with custom configuration
sra-pipeline batch \
  --sra-ids SRR123456,SRR123457,SRR123458 \
  --output-dir ./batch_results \
  --threads 8 \
  --config batch_config.env
```

**Note**: The batch command automatically merges VCF files from all samples into a single `merged_variants.vcf.gz` file, while keeping individual sample VCF files and cleaning up intermediate files.

### Workflow 4: ML Feature Table Generation

```bash
# Create ML-ready feature table from pipeline results
sra-pipeline create-ml-table \
  --input-dir ./batch_results \
  --output-file ml_features.csv \
  --format csv

# Create feature table in different formats
sra-pipeline create-ml-table \
  --input-dir ./batch_results \
  --output-file ml_features.parquet \
  --format parquet

# Create feature table with TSV format
sra-pipeline create-ml-table \
  --input-dir ./batch_results \
  --output-file ml_features.tsv \
  --format tsv
```

**Note**: The ML feature table combines features from all samples into a single tabular format suitable for machine learning training.

## ⚙️ Customization Options

### Performance Tuning

```bash
# Use more threads for faster processing
sra-pipeline run --sra-id SRR123456 --threads 16

# Reduce memory usage
export SRA_PIPELINE_MAX_MEMORY_GB=8
sra-pipeline run --sra-id SRR123456 --threads 4
```

### Quality Control

```bash
# Stricter quality filters
export SRA_PIPELINE_MIN_QUALITY_SCORE=30
export SRA_PIPELINE_MIN_COVERAGE=20
sra-pipeline run --sra-id SRR123456
```

### Logging

```bash
# Verbose logging for debugging
sra-pipeline run --sra-id SRR123456 --log-level DEBUG

# Save logs to file
sra-pipeline run --sra-id SRR123456 --log-file pipeline.log
```

## 🔍 Monitoring Progress

### Real-time Monitoring

The pipeline provides real-time progress updates:

```bash
sra-pipeline run --sra-id SRR123456
# Output shows progress bars and status updates
```

### Log Files

Check log files for detailed information:

```bash
# View pipeline logs
cat results/logs/pipeline.log

# View sample-specific logs
cat results/SRR123456/logs/processing.log
```

## 🐛 Basic Troubleshooting

### Common Issues

**"Command not found"**
```bash
# Reinstall pipeline
pip install -e .

# Check installation
which sra-pipeline
```

**"Tool not found"**
```bash
# Check if tools are installed
which bwa samtools bcftools

# Install missing tools
./scripts/install/install_tools.sh
```

**"Permission denied"**
```bash
# Check file permissions
ls -la output/

# Fix permissions
chmod 755 output/
```

### Getting Help

```bash
# Show help
sra-pipeline --help
sra-pipeline run --help

# Validate setup
sra-pipeline validate

# Test configuration
sra-pipeline run --sra-id SRR123456 --dry-run
```

## 📚 Next Steps

Now that you understand the basics:

1. **Learn about input data**: See [Input Data](input_data.md)
2. **Understand output format**: See [Output Format](output_format.md)
3. **Explore advanced options**: See [Command Line Interface](cli_reference.md)
4. **Try batch processing**: See [Batch Processing](batch_processing.md)

## 💡 Tips for Success

- **Start small**: Use smaller datasets for testing
- **Monitor resources**: Keep an eye on disk space and memory
- **Use dry runs**: Test configurations before full runs
- **Check logs**: Review log files for warnings or errors
- **Backup results**: Keep copies of important analyses

---

**Need help?** Check the [Troubleshooting](troubleshooting.md) section or create an issue on GitHub. 