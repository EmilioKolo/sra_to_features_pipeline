# Quick Start Guide

Get up and running with the SRA to Features Pipeline in minutes! This guide will walk you through your first analysis.

## 🚀 Prerequisites

Before starting, ensure you have:

- ✅ Python 3.8+ installed
- ✅ Bioinformatics tools installed (see [Installation Guide](installation.md))
- ✅ Pipeline installed: `pip install sra-to-features-pipeline`
- ✅ Internet connection for downloading SRA data

## 📋 Step 1: Validate Your Setup

First, let's make sure everything is working correctly:

```bash
# Check if the pipeline is installed
sra-pipeline --version

# Validate your setup
sra-pipeline validate
```

You should see output indicating that all required tools are available.

## 🎯 Step 2: Run Your First Analysis

### Option A: Using an SRA ID (Recommended for beginners)

```bash
# Run the pipeline with a test SRA ID
sra-pipeline run --sra-id SRR123456 --output-dir ./my_first_analysis
```

### Option B: Using FASTQ Files

If you already have FASTQ files:

```bash
# For single-end sequencing
sra-pipeline run --fastq sample.fastq.gz --output-dir ./my_first_analysis

# For paired-end sequencing
sra-pipeline run --fastq sample_R1.fastq.gz sample_R2.fastq.gz --output-dir ./my_first_analysis
```

## 📊 Step 3: Check Your Results

After the pipeline completes, explore your results:

```bash
# Navigate to your output directory
cd my_first_analysis

# List the contents
ls -la

# View the summary
cat SRR123456/summary.txt

# Examine the features (JSON format)
cat SRR123456/features.json
```

## 🔍 Understanding Your Output

The pipeline generates several files:

### `features.json`
Contains the extracted genomic features in structured JSON format:

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

### `summary.txt`
Human-readable summary of the analysis:

```
SRA to Features Pipeline Summary
========================================

total_variants: 15000
total_genes_with_variants: 500
total_cnv_regions: 25
mean_variants_per_bin: 12.5
mapping_rate: 95.0
mean_coverage: 30.0
```

## ⚙️ Customizing Your Analysis

### Using Multiple Threads

```bash
sra-pipeline run --sra-id SRR123456 --threads 8 --output-dir ./my_analysis
```

### Custom Configuration

Create a `.env` file in your working directory:

```ini
SRA_PIPELINE_THREADS=8
SRA_PIPELINE_MIN_QUALITY_SCORE=20
SRA_PIPELINE_MIN_COVERAGE=10
```

Then run:

```bash
sra-pipeline run --sra-id SRR123456 --output-dir ./my_analysis
```

### Dry Run (Testing)

Test your configuration without running the full pipeline:

```bash
sra-pipeline run --sra-id SRR123456 --dry-run
```

## 🐛 Troubleshooting Common Issues

### "Command not found: sra-pipeline"

```bash
# Reinstall the pipeline
pip install -e .

# Or install globally
pip install sra-to-features-pipeline
```

### "Tool not found" errors

```bash
# Check if tools are installed
which bwa samtools bcftools

# Install missing tools (see Installation Guide)
./scripts/install/install_tools.sh
```

### Permission errors

```bash
# Use user installation
pip install --user sra-to-features-pipeline

# Or create a virtual environment
python -m venv venv
source venv/bin/activate
pip install sra-to-features-pipeline
```

### Out of memory errors

```bash
# Reduce thread count
sra-pipeline run --sra-id SRR123456 --threads 2 --output-dir ./my_analysis

# Or increase system swap space
```

## 📚 Next Steps

Now that you've successfully run your first analysis:

1. **Explore the [User Guide](user_guide/)** for detailed usage instructions
2. **Check the [Configuration Guide](configuration.md)** for advanced settings
3. **Review the [Output Format](output_format.md)** to understand your results
4. **Try different SRA IDs** to explore various datasets

## 🎉 Congratulations!

You've successfully completed your first analysis with the SRA to Features Pipeline! 

The pipeline has extracted genomic features that can be used for:
- Machine learning model training
- Comparative genomics analysis
- Quality assessment of sequencing data
- Research and discovery

## 💡 Tips for Success

- **Start small**: Use smaller SRA datasets for testing
- **Monitor resources**: Keep an eye on disk space and memory usage
- **Use dry runs**: Test configurations before full runs
- **Check logs**: Review log files for any warnings or errors
- **Backup results**: Keep copies of important analysis results

---

**Need help?** Check our [Support Guide](support.md) or create an issue on GitHub. 