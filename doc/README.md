# SRA to Features Pipeline

A modern, production-ready bioinformatics pipeline for extracting genomic features from SRA (Sequence Read Archive) data for LLM training.

## 🚀 Features

- **Modern Architecture**: Built with Python 3.8+, type hints, and Pydantic for data validation
- **Comprehensive Logging**: Structured logging with performance monitoring
- **Error Handling**: Robust error handling with retry logic and graceful failures
- **Configuration Management**: Flexible configuration using environment variables or config files
- **Quality Control**: Built-in quality assessment and validation
- **Performance Monitoring**: Real-time performance metrics and resource usage tracking
- **Testing**: Comprehensive test suite with unit and integration tests
- **Documentation**: Complete API documentation and user guides
- **CLI Interface**: Modern command-line interface with rich output

## 📋 Requirements

### System Requirements
- Python 3.8 or higher
- 8GB RAM minimum (16GB+ recommended)
- 50GB+ disk space for reference genomes and intermediate files
- Unix-like operating system (Linux/macOS)

### External Tools
The pipeline requires the following bioinformatics tools to be installed and available in PATH:

- **BWA** (v0.7.17+) - Sequence alignment
- **SAMtools** (v1.10+) - SAM/BAM file manipulation
- **BCFtools** (v1.10+) - Variant calling
- **BEDtools** (v2.29+) - Genome arithmetic
- **FastQC** (v0.11+) - Quality control
- **fastq-dump** (SRA Toolkit) - SRA data download
- **tabix/bgzip** (HTSlib) - VCF compression and indexing
- **Java** (v8+) - For snpEff
- **CNVpytor** (v1.3+) - Copy number variation analysis

## 🛠️ Installation

### Option 1: Install from PyPI (Recommended)

```bash
pip install sra-to-features-pipeline
```

### Option 2: Install from Source

```bash
# Clone the repository
git clone https://github.com/your-org/sra-to-features-pipeline.git
cd sra-to-features-pipeline

# Install in development mode
pip install -e ".[dev]"
```

### Option 3: Using Conda

```bash
# Create a new conda environment
conda create -n sra-pipeline python=3.9
conda activate sra-pipeline

# Install the pipeline
pip install sra-to-features-pipeline

# Install bioinformatics tools
conda install -c bioconda bwa samtools bcftools bedtools fastqc sra-tools htslib
```

## ⚙️ Configuration

### Environment Variables

The pipeline can be configured using environment variables with the `SRA_PIPELINE_` prefix:

```bash
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

## 🚀 Quick Start

### 1. Validate Setup

First, validate that your system is properly configured:

```bash
sra-pipeline validate
```

### 2. Run with SRA ID

```bash
sra-pipeline run --sra-id SRR123456 --output-dir ./results --threads 8
```

### 3. Run with FASTQ Files

```bash
sra-pipeline run \
  --fastq sample_1.fastq.gz sample_2.fastq.gz \
  --output-dir ./results \
  --threads 8
```

### 4. Check Results

```bash
# View summary
cat ./results/SRR123456/summary.txt

# View features in JSON format
cat ./results/SRR123456/features.json
```

## 📖 Usage Examples

### Basic Usage

```bash
# Run pipeline with default settings
sra-pipeline run --sra-id SRR123456

# Run with custom output directory
sra-pipeline run --sra-id SRR123456 --output-dir /path/to/output

# Run with multiple threads
sra-pipeline run --sra-id SRR123456 --threads 16

# Run with custom configuration file
sra-pipeline run --sra-id SRR123456 --config my_config.env
```

### Advanced Usage

```bash
# Run with custom logging
sra-pipeline run \
  --sra-id SRR123456 \
  --log-level DEBUG \
  --log-file pipeline.log

# Run with FASTQ files (paired-end)
sra-pipeline run \
  --fastq sample_R1.fastq.gz sample_R2.fastq.gz \
  --output-dir ./results

# Run with FASTQ files (single-end)
sra-pipeline run \
  --fastq sample.fastq.gz \
  --output-dir ./results

# Dry run (validate without processing)
sra-pipeline run --sra-id SRR123456 --dry-run
```

### Batch Processing

```bash
# Process multiple SRA IDs
for sra_id in SRR123456 SRR123457 SRR123458; do
    sra-pipeline run --sra-id $sra_id --output-dir ./batch_results
done
```

## 📊 Output Format

The pipeline generates a structured output with the following components:

### Feature Set (JSON)

```json
{
  "sra_id": "SRR123456",
  "fragment_stats": {
    "mean": 150.0,
    "median": 150.0,
    "std": 25.0,
    "min": 50.0,
    "max": 300.0,
    "count": 1000000
  },
  "genomic_bins": [
    {
      "chromosome": "chr1",
      "start": 0,
      "end": 100000,
      "variant_count": 15
    }
  ],
  "gene_stats": [
    {
      "gene_name": "GENE1",
      "chromosome": "chr1",
      "start": 1000,
      "end": 5000,
      "total_variants": 50,
      "synonymous_variants": 30,
      "nonsynonymous_variants": 20,
      "dn_ds_ratio": 0.67
    }
  ],
  "cnv_regions": [
    {
      "chromosome": "chr1",
      "start": 10000,
      "end": 20000,
      "copy_number": 3.0,
      "confidence": 0.95,
      "type": "gain"
    }
  ],
  "quality_metrics": {
    "total_reads": 10000000,
    "mapped_reads": 9500000,
    "mapping_rate": 95.0,
    "mean_coverage": 30.0,
    "coverage_std": 5.0,
    "gc_content": 45.2,
    "duplication_rate": 12.5
  },
  "processing_time": 3600.5,
  "pipeline_version": "1.0.0"
}
```

### Summary Statistics

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

## 🧪 Testing

### Run Tests

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=src/sra_pipeline

# Run specific test categories
pytest -m unit
pytest -m integration
pytest -m "not slow"
```

### Test Configuration

```bash
# Test configuration validation
sra-pipeline validate --config test_config.env
```

## 📚 Documentation

### API Documentation

Generate API documentation:

```bash
sra-pipeline setup-docs --output-dir ./doc
```

### Development Documentation

- [Contributing Guidelines](CONTRIBUTING.md)
- [Development Setup](doc/development.md)
- [API Reference](doc/api/)
- [User Guide](doc/user_guide/)

## 🔧 Development

### Setup Development Environment

```bash
# Clone repository
git clone https://github.com/your-org/sra-to-features-pipeline.git
cd sra-to-features-pipeline

# Install development dependencies
pip install -e ".[dev]"

# Install pre-commit hooks
pre-commit install

# Run code formatting
black src/ test/

# Run linting
flake8 src/ test/

# Run type checking
mypy src/
```

### Project Structure

```
sra-to-features-pipeline/
├── src/sra_pipeline/          # Main package
│   ├── cli/                   # Command-line interface
│   ├── config/                # Configuration management
│   ├── core/                  # Core pipeline modules
│   │   ├── download.py        # SRA data download
│   │   ├── alignment.py       # BWA alignment
│   │   ├── variant_calling.py # Variant calling
│   │   ├── feature_extraction.py # Feature extraction
│   │   └── quality_control.py # Quality control
│   ├── models/                # Data models
│   └── utils/                 # Utilities
├── test/                      # Test suite
│   ├── unit/                  # Unit tests
│   └── integration/           # Integration tests
├── doc/                       # Documentation
├── scripts/                   # Utility scripts
├── requirements.txt           # Python dependencies
├── pyproject.toml            # Project configuration
└── README.md                 # This file
```

## 🤝 Contributing

We welcome contributions! Please see our [Contributing Guidelines](CONTRIBUTING.md) for details.

### Development Workflow

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests for new functionality
5. Ensure all tests pass
6. Submit a pull request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🆘 Support

### Getting Help

- **Documentation**: Check the [documentation](doc/) first
- **Issues**: Report bugs and request features on [GitHub Issues](https://github.com/your-org/sra-to-features-pipeline/issues)
- **Discussions**: Join our [GitHub Discussions](https://github.com/your-org/sra-to-features-pipeline/discussions)

### Common Issues

#### Installation Issues

```bash
# If you encounter permission issues
pip install --user sra-to-features-pipeline

# If you need to upgrade pip
python -m pip install --upgrade pip
```

#### Runtime Issues

```bash
# Check if all required tools are available
which bwa samtools bcftools bedtools fastqc fastq-dump

# Validate your configuration
sra-pipeline validate

# Check system resources
free -h
df -h
```

## 🔄 Version History

- **v1.0.0** - Initial release with modern architecture
- **v0.9.0** - Beta release with core functionality
- **v0.8.0** - Alpha release with basic pipeline

## 🙏 Acknowledgments

- The bioinformatics community for developing the tools this pipeline uses
- Contributors and users who provide feedback and improvements
- Funding agencies that support bioinformatics research

---

**Note**: This pipeline is designed for research use. Please ensure you have appropriate permissions and follow ethical guidelines when processing genomic data.

