# Input Data Guide

This guide explains the different types of input data that the SRA to Features Pipeline can process, including SRA IDs, FASTQ files, and data formats.

## 📊 Supported Input Types

The pipeline supports two main input types:

1. **SRA IDs** - Sequence Read Archive accession numbers
2. **FASTQ Files** - Raw sequencing data files

## 🔍 SRA (Sequence Read Archive) Data

### What are SRA IDs?

SRA IDs are unique identifiers for sequencing data stored in the Sequence Read Archive (SRA), a public repository maintained by NCBI.

### SRA ID Formats

| Format | Example | Description |
|--------|---------|-------------|
| **SRR** | `SRR123456` | Run accession (most common) |
| **SRX** | `SRX123456` | Experiment accession |
| **SRS** | `SRS123456` | Sample accession |
| **SRP** | `SRP123456` | Study accession |

### Finding SRA IDs

#### NCBI SRA Database
1. Visit [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra)
2. Search for your organism, condition, or study
3. Browse results and note the SRA IDs

#### Example Search
```
Search: "human cancer RNA-seq"
Results: SRR123456, SRR123457, SRR123458...
```

#### Programmatic Access
```bash
# Search SRA programmatically
curl "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=tax_eq(9606)&fields=run_accession&limit=10"
```

### SRA Data Types

| Data Type | Description | Common Use Cases |
|-----------|-------------|------------------|
| **WGS** | Whole Genome Sequencing | Variant discovery, genome assembly |
| **RNA-seq** | RNA Sequencing | Gene expression, transcriptomics |
| **ChIP-seq** | Chromatin Immunoprecipitation | Protein-DNA interactions |
| **ATAC-seq** | Assay for Transposase-Accessible Chromatin | Chromatin accessibility |
| **Bisulfite-seq** | Bisulfite Sequencing | DNA methylation |

### SRA Data Quality

#### Factors to Consider
- **Read Length**: Longer reads provide better alignment
- **Coverage**: Higher coverage improves variant detection
- **Quality Scores**: Higher quality scores indicate better data
- **Library Type**: Single-end vs paired-end sequencing

#### Checking SRA Metadata
```bash
# Get metadata for an SRA ID
curl "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=run_accession=SRR123456&fields=run_accession,fastq_ftp,read_count,base_count,instrument_platform&format=json"
```

## 📁 FASTQ Files

### FASTQ Format

FASTQ files contain raw sequencing data with quality scores:

```
@SRR123456.1 HWI-ST123:123:C0TUSCACXX:1:1101:1208:2450 length=100
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
+
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
```

### FASTQ File Structure

Each read consists of 4 lines:
1. **Header**: Starts with `@`, contains read identifier
2. **Sequence**: The actual DNA/RNA sequence
3. **Separator**: Usually just `+`
4. **Quality**: Quality scores for each base

### FASTQ File Types

#### Single-End Sequencing
```
sample.fastq.gz
```

#### Paired-End Sequencing
```
sample_R1.fastq.gz  # Forward reads
sample_R2.fastq.gz  # Reverse reads
```

### FASTQ File Compression

The pipeline supports compressed FASTQ files:

| Format | Extension | Compression |
|--------|-----------|-------------|
| **Gzipped** | `.fastq.gz` | Most common |
| **Bzipped** | `.fastq.bz2` | Higher compression |
| **Uncompressed** | `.fastq` | No compression |

### FASTQ File Naming Conventions

#### Common Patterns
```
# Illumina naming
sample_S1_L001_R1_001.fastq.gz
sample_S1_L001_R2_001.fastq.gz

# Simple naming
sample_R1.fastq.gz
sample_R2.fastq.gz

# SRA naming
SRR123456_1.fastq.gz
SRR123456_2.fastq.gz
```

## 🎯 Choosing Input Data

### For Variant Discovery
- **Recommended**: WGS data with high coverage (30x+)
- **Format**: Paired-end sequencing
- **Quality**: High quality scores (Q30+)

### For Gene Expression
- **Recommended**: RNA-seq data
- **Format**: Paired-end sequencing preferred
- **Quality**: Good quality scores (Q20+)

### For Quality Assessment
- **Recommended**: Any sequencing data
- **Format**: Single-end or paired-end
- **Quality**: Any quality level (pipeline will assess)

## 📋 Data Requirements

### Minimum Requirements

| Requirement | SRA Data | FASTQ Files |
|-------------|----------|-------------|
| **File Size** | Any size | > 1MB |
| **Read Length** | Any length | > 20bp |
| **Coverage** | Any coverage | Any coverage |
| **Quality** | Any quality | Any quality |

### Recommended Requirements

| Requirement | SRA Data | FASTQ Files |
|-------------|----------|-------------|
| **File Size** | > 100MB | > 100MB |
| **Read Length** | > 50bp | > 50bp |
| **Coverage** | > 10x | > 10x |
| **Quality** | Q20+ | Q20+ |

## 🔧 Preparing Input Data

### Downloading SRA Data

The pipeline automatically downloads SRA data, but you can also download manually:

```bash
# Using fastq-dump
fastq-dump --gzip --split-files SRR123456

# Using fasterq-dump (faster)
fasterq-dump SRR123456
gzip SRR123456_*.fastq
```

### Preparing FASTQ Files

#### Quality Control
```bash
# Run FastQC for quality assessment
fastqc sample_R1.fastq.gz sample_R2.fastq.gz

# Trim low-quality bases
trimmomatic PE sample_R1.fastq.gz sample_R2.fastq.gz \
    sample_R1_trimmed.fastq.gz sample_R1_unpaired.fastq.gz \
    sample_R2_trimmed.fastq.gz sample_R2_unpaired.fastq.gz \
    ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

#### File Organization
```bash
# Create organized directory structure
mkdir -p data/fastq
mv *.fastq.gz data/fastq/

# Check file integrity
md5sum data/fastq/*.fastq.gz > checksums.txt
```

## 📊 Data Validation

### Validating SRA IDs

```bash
# Check if SRA ID exists
curl -s "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=run_accession=SRR123456&fields=run_accession&format=json" | grep -q "SRR123456" && echo "Valid" || echo "Invalid"
```

### Validating FASTQ Files

```bash
# Check file format
head -n 4 sample.fastq.gz | gunzip

# Check file integrity
gunzip -t sample.fastq.gz && echo "Valid" || echo "Corrupted"

# Count reads
zcat sample.fastq.gz | echo $((`wc -l`/4))
```

## 🐛 Common Input Issues

### SRA Data Issues

**"SRA ID not found"**
```bash
# Check SRA ID format
echo $SRA_ID | grep -E "^SR[RXSP][0-9]+$"

# Verify SRA ID exists
curl "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=run_accession=$SRA_ID"
```

**"Download failed"**
```bash
# Check internet connection
ping -c 3 ftp.ncbi.nlm.nih.gov

# Try alternative download method
prefetch SRR123456
```

### FASTQ File Issues

**"Invalid FASTQ format"**
```bash
# Check file format
head -n 8 sample.fastq.gz | gunzip

# Fix common issues
sed 's/^@/@SRR123456./' sample.fastq.gz > fixed.fastq.gz
```

**"Paired-end mismatch"**
```bash
# Count reads in each file
zcat sample_R1.fastq.gz | echo $((`wc -l`/4))
zcat sample_R2.fastq.gz | echo $((`wc -l`/4))

# Files should have the same number of reads
```

## 📚 Best Practices

### Data Selection
- **Choose appropriate data type** for your analysis
- **Check data quality** before processing
- **Verify file integrity** after download
- **Use consistent naming** conventions

### Data Management
- **Organize files** in logical directory structures
- **Keep backups** of original data
- **Document data sources** and processing steps
- **Use checksums** to verify file integrity

### Performance Considerations
- **Use compressed files** to save disk space
- **Process smaller files** for testing
- **Monitor disk space** during processing
- **Use appropriate file systems** for large datasets

## 🔗 Related Topics

- **[Basic Usage](basic_usage.md)** - How to use input data with the pipeline
- **[Output Format](output_format.md)** - Understanding pipeline outputs
- **[Quality Control](quality_control.md)** - Assessing data quality
- **[Troubleshooting](troubleshooting.md)** - Solving input-related issues

---

**Need help with input data?** Check the [Troubleshooting](troubleshooting.md) section or create an issue on GitHub. 