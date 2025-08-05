# Output Format Guide

This guide explains the output format of the SRA to Features Pipeline, including the structure of results files and how to interpret the extracted features.

## 📊 Output Structure

The pipeline generates a structured output with the following organization:

```
output_directory/
├── sample_id/
│   ├── features.json      # Main results (JSON format)
│   ├── summary.txt        # Human-readable summary
│   └── logs/              # Processing logs
│       ├── pipeline.log   # Pipeline execution log
│       ├── download.log   # Download log
│       ├── alignment.log  # Alignment log
│       └── qc.log         # Quality control log
└── logs/                  # Global pipeline logs
    ├── pipeline.log       # Main pipeline log
    └── performance.log    # Performance metrics
```

## 📄 Main Output Files

### `features.json`

The primary output file containing all extracted features in structured JSON format:

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
      "variant_count": 15,
      "coverage_mean": 25.5,
      "coverage_std": 5.2
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
      "dn_ds_ratio": 0.67,
      "expression_level": 125.5
    }
  ],
  "cnv_regions": [
    {
      "chromosome": "chr1",
      "start": 10000,
      "end": 20000,
      "copy_number": 3.0,
      "confidence": 0.95,
      "type": "gain",
      "genes_affected": ["GENE1", "GENE2"]
    }
  ],
  "quality_metrics": {
    "total_reads": 10000000,
    "mapped_reads": 9500000,
    "mapping_rate": 95.0,
    "mean_coverage": 30.0,
    "coverage_std": 5.0,
    "gc_content": 45.2,
    "duplication_rate": 12.5,
    "quality_scores": {
      "q10": 95.5,
      "q20": 90.2,
      "q30": 85.1
    }
  },
  "metadata": {
    "pipeline_version": "1.0.0",
    "processing_date": "2024-01-15T10:30:00Z",
    "reference_genome": "GRCh38",
    "analysis_parameters": {
      "min_quality_score": 20,
      "min_coverage": 10,
      "bin_size_gvs": 100000,
      "bin_size_cnv": 100000
    }
  },
  "processing_time": 3600.5
}
```

### `summary.txt`

Human-readable summary of the analysis:

```
SRA to Features Pipeline Summary
========================================

Sample Information:
  SRA ID: SRR123456
  Processing Date: 2024-01-15 10:30:00
  Pipeline Version: 1.0.0
  Reference Genome: GRCh38

Quality Metrics:
  Total Reads: 10,000,000
  Mapped Reads: 9,500,000
  Mapping Rate: 95.0%
  Mean Coverage: 30.0x
  Coverage Std Dev: 5.0x
  GC Content: 45.2%
  Duplication Rate: 12.5%

Variant Analysis:
  Total Variants: 15,000
  Genes with Variants: 500
  Mean Variants per Gene: 30.0
  Synonymous Variants: 9,000
  Nonsynonymous Variants: 6,000
  Dn/Ds Ratio: 0.67

Copy Number Variation:
  Total CNV Regions: 25
  Gains: 15
  Losses: 10
  High Confidence (>0.9): 20

Fragment Analysis:
  Mean Fragment Length: 150.0 bp
  Median Fragment Length: 150.0 bp
  Fragment Length Std Dev: 25.0 bp
  Fragment Length Range: 50-300 bp

Performance:
  Processing Time: 1h 0m 5s
  Memory Usage: 8.5 GB
  CPU Usage: 85%
```

## 🔍 Feature Descriptions

### Fragment Statistics

Fragment length distribution statistics from sequencing data:

| Field | Type | Description |
|-------|------|-------------|
| `mean` | float | Mean fragment length in base pairs |
| `median` | float | Median fragment length in base pairs |
| `std` | float | Standard deviation of fragment lengths |
| `min` | float | Minimum fragment length |
| `max` | float | Maximum fragment length |
| `count` | int | Total number of fragments analyzed |

### Genomic Bins

Variant counts and coverage statistics for genomic regions:

| Field | Type | Description |
|-------|------|-------------|
| `chromosome` | string | Chromosome name |
| `start` | int | Start position (0-based) |
| `end` | int | End position (1-based) |
| `variant_count` | int | Number of variants in this bin |
| `coverage_mean` | float | Mean coverage in this bin |
| `coverage_std` | float | Coverage standard deviation |

### Gene Statistics

Gene-level variant and expression analysis:

| Field | Type | Description |
|-------|------|-------------|
| `gene_name` | string | Gene identifier |
| `chromosome` | string | Chromosome name |
| `start` | int | Gene start position |
| `end` | int | Gene end position |
| `total_variants` | int | Total variants in gene |
| `synonymous_variants` | int | Synonymous variant count |
| `nonsynonymous_variants` | int | Nonsynonymous variant count |
| `dn_ds_ratio` | float | Ratio of nonsynonymous to synonymous variants |
| `expression_level` | float | Gene expression level (if RNA-seq) |

### CNV Regions

Copy number variation analysis results:

| Field | Type | Description |
|-------|------|-------------|
| `chromosome` | string | Chromosome name |
| `start` | int | Region start position |
| `end` | int | Region end position |
| `copy_number` | float | Estimated copy number |
| `confidence` | float | Confidence score (0-1) |
| `type` | string | CNV type ("gain", "loss", "neutral") |
| `genes_affected` | list | List of genes in this region |

### Quality Metrics

Comprehensive quality assessment metrics:

| Field | Type | Description |
|-------|------|-------------|
| `total_reads` | int | Total number of sequencing reads |
| `mapped_reads` | int | Number of reads mapped to reference |
| `mapping_rate` | float | Percentage of reads mapped |
| `mean_coverage` | float | Mean coverage across genome |
| `coverage_std` | float | Coverage standard deviation |
| `gc_content` | float | GC content percentage |
| `duplication_rate` | float | Duplicate read percentage |
| `quality_scores` | object | Quality score distributions |

## 📈 Interpreting Results

### Quality Assessment

#### Good Quality Indicators
- **Mapping Rate**: > 90%
- **Mean Coverage**: > 20x for variant calling
- **Quality Scores**: Q30 > 80%
- **Duplication Rate**: < 20%

#### Warning Signs
- **Mapping Rate**: < 70%
- **Mean Coverage**: < 10x
- **Quality Scores**: Q30 < 60%
- **Duplication Rate**: > 50%

### Variant Analysis

#### Variant Distribution
- **Synonymous vs Nonsynonymous**: Normal Dn/Ds ratio ~1.0
- **Variant Density**: Varies by genomic region
- **Gene Impact**: Higher impact variants are less common

#### CNV Analysis
- **Confidence Scores**: > 0.9 indicates high confidence
- **Copy Number**: 2.0 is normal, >3.0 or <1.0 indicates CNV
- **Gene Impact**: Check affected genes for biological relevance

## 🔧 Working with Output Files

### Parsing JSON Output

#### Python Example
```python
import json

# Load results
with open('results/SRR123456/features.json', 'r') as f:
    data = json.load(f)

# Extract specific features
total_variants = sum(bin['variant_count'] for bin in data['genomic_bins'])
mapping_rate = data['quality_metrics']['mapping_rate']
mean_coverage = data['quality_metrics']['mean_coverage']

print(f"Total variants: {total_variants}")
print(f"Mapping rate: {mapping_rate}%")
print(f"Mean coverage: {mean_coverage}x")
```

#### R Example
```r
library(jsonlite)

# Load results
data <- fromJSON("results/SRR123456/features.json")

# Extract features
total_variants <- sum(data$genomic_bins$variant_count)
mapping_rate <- data$quality_metrics$mapping_rate
mean_coverage <- data$quality_metrics$mean_coverage

cat("Total variants:", total_variants, "\n")
cat("Mapping rate:", mapping_rate, "%\n")
cat("Mean coverage:", mean_coverage, "x\n")
```

### Batch Analysis

#### Processing Multiple Samples
```bash
# Combine results from multiple samples
for sample in results/*/features.json; do
    echo "Processing $sample..."
    # Your analysis code here
done
```

#### Creating Summary Tables
```python
import json
import pandas as pd

results = []
for sample_file in glob.glob("results/*/features.json"):
    with open(sample_file, 'r') as f:
        data = json.load(f)
    
    results.append({
        'sample_id': data['sra_id'],
        'total_variants': sum(bin['variant_count'] for bin in data['genomic_bins']),
        'mapping_rate': data['quality_metrics']['mapping_rate'],
        'mean_coverage': data['quality_metrics']['mean_coverage'],
        'processing_time': data['processing_time']
    })

df = pd.DataFrame(results)
df.to_csv('summary_table.csv', index=False)
```

## 📊 Visualization

### Quality Metrics Dashboard
```python
import matplotlib.pyplot as plt
import seaborn as sns

# Create quality metrics plot
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Mapping rate distribution
sns.histplot(data=df, x='mapping_rate', ax=axes[0,0])
axes[0,0].set_title('Mapping Rate Distribution')

# Coverage distribution
sns.histplot(data=df, x='mean_coverage', ax=axes[0,1])
axes[0,1].set_title('Coverage Distribution')

# Variant count distribution
sns.histplot(data=df, x='total_variants', ax=axes[1,0])
axes[1,0].set_title('Variant Count Distribution')

# Processing time vs coverage
sns.scatterplot(data=df, x='mean_coverage', y='processing_time', ax=axes[1,1])
axes[1,1].set_title('Processing Time vs Coverage')

plt.tight_layout()
plt.savefig('quality_dashboard.png', dpi=300, bbox_inches='tight')
```

## 🐛 Common Output Issues

### Missing Files
```bash
# Check if output files exist
ls -la results/SRR123456/

# Check for error logs
cat results/SRR123456/logs/pipeline.log | grep ERROR
```

### Invalid JSON
```bash
# Validate JSON format
python -m json.tool results/SRR123456/features.json > /dev/null && echo "Valid JSON" || echo "Invalid JSON"
```

### Incomplete Results
```bash
# Check processing time
grep "processing_time" results/SRR123456/features.json

# Check for timeout errors
grep -i "timeout\|error" results/SRR123456/logs/*.log
```

## 📚 Best Practices

### Data Management
- **Keep original outputs** for reproducibility
- **Create backups** of important results
- **Document analysis parameters** used
- **Version control** your analysis scripts

### Quality Control
- **Review quality metrics** before downstream analysis
- **Check for systematic biases** across samples
- **Validate results** with known datasets
- **Compare with published benchmarks**

### Performance Monitoring
- **Track processing times** for optimization
- **Monitor resource usage** during analysis
- **Log analysis parameters** for reproducibility
- **Document any issues** encountered

## 🔗 Related Topics

- **[Basic Usage](basic_usage.md)** - How to generate output files
- **[Input Data](input_data.md)** - Understanding input data formats
- **[Quality Control](quality_control.md)** - Interpreting quality metrics
- **[Troubleshooting](troubleshooting.md)** - Solving output-related issues

---

**Need help interpreting results?** Check the [Troubleshooting](troubleshooting.md) section or create an issue on GitHub. 