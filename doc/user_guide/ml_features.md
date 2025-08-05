# 🤖 ML Feature Tables

This guide explains how to create and use machine learning-ready feature tables from the SRA to Features Pipeline.

## 📊 Overview

The pipeline can generate tabular feature datasets suitable for machine learning training. These feature tables combine genomic, quality, and variant information from multiple samples into a structured format.

## 🚀 Quick Start

### Basic Usage

```bash
# Process multiple samples
sra-pipeline batch --sra-ids SRR123456,SRR123457,SRR123458 --output-dir ./batch_results

# Create ML feature table
sra-pipeline create-ml-table \
  --input-dir ./batch_results \
  --output-file ml_features.csv \
  --format csv
```

### Complete Workflow

```bash
# 1. Process multiple samples
sra-pipeline batch \
  --sra-ids SRR123456,SRR123457,SRR123458,SRR123459,SRR123460 \
  --output-dir ./batch_results \
  --threads 8

# 2. Create ML feature table
sra-pipeline create-ml-table \
  --input-dir ./batch_results \
  --output-file ml_features.csv \
  --format csv

# 3. Verify the output
head -5 ml_features.csv
```

## 📋 Feature Types

The ML feature table includes the following feature categories:

### 🔬 Quality Metrics Features
- `total_reads`: Total number of sequencing reads
- `mapped_reads`: Number of reads mapped to reference
- `mapping_rate`: Percentage of reads successfully mapped
- `mean_coverage`: Average coverage across the genome
- `coverage_std`: Standard deviation of coverage
- `gc_content`: GC content percentage
- `duplication_rate`: Duplicate read rate
- `quality_score_*`: Various quality scores

### 📏 Fragment Length Features
- `fragment_mean`: Mean fragment length
- `fragment_median`: Median fragment length
- `fragment_std`: Standard deviation of fragment lengths
- `fragment_min`: Minimum fragment length
- `fragment_max`: Maximum fragment length
- `fragment_count`: Number of fragments analyzed

### 🧬 Genomic Variant Features
- `total_genomic_variants`: Total number of variants
- `mean_variants_per_bin`: Average variants per genomic bin
- `std_variants_per_bin`: Standard deviation of variants per bin
- `max_variants_per_bin`: Maximum variants in any bin
- `min_variants_per_bin`: Minimum variants in any bin
- `bins_with_variants`: Number of bins containing variants
- `total_genomic_bins`: Total number of genomic bins
- `variant_density`: Variant density across genome

### 🧪 Gene-Level Features
- `total_genes_with_variants`: Number of genes with variants
- `total_gene_variants`: Total variants in genes
- `mean_variants_per_gene`: Average variants per gene
- `std_variants_per_gene`: Standard deviation of variants per gene
- `total_synonymous_variants`: Number of synonymous variants
- `total_nonsynonymous_variants`: Number of nonsynonymous variants
- `mean_dn_ds_ratio`: Average dN/dS ratio
- `std_dn_ds_ratio`: Standard deviation of dN/dS ratios
- `genes_with_high_variants`: Genes with >10 variants

### 🔄 Copy Number Variation Features
- `total_cnv_regions`: Total CNV regions detected
- `cnv_gains`: Number of copy number gains
- `cnv_losses`: Number of copy number losses
- `mean_copy_number`: Average copy number
- `std_copy_number`: Standard deviation of copy numbers
- `mean_cnv_confidence`: Average CNV confidence score
- `high_confidence_cnvs`: CNVs with confidence >0.9
- `cnv_burden`: Total CNV burden

## 📁 Output Formats

The pipeline supports multiple output formats:

### CSV Format
```bash
sra-pipeline create-ml-table \
  --input-dir ./batch_results \
  --output-file ml_features.csv \
  --format csv
```

**Advantages:**
- Human-readable
- Compatible with most tools
- Easy to inspect and edit

### TSV Format
```bash
sra-pipeline create-ml-table \
  --input-dir ./batch_results \
  --output-file ml_features.tsv \
  --format tsv
```

**Advantages:**
- Tab-separated values
- Good for data with commas
- Compatible with many bioinformatics tools

### Parquet Format
```bash
sra-pipeline create-ml-table \
  --input-dir ./batch_results \
  --output-file ml_features.parquet \
  --format parquet
```

**Advantages:**
- Compressed binary format
- Fast read/write
- Preserves data types
- Efficient for large datasets

### JSON Format
```bash
sra-pipeline create-ml-table \
  --input-dir ./batch_results \
  --output-file ml_features.json \
  --format json
```

**Advantages:**
- Human-readable
- Hierarchical structure
- Good for web applications

## 🔍 Working with Feature Tables

### Python Example

```python
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier

# Load feature table
df = pd.read_csv('ml_features.csv')

# Basic exploration
print(f"Dataset shape: {df.shape}")
print(f"Features: {df.columns.tolist()}")
print(f"Sample IDs: {df['sample_id'].tolist()}")

# Check for missing values
missing_values = df.isnull().sum()
print("Missing values:")
print(missing_values[missing_values > 0])

# Basic statistics
print("\nFeature statistics:")
print(df.describe())

# Prepare features for ML (exclude non-feature columns)
feature_columns = [col for col in df.columns 
                  if col not in ['sample_id', 'sample_name', 'pipeline_version']]
X = df[feature_columns]

# If you have labels (example)
# y = df['label']  # Replace with your actual label column
# X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)
```

### R Example

```r
library(tidyverse)
library(caret)

# Load feature table
df <- read_csv("ml_features.csv")

# Basic exploration
print(paste("Dataset shape:", nrow(df), "x", ncol(df)))
print("Features:")
print(colnames(df))

# Check for missing values
missing_summary <- df %>% 
  summarise_all(~sum(is.na(.))) %>%
  gather(key = "feature", value = "missing_count") %>%
  filter(missing_count > 0)

print("Missing values:")
print(missing_summary)

# Basic statistics
print("Feature statistics:")
print(summary(df))

# Prepare features for ML
feature_columns <- setdiff(colnames(df), 
                          c("sample_id", "sample_name", "pipeline_version"))
X <- df[, feature_columns]

# If you have labels (example)
# y <- df$label  # Replace with your actual label column
# train_index <- createDataPartition(y, p = 0.8, list = FALSE)
# X_train <- X[train_index, ]
# X_test <- X[-train_index, ]
```

## 🎯 ML Use Cases

### Classification Tasks
- **Disease Classification**: Predict disease status from genomic features
- **Tissue Classification**: Classify samples by tissue type
- **Quality Classification**: Identify high/low quality samples

### Regression Tasks
- **Age Prediction**: Predict biological age from genomic features
- **Expression Prediction**: Predict gene expression levels
- **Survival Prediction**: Predict survival outcomes

### Clustering Tasks
- **Sample Clustering**: Group similar samples
- **Feature Clustering**: Identify correlated features
- **Population Structure**: Identify genetic populations

## 🔧 Feature Engineering

### Creating Derived Features

```python
import pandas as pd
import numpy as np

# Load feature table
df = pd.read_csv('ml_features.csv')

# Create derived features
df['mapping_efficiency'] = df['mapped_reads'] / df['total_reads']
df['coverage_cv'] = df['coverage_std'] / df['mean_coverage']  # Coefficient of variation
df['variant_rate'] = df['total_genomic_variants'] / df['total_reads']
df['gene_variant_density'] = df['total_gene_variants'] / df['total_genes_with_variants']

# Create interaction features
df['quality_variant_score'] = df['mapping_rate'] * df['mean_coverage'] / 100

# Create categorical features
df['coverage_category'] = pd.cut(df['mean_coverage'], 
                                bins=[0, 10, 30, 100, np.inf], 
                                labels=['Low', 'Medium', 'High', 'Very High'])

# Create binary features
df['has_cnv'] = (df['total_cnv_regions'] > 0).astype(int)
df['high_variant_burden'] = (df['total_genomic_variants'] > df['total_genomic_variants'].median()).astype(int)
```

### Feature Selection

```python
from sklearn.feature_selection import SelectKBest, f_classif, mutual_info_classif
from sklearn.ensemble import RandomForestClassifier

# Load and prepare data
df = pd.read_csv('ml_features.csv')
feature_columns = [col for col in df.columns 
                  if col not in ['sample_id', 'sample_name', 'pipeline_version']]
X = df[feature_columns]

# Remove features with too many missing values
missing_threshold = 0.5
missing_ratio = X.isnull().sum() / len(X)
X_clean = X.loc[:, missing_ratio < missing_threshold]

# Fill remaining missing values
X_filled = X_clean.fillna(X_clean.median())

# Feature selection methods
# 1. Univariate selection
selector = SelectKBest(score_func=f_classif, k=20)
X_selected = selector.fit_transform(X_filled, y)  # y is your target variable
selected_features = X_filled.columns[selector.get_support()]

# 2. Feature importance from Random Forest
rf = RandomForestClassifier(n_estimators=100, random_state=42)
rf.fit(X_filled, y)
feature_importance = pd.DataFrame({
    'feature': X_filled.columns,
    'importance': rf.feature_importances_
}).sort_values('importance', ascending=False)

# Select top features
top_features = feature_importance.head(20)['feature'].tolist()
X_top = X_filled[top_features]
```

## 📈 Performance Monitoring

### Feature Table Quality Metrics

```python
def assess_feature_table_quality(df):
    """Assess the quality of a feature table."""
    
    quality_report = {
        'sample_count': len(df),
        'feature_count': len(df.columns),
        'missing_values': df.isnull().sum().sum(),
        'missing_percentage': (df.isnull().sum().sum() / (len(df) * len(df.columns))) * 100,
        'duplicate_samples': df['sample_id'].duplicated().sum(),
        'constant_features': len([col for col in df.columns if df[col].nunique() == 1]),
        'high_correlation_pairs': 0,  # Would need to calculate
    }
    
    # Calculate correlations
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    if len(numeric_cols) > 1:
        corr_matrix = df[numeric_cols].corr()
        high_corr_pairs = np.where(np.abs(corr_matrix) > 0.95)
        quality_report['high_correlation_pairs'] = len(set(zip(high_corr_pairs[0], high_corr_pairs[1]))) // 2
    
    return quality_report

# Assess quality
quality_report = assess_feature_table_quality(df)
print("Feature Table Quality Report:")
for metric, value in quality_report.items():
    print(f"  {metric}: {value}")
```

## 🐛 Troubleshooting

### Common Issues

**"No feature sets found"**
```bash
# Check if features.json files exist
find ./batch_results -name "features.json"

# Ensure pipeline completed successfully
ls -la ./batch_results/*/features.json
```

**"Missing values in feature table"**
```python
# Check which features have missing values
missing_summary = df.isnull().sum()
print(missing_summary[missing_summary > 0])

# Fill missing values
df_filled = df.fillna(df.median())  # For numeric features
```

**"Feature table is empty"**
```bash
# Check if pipeline generated features
ls -la ./batch_results/*/features.json

# Verify feature extraction completed
grep -r "Feature extraction completed" ./batch_results/logs/
```

### Performance Optimization

**For Large Datasets:**
```bash
# Use Parquet format for better performance
sra-pipeline create-ml-table \
  --input-dir ./batch_results \
  --output-file ml_features.parquet \
  --format parquet

# Use parallel processing for feature extraction
sra-pipeline batch --sra-ids $(cat sra_list.txt) --threads 16
```

## 📚 Additional Resources

- [Pandas Documentation](https://pandas.pydata.org/docs/)
- [Scikit-learn Feature Selection](https://scikit-learn.org/stable/modules/feature_selection.html)
- [Feature Engineering Best Practices](https://www.kaggle.com/learn/feature-engineering)
- [Genomic Feature Engineering](https://www.nature.com/articles/s41588-019-0418-7) 