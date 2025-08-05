# Troubleshooting Guide

This guide helps you identify and resolve common issues when using the SRA to Features Pipeline.

## 🔍 Quick Diagnosis

### Check Pipeline Status

```bash
# Check if pipeline is installed
which sra-pipeline

# Check pipeline version
sra-pipeline --version

# Validate your setup
sra-pipeline validate
```

### Check System Resources

```bash
# Check available memory
free -h

# Check disk space
df -h

# Check CPU cores
nproc

# Check if tools are available
which bwa samtools bcftools bedtools fastqc fastq-dump
```

## 🐛 Common Issues and Solutions

### Installation Issues

#### "Command not found: sra-pipeline"

**Problem**: The pipeline is not installed or not in your PATH.

**Solutions**:
```bash
# Reinstall the pipeline
pip install -e .

# Check installation
pip list | grep sra-pipeline

# Add to PATH if needed
export PATH="$HOME/.local/bin:$PATH"
```

#### "Permission denied" during installation

**Problem**: Insufficient permissions to install packages.

**Solutions**:
```bash
# Install with user permissions
pip install --user sra-to-features-pipeline

# Or use a virtual environment
python -m venv venv
source venv/bin/activate
pip install sra-to-features-pipeline
```

### Tool Installation Issues

#### "Tool not found" errors

**Problem**: Required bioinformatics tools are not installed.

**Solutions**:
```bash
# Run the automated installation script
./scripts/install/install_tools.sh

# Or install manually
# Ubuntu/Debian
sudo apt update && sudo apt install -y bwa samtools bcftools bedtools fastqc

# macOS
brew install bwa samtools bcftools bedtools fastqc

# Check installation
which bwa samtools bcftools bedtools fastqc
```

#### "Java not found" error

**Problem**: Java is required for some tools but not installed.

**Solutions**:
```bash
# Ubuntu/Debian
sudo apt install -y default-jre

# macOS
brew install openjdk@8

# Check Java installation
java -version
```

### Input Data Issues

#### "SRA ID not found" error

**Problem**: The SRA ID is invalid or doesn't exist.

**Solutions**:
```bash
# Check SRA ID format
echo $SRA_ID | grep -E "^SR[RXSP][0-9]+$"

# Verify SRA ID exists
curl "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=run_accession=$SRA_ID"

# Try a known valid SRA ID for testing
sra-pipeline run --sra-id SRR123456 --dry-run
```

#### "FASTQ file not found" error

**Problem**: FASTQ files are missing or have incorrect paths.

**Solutions**:
```bash
# Check if files exist
ls -la *.fastq.gz

# Check file permissions
ls -la sample.fastq.gz

# Fix file permissions if needed
chmod 644 *.fastq.gz

# Use absolute paths
sra-pipeline run --fastq /absolute/path/to/sample.fastq.gz
```

#### "Invalid FASTQ format" error

**Problem**: FASTQ files are corrupted or in wrong format.

**Solutions**:
```bash
# Check file format
head -n 8 sample.fastq.gz | gunzip

# Check file integrity
gunzip -t sample.fastq.gz && echo "Valid" || echo "Corrupted"

# Count reads to verify format
zcat sample.fastq.gz | echo $((`wc -l`/4))
```

### Configuration Issues

#### "Configuration file not found" error

**Problem**: Configuration file is missing or has wrong path.

**Solutions**:
```bash
# Check if .env file exists
ls -la .env

# Create default configuration
cat > .env << EOF
SRA_PIPELINE_THREADS=8
SRA_PIPELINE_MIN_QUALITY_SCORE=20
SRA_PIPELINE_MIN_COVERAGE=10
EOF

# Use environment variables instead
export SRA_PIPELINE_THREADS=8
export SRA_PIPELINE_MIN_QUALITY_SCORE=20
```

#### "Invalid configuration" error

**Problem**: Configuration values are invalid or out of range.

**Solutions**:
```bash
# Validate configuration
sra-pipeline validate --config .env

# Check configuration values
cat .env

# Use reasonable defaults
export SRA_PIPELINE_THREADS=4
export SRA_PIPELINE_MIN_QUALITY_SCORE=20
export SRA_PIPELINE_MIN_COVERAGE=10
```

### Performance Issues

#### "Out of memory" error

**Problem**: Insufficient memory for processing.

**Solutions**:
```bash
# Reduce memory usage
export SRA_PIPELINE_MAX_MEMORY_GB=8

# Reduce thread count
sra-pipeline run --sra-id SRR123456 --threads 4

# Check available memory
free -h

# Close other applications to free memory
```

#### "Process killed" or "Killed" message

**Problem**: System killed the process due to memory pressure.

**Solutions**:
```bash
# Reduce memory usage
export SRA_PIPELINE_MAX_MEMORY_GB=4

# Use fewer threads
sra-pipeline run --sra-id SRR123456 --threads 2

# Increase swap space (Linux)
sudo fallocate -l 4G /swapfile
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile
```

#### "Timeout" errors

**Problem**: Operations are taking too long and timing out.

**Solutions**:
```bash
# Increase timeout
export SRA_PIPELINE_TIMEOUT_SECONDS=7200

# Use more threads for faster processing
sra-pipeline run --sra-id SRR123456 --threads 16

# Check system performance
top
iostat
```

### Network Issues

#### "Download failed" error

**Problem**: Network issues preventing SRA data download.

**Solutions**:
```bash
# Check internet connection
ping -c 3 ftp.ncbi.nlm.nih.gov

# Try alternative download method
prefetch SRR123456

# Use a different network if available
# Check firewall settings
```

#### "Connection timeout" error

**Problem**: Slow or unreliable network connection.

**Solutions**:
```bash
# Increase timeout
export SRA_PIPELINE_TIMEOUT_SECONDS=7200

# Try downloading manually
fastq-dump --gzip SRR123456

# Use a more stable network connection
```

### Output Issues

#### "Permission denied" when writing output

**Problem**: Insufficient permissions to write to output directory.

**Solutions**:
```bash
# Check directory permissions
ls -la output/

# Fix permissions
chmod 755 output/

# Use a different output directory
sra-pipeline run --sra-id SRR123456 --output-dir ./my_results

# Run with appropriate permissions
sudo chown $USER:$USER output/
```

#### "Disk space full" error

**Problem**: Insufficient disk space for output files.

**Solutions**:
```bash
# Check disk space
df -h

# Clean up temporary files
rm -rf tmp/ *.tmp

# Use a different disk with more space
sra-pipeline run --sra-id SRR123456 --output-dir /path/to/large/disk/results
```

#### "Invalid JSON" in output

**Problem**: Output file is corrupted or incomplete.

**Solutions**:
```bash
# Validate JSON format
python -m json.tool results/SRR123456/features.json > /dev/null

# Check for incomplete writes
tail -n 5 results/SRR123456/features.json

# Re-run the analysis
sra-pipeline run --sra-id SRR123456 --output-dir ./new_results
```

## 🔧 Debugging Techniques

### Enable Debug Logging

```bash
# Run with debug logging
sra-pipeline run --sra-id SRR123456 --log-level DEBUG --log-file debug.log

# Check debug logs
cat debug.log | grep ERROR
cat debug.log | grep WARNING
```

### Validate Setup

```bash
# Comprehensive validation
sra-pipeline validate --verbose

# Check individual components
which bwa samtools bcftools bedtools fastqc fastq-dump
bwa version
samtools version
bcftools version
```

### Dry Run Testing

```bash
# Test configuration without running
sra-pipeline run --sra-id SRR123456 --dry-run

# Test with different parameters
sra-pipeline run --sra-id SRR123456 --threads 4 --dry-run
```

## 📊 Performance Optimization

### Memory Optimization

```bash
# For memory-constrained systems
export SRA_PIPELINE_MAX_MEMORY_GB=8
sra-pipeline run --sra-id SRR123456 --threads 4
```

### Speed Optimization

```bash
# For speed-optimized runs
export SRA_PIPELINE_MAX_MEMORY_GB=32
sra-pipeline run --sra-id SRR123456 --threads 16
```

### Quality Optimization

```bash
# For high-quality results
export SRA_PIPELINE_MIN_QUALITY_SCORE=30
export SRA_PIPELINE_MIN_COVERAGE=20
sra-pipeline run --sra-id SRR123456
```

## 🆘 Getting Help

### Before Asking for Help

1. **Check this troubleshooting guide**
2. **Run validation**: `sra-pipeline validate`
3. **Check logs**: Look for error messages
4. **Try dry run**: `sra-pipeline run --sra-id SRR123456 --dry-run`
5. **Check system resources**: Memory, disk space, CPU

### When Reporting Issues

Include the following information:

```bash
# System information
uname -a
python --version
sra-pipeline --version

# Error details
sra-pipeline run --sra-id SRR123456 --log-level DEBUG 2>&1 | tee error.log

# Configuration
cat .env
env | grep SRA_PIPELINE

# System resources
free -h
df -h
nproc
```

### Useful Commands for Diagnosis

```bash
# Check all tools
which bwa samtools bcftools bedtools fastqc fastq-dump java

# Check versions
bwa version
samtools version
bcftools version
fastqc --version
fastq-dump --version
java -version

# Check disk space
df -h

# Check memory
free -h

# Check CPU
nproc
lscpu

# Check network
ping -c 3 ftp.ncbi.nlm.nih.gov
```

## 🔗 Related Topics

- **[Basic Usage](basic_usage.md)** - Getting started guide
- **[Configuration Guide](../configuration.md)** - Configuration options
- **[CLI Reference](cli_reference.md)** - Command-line interface
- **[Performance Optimization](performance.md)** - Performance tuning

---

**Still having issues?** Create an issue on GitHub with the information above. 