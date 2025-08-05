# SRA to Features Pipeline

A modern, production-ready bioinformatics pipeline for extracting genomic features from SRA (Sequence Read Archive) data for LLM training.

## 🚀 Quick Start

```bash
# Install bioinformatics tools
./scripts/install/install_tools.sh

# Install pipeline (handles build_editable errors automatically)
./install_dev.sh

# Run pipeline
sra-pipeline run --sra-id SRR123456 --output-dir ./results
```

## 📖 Documentation

- **[User Guide](doc/user_guide/)** - Complete usage instructions and examples
- **[API Reference](doc/api/)** - Detailed API documentation
- **[Contributing Guidelines](doc/CONTRIBUTING.md)** - How to contribute to the project
- **[Changelog](doc/CHANGELOG.md)** - Version history and changes
- **[Restructuring Summary](doc/RESTRUCTURING_SUMMARY.md)** - Detailed overview of the project architecture

## 🎯 Features

- **Modern Architecture**: Type-safe, modular Python package
- **Flexible Configuration**: Environment variables or configuration files
- **Rich CLI**: Progress tracking and comprehensive help
- **Quality Control**: Built-in quality assessment and validation
- **Performance Monitoring**: Real-time metrics and resource tracking
- **Comprehensive Testing**: Unit and integration test coverage
- **Structured Output**: JSON-based feature extraction results

## 🛠️ Installation

### Prerequisites

- Python 3.8+
- Bioinformatics tools (BWA, SAMtools, BCFtools, etc.)

### Quick Installation

```bash
# Clone repository
git clone https://github.com/your-org/sra-to-features-pipeline.git
cd sra-to-features-pipeline

# Install bioinformatics tools
./scripts/install/install_tools.sh

# Install pipeline (handles build_editable errors automatically)
./install_dev.sh
```

## 📊 Output

The pipeline generates structured JSON output containing:

- **Fragment Statistics**: Read length distributions
- **Genomic Bins**: Variant counts per genomic regions
- **Gene Statistics**: Gene-level variant analysis
- **CNV Regions**: Copy number variation detection
- **Quality Metrics**: Data quality assessment

## 🔧 Development

```bash
# Install development dependencies
pip install -e ".[dev]"

# Run tests
pytest

# Code quality
black src/ test/
flake8 src/ test/
mypy src/
```

## 📄 License

This project is licensed under the MIT License - see [LICENSE](doc/LICENSE) for details.

## 🤝 Contributing

We welcome contributions! Please see our [Contributing Guidelines](doc/CONTRIBUTING.md) for details.

## 🆘 Support

- **Documentation**: Check the [documentation](doc/) first
- **Issues**: Report bugs on [GitHub Issues](https://github.com/your-org/sra-to-features-pipeline/issues)
- **Discussions**: Join our [GitHub Discussions](https://github.com/your-org/sra-to-features-pipeline/discussions)

---

**For detailed documentation, please visit the [doc/](doc/) directory.**

