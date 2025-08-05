# Changelog

All notable changes to the SRA to Features Pipeline project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Initial project structure and architecture
- Modern Python packaging with pyproject.toml
- Comprehensive configuration management with Pydantic
- Structured logging with structlog
- Performance monitoring and metrics
- Type hints throughout the codebase
- Comprehensive test framework
- Modern CLI interface with Click and Rich
- Data models with Pydantic validation
- Error handling and retry logic
- Quality control framework
- Documentation generation

### Changed
- Complete restructuring from original monolithic code
- Modernized code architecture with proper separation of concerns
- Improved error handling and logging
- Enhanced configuration management
- Better testing infrastructure

### Deprecated
- Original Docker-based approach (replaced with native installation)
- Old configuration system (replaced with Pydantic-based config)
- Legacy logging system (replaced with structured logging)

### Removed
- Docker dependency for basic usage
- Hardcoded paths and configurations
- Monolithic script structure

### Fixed
- Configuration validation issues
- Error handling gaps
- Logging inconsistencies
- Type safety issues

### Security
- Improved input validation
- Better error message handling
- Secure configuration management

## [1.0.0] - 2024-01-XX

### Added
- Initial release of the restructured pipeline
- Core pipeline functionality:
  - SRA data download with retry logic
  - BWA alignment with performance monitoring
  - Variant calling with bcftools
  - Feature extraction framework
  - Quality control assessment
- Modern CLI interface with rich output
- Comprehensive configuration system
- Structured logging and performance monitoring
- Type-safe data models
- Extensive test coverage
- Complete documentation

### Features
- **SRA Data Processing**: Download and process SRA data with robust error handling
- **Sequence Alignment**: BWA-based alignment with quality metrics
- **Variant Calling**: bcftools-based variant detection and filtering
- **Feature Extraction**: Genomic feature extraction for ML training
- **Quality Control**: Comprehensive quality assessment and validation
- **Performance Monitoring**: Real-time performance tracking and resource usage
- **Configuration Management**: Flexible configuration via environment variables or files
- **Error Handling**: Robust error handling with retry logic and graceful failures

### Technical Improvements
- **Modern Python**: Python 3.8+ with type hints and modern packaging
- **Data Validation**: Pydantic-based data validation and serialization
- **Structured Logging**: JSON-formatted logs with performance metrics
- **Testing**: Comprehensive unit and integration test suite
- **Documentation**: Complete API documentation and user guides
- **CLI**: Modern command-line interface with rich output and progress tracking

## [0.9.0] - 2023-XX-XX

### Added
- Beta release with core functionality
- Basic pipeline implementation
- Docker containerization
- Simple configuration system

### Changed
- Initial pipeline structure
- Basic error handling
- Simple logging system

## [0.8.0] - 2023-XX-XX

### Added
- Alpha release with basic pipeline
- Initial codebase structure
- Basic SRA processing functionality

---

## Version History Summary

### Major Versions

- **v1.0.0**: Complete restructuring with modern architecture
- **v0.9.0**: Beta release with core functionality
- **v0.8.0**: Alpha release with basic pipeline

### Key Architectural Changes

#### v1.0.0 (Current)
- **Modern Python Architecture**: Type hints, Pydantic, structured logging
- **Modular Design**: Separated concerns with clear module boundaries
- **Comprehensive Testing**: Unit and integration test coverage
- **Production Ready**: Error handling, monitoring, and documentation
- **Native Installation**: No Docker dependency for basic usage

#### v0.9.0 (Previous)
- **Docker-based**: Containerized approach
- **Monolithic**: Single script structure
- **Basic Configuration**: Simple config file approach
- **Limited Testing**: Basic test coverage

#### v0.8.0 (Original)
- **Prototype**: Initial implementation
- **Basic Functionality**: Core pipeline steps
- **Minimal Documentation**: Basic usage instructions

---

## Migration Guide

### From v0.9.0 to v1.0.0

#### Breaking Changes
- Configuration format changed from INI to environment variables/Pydantic
- CLI interface completely redesigned
- Output format changed to structured JSON
- Docker no longer required for basic usage

#### Migration Steps
1. **Installation**: Use `pip install sra-to-features-pipeline` instead of Docker
2. **Configuration**: Convert INI config to environment variables or .env file
3. **CLI**: Use new `sra-pipeline` command instead of `python main.py`
4. **Output**: Check new JSON output format in `features.json` files

#### Configuration Migration
```bash
# Old (v0.9.0)
# config.ini
[Paths]
BASE_DIR = /content

# New (v1.0.0)
# .env file
SRA_PIPELINE_BASE_DIR=/content
```

#### CLI Migration
```bash
# Old (v0.9.0)
python main.py --sra-id SRR123456 --output-dir ./output

# New (v1.0.0)
sra-pipeline run --sra-id SRR123456 --output-dir ./output
```

---

## Contributing

To add entries to this changelog:

1. **For Unreleased Changes**: Add entries under the `[Unreleased]` section
2. **For New Versions**: Move `[Unreleased]` entries to a new version section
3. **Format**: Follow the existing format with proper categories
4. **Descriptions**: Be clear and concise about what changed

### Categories
- **Added**: New features
- **Changed**: Changes in existing functionality
- **Deprecated**: Soon-to-be removed features
- **Removed**: Removed features
- **Fixed**: Bug fixes
- **Security**: Security-related changes 