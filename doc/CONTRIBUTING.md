# Contributing to SRA to Features Pipeline

Thank you for your interest in contributing to the SRA to Features Pipeline! This document provides guidelines and information for contributors.

## 🤝 How to Contribute

### Types of Contributions

We welcome various types of contributions:

- **Bug Reports**: Report bugs and issues
- **Feature Requests**: Suggest new features or improvements
- **Code Contributions**: Submit code changes and improvements
- **Documentation**: Improve or add documentation
- **Testing**: Add tests or improve test coverage
- **Performance**: Optimize code performance

### Getting Started

1. **Fork the Repository**: Click the "Fork" button on GitHub
2. **Clone Your Fork**: `git clone https://github.com/your-username/sra-to-features-pipeline.git`
3. **Create a Branch**: `git checkout -b feature/your-feature-name`
4. **Make Changes**: Implement your changes
5. **Test Your Changes**: Ensure all tests pass
6. **Submit a Pull Request**: Create a PR with a clear description

## 🛠️ Development Setup

### Prerequisites

- Python 3.8 or higher
- Git
- pip or conda

### Installation

```bash
# Clone your fork
git clone https://github.com/your-username/sra-to-features-pipeline.git
cd sra-to-features-pipeline

# Install in development mode
pip install -e ".[dev]"

# Install pre-commit hooks
pre-commit install
```

### Development Tools

The project uses several development tools:

- **Black**: Code formatting
- **Flake8**: Linting
- **MyPy**: Type checking
- **Pytest**: Testing
- **Pre-commit**: Git hooks

## 📝 Code Style

### Python Style Guide

We follow PEP 8 with some modifications:

- **Line Length**: 88 characters (Black default)
- **Type Hints**: Required for all functions
- **Docstrings**: Use Google-style docstrings
- **Imports**: Grouped and sorted (isort)

### Code Formatting

```bash
# Format code
black src/ test/

# Sort imports
isort src/ test/

# Check formatting
black --check src/ test/
```

### Linting

```bash
# Run linter
flake8 src/ test/

# Run type checker
mypy src/
```

## 🧪 Testing

### Running Tests

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=src/sra_pipeline

# Run specific test categories
pytest -m unit
pytest -m integration
pytest -m "not slow"

# Run tests in parallel
pytest -n auto
```

### Writing Tests

- Write tests for all new functionality
- Use descriptive test names
- Follow the Arrange-Act-Assert pattern
- Mock external dependencies
- Test both success and failure cases

### Test Structure

```python
def test_function_name():
    """Test description."""
    # Arrange
    input_data = "test"
    
    # Act
    result = function(input_data)
    
    # Assert
    assert result == "expected"
```

## 📚 Documentation

### Docstring Format

Use Google-style docstrings:

```python
def function_name(param1: str, param2: int) -> bool:
    """Short description of function.
    
    Longer description if needed.
    
    Args:
        param1: Description of param1
        param2: Description of param2
        
    Returns:
        Description of return value
        
    Raises:
        ValueError: When something goes wrong
        
    Example:
        >>> function_name("test", 42)
        True
    """
    pass
```

### API Documentation

- Document all public functions and classes
- Include type hints
- Provide examples where helpful
- Update documentation when changing APIs

## 🔄 Pull Request Process

### Before Submitting

1. **Ensure Tests Pass**: Run the full test suite
2. **Check Code Style**: Run formatting and linting tools
3. **Update Documentation**: Update docstrings and README if needed
4. **Add Tests**: Include tests for new functionality
5. **Update CHANGELOG**: Document your changes

### Pull Request Template

```markdown
## Description
Brief description of changes

## Type of Change
- [ ] Bug fix
- [ ] New feature
- [ ] Breaking change
- [ ] Documentation update

## Testing
- [ ] Unit tests pass
- [ ] Integration tests pass
- [ ] Manual testing completed

## Checklist
- [ ] Code follows style guidelines
- [ ] Self-review completed
- [ ] Documentation updated
- [ ] Tests added/updated
- [ ] CHANGELOG updated
```

### Review Process

1. **Automated Checks**: CI/CD pipeline runs tests and checks
2. **Code Review**: At least one maintainer reviews the PR
3. **Address Feedback**: Make requested changes
4. **Merge**: PR is merged after approval

## 🐛 Bug Reports

### Bug Report Template

```markdown
## Bug Description
Clear description of the bug

## Steps to Reproduce
1. Step 1
2. Step 2
3. Step 3

## Expected Behavior
What should happen

## Actual Behavior
What actually happens

## Environment
- OS: [e.g., Ubuntu 20.04]
- Python Version: [e.g., 3.9.7]
- Pipeline Version: [e.g., 1.0.0]

## Additional Information
Screenshots, logs, etc.
```

## 💡 Feature Requests

### Feature Request Template

```markdown
## Feature Description
Clear description of the feature

## Use Case
Why this feature is needed

## Proposed Solution
How you think it should work

## Alternatives Considered
Other approaches you considered

## Additional Information
Any other relevant information
```

## 📋 Issue Labels

We use the following labels to categorize issues:

- **bug**: Something isn't working
- **enhancement**: New feature or request
- **documentation**: Improvements or additions to documentation
- **good first issue**: Good for newcomers
- **help wanted**: Extra attention is needed
- **question**: Further information is requested
- **wontfix**: This will not be worked on

## 🏷️ Versioning

We follow [Semantic Versioning](https://semver.org/):

- **MAJOR**: Incompatible API changes
- **MINOR**: New functionality in a backwards-compatible manner
- **PATCH**: Backwards-compatible bug fixes

## 📞 Getting Help

### Communication Channels

- **GitHub Issues**: For bug reports and feature requests
- **GitHub Discussions**: For questions and general discussion
- **Email**: For sensitive issues or private discussions

### Code of Conduct

We are committed to providing a welcoming and inspiring community for all. Please read our [Code of Conduct](CODE_OF_CONDUCT.md) for details.

## 🙏 Recognition

Contributors will be recognized in:

- **README.md**: For significant contributions
- **CHANGELOG.md**: For all contributions
- **GitHub Contributors**: Automatic recognition

## 📄 License

By contributing, you agree that your contributions will be licensed under the same license as the project (MIT License).

---

Thank you for contributing to the SRA to Features Pipeline! 🎉 