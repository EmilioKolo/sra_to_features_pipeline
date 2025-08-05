#!/bin/bash

# Development Installation Script for SRA to Features Pipeline
# This script handles the "build_editable" error and provides multiple installation options

set -e  # Exit on any error

echo "🚀 SRA to Features Pipeline - Development Installation"
echo "======================================================"

# Check if we're in the right directory
if [ ! -f "pyproject.toml" ] || [ ! -f "setup.py" ]; then
    echo "❌ Error: Please run this script from the pipeline root directory"
    exit 1
fi

# Function to try installation method
try_install() {
    local method=$1
    local command=$2
    
    echo "🔄 Trying installation method: $method"
    echo "Command: $command"
    
    if eval "$command"; then
        echo "✅ Successfully installed using method: $method"
        return 0
    else
        echo "❌ Failed to install using method: $method"
        return 1
    fi
}

# Check Python version
echo "🐍 Checking Python version..."
python_version=$(python3 --version 2>&1 | cut -d' ' -f2 | cut -d'.' -f1,2)
echo "Python version: $python_version"

# Check if pip is available
if ! command -v pip &> /dev/null; then
    echo "❌ Error: pip is not installed. Please install pip first."
    exit 1
fi

# Upgrade pip and setuptools
echo "⬆️  Upgrading pip and setuptools..."
pip install --upgrade pip setuptools wheel

# Try different installation methods
echo ""
echo "🔧 Attempting installation..."

# Method 1: Standard editable installation
if try_install "Standard editable" "pip install -e ."; then
    echo "🎉 Installation completed successfully!"
    exit 0
fi

# Method 2: Editable installation with no build isolation
if try_install "Editable with no build isolation" "pip install -e . --no-build-isolation"; then
    echo "🎉 Installation completed successfully!"
    exit 0
fi

# Method 3: Use setup.py directly
if try_install "setup.py develop" "python setup.py develop"; then
    echo "🎉 Installation completed successfully!"
    exit 0
fi

# Method 4: Regular installation (not editable)
if try_install "Regular installation" "pip install ."; then
    echo "🎉 Installation completed successfully!"
    echo "⚠️  Note: This is not an editable installation. Changes to source code won't be reflected."
    exit 0
fi

# Method 5: Try with conda if available
if command -v conda &> /dev/null; then
    echo "🐍 Conda detected, trying conda installation..."
    if try_install "Conda installation" "conda install -c conda-forge pip setuptools && pip install -e ."; then
        echo "🎉 Installation completed successfully!"
        exit 0
    fi
fi

# If all methods fail
echo ""
echo "❌ All installation methods failed!"
echo ""
echo "🔧 Manual troubleshooting steps:"
echo "1. Check your Python version (requires 3.8+)"
echo "2. Ensure you have write permissions to the current directory"
echo "3. Try creating a virtual environment:"
echo "   python3 -m venv venv"
echo "   source venv/bin/activate"
echo "   pip install -e ."
echo "4. Check if all dependencies are available:"
echo "   pip install numpy pandas pydantic click rich"
echo "5. Try installing dependencies first:"
echo "   pip install -r requirements.txt"
echo "   pip install -e ."
echo ""
echo "📚 For more help, check the installation documentation:"
echo "   doc/installation.md"

exit 1 