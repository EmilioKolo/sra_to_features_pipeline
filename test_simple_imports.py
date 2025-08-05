#!/usr/bin/env python3
"""
Simple test for import functionality.
"""

import sys
from pathlib import Path

# Add src to path for testing
sys.path.insert(0, str(Path(__file__).parent / "src"))

def test_imports():
    """Test that all modules can be imported correctly."""
    
    print("🧪 Testing imports...")
    
    try:
        # Test feature_extraction import
        from sra_pipeline.core import feature_extraction
        print("✅ feature_extraction imported successfully")
    except ImportError as e:
        print(f"❌ Failed to import feature_extraction: {e}")
        return False
    
    try:
        # Test models import
        from sra_pipeline.models.features import FeatureSet
        print("✅ FeatureSet imported successfully")
    except ImportError as e:
        print(f"❌ Failed to import FeatureSet: {e}")
        return False
    
    try:
        # Test config import
        from sra_pipeline.config.settings import PipelineConfig
        print("✅ PipelineConfig imported successfully")
    except ImportError as e:
        print(f"❌ Failed to import PipelineConfig: {e}")
        return False
    
    try:
        # Test core modules import
        from sra_pipeline.core import pipeline, download, alignment, variant_calling
        print("✅ Core modules imported successfully")
    except ImportError as e:
        print(f"❌ Failed to import core modules: {e}")
        return False
    
    try:
        # Test utils import
        from sra_pipeline.utils import setup_logging
        print("✅ Utils imported successfully")
    except ImportError as e:
        print(f"❌ Failed to import utils: {e}")
        return False
    
    try:
        # Test ML features import
        from sra_pipeline.utils.ml_features import MLFeatureTable
        print("✅ ML features imported successfully")
    except ImportError as e:
        print(f"❌ Failed to import ML features: {e}")
        return False
    
    try:
        # Test CLI import
        from sra_pipeline.cli import main
        print("✅ CLI imported successfully")
    except ImportError as e:
        print(f"❌ Failed to import CLI: {e}")
        return False
    
    print("🎉 All import tests passed!")
    return True

if __name__ == "__main__":
    success = test_imports()
    sys.exit(0 if success else 1) 