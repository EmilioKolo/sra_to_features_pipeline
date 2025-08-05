#!/usr/bin/env python3
"""
Tests for import functionality.
"""

import pytest
import sys
from pathlib import Path

# Add src to path for testing
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))


def test_feature_extraction_imports():
    """Test that feature_extraction.py can be imported correctly."""
    try:
        from sra_pipeline.core import feature_extraction
        assert feature_extraction is not None
        print("✅ feature_extraction imported successfully")
    except ImportError as e:
        pytest.fail(f"Failed to import feature_extraction: {e}")


def test_models_imports():
    """Test that models can be imported correctly."""
    try:
        from sra_pipeline.models.features import (
            FeatureSet,
            QualityMetrics,
            FragmentLengthStats,
            GenomicBin,
            GeneVariantStats,
            CNVRegion
        )
        assert FeatureSet is not None
        assert QualityMetrics is not None
        assert FragmentLengthStats is not None
        assert GenomicBin is not None
        assert GeneVariantStats is not None
        assert CNVRegion is not None
        print("✅ All feature models imported successfully")
    except ImportError as e:
        pytest.fail(f"Failed to import feature models: {e}")


def test_config_imports():
    """Test that config can be imported correctly."""
    try:
        from sra_pipeline.config.settings import PipelineConfig
        assert PipelineConfig is not None
        print("✅ PipelineConfig imported successfully")
    except ImportError as e:
        pytest.fail(f"Failed to import PipelineConfig: {e}")


def test_core_imports():
    """Test that core modules can be imported correctly."""
    try:
        from sra_pipeline.core import (
            pipeline,
            download,
            alignment,
            variant_calling,
            feature_extraction,
            quality_control
        )
        assert pipeline is not None
        assert download is not None
        assert alignment is not None
        assert variant_calling is not None
        assert feature_extraction is not None
        assert quality_control is not None
        print("✅ All core modules imported successfully")
    except ImportError as e:
        pytest.fail(f"Failed to import core modules: {e}")


def test_utils_imports():
    """Test that utils can be imported correctly."""
    try:
        from sra_pipeline.utils import (
            setup_logging,
            PipelineLogger,
            PerformanceMonitor,
            log_command,
            log_error
        )
        assert setup_logging is not None
        assert PipelineLogger is not None
        assert PerformanceMonitor is not None
        assert log_command is not None
        assert log_error is not None
        print("✅ All utils imported successfully")
    except ImportError as e:
        pytest.fail(f"Failed to import utils: {e}")


def test_ml_features_imports():
    """Test that ML features can be imported correctly."""
    try:
        from sra_pipeline.utils.ml_features import (
            MLFeatureTable,
            create_ml_feature_table_from_directory
        )
        assert MLFeatureTable is not None
        assert create_ml_feature_table_from_directory is not None
        print("✅ ML features imported successfully")
    except ImportError as e:
        pytest.fail(f"Failed to import ML features: {e}")


def test_cli_imports():
    """Test that CLI can be imported correctly."""
    try:
        from sra_pipeline.cli import main
        assert main is not None
        print("✅ CLI imported successfully")
    except ImportError as e:
        pytest.fail(f"Failed to import CLI: {e}")


if __name__ == "__main__":
    print("🧪 Testing all imports...")
    
    test_feature_extraction_imports()
    test_models_imports()
    test_config_imports()
    test_core_imports()
    test_utils_imports()
    test_ml_features_imports()
    test_cli_imports()
    
    print("🎉 All import tests passed!") 