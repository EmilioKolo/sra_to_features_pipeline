#!/usr/bin/env python3
"""
Tests for ML feature table functionality.
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import tempfile
import json
import structlog

from src.sra_pipeline.utils.ml_features import MLFeatureTable, create_ml_feature_table_from_directory
from src.sra_pipeline.models.features import FeatureSet, QualityMetrics, FragmentLengthStats


@pytest.fixture
def sample_feature_set():
    """Create a sample feature set for testing."""
    quality_metrics = QualityMetrics(
        total_reads=1000000,
        mapped_reads=950000,
        mapping_rate=95.0,
        mean_coverage=30.0,
        coverage_std=5.0,
        gc_content=45.2,
        duplication_rate=12.5,
        quality_scores={"q10": 95.5, "q20": 90.2, "q30": 85.1}
    )
    
    fragment_stats = FragmentLengthStats(
        mean=150.0,
        median=150.0,
        std=25.0,
        min=50.0,
        max=300.0,
        count=1000000
    )
    
    return FeatureSet(
        sra_id="SRR123456",
        sample_name="test_sample",
        fragment_stats=fragment_stats,
        genomic_bins=[],
        gene_stats=[],
        cnv_regions=[],
        quality_metrics=quality_metrics,
        metadata={"test": "data"},
        processing_time=3600.5,
        pipeline_version="1.0.0"
    )


@pytest.fixture
def logger():
    """Create a test logger."""
    return structlog.get_logger()


def test_ml_feature_table_initialization(logger):
    """Test MLFeatureTable initialization."""
    ml_table = MLFeatureTable(logger)
    assert len(ml_table.feature_sets) == 0


def test_add_feature_set(logger, sample_feature_set):
    """Test adding feature sets."""
    ml_table = MLFeatureTable(logger)
    ml_table.add_feature_set(sample_feature_set)
    
    assert len(ml_table.feature_sets) == 1
    assert ml_table.feature_sets[0].sra_id == "SRR123456"


def test_create_sample_features_table(logger, sample_feature_set):
    """Test creating sample features table."""
    ml_table = MLFeatureTable(logger)
    ml_table.add_feature_set(sample_feature_set)
    
    df = ml_table.create_sample_features_table()
    
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 1  # One sample
    assert df.iloc[0]['sample_id'] == "SRR123456"
    assert df.iloc[0]['total_reads'] == 1000000
    assert df.iloc[0]['mapping_rate'] == 95.0
    assert df.iloc[0]['fragment_mean'] == 150.0


def test_save_feature_table_csv(logger, sample_feature_set, tmp_path):
    """Test saving feature table as CSV."""
    ml_table = MLFeatureTable(logger)
    ml_table.add_feature_set(sample_feature_set)
    
    output_file = tmp_path / "test_features.csv"
    ml_table.save_feature_table(output_file, format='csv')
    
    assert output_file.exists()
    
    # Verify the file can be read back
    df = pd.read_csv(output_file)
    assert len(df) == 1
    assert df.iloc[0]['sample_id'] == "SRR123456"


def test_save_feature_table_parquet(logger, sample_feature_set, tmp_path):
    """Test saving feature table as Parquet."""
    ml_table = MLFeatureTable(logger)
    ml_table.add_feature_set(sample_feature_set)
    
    output_file = tmp_path / "test_features.parquet"
    ml_table.save_feature_table(output_file, format='parquet')
    
    assert output_file.exists()
    
    # Verify the file can be read back
    df = pd.read_parquet(output_file)
    assert len(df) == 1
    assert df.iloc[0]['sample_id'] == "SRR123456"


def test_get_feature_summary(logger, sample_feature_set):
    """Test getting feature summary."""
    ml_table = MLFeatureTable(logger)
    ml_table.add_feature_set(sample_feature_set)
    
    summary = ml_table.get_feature_summary()
    
    assert summary['sample_count'] == 1
    assert summary['feature_count'] > 0
    assert summary['missing_values'] >= 0
    assert 'sample_id' in summary['feature_names']


def test_create_ml_feature_table_from_directory(logger, sample_feature_set, tmp_path):
    """Test creating ML feature table from directory."""
    # Create a temporary directory structure
    sample_dir = tmp_path / "SRR123456"
    sample_dir.mkdir()
    
    # Save feature set as JSON
    feature_file = sample_dir / "features.json"
    with open(feature_file, 'w') as f:
        json.dump(sample_feature_set.dict(), f, indent=2)
    
    # Create ML feature table
    output_file = tmp_path / "ml_features.csv"
    df = create_ml_feature_table_from_directory(
        input_directory=tmp_path,
        output_path=output_file,
        format='csv',
        logger=logger
    )
    
    assert output_file.exists()
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 1
    assert df.iloc[0]['sample_id'] == "SRR123456"


def test_empty_feature_sets(logger):
    """Test handling empty feature sets."""
    ml_table = MLFeatureTable(logger)
    
    with pytest.raises(ValueError, match="No feature sets available"):
        ml_table.create_sample_features_table()


def test_invalid_format(logger, sample_feature_set, tmp_path):
    """Test handling invalid format."""
    ml_table = MLFeatureTable(logger)
    ml_table.add_feature_set(sample_feature_set)
    
    output_file = tmp_path / "test_features.txt"
    
    with pytest.raises(ValueError, match="Unsupported format"):
        ml_table.save_feature_table(output_file, format='invalid')


if __name__ == "__main__":
    pytest.main([__file__]) 