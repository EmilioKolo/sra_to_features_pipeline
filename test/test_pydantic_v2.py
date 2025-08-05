#!/usr/bin/env python3
"""
Tests for Pydantic v2 compatibility.
"""

import pytest
from pathlib import Path
import structlog

from src.sra_pipeline.config.settings import PipelineConfig
from src.sra_pipeline.models.features import (
    FeatureSet, 
    QualityMetrics, 
    FragmentLengthStats,
    GenomicBin,
    GeneVariantStats,
    CNVRegion
)


def test_pipeline_config_creation():
    """Test that PipelineConfig can be created with Pydantic v2."""
    # Test with minimal required fields
    config = PipelineConfig(
        reference_fasta=Path("/tmp/test.fasta"),
        reference_gff=Path("/tmp/test.gff"),
        bed_genes=Path("/tmp/genes.bed"),
        genome_sizes=Path("/tmp/genome.sizes"),
        kraken_db=Path("/tmp/kraken"),
        snpeff_dir=Path("/tmp/snpeff")
    )
    
    assert config.reference_fasta == Path("/tmp/test.fasta")
    assert config.threads == 1  # Default value
    assert config.bin_size_gvs == 100000  # Default value


def test_pipeline_config_validation():
    """Test that PipelineConfig validation works correctly."""
    # Test invalid thread count
    with pytest.raises(ValueError, match="Threads must be positive"):
        PipelineConfig(
            reference_fasta=Path("/tmp/test.fasta"),
            reference_gff=Path("/tmp/test.gff"),
            bed_genes=Path("/tmp/genes.bed"),
            genome_sizes=Path("/tmp/genome.sizes"),
            kraken_db=Path("/tmp/kraken"),
            snpeff_dir=Path("/tmp/snpeff"),
            threads=0
        )


def test_feature_models_creation():
    """Test that feature models can be created with Pydantic v2."""
    # Test QualityMetrics
    qm = QualityMetrics(
        total_reads=1000000,
        mapped_reads=950000,
        mapping_rate=95.0,
        mean_coverage=30.0,
        coverage_std=5.0,
        gc_content=45.2,
        duplication_rate=12.5,
        quality_scores={"q10": 95.5, "q20": 90.2, "q30": 85.1}
    )
    
    assert qm.total_reads == 1000000
    assert qm.mapping_rate == 95.0
    
    # Test FragmentLengthStats
    fs = FragmentLengthStats(
        mean=150.0,
        median=150.0,
        std=25.0,
        min=50.0,
        max=300.0,
        count=1000000
    )
    
    assert fs.mean == 150.0
    assert fs.count == 1000000
    
    # Test GenomicBin
    gb = GenomicBin(
        chromosome="chr1",
        start=0,
        end=100000,
        variant_count=10
    )
    
    assert gb.chromosome == "chr1"
    assert gb.variant_count == 10
    
    # Test GeneVariantStats
    gvs = GeneVariantStats(
        gene_name="TEST_GENE",
        chromosome="chr1",
        start=1000,
        end=5000,
        total_variants=50,
        synonymous_variants=30,
        nonsynonymous_variants=20,
        dn_ds_ratio=0.67
    )
    
    assert gvs.gene_name == "TEST_GENE"
    assert gvs.dn_ds_ratio == 0.67
    
    # Test CNVRegion
    cnv = CNVRegion(
        chromosome="chr1",
        start=10000,
        end=20000,
        copy_number=3.0,
        confidence=0.95,
        type="gain"
    )
    
    assert cnv.chromosome == "chr1"
    assert cnv.copy_number == 3.0
    assert cnv.type == "gain"


def test_feature_set_creation():
    """Test that FeatureSet can be created with Pydantic v2."""
    qm = QualityMetrics(
        total_reads=1000000,
        mapped_reads=950000,
        mapping_rate=95.0,
        mean_coverage=30.0,
        coverage_std=5.0,
        gc_content=45.2,
        duplication_rate=12.5,
        quality_scores={"q10": 95.5, "q20": 90.2, "q30": 85.1}
    )
    
    fs = FragmentLengthStats(
        mean=150.0,
        median=150.0,
        std=25.0,
        min=50.0,
        max=300.0,
        count=1000000
    )
    
    feature_set = FeatureSet(
        sra_id="SRR123456",
        sample_name="test_sample",
        fragment_stats=fs,
        genomic_bins=[],
        gene_stats=[],
        cnv_regions=[],
        quality_metrics=qm,
        metadata={"test": "data"},
        processing_time=3600.5,
        pipeline_version="1.0.0"
    )
    
    assert feature_set.sra_id == "SRR123456"
    assert feature_set.quality_metrics.mapping_rate == 95.0
    assert feature_set.fragment_stats.mean == 150.0


def test_feature_set_serialization():
    """Test that FeatureSet serialization works with Pydantic v2."""
    qm = QualityMetrics(
        total_reads=1000000,
        mapped_reads=950000,
        mapping_rate=95.0,
        mean_coverage=30.0,
        coverage_std=5.0,
        gc_content=45.2,
        duplication_rate=12.5,
        quality_scores={"q10": 95.5, "q20": 90.2, "q30": 85.1}
    )
    
    feature_set = FeatureSet(
        sra_id="SRR123456",
        sample_name="test_sample",
        fragment_stats=None,
        genomic_bins=[],
        gene_stats=[],
        cnv_regions=[],
        quality_metrics=qm,
        metadata={"test": "data"},
        processing_time=3600.5,
        pipeline_version="1.0.0"
    )
    
    # Test to_dict method
    data_dict = feature_set.to_dict()
    assert isinstance(data_dict, dict)
    assert data_dict["sra_id"] == "SRR123456"
    assert data_dict["quality_metrics"]["mapping_rate"] == 95.0
    
    # Test to_json method
    json_str = feature_set.to_json()
    assert isinstance(json_str, str)
    assert "SRR123456" in json_str
    assert "95.0" in json_str


def test_validation_errors():
    """Test that validation errors work correctly with Pydantic v2."""
    # Test negative values
    with pytest.raises(ValueError, match="Fragment length statistics must be non-negative"):
        FragmentLengthStats(
            mean=-150.0,
            median=150.0,
            std=25.0,
            min=50.0,
            max=300.0,
            count=1000000
        )
    
    # Test invalid CNV type
    with pytest.raises(ValueError, match="CNV type must be 'gain' or 'loss'"):
        CNVRegion(
            chromosome="chr1",
            start=10000,
            end=20000,
            copy_number=3.0,
            confidence=0.95,
            type="invalid"
        )
    
    # Test invalid confidence
    with pytest.raises(ValueError, match="Confidence must be between 0 and 1"):
        CNVRegion(
            chromosome="chr1",
            start=10000,
            end=20000,
            copy_number=3.0,
            confidence=1.5,
            type="gain"
        )


if __name__ == "__main__":
    pytest.main([__file__]) 