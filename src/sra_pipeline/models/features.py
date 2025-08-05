"""
Data models for features extracted by the SRA to Features Pipeline.
"""

from typing import Dict, List, Optional, Union
from pydantic import BaseModel, Field, validator
import numpy as np


class FragmentLengthStats(BaseModel):
    """Statistics for fragment lengths in paired-end sequencing."""
    
    mean: float = Field(description="Mean fragment length")
    median: float = Field(description="Median fragment length")
    std: float = Field(description="Standard deviation of fragment lengths")
    min: float = Field(description="Minimum fragment length")
    max: float = Field(description="Maximum fragment length")
    count: int = Field(description="Number of fragments analyzed")
    
    @validator('mean', 'median', 'std', 'min', 'max')
    def validate_positive(cls, v):
        """Validate that length statistics are positive."""
        if v < 0:
            raise ValueError("Fragment length statistics must be non-negative")
        return v
    
    @validator('count')
    def validate_count(cls, v):
        """Validate that count is positive."""
        if v <= 0:
            raise ValueError("Fragment count must be positive")
        return v


class GenomicBin(BaseModel):
    """Represents a genomic bin with variant counts."""
    
    chromosome: str = Field(description="Chromosome name")
    start: int = Field(description="Start position (0-based)")
    end: int = Field(description="End position (1-based)")
    variant_count: int = Field(description="Number of variants in this bin")
    
    @validator('start', 'end')
    def validate_positions(cls, v):
        """Validate that positions are non-negative."""
        if v < 0:
            raise ValueError("Genomic positions must be non-negative")
        return v
    
    @validator('end')
    def validate_end_after_start(cls, v, values):
        """Validate that end position is after start position."""
        if 'start' in values and v <= values['start']:
            raise ValueError("End position must be after start position")
        return v
    
    @validator('variant_count')
    def validate_variant_count(cls, v):
        """Validate that variant count is non-negative."""
        if v < 0:
            raise ValueError("Variant count must be non-negative")
        return v


class GeneVariantStats(BaseModel):
    """Variant statistics for a specific gene."""
    
    gene_name: str = Field(description="Gene name")
    chromosome: str = Field(description="Chromosome name")
    start: int = Field(description="Gene start position")
    end: int = Field(description="Gene end position")
    total_variants: int = Field(description="Total number of variants")
    synonymous_variants: int = Field(description="Number of synonymous variants")
    nonsynonymous_variants: int = Field(description="Number of nonsynonymous variants")
    dn_ds_ratio: Optional[float] = Field(
        default=None, 
        description="dN/dS ratio (nonsynonymous/synonymous substitution rate)"
    )
    
    @validator('total_variants', 'synonymous_variants', 'nonsynonymous_variants')
    def validate_variant_counts(cls, v):
        """Validate that variant counts are non-negative."""
        if v < 0:
            raise ValueError("Variant counts must be non-negative")
        return v
    
    @validator('dn_ds_ratio')
    def validate_dn_ds_ratio(cls, v):
        """Validate that dN/dS ratio is positive if provided."""
        if v is not None and v < 0:
            raise ValueError("dN/dS ratio must be non-negative")
        return v


class CNVRegion(BaseModel):
    """Copy number variation region."""
    
    chromosome: str = Field(description="Chromosome name")
    start: int = Field(description="Start position")
    end: int = Field(description="End position")
    copy_number: float = Field(description="Copy number value")
    confidence: float = Field(description="Confidence score")
    type: str = Field(description="CNV type (gain/loss)")
    
    @validator('copy_number')
    def validate_copy_number(cls, v):
        """Validate that copy number is positive."""
        if v <= 0:
            raise ValueError("Copy number must be positive")
        return v
    
    @validator('confidence')
    def validate_confidence(cls, v):
        """Validate that confidence is between 0 and 1."""
        if not 0 <= v <= 1:
            raise ValueError("Confidence must be between 0 and 1")
        return v
    
    @validator('type')
    def validate_type(cls, v):
        """Validate that CNV type is valid."""
        if v not in ['gain', 'loss']:
            raise ValueError("CNV type must be 'gain' or 'loss'")
        return v


class QualityMetrics(BaseModel):
    """Quality control metrics for the pipeline."""
    
    total_reads: int = Field(description="Total number of reads")
    mapped_reads: int = Field(description="Number of mapped reads")
    mapping_rate: float = Field(description="Mapping rate percentage")
    mean_coverage: float = Field(description="Mean coverage across genome")
    coverage_std: float = Field(description="Standard deviation of coverage")
    gc_content: float = Field(description="GC content percentage")
    duplication_rate: float = Field(description="Duplication rate percentage")
    quality_scores: Dict[str, float] = Field(
        default_factory=dict,
        description="Quality scores for different metrics"
    )
    
    @validator('mapping_rate', 'gc_content', 'duplication_rate')
    def validate_percentages(cls, v):
        """Validate that percentages are between 0 and 100."""
        if not 0 <= v <= 100:
            raise ValueError("Percentage values must be between 0 and 100")
        return v
    
    @validator('total_reads', 'mapped_reads')
    def validate_read_counts(cls, v):
        """Validate that read counts are non-negative."""
        if v < 0:
            raise ValueError("Read counts must be non-negative")
        return v
    
    @validator('mean_coverage', 'coverage_std')
    def validate_coverage(cls, v):
        """Validate that coverage values are non-negative."""
        if v < 0:
            raise ValueError("Coverage values must be non-negative")
        return v


class FeatureSet(BaseModel):
    """Complete set of features extracted from SRA data."""
    
    # Basic information
    sra_id: str = Field(description="SRA accession ID")
    sample_name: Optional[str] = Field(default=None, description="Sample name")
    
    # Fragment length statistics (for paired-end data)
    fragment_stats: Optional[FragmentLengthStats] = Field(
        default=None, 
        description="Fragment length statistics"
    )
    
    # Genomic variant bins
    genomic_bins: List[GenomicBin] = Field(
        default_factory=list,
        description="Genomic bins with variant counts"
    )
    
    # Gene-level variant statistics
    gene_stats: List[GeneVariantStats] = Field(
        default_factory=list,
        description="Gene-level variant statistics"
    )
    
    # Copy number variations
    cnv_regions: List[CNVRegion] = Field(
        default_factory=list,
        description="Copy number variation regions"
    )
    
    # Quality metrics
    quality_metrics: QualityMetrics = Field(description="Quality control metrics")
    
    # Additional metadata
    metadata: Dict[str, Union[str, int, float, bool]] = Field(
        default_factory=dict,
        description="Additional metadata"
    )
    
    # Processing information
    processing_time: float = Field(description="Total processing time in seconds")
    pipeline_version: str = Field(description="Pipeline version used")
    
    @validator('processing_time')
    def validate_processing_time(cls, v):
        """Validate that processing time is positive."""
        if v < 0:
            raise ValueError("Processing time must be non-negative")
        return v
    
    def to_dict(self) -> Dict:
        """Convert feature set to dictionary format."""
        return self.dict()
    
    def to_json(self) -> str:
        """Convert feature set to JSON string."""
        return self.json(indent=2)
    
    def get_summary_stats(self) -> Dict[str, Union[int, float]]:
        """Get summary statistics for the feature set."""
        return {
            "total_variants": sum(bin.variant_count for bin in self.genomic_bins),
            "total_genes_with_variants": len(self.gene_stats),
            "total_cnv_regions": len(self.cnv_regions),
            "mean_variants_per_bin": np.mean([bin.variant_count for bin in self.genomic_bins]) if self.genomic_bins else 0,
            "mapping_rate": self.quality_metrics.mapping_rate,
            "mean_coverage": self.quality_metrics.mean_coverage,
        }
    
    def filter_by_quality(self, min_mapping_rate: float = 80.0, min_coverage: float = 10.0) -> bool:
        """Check if the sample meets quality thresholds."""
        return (
            self.quality_metrics.mapping_rate >= min_mapping_rate and
            self.quality_metrics.mean_coverage >= min_coverage
        ) 