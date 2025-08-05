"""
Configuration settings for the SRA to Features Pipeline.
"""

import os
from pathlib import Path
from typing import Optional, List
from pydantic import Field, field_validator
from pydantic_settings import BaseSettings


class PipelineConfig(BaseSettings):
    """Configuration for the SRA to Features Pipeline."""
    
    # Base paths
    base_dir: Path = Field(default=Path("/content"), description="Base directory for pipeline")
    output_dir: Path = Field(default=Path("/content/data/output"), description="Output directory")
    
    # Reference files
    reference_fasta: Path = Field(description="Reference genome FASTA file")
    reference_gff: Path = Field(description="Reference genome GFF file")
    bed_file: Path = Field(default=Path("regions.bed"), description="BED file with regions of interest")
    bed_genes: Path = Field(description="BED file with gene regions")
    genome_sizes: Path = Field(description="Genome sizes file")
    
    # External tools
    kraken_db: Path = Field(description="Kraken2 database directory")
    snpeff_dir: Path = Field(description="snpEff installation directory")
    
    # Parameters
    genome_name: str = Field(default="custom_human", description="Genome name for snpEff")
    bin_size_gvs: int = Field(default=100000, description="Bin size for genomic variants")
    bin_size_cnv: int = Field(default=100000, description="Bin size for CNV analysis")
    threads: int = Field(default=1, description="Number of threads to use")
    
    # URLs for downloads
    fasta_url: str = Field(
        default="https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
        description="URL for reference FASTA"
    )
    gff_url: str = Field(
        default="https://ftp.ensembl.org/pub/release-109/gff3/homo_sapiens/Homo_sapiens.GRCh38.109.gff3.gz",
        description="URL for reference GFF"
    )
    snpeff_url: str = Field(
        default="https://sourceforge.net/projects/snpeff/files/snpEff_v4_3t_core.zip",
        description="URL for snpEff download"
    )
    kraken2_db_url: str = Field(
        default="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20201202.tar.gz",
        description="URL for Kraken2 database"
    )
    
    # Logging
    log_level: str = Field(default="INFO", description="Logging level")
    log_file: Optional[Path] = Field(default=None, description="Log file path")
    
    # Performance
    max_memory_gb: int = Field(default=8, description="Maximum memory usage in GB")
    timeout_seconds: int = Field(default=3600, description="Timeout for operations in seconds")
    
    # Quality control
    min_quality_score: int = Field(default=20, description="Minimum quality score for variants")
    min_coverage: int = Field(default=10, description="Minimum coverage for variant calling")
    
    @field_validator('base_dir', 'output_dir', 'reference_fasta', 'reference_gff', 
                    'bed_file', 'bed_genes', 'genome_sizes', 'kraken_db', 'snpeff_dir')
    @classmethod
    def validate_paths(cls, v):
        """Validate that paths are absolute or can be resolved."""
        if isinstance(v, str):
            v = Path(v)
        return v
    
    @field_validator('threads')
    @classmethod
    def validate_threads(cls, v):
        """Validate thread count is positive."""
        if v <= 0:
            raise ValueError("Threads must be positive")
        return v
    
    @field_validator('bin_size_gvs', 'bin_size_cnv')
    @classmethod
    def validate_bin_sizes(cls, v):
        """Validate bin sizes are positive."""
        if v <= 0:
            raise ValueError("Bin sizes must be positive")
        return v
    
    @field_validator('max_memory_gb')
    @classmethod
    def validate_memory(cls, v):
        """Validate memory usage is reasonable."""
        if v <= 0 or v > 128:
            raise ValueError("Memory usage must be between 1 and 128 GB")
        return v
    
    model_config = {
        "env_prefix": "SRA_PIPELINE_",
        "case_sensitive": False,
        "env_file": ".env"
    }
    
    def get_data_dir(self) -> Path:
        """Get the data directory path."""
        return self.base_dir / "data"
    
    def get_install_dir(self) -> Path:
        """Get the installation directory path."""
        return self.base_dir / "install"
    
    def get_bin_dir(self) -> Path:
        """Get the binary directory path."""
        return self.base_dir / "bin"
    
    def get_log_dir(self) -> Path:
        """Get the log directory path."""
        return self.output_dir / "logs"
    
    def get_tmp_dir(self) -> Path:
        """Get the temporary directory path."""
        return self.get_data_dir() / "tmp"
    
    def ensure_directories(self) -> None:
        """Ensure all required directories exist."""
        directories = [
            self.get_data_dir(),
            self.get_install_dir(),
            self.get_bin_dir(),
            self.get_log_dir(),
            self.get_tmp_dir(),
            self.output_dir,
        ]
        
        for directory in directories:
            directory.mkdir(parents=True, exist_ok=True)
    
    def validate_setup(self) -> List[str]:
        """Validate that the pipeline is properly set up."""
        errors = []
        
        # Check required files exist
        required_files = [
            self.reference_fasta,
            self.reference_gff,
            self.bed_genes,
            self.genome_sizes,
        ]
        
        for file_path in required_files:
            if not file_path.exists():
                errors.append(f"Required file not found: {file_path}")
        
        # Check required directories exist
        required_dirs = [
            self.kraken_db,
            self.snpeff_dir,
        ]
        
        for dir_path in required_dirs:
            if not dir_path.exists():
                errors.append(f"Required directory not found: {dir_path}")
        
        # Check external tools are available
        required_tools = [
            "bwa", "samtools", "bcftools", "bedtools", 
            "fastq-dump", "tabix", "bgzip", "java"
        ]
        
        for tool in required_tools:
            if not self._check_tool_available(tool):
                errors.append(f"Required tool not found: {tool}")
        
        return errors
    
    def _check_tool_available(self, tool: str) -> bool:
        """Check if a tool is available in PATH."""
        import shutil
        return shutil.which(tool) is not None 