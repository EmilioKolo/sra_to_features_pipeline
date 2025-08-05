"""
Main pipeline class for the SRA to Features Pipeline.
"""

import time
from pathlib import Path
from typing import List, Optional, Dict, Any
import structlog

from ..config.settings import PipelineConfig
from ..models.features import FeatureSet, QualityMetrics
from ..utils import PipelineLogger, PerformanceMonitor, log_command, log_error
from . import download, alignment, variant_calling, feature_extraction, quality_control


class Pipeline:
    """Main pipeline class for processing SRA data and extracting features."""
    
    def __init__(self, config: PipelineConfig, logger: structlog.BoundLogger):
        """
        Initialize the pipeline.
        
        Args:
            config: Pipeline configuration
            logger: Structured logger instance
        """
        self.config = config
        self.logger = logger
        self.monitor = PerformanceMonitor(logger)
        
        # Ensure directories exist
        self.config.ensure_directories()
        
        # Validate setup
        errors = self.config.validate_setup()
        if errors:
            error_msg = "Pipeline setup validation failed:\n" + "\n".join(f"  - {e}" for e in errors)
            raise RuntimeError(error_msg)
        
        self.logger.info("Pipeline initialized successfully", config_summary=self._get_config_summary())
    
    def run_sra(self, sra_id: str) -> FeatureSet:
        """
        Run the pipeline on an SRA ID.
        
        Args:
            sra_id: SRA accession ID
            
        Returns:
            FeatureSet containing extracted features
        """
        with PipelineLogger(self.logger, f"pipeline_sra_{sra_id}") as plog:
            plog.add_context(sra_id=sra_id)
            
            # Download FASTQ files
            fastq_files = self._download_sra_data(sra_id)
            
            # Run pipeline on FASTQ files
            return self._run_pipeline(sra_id, fastq_files)
    
    def run_fastq(self, fastq_files: List[Path]) -> FeatureSet:
        """
        Run the pipeline on FASTQ files.
        
        Args:
            fastq_files: List of FASTQ file paths
            
        Returns:
            FeatureSet containing extracted features
        """
        # Generate sample ID from FASTQ files
        sample_id = self._generate_sample_id(fastq_files)
        
        with PipelineLogger(self.logger, f"pipeline_fastq_{sample_id}") as plog:
            plog.add_context(sample_id=sample_id, fastq_files=[str(f) for f in fastq_files])
            
            return self._run_pipeline(sample_id, fastq_files)
    
    def _download_sra_data(self, sra_id: str) -> List[Path]:
        """Download FASTQ files for an SRA ID."""
        with PipelineLogger(self.logger, f"download_sra_{sra_id}") as plog:
            plog.add_context(sra_id=sra_id)
            
            try:
                fastq_files = download.download_sra_data(
                    sra_id=sra_id,
                    output_dir=self.config.get_data_dir() / "fastq",
                    logger=self.logger
                )
                
                plog.log_progress(f"Downloaded {len(fastq_files)} FASTQ files")
                return fastq_files
                
            except Exception as e:
                log_error(self.logger, e, context={"sra_id": sra_id, "operation": "download"})
                raise
    
    def _run_pipeline(self, sample_id: str, fastq_files: List[Path]) -> FeatureSet:
        """Run the complete pipeline workflow."""
        start_time = time.time()
        
        try:
            # Step 1: Quality Control
            qc_results = self._run_quality_control(fastq_files)
            
            # Step 2: Alignment
            bam_file = self._run_alignment(sample_id, fastq_files)
            
            # Step 3: Variant Calling
            vcf_file = self._run_variant_calling(sample_id, bam_file)
            
            # Step 4: Feature Extraction
            features = self._extract_features(sample_id, bam_file, vcf_file)
            
            # Step 5: Create FeatureSet
            feature_set = self._create_feature_set(
                sample_id=sample_id,
                features=features,
                qc_results=qc_results,
                processing_time=time.time() - start_time
            )
            
            # Step 6: Save results
            self._save_results(sample_id, feature_set)
            
            return feature_set
            
        except Exception as e:
            log_error(self.logger, e, context={"sample_id": sample_id, "operation": "pipeline"})
            raise
    
    def _run_quality_control(self, fastq_files: List[Path]) -> Dict[str, Any]:
        """Run quality control on FASTQ files."""
        with PipelineLogger(self.logger, "quality_control") as plog:
            plog.add_context(fastq_files=[str(f) for f in fastq_files])
            
            try:
                qc_results = quality_control.run_quality_control(
                    fastq_files=fastq_files,
                    output_dir=self.config.get_tmp_dir(),
                    logger=self.logger
                )
                
                plog.log_progress("Quality control completed")
                return qc_results
                
            except Exception as e:
                log_error(self.logger, e, context={"operation": "quality_control"})
                raise
    
    def _run_alignment(self, sample_id: str, fastq_files: List[Path]) -> Path:
        """Run alignment of FASTQ files to reference genome."""
        with PipelineLogger(self.logger, f"alignment_{sample_id}") as plog:
            plog.add_context(sample_id=sample_id, fastq_files=[str(f) for f in fastq_files])
            
            try:
                bam_file = alignment.run_alignment(
                    sample_id=sample_id,
                    fastq_files=fastq_files,
                    reference_fasta=self.config.reference_fasta,
                    output_dir=self.config.get_data_dir() / "bam",
                    threads=self.config.threads,
                    logger=self.logger
                )
                
                plog.log_progress(f"Alignment completed: {bam_file}")
                return bam_file
                
            except Exception as e:
                log_error(self.logger, e, context={"sample_id": sample_id, "operation": "alignment"})
                raise
    
    def _run_variant_calling(self, sample_id: str, bam_file: Path) -> Path:
        """Run variant calling on aligned BAM file."""
        with PipelineLogger(self.logger, f"variant_calling_{sample_id}") as plog:
            plog.add_context(sample_id=sample_id, bam_file=str(bam_file))
            
            try:
                vcf_file = variant_calling.run_variant_calling(
                    sample_id=sample_id,
                    bam_file=bam_file,
                    reference_fasta=self.config.reference_fasta,
                    output_dir=self.config.get_data_dir() / "vcf",
                    min_quality_score=self.config.min_quality_score,
                    min_coverage=self.config.min_coverage,
                    logger=self.logger
                )
                
                plog.log_progress(f"Variant calling completed: {vcf_file}")
                return vcf_file
                
            except Exception as e:
                log_error(self.logger, e, context={"sample_id": sample_id, "operation": "variant_calling"})
                raise
    
    def _extract_features(self, sample_id: str, bam_file: Path, vcf_file: Path) -> Dict[str, Any]:
        """Extract features from aligned and variant-called data."""
        with PipelineLogger(self.logger, f"feature_extraction_{sample_id}") as plog:
            plog.add_context(sample_id=sample_id, bam_file=str(bam_file), vcf_file=str(vcf_file))
            
            try:
                features = feature_extraction.extract_features(
                    sample_id=sample_id,
                    bam_file=bam_file,
                    vcf_file=vcf_file,
                    reference_gff=self.config.reference_gff,
                    bed_genes=self.config.bed_genes,
                    genome_sizes=self.config.genome_sizes,
                    bin_size_gvs=self.config.bin_size_gvs,
                    bin_size_cnv=self.config.bin_size_cnv,
                    output_dir=self.config.get_tmp_dir(),
                    logger=self.logger
                )
                
                plog.log_progress("Feature extraction completed")
                return features
                
            except Exception as e:
                log_error(self.logger, e, context={"sample_id": sample_id, "operation": "feature_extraction"})
                raise
    
    def _create_feature_set(
        self, 
        sample_id: str, 
        features: Dict[str, Any], 
        qc_results: Dict[str, Any],
        processing_time: float
    ) -> FeatureSet:
        """Create a FeatureSet from extracted features and QC results."""
        try:
            # Create quality metrics
            quality_metrics = QualityMetrics(
                total_reads=qc_results.get("total_reads", 0),
                mapped_reads=qc_results.get("mapped_reads", 0),
                mapping_rate=qc_results.get("mapping_rate", 0.0),
                mean_coverage=qc_results.get("mean_coverage", 0.0),
                coverage_std=qc_results.get("coverage_std", 0.0),
                gc_content=qc_results.get("gc_content", 0.0),
                duplication_rate=qc_results.get("duplication_rate", 0.0),
                quality_scores=qc_results.get("quality_scores", {})
            )
            
            # Create feature set
            feature_set = FeatureSet(
                sra_id=sample_id,
                fragment_stats=features.get("fragment_stats"),
                genomic_bins=features.get("genomic_bins", []),
                gene_stats=features.get("gene_stats", []),
                cnv_regions=features.get("cnv_regions", []),
                quality_metrics=quality_metrics,
                metadata=features.get("metadata", {}),
                processing_time=processing_time,
                pipeline_version="1.0.0"
            )
            
            return feature_set
            
        except Exception as e:
            log_error(self.logger, e, context={"sample_id": sample_id, "operation": "create_feature_set"})
            raise
    
    def _save_results(self, sample_id: str, feature_set: FeatureSet):
        """Save pipeline results."""
        with PipelineLogger(self.logger, f"save_results_{sample_id}") as plog:
            plog.add_context(sample_id=sample_id)
            
            try:
                # Create output directory for this sample
                sample_output_dir = self.config.output_dir / sample_id
                sample_output_dir.mkdir(parents=True, exist_ok=True)
                
                # Save feature set as JSON
                feature_file = sample_output_dir / "features.json"
                with open(feature_file, "w") as f:
                    f.write(feature_set.to_json())
                
                # Save summary statistics
                summary_file = sample_output_dir / "summary.txt"
                summary = feature_set.get_summary_stats()
                with open(summary_file, "w") as f:
                    f.write("SRA to Features Pipeline Summary\n")
                    f.write("=" * 40 + "\n\n")
                    for key, value in summary.items():
                        f.write(f"{key}: {value}\n")
                
                plog.log_progress(f"Results saved to {sample_output_dir}")
                
            except Exception as e:
                log_error(self.logger, e, context={"sample_id": sample_id, "operation": "save_results"})
                raise
    
    def _generate_sample_id(self, fastq_files: List[Path]) -> str:
        """Generate a sample ID from FASTQ file names."""
        if len(fastq_files) == 1:
            return fastq_files[0].stem.replace("_1", "").replace("_2", "")
        else:
            # For paired-end, use common prefix
            names = [f.stem for f in fastq_files]
            common_prefix = ""
            for chars in zip(*names):
                if len(set(chars)) == 1:
                    common_prefix += chars[0]
                else:
                    break
            return common_prefix.rstrip("_")
    
    def _get_config_summary(self) -> Dict[str, Any]:
        """Get a summary of the configuration for logging."""
        return {
            "base_dir": str(self.config.base_dir),
            "output_dir": str(self.config.output_dir),
            "threads": self.config.threads,
            "bin_size_gvs": self.config.bin_size_gvs,
            "bin_size_cnv": self.config.bin_size_cnv,
        } 