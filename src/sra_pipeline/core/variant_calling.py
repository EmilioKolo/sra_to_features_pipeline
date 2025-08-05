"""
Variant calling functionality using bcftools.
"""

import subprocess
from pathlib import Path
import structlog

from ..utils import log_command


def run_variant_calling(
    sample_id: str,
    bam_file: Path,
    reference_fasta: Path,
    output_dir: Path,
    min_quality_score: int,
    min_coverage: int,
    logger: structlog.BoundLogger
) -> Path:
    """
    Run variant calling using bcftools.
    
    Args:
        sample_id: Sample identifier
        bam_file: Input BAM file
        reference_fasta: Reference genome FASTA file
        output_dir: Output directory for VCF files
        min_quality_score: Minimum quality score for variants
        min_coverage: Minimum coverage for variant calling
        logger: Logger instance
        
    Returns:
        Path to the final VCF file
        
    Raises:
        RuntimeError: If variant calling fails
    """
    logger.info("Starting variant calling", 
                sample_id=sample_id, 
                bam_file=str(bam_file),
                reference=str(reference_fasta),
                min_quality=min_quality_score,
                min_coverage=min_coverage)
    
    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Generate pileup
    bcf_file = _run_mpileup(sample_id, bam_file, reference_fasta, output_dir, logger)
    
    # Step 2: Call variants
    vcf_file = _call_variants(sample_id, bcf_file, output_dir, min_quality_score, logger)
    
    # Step 3: Filter variants
    filtered_vcf_file = _filter_variants(sample_id, vcf_file, output_dir, min_quality_score, min_coverage, logger)
    
    # Step 4: Compress and index VCF
    final_vcf_file = _compress_and_index_vcf(sample_id, filtered_vcf_file, logger)
    
    logger.info("Variant calling completed", 
                sample_id=sample_id, 
                output_vcf=str(final_vcf_file))
    
    return final_vcf_file


def _run_mpileup(
    sample_id: str,
    bam_file: Path,
    reference_fasta: Path,
    output_dir: Path,
    logger: structlog.BoundLogger
) -> Path:
    """Run bcftools mpileup to generate pileup."""
    bcf_file = output_dir / f"{sample_id}.bcf"
    
    cmd = [
        "bcftools", "mpileup",
        "-f", str(reference_fasta),
        str(bam_file),
        "-o", str(bcf_file)
    ]
    
    log_command(logger, " ".join(cmd), sample_id=sample_id)
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
            timeout=3600  # 1 hour timeout
        )
        
        logger.info("Mpileup completed", 
                    sample_id=sample_id, 
                    bcf_file=str(bcf_file))
        
        return bcf_file
        
    except subprocess.TimeoutExpired:
        raise RuntimeError(f"Mpileup timed out for sample: {sample_id}")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Mpileup failed for sample {sample_id}: {e.stderr}")


def _call_variants(
    sample_id: str,
    bcf_file: Path,
    output_dir: Path,
    min_quality_score: int,
    logger: structlog.BoundLogger
) -> Path:
    """Call variants using bcftools call."""
    vcf_file = output_dir / f"{sample_id}.vcf"
    
    cmd = [
        "bcftools", "call",
        "-mv",  # Multi-allelic variant caller
        "-o", str(vcf_file),
        str(bcf_file)
    ]
    
    log_command(logger, " ".join(cmd), sample_id=sample_id)
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
            timeout=1800  # 30 minute timeout
        )
        
        logger.info("Variant calling completed", 
                    sample_id=sample_id, 
                    vcf_file=str(vcf_file))
        
        return vcf_file
        
    except subprocess.TimeoutExpired:
        raise RuntimeError(f"Variant calling timed out for sample: {sample_id}")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Variant calling failed for sample {sample_id}: {e.stderr}")


def _filter_variants(
    sample_id: str,
    vcf_file: Path,
    output_dir: Path,
    min_quality_score: int,
    min_coverage: int,
    logger: structlog.BoundLogger
) -> Path:
    """Filter variants based on quality criteria."""
    filtered_vcf_file = output_dir / f"{sample_id}_filtered.vcf"
    
    # Create filter expression
    filter_expr = f"QUAL >= {min_quality_score} && INFO/DP >= {min_coverage}"
    
    cmd = [
        "bcftools", "filter",
        "-i", filter_expr,
        "-o", str(filtered_vcf_file),
        str(vcf_file)
    ]
    
    log_command(logger, " ".join(cmd), sample_id=sample_id, filter_expr=filter_expr)
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
            timeout=600  # 10 minute timeout
        )
        
        logger.info("Variant filtering completed", 
                    sample_id=sample_id, 
                    filtered_vcf=str(filtered_vcf_file))
        
        return filtered_vcf_file
        
    except subprocess.TimeoutExpired:
        raise RuntimeError(f"Variant filtering timed out for sample: {sample_id}")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Variant filtering failed for sample {sample_id}: {e.stderr}")


def _compress_and_index_vcf(
    sample_id: str,
    vcf_file: Path,
    logger: structlog.BoundLogger
) -> Path:
    """Compress and index VCF file."""
    compressed_vcf_file = vcf_file.with_suffix(".vcf.gz")
    
    # Compress VCF
    cmd = ["bgzip", "-c", str(vcf_file), ">", str(compressed_vcf_file)]
    log_command(logger, " ".join(cmd), sample_id=sample_id)
    
    try:
        result = subprocess.run(
            " ".join(cmd),
            shell=True,
            capture_output=True,
            text=True,
            check=True,
            timeout=300  # 5 minute timeout
        )
        
        # Index VCF
        cmd = ["tabix", "-p", "vcf", str(compressed_vcf_file)]
        log_command(logger, " ".join(cmd), sample_id=sample_id)
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
            timeout=300  # 5 minute timeout
        )
        
        logger.info("VCF compression and indexing completed", 
                    sample_id=sample_id, 
                    compressed_vcf=str(compressed_vcf_file))
        
        return compressed_vcf_file
        
    except subprocess.TimeoutExpired:
        raise RuntimeError(f"VCF compression/indexing timed out for sample: {sample_id}")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"VCF compression/indexing failed for sample {sample_id}: {e.stderr}") 