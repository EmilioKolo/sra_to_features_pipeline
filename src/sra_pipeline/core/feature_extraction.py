"""
Feature extraction functionality for genomic data.
"""

from pathlib import Path
from typing import Dict, Any, List
from ..utils import log_command
import statistics
import structlog
import subprocess

from ..models.features import (
    FragmentLengthStats, 
    GenomicBin, 
    GeneVariantStats, 
    CNVRegion
)


def extract_features(
    sample_id: str,
    bam_file: Path,
    vcf_file: Path,
    reference_gff: Path,
    bed_genes: Path,
    genome_sizes: Path,
    bin_size_gvs: int,
    bin_size_cnv: int,
    output_dir: Path,
    logger: structlog.BoundLogger
) -> Dict[str, Any]:
    """
    Extract genomic features from aligned and variant-called data.
    
    Args:
        sample_id: Sample identifier
        bam_file: Aligned BAM file
        vcf_file: Variant VCF file
        reference_gff: Reference genome GFF file
        bed_genes: BED file with gene regions
        genome_sizes: Genome sizes file
        bin_size_gvs: Bin size for genomic variants
        bin_size_cnv: Bin size for CNV analysis
        output_dir: Output directory for intermediate files
        logger: Logger instance
        
    Returns:
        Dictionary containing extracted features
    """
    logger.info("Starting feature extraction", 
                sample_id=sample_id, 
                bam_file=str(bam_file),
                vcf_file=str(vcf_file))
    
    features = {}
    
    # Extract fragment length statistics (for paired-end data)
    features["fragment_stats"] = _extract_fragment_lengths(sample_id, bam_file, logger)
    
    # Extract genomic variant bins
    features["genomic_bins"] = _extract_genomic_bins(vcf_file, genome_sizes, bin_size_gvs, logger)
    
    # Extract gene-level variant statistics
    features["gene_stats"] = _extract_gene_variants(vcf_file, bed_genes, logger)
    
    # Extract copy number variations
    features["cnv_regions"] = _extract_cnv_regions(sample_id, bam_file, bin_size_cnv, logger)
    
    # Add metadata
    features["metadata"] = {
        "sample_id": sample_id,
        "bin_size_gvs": bin_size_gvs,
        "bin_size_cnv": bin_size_cnv,
        "pipeline_version": "1.0.0"
    }
    
    logger.info("Feature extraction completed", 
                sample_id=sample_id,
                features_extracted=list(features.keys()))
    
    return features


def _extract_fragment_lengths(
    sample_id: str,
    bam_file: Path,
    logger: structlog.BoundLogger,
    min_size: int = 0,
    max_size: int = 1000
) -> FragmentLengthStats:
    """Extract fragment length statistics from BAM file."""
    logger.info("Extracting fragment length statistics", 
                sample_id=sample_id)

    # First command: samtools view
    cmd1 = ["samtools", "view", "-f", "0x2", str(bam_file)]

    log_command(logger, " ".join(cmd1), sample_id=sample_id)

    p1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)

    # Second command: awk
    awk_command = f'{{if ($9>{min_size} && $9<{max_size}) print $9}}'
    cmd2 = ['awk', awk_command]

    log_command(logger, " ".join(cmd2), sample_id=sample_id)

    p2 = subprocess.Popen(cmd2, stdin=p1.stdout, stdout=subprocess.PIPE)
    p1.stdout.close()
    
    # Get the output of the second command
    output, _ = p2.communicate()

    if not _ is None:
        _decoded = _.decode('utf-8')
        logger.warning(f'_: {_decoded}', sample_id=sample_id)

    # Parse the output to get fragment lengths
    fragment_lengths = output.decode('utf-8').strip().split('\n')
    # Convert to integers
    fragment_lengths = list(map(int, fragment_lengths))

    # Check if there are any fragment lengths
    if not fragment_lengths:
        logger.warning("No fragment lengths found", sample_id=sample_id)
        return FragmentLengthStats(
            mean=0.0,
            median=0.0,
            std=0.0,
            min=0.0,
            max=0.0,
            count=0
        )
    
    logger.info("Fragment lengths extracted",
                sample_id=sample_id,
                bam_file=str(bam_file))
    # Calculate and return required stats
    return FragmentLengthStats(
        mean=statistics.mean(fragment_lengths),
        median=statistics.median(fragment_lengths),
        std=statistics.stdev(fragment_lengths),
        min=min(fragment_lengths),
        max=max(fragment_lengths),
        count=len(fragment_lengths)
    )


def _extract_genomic_bins(
    vcf_file: Path,
    genome_sizes: Path,
    bin_size: int,
    logger: structlog.BoundLogger
) -> List[GenomicBin]:
    """Extract variant counts per genomic bin."""
    logger.info("Extracting genomic variant bins", 
                vcf_file=str(vcf_file), 
                bin_size=bin_size)
    
    # This is a placeholder implementation
    # In a real implementation, you would:
    # 1. Read the genome sizes file
    # 2. Create bins based on bin_size
    # 3. Count variants in each bin using bedtools
    
    # For now, return dummy data
    bins = []
    chromosomes = ["chr1", "chr2", "chr3", "chrX", "chrY"]
    
    for chrom in chromosomes:
        for i in range(0, 1000000, bin_size):
            bins.append(GenomicBin(
                chromosome=chrom,
                start=i,
                end=i + bin_size,
                variant_count=10  # Dummy count
            ))
    
    return bins


def _extract_gene_variants(
    vcf_file: Path,
    bed_genes: Path,
    logger: structlog.BoundLogger
) -> List[GeneVariantStats]:
    """Extract gene-level variant statistics."""
    logger.info("Extracting gene-level variant statistics", 
                vcf_file=str(vcf_file), 
                bed_genes=str(bed_genes))
    
    # This is a placeholder implementation
    # In a real implementation, you would:
    # 1. Read the BED genes file
    # 2. Use bedtools to intersect variants with genes
    # 3. Calculate synonymous/nonsynonymous ratios
    
    # For now, return dummy data
    genes = []
    gene_names = ["GENE1", "GENE2", "GENE3", "GENE4", "GENE5"]
    
    for i, gene_name in enumerate(gene_names):
        genes.append(GeneVariantStats(
            gene_name=gene_name,
            chromosome=f"chr{i+1}",
            start=i * 10000,
            end=(i + 1) * 10000,
            total_variants=50,
            synonymous_variants=30,
            nonsynonymous_variants=20,
            dn_ds_ratio=0.67
        ))
    
    return genes


def _extract_cnv_regions(
    sample_id: str,
    bam_file: Path,
    bin_size: int,
    logger: structlog.BoundLogger
) -> List[CNVRegion]:
    """Extract copy number variation regions."""
    logger.info("Extracting CNV regions", 
                sample_id=sample_id, 
                bam_file=str(bam_file),
                bin_size=bin_size)
    
    # This is a placeholder implementation
    # In a real implementation, you would:
    # 1. Use CNVpytor or similar tool
    # 2. Analyze read depth patterns
    # 3. Identify copy number variations
    
    # For now, return dummy data
    cnv_regions = []
    chromosomes = ["chr1", "chr2", "chr3"]
    
    for i, chrom in enumerate(chromosomes):
        cnv_regions.append(CNVRegion(
            chromosome=chrom,
            start=i * 50000,
            end=(i + 1) * 50000,
            copy_number=3.0,  # Amplification
            confidence=0.95,
            type="gain"
        ))
    
    return cnv_regions 