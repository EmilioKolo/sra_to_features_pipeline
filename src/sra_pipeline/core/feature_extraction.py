"""
Feature extraction functionality for genomic data.
"""

from pathlib import Path
from typing import Dict, Any, List
from ..utils import log_command
import gzip
import pandas as pd
import statistics
import structlog
import subprocess
import sys
import tempfile

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
    snpeff_file: Path,
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
    features["genomic_bins"] = _extract_genomic_bins(sample_id, vcf_file, genome_sizes, bin_size_gvs, logger)
    
    # Extract gene-level variant statistics
    features["gene_stats"] = _extract_gene_variants(sample_id, snpeff_file, bed_genes, reference_gff, logger)
    
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
    sample_id: str,
    vcf_file: Path,
    genome_sizes: Path,
    bin_size: int,
    logger: structlog.BoundLogger
) -> List[GenomicBin]:
    """Extract variant counts per genomic bin."""
    logger.info("Extracting genomic variant bins", 
                vcf_file=str(vcf_file), 
                bin_size=bin_size)
    # Initialize the output list
    bins = []
    # Create a temporary directory to store intermediate BED files
    with tempfile.TemporaryDirectory() as tmpdir:
        bins_bed = Path(tmpdir) / "genome_bins.bed"
        variants_bed = Path(tmpdir) / "variants.bed"
        intersected = Path(tmpdir) / "intersected.bed"

        # Generate genome bins using bedtools makewindows
        cmd_makewindows = ['bedtools', 'makewindows', '-g', 
                           str(genome_sizes), '-w', str(bin_size)]
        log_command(logger, " ".join(cmd_makewindows),
                    sample_id=sample_id)
        try:
            with open(bins_bed, "w") as outfile:
                subprocess.run(
                    cmd_makewindows,
                    stdout=outfile,
                    check=True,
                    text=True
                )
            logger.info(f"Genome bins saved to {bins_bed}", 
                        sample_id=sample_id)
        except FileNotFoundError:
            logger.error("bedtools command not found. Please ensure bedtools is installed and in your PATH.",
                         sample_id=sample_id)
            return []
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to generate genome bins using bedtools: {e}",
                         sample_id=sample_id)
            logger.error(f"Stderr: {e.stderr}",
                         sample_id=sample_id)
            return []
        except Exception as e:
            logger.error(f"An unexpected error occurred during bedtools makewindows: {e}",
                         sample_id=sample_id)
            return []

        # Convert VCF to BED format
        logger.info(f"Converting VCF file {vcf_file} to BED format...",
                    sample_id=sample_id)
        try:
            with gzip.open(vcf_file, 'rt') as vcf_in, \
                 open(variants_bed, "w") as bed_out:
                for line in vcf_in:
                    if line.startswith("#"):
                        continue
                    fields = line.strip().split("\t")
                    chrom = fields[0]
                    start = int(fields[1]) - 1 # BED is 0-based
                    ref = fields[3]
                    end = start + max(len(ref), 1) # Minimum 1 base length for variants
                    bed_out.write(f"{chrom}\t{start}\t{end}\n")
            logger.info(f"Variants converted to BED and saved to {variants_bed}",
                        sample_id=sample_id)
        except FileNotFoundError:
            logger.error(f"VCF file not found at {vcf_file}.",
                         sample_id=sample_id)
            return []
        except gzip.BadGzipFile:
            logger.error(f"VCF file {vcf_file} is not a valid gzipped file.",
                         sample_id=sample_id)
            return []
        except Exception as e:
            logger.error(f"Failed to convert VCF to BED: {e}",
                         sample_id=sample_id)
            return []

        # Count variants per bin using bedtools intersect
        cmd_intersect = ["bedtools", "intersect", "-a", str(bins_bed), 
                         "-b", str(variants_bed), "-c"]
        log_command(logger, " ".join(cmd_intersect),
                    sample_id=sample_id)
        try:
            with open(intersected, "w") as outfile:
                subprocess.run(
                    cmd_intersect,
                    stdout=outfile,
                    check=True,
                    text=True
                )
            logger.info(f"Intersection results saved to {intersected}",
                        sample_id=sample_id)
        except FileNotFoundError:
            logger.error("bedtools command not found. Please ensure bedtools is installed and in your PATH.",
                         sample_id=sample_id)
            return []
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to intersect variants with genome bins: {e}",
                         sample_id=sample_id)
            logger.error(f"Stderr: {e.stderr}",
                         sample_id=sample_id)
            return []
        except Exception as e:
            logger.error(f"An unexpected error occurred during bedtools intersect: {e}",
                         sample_id=sample_id)
            return []

        # Parse intersection output into a list of dictionaries
        logger.info("Parsing intersection output...", sample_id=sample_id)
        try:
            df = pd.read_csv(
                intersected,
                sep="\t",
                header=None,
                names=["chrom", "start", "end", "variant_count"]
            )
        except pd.errors.EmptyDataError:
            logger.warning("Intersection result file is empty. No variants found in bins.",
                           sample_id=sample_id)
            return []
        except Exception as e:
            logger.error(f"Failed to read intersection result: {e}",
                         sample_id=sample_id)
            return []

        if df.empty:
            logger.warning("No intersected variants found after parsing.",
                           sample_id=sample_id)
        else:
            # Load the output to the bins list
            for _, row in df.iterrows():
                bins.append(GenomicBin(
                    chromosome=row['chrom'],
                    start=row['start'],
                    end=row['end'],
                    variant_count=row["variant_count"]
                ))
            logger.info(f"Successfully processed {len(bins)} genomic bins.",
                        sample_id=sample_id)
    return bins


def _extract_gene_variants(
    sample_id: str,
    vcf_file: Path,
    bed_genes: Path,
    reference_gff: Path,
    logger: structlog.BoundLogger
) -> List[GeneVariantStats]:
    """Extract gene-level variant statistics."""
    logger.info("Extracting gene-level variant statistics",
                sample_id=sample_id,
                vcf_file=str(vcf_file), 
                bed_genes=str(bed_genes))
    
    # Initialize the variable that is returned
    genes = []
    # Get gene names and positions from the gff file
    gene_data = []
    with open(reference_gff, 'r') as gff_file:
        for line in gff_file:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) > 8 and "gene" in parts[2]:
                attr = parts[8]
                gene_name = attr.split("gene_name=")[-1].split(";")[0]
                gene_data.append({
                    "gene_name": gene_name,
                    "chr": parts[0],
                    "start": int(parts[3]),
                    "end": int(parts[4])
                })
    
    # Create a BED data string in memory
    bed_data: str = ""
    for gene in gene_data:
        bed_data += f"{gene['chr']}\t{gene['start']}\t" + \
                    f"{gene['end']}\t{gene['gene_name']}\n"
    
    # Run bcftools intersect and capture the output
    cmd = ["bcftools", "view", "-R", "-", str(vcf_file)]
    
    log_command(logger, " ".join(cmd), vcf_file=str(vcf_file))

    try:
        result = subprocess.run(
            cmd,
            input=bed_data,
            capture_output=True,
            text=True,
            timeout=600  # 10 minute timeout
        )

    except subprocess.TimeoutExpired:
        raise RuntimeError(f"bcftools view timed out for sample: {sample_id}")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"bcftools view failed for sample {sample_id}: {e.stderr}")

    # Parse the bcftools output to count variants per gene
    varcount_all = {gene['gene_name']: 0 for gene in gene_data}
    varcount_syn = {gene['gene_name']: 0 for gene in gene_data}
    varcount_nonsyn = {gene['gene_name']: 0 for gene in gene_data}
    varcount_indel = {gene['gene_name']: 0 for gene in gene_data}
    # Start a set of unhandled genetic effects
    unhandled_effects = set()
    n_unhandled = 0
    if result.returncode == 0:
        for line in result.stdout.splitlines():
            if line.startswith('#'):
                continue
            parts = line.split('\t')
            chr = parts[0]
            start = int(parts[1])
            # Determine which gene the variant belongs to
            for gene in gene_data:
                if gene['chr'] == chr and \
                    gene['start'] <= start <= gene['end']:
                    # Count the variant to varcount_all
                    varcount_all[gene['gene_name']] += 1
                    # Define if the variant is synonymous or nonsynonymous
                    info_field = parts[7].split(';')
                    if "INDEL"==info_field[0].upper():
                        varcount_indel[gene['gene_name']] += 1
                    else:
                        ann_field = _extract_ann(info_field)
                        l_ann = ann_field.split(',')
                        nonsynonymous = False
                        synonymous = False
                        for ann in l_ann:
                            if 'missense_variant' in ann.lower() or \
                               'nonsense_variant' in ann.lower() or \
                                'stop_gained' in ann.lower() or \
                                'stop_lost' in ann.lower():
                                varcount_nonsyn[gene['gene_name']] += 1
                                nonsynonymous = True
                            if 'synonymous_variant' in ann.lower():
                                varcount_syn[gene['gene_name']] += 1
                                synonymous = True
                            if (not synonymous and not nonsynonymous) and \
                                not 'upstream' in ann.lower() and \
                                not 'downstream' in ann.lower() and \
                                not 'intron' in ann.lower() and \
                                not 'intergenic' in ann.lower() and \
                                not 'non_coding' in ann.lower() and \
                                not 'intragenic' in ann.lower() and \
                                not 'utr_variant' in ann.lower():
                                curr_eff = ann.split('|')[1] if '|' in ann else ann
                                unhandled_effects.add(curr_eff)
                                n_unhandled += 1
    # Check if there were unhandled effects
    if unhandled_effects:
        logger.warning(
            f"There were {n_unhandled} unhandled effects of the following types: {unhandled_effects}",
            sample_id=sample_id)
    else:
        logger.warning(
            f"Bad returncode: {result.returncode}\n{result.stderr}",
            sample_id=sample_id
        )
        return []
    # Define values per gene
    for gene in gene_data:
        n_var = varcount_all.get(gene['gene_name'], 0)
        n_syn = varcount_syn.get(gene['gene_name'], 0)
        n_nonsyn = varcount_nonsyn.get(gene['gene_name'], 0)
        n_indel = varcount_indel.get(gene['gene_name'], 0)
        dn_ds_ratio = (n_nonsyn / n_syn) if n_syn > 0 else float('inf')
        genes.append(
            GeneVariantStats(
                gene_name=gene['gene_name'],
                chromosome=gene['chr'],
                start=gene['start'],
                end=gene['end'],
                total_variants=n_var,
                synonymous_variants=n_syn,
                nonsynonymous_variants=n_nonsyn,
                dn_ds_ratio=dn_ds_ratio,
                indel_variants=n_indel
            )
        )
    
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
    # Initialize the output list
    cnv_regions = []

    # Create a temporary directory to store intermediate cnvpytor files
    with tempfile.TemporaryDirectory() as tmpdir:
        logger.info("Defining CNVpytor root file",
                    sample_id=sample_id)
        # Define the root file for CNVpytor
        ROOT_FILE = Path(tmpdir) / f"{sample_id}.pytor"
        
        # Remove existing root file if it exists
        if ROOT_FILE.exists():
            ROOT_FILE.unlink()
        
        logger.info(f"CNVpytor root file will be: {ROOT_FILE}", 
                    sample_id=sample_id)
        # Check that bam file exists
        if not bam_file.exists():
            logger.error(f"BAM file not found at {bam_file}",
                         sample_id=sample_id)
            sys.exit(1)
        
        # Wrap CNVpytor scripts in try/except
        try:
            # Run CNVpytor scripts in order
            cnv_call_file = _run_cnvpytor(
                sample_id,
                bam_file,
                ROOT_FILE,
                bin_size,
                Path(tmpdir),
                logger
            )

            # Read CNVpytor output and extract CNV features
            cnv_regions = _read_cnvpytor_out(
                sample_id,
                cnv_call_file,
                logger
            )
        except Exception as e:
            logger.warning(f'CNVpytor failed for sample with error: {e}',
                           sample_id=sample_id)
            cnv_regions = []
    return cnv_regions

def _extract_ann(info_field:list[str]) -> str:
    """Extracts the ANN field from a snpEff vcf into a string."""
    # Go through info_field
    for field in info_field:
        if field.startswith('ANN='):
            return field.split('ANN=')[1]
    return ""

def _run_cnvpytor(
    sample_id: str,
    bam_file: Path,
    root_file: Path,
    bin_size: int,
    output_dir: Path,
    logger: structlog.BoundLogger
) -> Path:
    """Runs several cnvpytor scripts in order.
    Returns the path to the cnv call file."""

    # Define output file
    cnv_call_file = output_dir / f"cnv_calls_{sample_id}.txt"

    logger.info("Processing Read Depth data...",
                sample_id=sample_id)

    # Start by processing read depth data
    cmd1 = ['cnvpytor', '-root', str(root_file), '-rd', str(bam_file)]

    log_command(logger, " ".join(cmd1), sample_id=sample_id)

    try:
        result = subprocess.run(cmd1, capture_output=True, check=True,
                                timeout=600) # 10 minute timeout
    except subprocess.TimeoutExpired:
        raise RuntimeError(f"cnvpytor read depth timed out for sample: {sample_id}")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"cnvpytor read depth failed for sample {sample_id}: {e.stderr}")

    logger.info("Read Depth processing complete.",
                sample_id=sample_id)

    logger.info("Generating histograms and partitioning data...",
                sample_id=sample_id)

    # Create histograms for RD
    cmd2 = ['cnvpytor', '-root', str(root_file), '-his', str(bin_size)]

    log_command(logger, " ".join(cmd2), sample_id=sample_id)

    try:
        result = subprocess.run(cmd2, capture_output=True, check=True,
                                timeout=600) # 10 minute timeout
    except subprocess.TimeoutExpired:
        raise RuntimeError(f"cnvpytor histograms timed out for sample: {sample_id}")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"cnvpytor histograms failed for sample {sample_id}: {e.stderr}")

    logger.info("Histograms created.",
                sample_id=sample_id)
    
    logger.info("Starting with partitioning...",
                sample_id=sample_id)
    # Partition data for CNV calling
    cmd3 = ['cnvpytor', '-root', str(root_file),
            '-partition', str(bin_size)]
    
    log_command(logger, " ".join(cmd3), sample_id=sample_id)

    try:
        result = subprocess.run(cmd3, capture_output=True, check=True,
                                timeout=600) # 10 minute timeout
    except subprocess.TimeoutExpired:
        raise RuntimeError(f"cnvpytor partitioning timed out for sample: {sample_id}")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"cnvpytor partitioning failed for sample {sample_id}: {e}")

    logger.info("Partitioning complete.",
                sample_id=sample_id)

    # Call CNVs
    cmd4 = ['cnvpytor', '-root', str(root_file),
            '-call', str(bin_size), '>', str(cnv_call_file)]
    
    log_command(logger, " ".join(cmd4), sample_id=sample_id)

    try:
        result = subprocess.run(
            " ".join(cmd4),
            shell=True,
            capture_output=True,
            text=True,
            check=True,
            timeout=1800 # 30 minute timeout
        )
    except subprocess.TimeoutExpired:
        raise RuntimeError(f"cnvpytor CNV calling timed out for sample: {sample_id}")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"cnvpytor partitioning failed for sample {sample_id}: {e}")
    
    logger.info("CNV calling complete.",
                sample_id=sample_id)
    
    return cnv_call_file

def _read_cnvpytor_out(
    sample_id: str,
    cnv_call_file: Path,
    logger: structlog.BoundLogger
) -> list[CNVRegion]:
    """Extracts CNV features from a cnvpytor call output file."""
    # Define output list
    cnv_regions = []
    # Extract calls from cnv_call_file
    calls = []
    with open(cnv_call_file, 'r') as f:
        for line in f.readlines():
            calls.append(line.rstrip('\n').split('\t'))
    # Check CNV calls
    for call in calls:
        # Define type from call[0]
        if call[0] == 'deletion':
            cnv_type = 'loss'
        elif call[0] == 'duplication':
            cnv_type = 'gain'
        else:
            logger.error(f"Unknown CNV type: {call[0]} in sample {sample_id}. Skipping this call.",
                         sample_id=sample_id)
            continue
        # Define chromosome, start and end from call[1]
        chrom: str = call[1].split(':')[0]
        start_end: list = call[1].split(':')[1].split('-')
        start: int = int(start_end[0])
        end: int = int(start_end[1])
        # Define copy number and confidence
        copy_number: float = float(call[3]) # read depth normalized to average depth
        if float(call[6]) <= 1.0:
            confidence: float = 1.0 - float(call[6]) # 1 - p-value (e-val3)
        else:
            logger.warning(f"P-value {call[6]} is greater than 1.0 for {call}. Setting confidence to 0.0.",
                           sample_id=sample_id)
            confidence: float = 0.0
        cnv_regions.append(CNVRegion(
            chromosome=chrom,
            start=start,
            end=end,
            copy_number=copy_number,
            confidence=confidence,
            type=cnv_type
        ))
    return cnv_regions