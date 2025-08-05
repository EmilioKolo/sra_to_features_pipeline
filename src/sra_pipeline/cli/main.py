"""
Main command-line interface for the SRA to Features Pipeline.
"""

import sys
from pathlib import Path
from typing import Optional, List
import click
from rich.console import Console
from rich.table import Table
from rich.progress import Progress, SpinnerColumn, TextColumn

from ..config.settings import PipelineConfig
from ..core.pipeline import Pipeline
from ..utils import setup_logging, PipelineLogger, PerformanceMonitor
from .. import __version__


console = Console()


@click.group()
@click.version_option(version=__version__, prog_name="SRA Pipeline")
def cli():
    """SRA to Features Pipeline - Extract features from SRA data for LLM training."""
    pass


@cli.command()
@click.option(
    "--sra-id",
    help="SRA Run Accession ID (e.g., SRR123456)",
    type=str,
)
@click.option(
    "--fastq",
    multiple=True,
    help="FASTQ file paths (can specify multiple for paired-end)",
    type=click.Path(exists=True, path_type=Path),
)
@click.option(
    "--output-dir",
    default="./output",
    help="Output directory for results",
    type=click.Path(path_type=Path),
)
@click.option(
    "--config",
    help="Configuration file path",
    type=click.Path(exists=True, path_type=Path),
)
@click.option(
    "--threads",
    default=1,
    help="Number of threads to use",
    type=int,
)
@click.option(
    "--log-level",
    default="INFO",
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR"]),
    help="Logging level",
)
@click.option(
    "--log-file",
    help="Log file path",
    type=click.Path(path_type=Path),
)
@click.option(
    "--dry-run",
    is_flag=True,
    help="Perform a dry run without executing the pipeline",
)
def run(
    sra_id: Optional[str],
    fastq: List[Path],
    output_dir: Path,
    config: Optional[Path],
    threads: int,
    log_level: str,
    log_file: Optional[Path],
    dry_run: bool,
):
    """Run the SRA to Features Pipeline."""
    
    # Validate input
    if not sra_id and not fastq:
        console.print("[red]Error: Must provide either --sra-id or --fastq[/red]")
        sys.exit(1)
    
    if sra_id and fastq:
        console.print("[red]Error: Cannot specify both --sra-id and --fastq[/red]")
        sys.exit(1)
    
    if len(fastq) > 2:
        console.print("[red]Error: Maximum 2 FASTQ files allowed for paired-end sequencing[/red]")
        sys.exit(1)
    
    # Load configuration
    try:
        if config:
            pipeline_config = PipelineConfig(_env_file=config)
        else:
            pipeline_config = PipelineConfig()
        
        # Override config with CLI options
        pipeline_config.output_dir = output_dir
        pipeline_config.threads = threads
        pipeline_config.log_level = log_level
        pipeline_config.log_file = log_file
        
    except Exception as e:
        console.print(f"[red]Error loading configuration: {e}[/red]")
        sys.exit(1)
    
    # Setup logging
    logger = setup_logging(
        log_level=pipeline_config.log_level,
        log_file=pipeline_config.log_file,
        log_format="console"
    )
    
    # Display pipeline information
    with console.status("[bold green]Initializing pipeline..."):
        console.print(f"[bold blue]SRA to Features Pipeline v{__version__}[/bold blue]")
        console.print(f"Output directory: {output_dir}")
        console.print(f"Threads: {threads}")
        console.print(f"Log level: {log_level}")
        
        if sra_id:
            console.print(f"SRA ID: {sra_id}")
        else:
            console.print(f"FASTQ files: {[str(f) for f in fastq]}")
    
    if dry_run:
        console.print("[yellow]Dry run mode - no actual processing will occur[/yellow]")
        return
    
    # Validate required configuration
    missing_fields = []
    if pipeline_config.reference_fasta is None:
        missing_fields.append("reference_fasta")
    if pipeline_config.reference_gff is None:
        missing_fields.append("reference_gff")
    if pipeline_config.bed_genes is None:
        missing_fields.append("bed_genes")
    if pipeline_config.genome_sizes is None:
        missing_fields.append("genome_sizes")
    
    if missing_fields:
        console.print(f"[red]Error: Missing required configuration fields: {', '.join(missing_fields)}[/red]")
        console.print("[yellow]Please provide these fields via:[/yellow]")
        console.print("  - Environment variables (SRA_PIPELINE_REFERENCE_FASTA, etc.)")
        console.print("  - Configuration file (--config option)")
        console.print("  - .env file in the current directory")
        console.print("\n[yellow]Example configuration:[/yellow]")
        console.print("  export SRA_PIPELINE_REFERENCE_FASTA=/path/to/reference.fasta")
        console.print("  export SRA_PIPELINE_REFERENCE_GFF=/path/to/reference.gff")
        console.print("  export SRA_PIPELINE_BED_GENES=/path/to/genes.bed")
        console.print("  export SRA_PIPELINE_GENOME_SIZES=/path/to/genome.sizes")
        sys.exit(1)
    
    # Initialize pipeline
    try:
        pipeline = Pipeline(pipeline_config, logger)
    except Exception as e:
        console.print(f"[red]Error initializing pipeline: {e}[/red]")
        logger.error("Pipeline initialization failed", error=str(e))
        sys.exit(1)
    
    # Run pipeline
    try:
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=console,
        ) as progress:
            task = progress.add_task("Running pipeline...", total=None)
            
            if sra_id:
                feature_set = pipeline.run_sra(sra_id)
            else:
                feature_set = pipeline.run_fastq(fastq)
            
            progress.update(task, description="Pipeline completed successfully!")
        
        # Display results
        display_results(feature_set)
        
    except Exception as e:
        console.print(f"[red]Pipeline failed: {e}[/red]")
        logger.error("Pipeline execution failed", error=str(e))
        sys.exit(1)


@cli.command()
@click.option(
    "--sra-ids",
    required=True,
    help="Comma-separated list of SRA Run Accession IDs (e.g., SRR123456,SRR123457)",
    type=str,
)
@click.option(
    "--output-dir",
    default="./output",
    help="Output directory for results",
    type=click.Path(path_type=Path),
)
@click.option(
    "--config",
    help="Configuration file path",
    type=click.Path(exists=True, path_type=Path),
)
@click.option(
    "--threads",
    default=1,
    help="Number of threads to use",
    type=int,
)
@click.option(
    "--log-level",
    default="INFO",
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR"]),
    help="Logging level",
)
@click.option(
    "--log-file",
    help="Log file path",
    type=click.Path(path_type=Path),
)
@click.option(
    "--no-merge-vcfs",
    is_flag=True,
    help="Disable VCF merging (keep individual VCF files only)",
)
def batch(
    sra_ids: str,
    output_dir: Path,
    config: Optional[Path],
    threads: int,
    log_level: str,
    log_file: Optional[Path],
    no_merge_vcfs: bool,
):
    """Run the SRA to Features Pipeline on multiple SRA IDs with VCF merging."""
    
    # Parse SRA IDs
    sra_id_list = [s.strip() for s in sra_ids.split(",") if s.strip()]
    if not sra_id_list:
        console.print("[red]Error: No valid SRA IDs provided[/red]")
        sys.exit(1)
    
    console.print(f"[green]Processing {len(sra_id_list)} SRA IDs: {', '.join(sra_id_list)}[/green]")
    
    # Load configuration
    try:
        if config:
            pipeline_config = PipelineConfig(_env_file=config)
        else:
            pipeline_config = PipelineConfig()
        
        # Override config with CLI options
        pipeline_config.output_dir = output_dir
        pipeline_config.threads = threads
        pipeline_config.log_level = log_level
        pipeline_config.log_file = log_file
        
    except Exception as e:
        console.print(f"[red]Error loading configuration: {e}[/red]")
        sys.exit(1)
    
    # Setup logging
    logger = setup_logging(
        log_level=pipeline_config.log_level,
        log_file=pipeline_config.log_file,
    )
    
    try:
        # Initialize pipeline
        pipeline = Pipeline(pipeline_config, logger)
        
        # Run batch processing
        merge_vcfs = not no_merge_vcfs
        feature_sets = pipeline.run_batch(sra_id_list, merge_vcfs=merge_vcfs)
        
        # Display results
        console.print(f"\n[green]✓ Batch processing completed successfully![/green]")
        console.print(f"Processed {len(feature_sets)} samples")
        
        if merge_vcfs:
            merged_vcf = output_dir / "merged_variants.vcf.gz"
            if merged_vcf.exists():
                console.print(f"Merged VCF file: {merged_vcf}")
        
        # Display summary for each sample
        for feature_set in feature_sets:
            display_results(feature_set)
        
    except Exception as e:
        console.print(f"[red]Error during batch processing: {e}[/red]")
        sys.exit(1)


@cli.command()
@click.option(
    "--config",
    help="Configuration file path",
    type=click.Path(exists=True, path_type=Path),
)
def validate(config: Optional[Path]):
    """Validate pipeline configuration and dependencies."""
    
    try:
        if config:
            pipeline_config = PipelineConfig(_env_file=config)
        else:
            pipeline_config = PipelineConfig()
        
        console.print("[bold blue]Validating pipeline configuration...[/bold blue]")
        
        # Check configuration
        errors = pipeline_config.validate_setup()
        
        if errors:
            console.print("[red]Configuration validation failed:[/red]")
            for error in errors:
                console.print(f"  [red]• {error}[/red]")
            sys.exit(1)
        else:
            console.print("[green]✓ Configuration validation passed[/green]")
        
        # Display configuration summary
        display_config_summary(pipeline_config)
        
    except Exception as e:
        console.print(f"[red]Validation failed: {e}[/red]")
        sys.exit(1)


@cli.command()
@click.option(
    "--input-dir",
    required=True,
    help="Input directory containing pipeline results",
    type=click.Path(exists=True, path_type=Path),
)
@click.option(
    "--output-file",
    required=True,
    help="Output file path for ML feature table",
    type=click.Path(path_type=Path),
)
@click.option(
    "--format",
    default="csv",
    type=click.Choice(["csv", "tsv", "parquet", "json"]),
    help="Output format for feature table",
)
@click.option(
    "--log-level",
    default="INFO",
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR"]),
    help="Logging level",
)
def create_ml_table(
    input_dir: Path,
    output_file: Path,
    format: str,
    log_level: str,
):
    """Create ML-ready feature table from pipeline results."""
    
    # Setup logging
    logger = setup_logging(log_level=log_level)
    
    try:
        from ..utils.ml_features import create_ml_feature_table_from_directory
        
        # Create ML feature table
        df = create_ml_feature_table_from_directory(
            input_directory=input_dir,
            output_path=output_file,
            logger=logger,
            format=format
        )
        
        # Display summary
        console.print(f"\n[green]✓ ML feature table created successfully![/green]")
        console.print(f"Output file: {output_file}")
        console.print(f"Format: {format}")
        console.print(f"Shape: {df.shape[0]} samples × {df.shape[1]} features")
        
        # Show feature summary
        summary = {
            "Sample Count": df.shape[0],
            "Feature Count": df.shape[1],
            "Numeric Features": len(df.select_dtypes(include=['number']).columns),
            "Categorical Features": len(df.select_dtypes(include=['object']).columns),
            "Missing Values": df.isnull().sum().sum(),
        }
        
        table = Table(title="Feature Table Summary")
        table.add_column("Metric", style="cyan")
        table.add_column("Value", style="magenta")
        
        for metric, value in summary.items():
            table.add_row(metric, str(value))
        
        console.print(table)
        
        # Show sample of features
        console.print(f"\n[bold]Sample Features:[/bold]")
        feature_sample = df.columns[:10].tolist()
        console.print(", ".join(feature_sample))
        if len(df.columns) > 10:
            console.print(f"... and {len(df.columns) - 10} more features")
        
    except Exception as e:
        console.print(f"[red]Error creating ML feature table: {e}[/red]")
        sys.exit(1)


@cli.command()
@click.option(
    "--output-dir",
    default="./doc",
    help="Output directory for documentation",
    type=click.Path(path_type=Path),
)
def setup_docs(output_dir: Path):
    """Generate documentation for the pipeline."""
    
    console.print("[bold blue]Generating documentation...[/bold blue]")
    
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Generate basic documentation
        generate_documentation(output_dir)
        
        console.print(f"[green]Documentation generated in: {output_dir}[/green]")
        
    except Exception as e:
        console.print(f"[red]Documentation generation failed: {e}[/red]")
        sys.exit(1)


def display_results(feature_set):
    """Display pipeline results in a formatted table."""
    
    console.print("\n[bold green]Pipeline Results[/bold green]")
    
    # Summary statistics
    summary = feature_set.get_summary_stats()
    
    table = Table(title="Feature Summary")
    table.add_column("Metric", style="cyan")
    table.add_column("Value", style="magenta")
    
    table.add_row("Total Variants", str(summary["total_variants"]))
    table.add_row("Genes with Variants", str(summary["total_genes_with_variants"]))
    table.add_row("CNV Regions", str(summary["total_cnv_regions"]))
    table.add_row("Mean Variants per Bin", f"{summary['mean_variants_per_bin']:.2f}")
    table.add_row("Mapping Rate", f"{summary['mapping_rate']:.1f}%")
    table.add_row("Mean Coverage", f"{summary['mean_coverage']:.1f}x")
    
    console.print(table)
    
    # Quality assessment
    if feature_set.filter_by_quality():
        console.print("[green]✓ Sample passes quality thresholds[/green]")
    else:
        console.print("[yellow]⚠ Sample does not meet quality thresholds[/yellow]")


def display_config_summary(config: PipelineConfig):
    """Display configuration summary."""
    
    table = Table(title="Configuration Summary")
    table.add_column("Setting", style="cyan")
    table.add_column("Value", style="magenta")
    
    table.add_row("Base Directory", str(config.base_dir))
    table.add_row("Output Directory", str(config.output_dir))
    table.add_row("Reference FASTA", str(config.reference_fasta) if config.reference_fasta else "Not set")
    table.add_row("Reference GFF", str(config.reference_gff) if config.reference_gff else "Not set")
    table.add_row("Threads", str(config.threads))
    table.add_row("Bin Size (GVS)", str(config.bin_size_gvs))
    table.add_row("Bin Size (CNV)", str(config.bin_size_cnv))
    
    console.print(table)


def generate_documentation(output_dir: Path):
    """Generate basic documentation."""
    
    # Create README
    readme_content = f"""# SRA to Features Pipeline Documentation

## Overview
This pipeline extracts genomic features from SRA data for LLM training.

## Installation
```bash
pip install sra-to-features-pipeline
```

## Usage
```bash
# Run with SRA ID
sra-pipeline run --sra-id SRR123456 --output-dir ./results

# Run with FASTQ files
sra-pipeline run --fastq sample_1.fastq.gz sample_2.fastq.gz --output-dir ./results

# Validate configuration
sra-pipeline validate

# Generate documentation
sra-pipeline setup-docs
```

## Configuration
The pipeline uses environment variables or a configuration file for settings.

## Output
The pipeline generates:
- Feature set in JSON format
- Quality metrics
- Log files
- Intermediate files (optional)

For more information, see the full documentation.
"""
    
    with open(output_dir / "README.md", "w") as f:
        f.write(readme_content)


def main():
    """Main entry point for the CLI."""
    cli()


if __name__ == "__main__":
    main() 