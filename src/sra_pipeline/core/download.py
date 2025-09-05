"""
SRA data download functionality.
"""

import subprocess
import time
import random
from pathlib import Path
from typing import List, Dict, Any, Optional
import shutil
import structlog
import requests

from ..utils import log_command, log_error


def download_sra_data(
    sra_id: str,
    output_dir: Path,
    logger: structlog.BoundLogger,
    max_retries: int = 5
) -> List[Path]:
    """
    Download FASTQ files for an SRA ID.
    
    Args:
        sra_id: SRA accession ID
        output_dir: Directory to save FASTQ files
        logger: Logger instance
        max_retries: Maximum number of retry attempts
        
    Returns:
        List of downloaded FASTQ file paths
        
    Raises:
        RuntimeError: If download fails after all retries
    """
    logger.info("Starting SRA data download", sra_id=sra_id, output_dir=str(output_dir))
    
    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Get SRA metadata
    sra_metadata = _get_sra_metadata(sra_id, logger)
    if not sra_metadata:
        raise RuntimeError(f"Failed to retrieve metadata for SRA ID: {sra_id}")
    
    # Determine FASTQ file names
    fastq_files = _determine_fastq_files(sra_id, sra_metadata, output_dir)

    # Download FASTQ files
    for attempt in range(max_retries):
        try:
            logger.info(f"Starting download attempt {attempt+1}.", 
                        sra_id=str(sra_id), attempt=attempt+1)
            
            _download_fastq_files(sra_id, output_dir, logger)
            break
            
        except Exception as e:
            logger.warning(f"Download attempt {attempt+1} failed", 
                            sra_id=str(sra_id), error=str(e))
            
            if attempt < max_retries - 1:
                # Wait before retry with exponential backoff
                wait_time = (2 ** attempt) + random.uniform(0, 1)
                logger.info(
                    f"Waiting {wait_time:.1f} seconds before retry"
                )
                time.sleep(wait_time)
            else:
                err = f"Failed to download {sra_id} after {max_retries} attempts"
                raise RuntimeError(err)

    logger.info("SRA data download completed", 
                sra_id=sra_id, downloaded_files=fastq_files)
    
    return fastq_files


def _get_sra_metadata(sra_id: str, logger: structlog.BoundLogger) -> Optional[Dict[str, Any]]:
    """Get metadata for an SRA ID from the ENA API."""
    url = "https://www.ebi.ac.uk/ena/portal/api/search"
    params = {
        "result": "read_run",
        "query": f"run_accession={sra_id}",
        "fields": "run_accession,fastq_ftp,sra_ftp,experiment_accession,sample_accession,study_accession,library_name,library_strategy,library_source,library_selection,instrument_platform,instrument_model,base_count,read_count,scientific_name,tax_id",
        "format": "json",
        "limit": 1
    }
    
    try:
        logger.info("Fetching SRA metadata", sra_id=sra_id, url=url)
        response = requests.get(url, params=params, timeout=30)
        response.raise_for_status()
        
        data = response.json()
        if not data:
            logger.error("No metadata found for SRA ID", sra_id=sra_id)
            return None
        
        metadata = data[0]
        logger.info("SRA metadata retrieved successfully", 
                    sra_id=sra_id, 
                    fastq_ftp=metadata.get("fastq_ftp"),
                    read_count=metadata.get("read_count"))
        
        return metadata
        
    except requests.exceptions.RequestException as e:
        logger.error("Failed to fetch SRA metadata", sra_id=sra_id, error=str(e))
        return None


def _determine_fastq_files(
    sra_id: str, 
    metadata: Dict[str, Any], 
    output_dir: Path
) -> List[Path]:
    """Determine FASTQ file names from SRA metadata."""
    fastq_ftp = metadata.get("fastq_ftp", "")
    
    if not fastq_ftp:
        raise ValueError(f"No FASTQ FTP links found for SRA ID: {sra_id}")
    
    ftp_links = fastq_ftp.split(";")
    fastq_files = []
    
    for i, ftp_link in enumerate(ftp_links):
        if not ftp_link.strip():
            continue
            
        # Extract filename from FTP link
        filename = ftp_link.split("/")[-1]
        
        # For paired-end data, ensure proper naming
        if len(ftp_links) == 2:
            if i == 0:
                filename = f"{sra_id}_1.fastq.gz"
            else:
                filename = f"{sra_id}_2.fastq.gz"
        else:
            filename = f"{sra_id}.fastq.gz"
        
        fastq_files.append(output_dir / filename)
    
    return fastq_files


def _download_fastq_files(
    sra_id: str,
    output_dir: Path,
    logger: structlog.BoundLogger
) -> None:
    """
    Download and convert an SRA file to FASTQ using prefetch 
    and fasterq-dump.
    """
    # Define the base directory for SRA download
    sra_dir = output_dir / sra_id

    # Prefetch the SRA archive
    logger.info("Starting prefetch", sra_id=sra_id)
    prefetch_cmd = [
        "prefetch",
        sra_id,
        "--output-directory", 
        str(output_dir)
    ]
    try:
        subprocess.run(
            prefetch_cmd,
            capture_output=True,
            text=True,
            check=True,
            #timeout=3600  # 1 hour timeout
        )
        logger.info("Prefetch completed", sra_id=sra_id, 
                    sra_dir=str(sra_dir))
    except subprocess.TimeoutExpired:
        raise RuntimeError(f"Download timed out for SRA ID: {sra_id}")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"Prefetch failed for SRA ID {sra_id}: {e.stderr}"
        )

    # Convert the SRA archive to FASTQ using fasterq-dump
    logger.info("Starting fasterq-dump", sra_id=sra_id)
    fasterq_dump_cmd = [
        "fasterq-dump",
        sra_id,
        "--threads", "8",
        "--outdir", str(output_dir),
        "--split-files" # Flag is ignored for single-end files
    ]
    try:
        subprocess.run(
            fasterq_dump_cmd,
            capture_output=True,
            text=True,
            check=True,
            #timeout=3600  # 1 hour timeout
        )
        logger.info("Fasterq-dump completed", sra_id=sra_id)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"Fasterq-dump failed for SRA ID {sra_id}: {e.stderr}"
        )

    # Clean up the SRA archive folder
    if sra_dir.exists():
        logger.info("Starting cleanup of SRA directory", sra_id=sra_id)
        try:
            shutil.rmtree(sra_dir)
            logger.info("SRA directory removed", sra_id=sra_id)
        except OSError as e:
            logger.error(f"Error removing SRA directory for {sra_id}: {e}", 
                         sra_dir=sra_dir)
    
    # Gzip the output files
    logger.info("Starting gzipping of FASTQ files", sra_id=sra_id)
    try:
        fastq_files = list(output_dir.glob(f"{sra_id}*.fastq"))
        if not fastq_files:
            logger.warning("No FASTQ files found to gzip", sra_id=sra_id)
        
        for fastq_file in fastq_files:
            gzip_cmd = ["gzip", "--force", str(fastq_file)]
            subprocess.run(gzip_cmd, check=True)
            logger.info(f"Gzipped file: {fastq_file.name}")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Gzipping failed for SRA ID {sra_id}: {e.stderr}")

    return None


def validate_fastq_files(fastq_files: List[Path], logger: structlog.BoundLogger) -> bool:
    """Validate that downloaded FASTQ files are complete and readable."""
    for fastq_file in fastq_files:
        if not fastq_file.exists():
            logger.error("FASTQ file does not exist", file=str(fastq_file))
            return False
        
        if fastq_file.stat().st_size == 0:
            logger.error("FASTQ file is empty", file=str(fastq_file))
            return False
        
        # Check if file is gzipped and readable
        try:
            import gzip
            with gzip.open(fastq_file, 'rt') as f:
                # Read first few lines to check format
                for i, line in enumerate(f):
                    if i >= 4:  # Check first read
                        break
                    if i % 4 == 0 and not line.startswith('@'):
                        logger.error("Invalid FASTQ format", file=str(fastq_file))
                        return False
        except Exception as e:
            logger.error("Failed to read FASTQ file", file=str(fastq_file), error=str(e))
            return False
    
    logger.info("FASTQ file validation passed", files=[str(f) for f in fastq_files])
    return True 