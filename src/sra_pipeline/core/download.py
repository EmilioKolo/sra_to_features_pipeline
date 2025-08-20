"""
SRA data download functionality.
"""

import os
import subprocess
import time
import random
from pathlib import Path
from typing import List, Dict, Any, Optional
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
    downloaded_files = []
    for i, fastq_file in enumerate(fastq_files):
        success = False
        for attempt in range(max_retries):
            try:
                logger.info(f"Downloading FASTQ file {i+1}/{len(fastq_files)}", 
                           file=str(fastq_file), attempt=attempt+1)
                
                _download_fastq_file(sra_id, fastq_file, output_dir, logger)
                downloaded_files.append(fastq_file)
                success = True
                break
                
            except Exception as e:
                logger.warning(f"Download attempt {attempt+1} failed", 
                              file=str(fastq_file), error=str(e))
                
                if attempt < max_retries - 1:
                    # Wait before retry with exponential backoff
                    wait_time = (2 ** attempt) + random.uniform(0, 1)
                    logger.info(f"Waiting {wait_time:.1f} seconds before retry")
                    time.sleep(wait_time)
                else:
                    raise RuntimeError(f"Failed to download {fastq_file} after {max_retries} attempts")
    
    logger.info("SRA data download completed", 
                sra_id=sra_id, downloaded_files=[str(f) for f in downloaded_files])
    
    return downloaded_files


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


def _download_fastq_file(
    sra_id: str, 
    fastq_file: Path, 
    output_dir: Path, 
    logger: structlog.BoundLogger
):
    """Download a single FASTQ file using fastq-dump."""
    # Use fastq-dump to download and convert to FASTQ format
    cmd = [
        "fastq-dump",
        "--gzip",
        "--outdir", str(output_dir),
        sra_id
    ]
    
    # For paired-end data, add split-files flag
    if "_1.fastq.gz" in str(fastq_file) or "_2.fastq.gz" in str(fastq_file):
        cmd.insert(1, "--split-files")
    
    log_command(logger, " ".join(cmd), sra_id=sra_id, output_file=str(fastq_file))
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
            #timeout=3600  # 1 hour timeout
        )
        
        logger.info("FASTQ download completed", 
                    sra_id=sra_id, 
                    output_file=str(fastq_file),
                    stdout=result.stdout[-200:] if result.stdout else None)
        
    except subprocess.TimeoutExpired:
        raise RuntimeError(f"Download timed out for SRA ID: {sra_id}")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Download failed for SRA ID {sra_id}: {e.stderr}")


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