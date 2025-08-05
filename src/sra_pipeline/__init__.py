"""
SRA to Features Pipeline

A bioinformatics pipeline for extracting features from SRA data for LLM training.
"""

__version__ = "1.0.0"
__author__ = "Bioinformatics Team"
__email__ = "team@example.com"

from .core.pipeline import Pipeline
from .config.settings import PipelineConfig
from .models.features import FeatureSet

__all__ = ["Pipeline", "PipelineConfig", "FeatureSet"] 