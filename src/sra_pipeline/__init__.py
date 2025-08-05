"""
SRA to Features Pipeline

A bioinformatics pipeline for extracting features from SRA data for LLM training.
"""

__version__ = "1.0.0"
__author__ = "Bioinformatics Team"
__email__ = "team@example.com"

# Lazy imports to avoid dependency issues
def get_pipeline():
    """Get the Pipeline class."""
    from .core.pipeline import Pipeline
    return Pipeline

def get_pipeline_config():
    """Get the PipelineConfig class."""
    from .config.settings import PipelineConfig
    return PipelineConfig

def get_feature_set():
    """Get the FeatureSet class."""
    from .models.features import FeatureSet
    return FeatureSet

__all__ = ["get_pipeline", "get_pipeline_config", "get_feature_set"] 