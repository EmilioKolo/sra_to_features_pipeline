"""
Configuration management for the SRA to Features Pipeline.
"""

# Lazy import to avoid dependency issues
def get_pipeline_config():
    """Get the PipelineConfig class."""
    from .settings import PipelineConfig
    return PipelineConfig

__all__ = ["get_pipeline_config"] 