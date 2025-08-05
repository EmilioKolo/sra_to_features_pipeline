"""
Core pipeline modules for the SRA to Features Pipeline.
"""

# Lazy imports to avoid dependency issues
def get_pipeline():
    """Get the Pipeline class."""
    from .pipeline import Pipeline
    return Pipeline

__all__ = ["get_pipeline"] 