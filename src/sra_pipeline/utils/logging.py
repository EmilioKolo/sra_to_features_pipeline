"""
Logging utilities for the SRA to Features Pipeline.
"""

import sys
import time
from pathlib import Path
from typing import Any, Dict, Optional
import structlog
from structlog.stdlib import LoggerFactory


def setup_logging(
    log_level: str = "INFO",
    log_file: Optional[Path] = None,
    log_format: str = "json"
) -> structlog.BoundLogger:
    """
    Set up structured logging for the pipeline.
    
    Args:
        log_level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        log_file: Optional log file path
        log_format: Log format ("json" or "console")
    
    Returns:
        Configured logger instance
    """
    # Configure structlog
    structlog.configure(
        processors=[
            structlog.stdlib.filter_by_level,
            structlog.stdlib.add_logger_name,
            structlog.stdlib.add_log_level,
            structlog.stdlib.PositionalArgumentsFormatter(),
            structlog.processors.TimeStamper(fmt="iso"),
            structlog.processors.StackInfoRenderer(),
            structlog.processors.format_exc_info,
            structlog.processors.UnicodeDecoder(),
            structlog.processors.JSONRenderer() if log_format == "json" else structlog.dev.ConsoleRenderer(),
        ],
        context_class=dict,
        logger_factory=LoggerFactory(),
        wrapper_class=structlog.stdlib.BoundLogger,
        cache_logger_on_first_use=True,
    )
    
    # Configure standard library logging
    import logging
    logging.basicConfig(
        format="%(message)s",
        stream=sys.stdout if log_file is None else open(log_file, "w"),
        level=getattr(logging, log_level.upper()),
    )
    
    return structlog.get_logger()


class PipelineLogger:
    """Context manager for pipeline logging with performance tracking."""
    
    def __init__(self, logger: structlog.BoundLogger, operation: str):
        """
        Initialize the pipeline logger.
        
        Args:
            logger: Structured logger instance
            operation: Name of the operation being logged
        """
        self.logger = logger
        self.operation = operation
        self.start_time = None
        self.context: Dict[str, Any] = {}
    
    def __enter__(self):
        """Enter the logging context."""
        self.start_time = time.time()
        self.logger.info(
            f"Starting {self.operation}",
            operation=self.operation,
            **self.context
        )
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Exit the logging context."""
        duration = time.time() - self.start_time
        
        if exc_type is None:
            self.logger.info(
                f"Completed {self.operation}",
                operation=self.operation,
                duration_seconds=duration,
                status="success",
                **self.context
            )
        else:
            self.logger.error(
                f"Failed {self.operation}",
                operation=self.operation,
                duration_seconds=duration,
                status="error",
                error_type=exc_type.__name__,
                error_message=str(exc_val),
                **self.context
            )
    
    def add_context(self, **kwargs):
        """Add context information to the logger."""
        self.context.update(kwargs)
        return self
    
    def log_progress(self, message: str, **kwargs):
        """Log progress information."""
        self.logger.info(
            message,
            operation=self.operation,
            **{**self.context, **kwargs}
        )


class PerformanceMonitor:
    """Monitor and log performance metrics."""
    
    def __init__(self, logger: structlog.BoundLogger):
        """
        Initialize the performance monitor.
        
        Args:
            logger: Structured logger instance
        """
        self.logger = logger
        self.metrics: Dict[str, float] = {}
        self.start_times: Dict[str, float] = {}
    
    def start_timer(self, name: str):
        """Start a timer for a named operation."""
        self.start_times[name] = time.time()
    
    def stop_timer(self, name: str) -> float:
        """Stop a timer and return the duration."""
        if name not in self.start_times:
            raise ValueError(f"Timer '{name}' was not started")
        
        duration = time.time() - self.start_times[name]
        self.metrics[name] = duration
        
        self.logger.info(
            f"Operation '{name}' completed",
            operation=name,
            duration_seconds=duration
        )
        
        del self.start_times[name]
        return duration
    
    def log_memory_usage(self):
        """Log current memory usage."""
        try:
            import psutil
            process = psutil.Process()
            memory_info = process.memory_info()
            
            self.logger.info(
                "Memory usage",
                memory_rss_mb=memory_info.rss / 1024 / 1024,
                memory_vms_mb=memory_info.vms / 1024 / 1024,
                memory_percent=process.memory_percent()
            )
        except ImportError:
            self.logger.warning("psutil not available, cannot log memory usage")
    
    def log_system_info(self):
        """Log system information."""
        try:
            import psutil
            
            self.logger.info(
                "System information",
                cpu_count=psutil.cpu_count(),
                memory_total_gb=psutil.virtual_memory().total / 1024 / 1024 / 1024,
                disk_usage_percent=psutil.disk_usage('/').percent
            )
        except ImportError:
            self.logger.warning("psutil not available, cannot log system info")
    
    def get_summary(self) -> Dict[str, float]:
        """Get a summary of all recorded metrics."""
        return self.metrics.copy()


def log_command(logger: structlog.BoundLogger, command: str, **kwargs):
    """Log a command being executed."""
    logger.info(
        "Executing command",
        command=command,
        **kwargs
    )


def log_file_operation(logger: structlog.BoundLogger, operation: str, file_path: Path, **kwargs):
    """Log a file operation."""
    logger.info(
        f"File {operation}",
        operation=operation,
        file_path=str(file_path),
        file_size_mb=file_path.stat().st_size / 1024 / 1024 if file_path.exists() else 0,
        **kwargs
    )


def log_error(logger: structlog.BoundLogger, error: Exception, context: Optional[Dict[str, Any]] = None):
    """Log an error with context."""
    logger.error(
        "Pipeline error",
        error_type=type(error).__name__,
        error_message=str(error),
        context=context or {}
    ) 