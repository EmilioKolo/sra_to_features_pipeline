"""
Logging utilities for the SRA to Features Pipeline.
"""


from pathlib import Path
from structlog.stdlib import LoggerFactory
from typing import Any, Dict, Callable, Optional

import asyncio
import csv
import datetime
import psutil
import structlog
import sys
import time

try:
    import pynvml
    pynvml.nvmlInit()
    GPU_AVAILABLE = True
except:
    GPU_AVAILABLE = False


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
    
    # Make sure the log_file exists if it is not None
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
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
    
    def __init__(
        self,
        logger: structlog.BoundLogger,
        csv_path: str = "performance_log.csv"
    ):
        """
        Initialize the performance monitor.
        
        Args:
            logger: Structured logger instance
        """
        self.logger = logger
        self.metrics: Dict[str, float] = {}
        self.start_times: Dict[str, float] = {}

        # Async monitor control
        self._stop_event = asyncio.Event()
        self._monitor_task: Optional[asyncio.Task] = None
        self.current_section: str = "init"

        # CSV logging
        self.csv_path = csv_path

        # peaks
        self.peak_ram_mb = 0.0
        self.peak_cpu = 0.0
        self.peak_gpu_mem_mb = 0.0
        self.peak_gpu_util = 0.0

        self.peak_ram_mb_section = 0.0
        self.peak_cpu_section = 0.0
        self.peak_gpu_mem_mb_section = 0.0
        self.peak_gpu_util_section = 0.0
    
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
    
    def init_csv(self):
        header = [
            "timestamp",
            "section",
            "cpu_percent",
            "ram_percent",
            "proc_ram_mb",
            "gpu_mem_mb",
            "gpu_util"
        ]
        with open(self.csv_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(header)
    
    async def _monitor_loop(self, interval: float):
        """Run until stop is requested."""
        process = psutil.Process()

        while not self._stop_event.is_set():

            timestamp = datetime.datetime.now().isoformat()
            cpu_total = psutil.cpu_percent(interval=None)
            ram_total = psutil.virtual_memory().percent
            proc_ram = process.memory_info().rss / 1024**2  # MB

            gpu_mem = None
            gpu_util = None
            if GPU_AVAILABLE:
                handle = pynvml.nvmlDeviceGetHandleByIndex(0)
                mem_info = pynvml.nvmlDeviceGetMemoryInfo(handle)
                util = pynvml.nvmlDeviceGetUtilizationRates(handle)
                gpu_mem = mem_info.used / 1024**2
                gpu_util = util.gpu

            # Update peaks
            self.peak_cpu = max(self.peak_cpu, cpu_total)
            self.peak_cpu_section = max(self.peak_cpu_section, cpu_total)
            self.peak_ram_mb = max(self.peak_ram_mb, proc_ram)
            self.peak_ram_mb_section = max(self.peak_ram_mb_section, proc_ram)
            if GPU_AVAILABLE and gpu_mem is not None:
                self.peak_gpu_mem_mb = max(self.peak_gpu_mem_mb, gpu_mem)
                self.peak_gpu_mem_mb_section = max(self.peak_gpu_mem_mb_section, gpu_mem)
                self.peak_gpu_util = max(self.peak_gpu_util, gpu_util)
                self.peak_gpu_util_section = max(self.peak_gpu_util_section, gpu_util)

            # Write CSV row
            with open(self.csv_path, "a", newline="") as f:
                writer = csv.writer(f)
                writer.writerow([
                    timestamp,
                    self.current_section,
                    cpu_total,
                    ram_total,
                    proc_ram,
                    gpu_mem,
                    gpu_util
                ])

            await asyncio.sleep(interval)
    
    def start_monitoring(
        self,
        interval: float = 1.0,
        section: str = "global"
    ):
        """Start asynchronous background monitoring."""
        self._stop_event.clear()
        self.current_section = section
        self._monitor_task = asyncio.create_task(self._monitor_loop(interval))

    async def stop_monitoring(self):
        """Stop async monitor and wait for its task to finish."""
        self._stop_event.set()
        if self._monitor_task:
            await self._monitor_task
    
    async def wrap_section(self, section_name: str, coro):
        """
        Run an async section while monitoring and tagging metrics
        with its section name.
        """
        self.current_section = section_name
        return await coro

    def report_peaks(self):
        self.logger.info(
            "Peak resource usage",
            peak_cpu_percent=self.peak_cpu,
            peak_ram_mb=self.peak_ram_mb,
            peak_gpu_mem_mb=self.peak_gpu_mem_mb,
            peak_gpu_util=self.peak_gpu_util
        )
        if self.current_section != "global":
            self.logger.info(
                f"Peak resource usage ({self.current_section})",
                peak_cpu_percent_section=self.peak_cpu_section,
                peak_ram_mb_section=self.peak_ram_mb_section,
                peak_gpu_mem_mb_section=self.peak_gpu_mem_mb_section,
                peak_gpu_util_section=self.peak_gpu_util_section
            )
        # Reset section peaks
        self.reset_section_peaks()
    
    def reset_section_peaks(self):
        """Reset peak metrics for the current section."""
        self.peak_ram_mb_section = 0.0
        self.peak_cpu_section = 0.0
        self.peak_gpu_mem_mb_section = 0.0
        self.peak_gpu_util_section = 0.0

    async def section(self, name: str, func: Callable, *args, **kwargs):
        # Reset section peaks
        self.reset_section_peaks()
        # Start the section timer
        self.current_section = name
        self.start_timer(name)

        result = await asyncio.to_thread(func, *args, **kwargs)

        self.report_peaks()

        self.stop_timer(name)
        return result


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