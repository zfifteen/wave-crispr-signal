"""
Data management and registry for wave-crispr-signal.

This package provides dataset loading, caching, and provenance tracking.
"""

from .registry import DatasetRegistry, load_dataset, register_dataset

__all__ = [
    "DatasetRegistry",
    "load_dataset",
    "register_dataset",
]
