"""
CRISPR Applications Package

This package provides CRISPR guide design and analysis tools using
signal-theoretic DNA analysis methods.
"""

from .crispr_guide_designer import CRISPRGuideDesigner
from .wave_crispr_metrics import WaveCRISPRMetrics
from .crispr_visualization import CRISPRVisualizer

__version__ = "1.0.0"
__author__ = "Wave CRISPR Signal Team"

__all__ = [
    'CRISPRGuideDesigner',
    'WaveCRISPRMetrics', 
    'CRISPRVisualizer'
]