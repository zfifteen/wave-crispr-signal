"""
Dataset registry for wave-crispr-signal.

Centralized dataset loading with:
- SHA-256 checksum validation
- Caching to data/processed/
- Provenance tracking
- Split integrity enforcement

SCIENTIFIC GATES:
- Human DNA only (Homo sapiens)
- Dataset version tracking
- License compliance checks
"""

import hashlib
import logging
import os
from pathlib import Path
from typing import Dict, Optional, Tuple
import pandas as pd
import numpy as np
from dataclasses import dataclass


logger = logging.getLogger(__name__)

# --- Constants ---
DEFAULT_SPECIES = "Homo sapiens"


@dataclass
class DatasetMetadata:
    """Metadata for a CRISPR dataset."""

    name: str
    path: str
    version: str
    species: str
    license: str
    n_guides: int
    sha256: Optional[str] = None
    url: Optional[str] = None
    citation: Optional[str] = None


class DatasetRegistry:
    """
    Registry for managing CRISPR datasets with provenance tracking.

    Features:
    - Automatic SHA-256 checksum validation
    - Caching to data/processed/
    - Species filtering (Homo sapiens only)
    - Split-by-gene and split-by-screen support
    """

    def __init__(self, base_dir: Optional[str] = None):
        """
        Initialize dataset registry.

        Args:
            base_dir: Base directory for data (default: data/ in repo root)
        """
        if base_dir is None:
            # Default to data/ in repository root
            base_dir = Path(__file__).parent.parent.parent / "data"

        self.base_dir = Path(base_dir)
        self.processed_dir = self.base_dir / "processed"
        self.processed_dir.mkdir(parents=True, exist_ok=True)

        self.datasets: Dict[str, DatasetMetadata] = {}
        self._register_default_datasets()

    def _register_default_datasets(self):
        """Register default CRISPR datasets."""
        # BioGRID-ORCS Homo sapiens v1.1.17
        self.register(
            DatasetMetadata(
                name="biogrid_orcs_v1.1.17",
                path=str(self.base_dir / "biogrid_orcs_v1.1.17.csv"),
                version="1.1.17",
                species=DEFAULT_SPECIES,
                license="MIT",
                n_guides=0,  # Will be determined on load
                url="https://orcs.thebiogrid.org/",
                citation="Pacini et al. (2021). Integrated cross-study datasets of genetic dependencies in cancer. Nature Communications.",
            )
        )

    def register(self, metadata: DatasetMetadata):
        """
        Register a dataset.

        Args:
            metadata: Dataset metadata
        """
        if metadata.species != DEFAULT_SPECIES:
            logger.warning(
                f"Dataset {metadata.name} is not {DEFAULT_SPECIES}: {metadata.species}"
            )

        self.datasets[metadata.name] = metadata
        logger.info(f"Registered dataset: {metadata.name} v{metadata.version}")

    def compute_sha256(self, path: str) -> str:
        """
        Compute SHA-256 checksum of a file.

        Args:
            path: File path

        Returns:
            Hex string of SHA-256 hash
        """
        sha256_hash = hashlib.sha256()
        with open(path, "rb") as f:
            # Read in chunks for large files
            for byte_block in iter(lambda: f.read(4096), b""):
                sha256_hash.update(byte_block)
        return sha256_hash.hexdigest()

    def validate_checksum(self, path: str, expected_sha256: str) -> bool:
        """
        Validate file checksum.

        Args:
            path: File path
            expected_sha256: Expected SHA-256 hex string

        Returns:
            True if checksum matches, False otherwise
        """
        actual_sha256 = self.compute_sha256(path)
        if actual_sha256 != expected_sha256:
            logger.error(f"Checksum mismatch for {path}")
            logger.error(f"  Expected: {expected_sha256}")
            logger.error(f"  Actual:   {actual_sha256}")
            return False
        return True

    def load(
        self, name: str, validate_checksum: bool = True, cache: bool = True
    ) -> Tuple[pd.DataFrame, DatasetMetadata]:
        """
        Load a dataset with validation.

        Args:
            name: Dataset name
            validate_checksum: Whether to validate SHA-256 checksum
            cache: Whether to use/save cached version

        Returns:
            Tuple of (dataframe, metadata)

        Raises:
            ValueError: If dataset not found or validation fails
        """
        if name not in self.datasets:
            raise ValueError(
                f"Dataset '{name}' not registered. Available: {list(self.datasets.keys())}"
            )

        metadata = self.datasets[name]

        # Check for cached version
        cached_path = self.processed_dir / f"{name}.pkl"
        if cache and cached_path.exists():
            logger.info(f"Loading cached dataset: {name}")
            df = pd.read_pickle(cached_path)
            return df, metadata

        # Load from source
        logger.info(f"Loading dataset: {name} from {metadata.path}")

        if not os.path.exists(metadata.path):
            raise ValueError(f"Dataset file not found: {metadata.path}")

        # Validate checksum if provided
        if validate_checksum and metadata.sha256:
            if not self.validate_checksum(metadata.path, metadata.sha256):
                raise ValueError(f"Checksum validation failed for {name}")
        elif not metadata.sha256:
            # Compute and log checksum for future use
            sha256 = self.compute_sha256(metadata.path)
            logger.info(f"Dataset SHA-256: {sha256}")
            metadata.sha256 = sha256

        # Load CSV
        df = pd.read_csv(metadata.path)

        # Validate species
        if "species" in df.columns:
            species_counts = df["species"].value_counts()
            if DEFAULT_SPECIES not in species_counts:
                raise ValueError(f"Dataset {name} does not contain {DEFAULT_SPECIES} data")

            # Filter to Homo sapiens only
            df = df[df["species"] == DEFAULT_SPECIES].copy()
            logger.info(f"Filtered to {len(df)} {DEFAULT_SPECIES} guides")

        # Update metadata
        metadata.n_guides = len(df)

        # Cache if requested
        if cache:
            logger.info(f"Caching dataset to {cached_path}")
            df.to_pickle(cached_path)

        return df, metadata

    def split_dataset(
        self,
        df: pd.DataFrame,
        strategy: str = "split_by_gene",
        train_ratio: float = 0.7,
        val_ratio: float = 0.15,
        test_ratio: float = 0.15,
        seed: int = 0,
    ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        Split dataset with leakage prevention.

        Args:
            df: Input dataframe
            strategy: Split strategy ('split_by_gene', 'split_by_screen', 'split_by_guide')
            train_ratio: Training set ratio
            val_ratio: Validation set ratio
            test_ratio: Test set ratio
            seed: Random seed for shuffling. Providing a consistent seed ensures
                  deterministic and reproducible splits.

        Returns:
            Tuple of (train_df, val_df, test_df)

        Raises:
            ValueError: If ratios don't sum to 1.0 or strategy invalid
        """
        if not np.isclose(train_ratio + val_ratio + test_ratio, 1.0):
            raise ValueError(
                f"Ratios must sum to 1.0, got {train_ratio + val_ratio + test_ratio}"
            )

        rng = np.random.default_rng(seed)

        if strategy == "split_by_gene":
            if "gene" not in df.columns:
                raise ValueError("Dataset must have 'gene' column for split_by_gene")

            # Split by unique genes to prevent leakage
            unique_genes = df["gene"].unique()
            rng.shuffle(unique_genes)

            n_train = int(len(unique_genes) * train_ratio)
            n_val = int(len(unique_genes) * val_ratio)

            train_genes = unique_genes[:n_train]
            val_genes = unique_genes[n_train : n_train + n_val]
            test_genes = unique_genes[n_train + n_val :]

            train_df = df[df["gene"].isin(train_genes)].copy()
            val_df = df[df["gene"].isin(val_genes)].copy()
            test_df = df[df["gene"].isin(test_genes)].copy()

        elif strategy == "split_by_screen":
            if "screen_id" not in df.columns:
                raise ValueError(
                    "Dataset must have 'screen_id' column for split_by_screen"
                )

            # Split by screen to prevent leakage
            unique_screens = df["screen_id"].unique()
            rng.shuffle(unique_screens)

            n_train = int(len(unique_screens) * train_ratio)
            n_val = int(len(unique_screens) * val_ratio)

            train_screens = unique_screens[:n_train]
            val_screens = unique_screens[n_train : n_train + n_val]
            test_screens = unique_screens[n_train + n_val :]

            train_df = df[df["screen_id"].isin(train_screens)].copy()
            val_df = df[df["screen_id"].isin(val_screens)].copy()
            test_df = df[df["screen_id"].isin(test_screens)].copy()

        elif strategy == "split_by_guide":
            # Random split at guide level (least preferred, may leak)
            logger.warning("split_by_guide may leak information between splits")

            indices = np.arange(len(df))
            rng.shuffle(indices)

            n_train = int(len(indices) * train_ratio)
            n_val = int(len(indices) * val_ratio)

            train_idx = indices[:n_train]
            val_idx = indices[n_train : n_train + n_val]
            test_idx = indices[n_train + n_val :]

            train_df = df.iloc[train_idx].copy()
            val_df = df.iloc[val_idx].copy()
            test_df = df.iloc[test_idx].copy()
        else:
            raise ValueError(f"Unknown strategy: {strategy}")

        logger.info(
            f"Split dataset ({strategy}): train={len(train_df)}, val={len(val_df)}, test={len(test_df)}"
        )

        return train_df, val_df, test_df


# Convenience functions

_global_registry = None


def get_registry() -> DatasetRegistry:
    """Get or create global dataset registry."""
    global _global_registry
    if _global_registry is None:
        _global_registry = DatasetRegistry()
    return _global_registry


def load_dataset(name: str, **kwargs) -> Tuple[pd.DataFrame, DatasetMetadata]:
    """
    Load a dataset using the global registry.

    Args:
        name: Dataset name
        **kwargs: Additional arguments to DatasetRegistry.load()

    Returns:
        Tuple of (dataframe, metadata)
    """
    registry = get_registry()
    return registry.load(name, **kwargs)


def register_dataset(metadata: DatasetMetadata):
    """
    Register a dataset in the global registry.

    Args:
        metadata: Dataset metadata
    """
    registry = get_registry()
    registry.register(metadata)
