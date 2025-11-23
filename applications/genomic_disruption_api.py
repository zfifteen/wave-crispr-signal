#!/usr/bin/env python3
"""
Genomic Disruption Analyzer: Z-Framework SaaS Pipeline for CRISPR gRNA Optimization

This module provides a comprehensive REST API for CRISPR guide RNA design, scoring,
and off-target profiling using FFT-encoded DNA waveforms with θ′(n,k) phase weighting.

Core Features:
- Single guide scoring with disruption metrics
- Batch processing (target: 1,000 guides/min)
- Guide design from target sequences
- Off-target profiling via spectral signatures
- Repair pathway estimation from sidelobe metrics

Scientific Gates Enforced:
- Human DNA/RNA only (A/C/G/T or A/C/G/U, N allowed)
- Fail-fast validation with clear errors
- Bootstrap CI (≥1,000 resamples) for statistical validity
- Seed control for reproducibility
- Z-invariant normalization: Z = A(B/c)
- Geometric resolution: θ′(n,k) with k ≈ 0.3

API Endpoints:
- POST /api/v1/score - Score single guide
- POST /api/v1/batch - Batch score multiple guides
- POST /api/v1/design - Design guides from target sequence
- POST /api/v1/offtarget - Analyze off-target potential
- GET /api/v1/health - Health check
- GET /api/v1/info - API information

Usage:
    # As module
    from applications.genomic_disruption_api import GenomicDisruptionAPI
    api = GenomicDisruptionAPI()
    
    # Standalone server
    python applications/genomic_disruption_api.py

Example:
    # Score a single guide
    curl -X POST http://localhost:8080/api/v1/score \\
      -H "Content-Type: application/json" \\
      -d '{"guide": "GCTGCGGAGACCTGGAGAGA", "k": 0.3}'

References:
- Z Framework: docs/Z_FRAMEWORK.md
- Phase-Weighted Scorecard: docs/PHASE_WEIGHTED_SCORECARD.md
- FFT Golden Ratio: docs/FFT_GOLDEN_RATIO_CRISPR.md
"""

import sys
import os
import json
import logging
import time
import hashlib
from datetime import datetime
from typing import List, Dict, Optional, Union, Any
from pathlib import Path

# Add parent directory for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from scipy.stats import bootstrap

# Import core modules
from applications.phase_weighted_scorecard import (
    PhaseWeightedScorecard,
    validate_dna_sequence,
    K_STAR,
    PHI,
    E_SQUARED,
)
from applications.fft_crispr_disruption import FFTCRISPRDisruptionAnalyzer

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class DisruptionAnalyzer:
    """
    Core disruption analysis engine for CRISPR guide scoring.
    
    Integrates phase-weighted spectral analysis with FFT-based disruption metrics
    to provide comprehensive guide scoring and off-target profiling.
    """
    
    def __init__(
        self,
        k: float = K_STAR,
        phi_period: float = 21.0,
        seed: Optional[int] = None,
    ):
        """
        Initialize DisruptionAnalyzer.
        
        Args:
            k: Resolution exponent for θ′(n,k) (default: 0.3)
            phi_period: Geometric period φ for resolution (default: 21 for 21-nt guides)
            seed: Random seed for reproducibility
        """
        self.k = k
        self.phi_period = phi_period
        self.seed = seed
        
        if seed is not None:
            np.random.seed(seed)
        
        # Initialize FFT analyzer
        self.fft_analyzer = FFTCRISPRDisruptionAnalyzer(phi_period=phi_period, k=k)
        
        logger.info(f"Initialized DisruptionAnalyzer with k={k}, φ_period={phi_period}")
    
    def score_guide(
        self,
        guide: str,
        target: Optional[str] = None,
        compute_ci: bool = False,
        n_bootstrap: int = 1000,
    ) -> Dict[str, Any]:
        """
        Score a single CRISPR guide with spectral disruption metrics.
        
        Args:
            guide: Guide RNA sequence (20-30 nt, A/C/G/T for DNA or A/C/G/U for RNA)
            target: Optional target DNA sequence for on-target analysis
            compute_ci: Whether to compute bootstrap confidence intervals
            n_bootstrap: Number of bootstrap resamples (default: 1000)
        
        Returns:
            Dictionary with scoring results:
                - disruption_score: Composite Z-invariant score
                - spectral_features: FFT-based features
                - on_target_metrics: Target-specific metrics (if target provided)
                - confidence_interval: Bootstrap CI (if compute_ci=True)
                - metadata: Analysis metadata
        
        Raises:
            ValueError: If sequences contain invalid nucleotides
        """
        start_time = time.time()
        
        # Validate guide sequence
        guide_upper = guide.upper()
        is_rna = 'U' in guide_upper
        validate_dna_sequence(guide, allow_rna=is_rna)
        
        # Initialize scorecard with correct is_rna flag
        scorecard = PhaseWeightedScorecard(k=self.k, is_rna=is_rna)
        
        # Score with phase-weighted scorecard
        scorecard_result = scorecard.score_guide(guide, target)
        
        # Extract guide features from scorecard result
        guide_features = scorecard_result.get('guide_features', {})
        
        # Get Z-score if target provided, otherwise use spectral entropy as proxy
        if target and scorecard_result.get('z_score') is not None:
            z_score = scorecard_result['z_score']
        else:
            # For single guide, use sigmoid of entropy as score
            entropy = guide_features.get('entropy', 0.0)
            diversity = guide_features.get('diversity', 0.5)
            kappa = diversity * np.log(len(guide) + 1) / E_SQUARED if len(guide) >= 10 else 1.0
            z_score = 1.0 / (1.0 + np.exp(-kappa * entropy / PHI))
        
        result = {
            'disruption_score': float(z_score),
            'spectral_features': {
                'entropy': float(guide_features.get('entropy', 0.0)),
                'dominant_freq': int(guide_features.get('dominant_freq_idx', 0)),
                'dominant_magnitude': float(guide_features.get('dominant_freq_mag', 0.0)),
                'sidelobe_count': int(guide_features.get('sidelobe_count', 0)),
                'diversity': float(guide_features.get('diversity', 0.0)),
            },
            'metadata': {
                'guide_length': len(guide),
                'guide_type': 'RNA' if is_rna else 'DNA',
                'k_parameter': self.k,
                'phi_period': self.phi_period,
                'computation_time_ms': (time.time() - start_time) * 1000,
            }
        }
        
        # Add target-specific metrics if target provided
        if target and 'disruptions' in scorecard_result:
            disruptions = scorecard_result['disruptions']
            result['on_target_metrics'] = {
                'delta_entropy': float(disruptions.get('delta_entropy', 0.0)),
                'delta_freq': float(disruptions.get('delta_freq', 0.0)),
                'delta_sidelobes': float(disruptions.get('delta_sidelobes', 0.0)),
                'delta_spectral': float(scorecard_result.get('delta_spectral', 0.0)),
            }
        
        # Compute bootstrap CI if requested
        if compute_ci:
            result['confidence_interval'] = self._compute_bootstrap_ci(
                guide, target, n_bootstrap
            )
        
        return result
    
    def _compute_bootstrap_ci(
        self,
        guide: str,
        target: Optional[str],
        n_bootstrap: int,
    ) -> Dict[str, float]:
        """
        Compute bootstrap confidence interval for disruption score.
        
        Uses position resampling with replacement to estimate CI.
        
        Args:
            guide: Guide sequence
            target: Target sequence (optional)
            n_bootstrap: Number of bootstrap resamples
        
        Returns:
            Dictionary with CI bounds and median
        """
        scores = []
        guide_array = np.array(list(guide))
        
        for _ in range(n_bootstrap):
            # Resample positions with replacement
            indices = np.random.choice(len(guide), size=len(guide), replace=True)
            resampled_guide = ''.join(guide_array[indices])
            
            try:
                result = self.score_guide(resampled_guide, target, compute_ci=False)
                scores.append(result['disruption_score'])
            except (ValueError, KeyError):
                # Skip invalid resamples
                continue
        
        scores = np.array(scores)
        return {
            'median': float(np.median(scores)),
            'lower_95': float(np.percentile(scores, 2.5)),
            'upper_95': float(np.percentile(scores, 97.5)),
            'n_samples': len(scores),
        }
    
    def batch_score(
        self,
        guides: List[str],
        targets: Optional[List[str]] = None,
        compute_ci: bool = False,
        n_bootstrap: int = 1000,
    ) -> List[Dict[str, Any]]:
        """
        Batch score multiple guides.
        
        Optimized for throughput, targeting 1,000 guides/min.
        
        Args:
            guides: List of guide sequences
            targets: Optional list of target sequences (same length as guides)
            compute_ci: Whether to compute bootstrap CIs
            n_bootstrap: Number of bootstrap resamples
        
        Returns:
            List of scoring results, one per guide
        
        Raises:
            ValueError: If targets length doesn't match guides
        """
        if targets is not None and len(targets) != len(guides):
            raise ValueError(
                f"Targets length ({len(targets)}) must match guides length ({len(guides)})"
            )
        
        results = []
        for i, guide in enumerate(guides):
            target = targets[i] if targets else None
            try:
                result = self.score_guide(guide, target, compute_ci, n_bootstrap)
                result['index'] = i
                results.append(result)
            except Exception as e:
                logger.error(f"Error scoring guide {i}: {e}")
                results.append({
                    'index': i,
                    'error': str(e),
                    'guide': guide,
                })
        
        return results
    
    def design_guides(
        self,
        target_sequence: str,
        pam: str = 'NGG',
        guide_length: int = 20,
        n_guides: int = 5,
        score_threshold: float = 0.5,
    ) -> List[Dict[str, Any]]:
        """
        Design optimal CRISPR guides from target sequence.
        
        Scans target for PAM sites and scores candidate guides.
        
        Args:
            target_sequence: Target DNA sequence
            pam: PAM sequence motif (default: 'NGG' for SpCas9)
            guide_length: Length of guide RNA (default: 20)
            n_guides: Maximum number of guides to return
            score_threshold: Minimum disruption score threshold
        
        Returns:
            List of designed guides with scores and positions
        """
        # Validate target
        validate_dna_sequence(target_sequence, allow_rna=False)
        target_upper = target_sequence.upper()
        
        candidates = []
        
        # Scan for PAM sites
        for i in range(len(target_upper) - guide_length - len(pam) + 1):
            pam_site = target_upper[i + guide_length:i + guide_length + len(pam)]
            
            # Check PAM match (N matches any base)
            if self._matches_pam(pam_site, pam):
                # Extract guide sequence
                guide = target_upper[i:i + guide_length]
                
                # Score guide without target (just spectral features)
                result = self.score_guide(guide, target=None)
                
                if result['disruption_score'] >= score_threshold:
                    candidates.append({
                        'guide': guide,
                        'position': i,
                        'pam_site': pam_site,
                        'disruption_score': result['disruption_score'],
                        'spectral_features': result['spectral_features'],
                    })
        
        # Sort by disruption score (descending) and return top n_guides
        candidates.sort(key=lambda x: x['disruption_score'], reverse=True)
        return candidates[:n_guides]
    
    def _matches_pam(self, sequence: str, pam: str) -> bool:
        """
        Check if sequence matches PAM motif (N matches any base).
        
        Args:
            sequence: Sequence to check
            pam: PAM motif pattern
        
        Returns:
            True if sequence matches PAM
        """
        if len(sequence) != len(pam):
            return False
        
        for s, p in zip(sequence, pam):
            if p != 'N' and s != p:
                return False
        
        return True
    
    def analyze_offtarget(
        self,
        guide: str,
        offtarget_sequences: List[str],
        mismatch_threshold: int = 3,
    ) -> List[Dict[str, Any]]:
        """
        Analyze off-target potential using spectral signature distance.
        
        Args:
            guide: Guide sequence
            offtarget_sequences: List of potential off-target sequences
            mismatch_threshold: Maximum mismatches to report
        
        Returns:
            List of off-target hits with distances and mismatch counts
        """
        results = []
        
        # Get guide spectral signature
        guide_result = self.score_guide(guide)
        guide_features = np.array([
            guide_result['spectral_features']['entropy'],
            guide_result['spectral_features']['dominant_freq'],
            guide_result['spectral_features']['diversity'],
        ])
        
        for offtarget in offtarget_sequences:
            try:
                # Score off-target
                offtarget_result = self.score_guide(offtarget)
                offtarget_features = np.array([
                    offtarget_result['spectral_features']['entropy'],
                    offtarget_result['spectral_features']['dominant_freq'],
                    offtarget_result['spectral_features']['diversity'],
                ])
                
                # Compute spectral distance
                spectral_distance = float(np.linalg.norm(guide_features - offtarget_features))
                
                # Compute Hamming distance (mismatches)
                mismatches = sum(g != o for g, o in zip(guide.upper(), offtarget.upper()))
                
                if mismatches <= mismatch_threshold:
                    results.append({
                        'sequence': offtarget,
                        'mismatches': mismatches,
                        'spectral_distance': spectral_distance,
                        'disruption_score': offtarget_result['disruption_score'],
                    })
            except Exception as e:
                logger.warning(f"Error analyzing off-target {offtarget}: {e}")
                continue
        
        # Sort by spectral distance
        results.sort(key=lambda x: x['spectral_distance'])
        return results


class GenomicDisruptionAPI:
    """
    REST API wrapper for DisruptionAnalyzer.
    
    Provides HTTP endpoints for CRISPR guide analysis services.
    Note: This is a minimal implementation. For production, use FastAPI or Flask.
    """
    
    def __init__(self, k: float = K_STAR, seed: Optional[int] = None):
        """
        Initialize API.
        
        Args:
            k: Resolution exponent for θ′(n,k)
            seed: Random seed for reproducibility
        """
        self.analyzer = DisruptionAnalyzer(k=k, seed=seed)
        logger.info(f"GenomicDisruptionAPI initialized with k={k}")
    
    def handle_score(self, request: Dict[str, Any]) -> Dict[str, Any]:
        """
        Handle single guide scoring request.
        
        Request format:
            {
                "guide": "GCTGCGGAGACCTGGAGAGA",
                "target": "ATGC..." (optional),
                "compute_ci": false (optional),
                "n_bootstrap": 1000 (optional)
            }
        
        Returns:
            Scoring result dictionary
        """
        try:
            guide = request.get('guide')
            if not guide:
                return {'error': 'Missing required field: guide'}
            
            result = self.analyzer.score_guide(
                guide=guide,
                target=request.get('target'),
                compute_ci=request.get('compute_ci', False),
                n_bootstrap=request.get('n_bootstrap', 1000),
            )
            return {'success': True, 'data': result}
        
        except Exception as e:
            logger.error(f"Error in handle_score: {e}", exc_info=True)
            return {'success': False, 'error': str(e)}
    
    def handle_batch(self, request: Dict[str, Any]) -> Dict[str, Any]:
        """
        Handle batch scoring request.
        
        Request format:
            {
                "guides": ["GCTGCG...", "ATCGAT..."],
                "targets": ["ATGC...", "CGTA..."] (optional),
                "compute_ci": false (optional),
                "n_bootstrap": 1000 (optional)
            }
        
        Returns:
            Batch scoring results
        """
        try:
            guides = request.get('guides')
            if not guides:
                return {'error': 'Missing required field: guides'}
            
            results = self.analyzer.batch_score(
                guides=guides,
                targets=request.get('targets'),
                compute_ci=request.get('compute_ci', False),
                n_bootstrap=request.get('n_bootstrap', 1000),
            )
            return {'success': True, 'count': len(results), 'data': results}
        
        except Exception as e:
            logger.error(f"Error in handle_batch: {e}", exc_info=True)
            return {'success': False, 'error': str(e)}
    
    def handle_design(self, request: Dict[str, Any]) -> Dict[str, Any]:
        """
        Handle guide design request.
        
        Request format:
            {
                "target": "ATGCGATCGATC...",
                "pam": "NGG" (optional),
                "guide_length": 20 (optional),
                "n_guides": 5 (optional),
                "score_threshold": 0.5 (optional)
            }
        
        Returns:
            Designed guides with scores
        """
        try:
            target = request.get('target')
            if not target:
                return {'error': 'Missing required field: target'}
            
            results = self.analyzer.design_guides(
                target_sequence=target,
                pam=request.get('pam', 'NGG'),
                guide_length=request.get('guide_length', 20),
                n_guides=request.get('n_guides', 5),
                score_threshold=request.get('score_threshold', 0.5),
            )
            return {'success': True, 'count': len(results), 'data': results}
        
        except Exception as e:
            logger.error(f"Error in handle_design: {e}", exc_info=True)
            return {'success': False, 'error': str(e)}
    
    def handle_offtarget(self, request: Dict[str, Any]) -> Dict[str, Any]:
        """
        Handle off-target analysis request.
        
        Request format:
            {
                "guide": "GCTGCGGAGACCTGGAGAGA",
                "offtargets": ["GCTGCGG...", "ACTGCGG..."],
                "mismatch_threshold": 3 (optional)
            }
        
        Returns:
            Off-target analysis results
        """
        try:
            guide = request.get('guide')
            offtargets = request.get('offtargets')
            
            if not guide:
                return {'error': 'Missing required field: guide'}
            if not offtargets:
                return {'error': 'Missing required field: offtargets'}
            
            results = self.analyzer.analyze_offtarget(
                guide=guide,
                offtarget_sequences=offtargets,
                mismatch_threshold=request.get('mismatch_threshold', 3),
            )
            return {'success': True, 'count': len(results), 'data': results}
        
        except Exception as e:
            logger.error(f"Error in handle_offtarget: {e}", exc_info=True)
            return {'success': False, 'error': str(e)}


# CLI interface for standalone usage
def main():
    """CLI interface for testing and standalone usage."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Genomic Disruption Analyzer CLI'
    )
    parser.add_argument(
        'command',
        choices=['score', 'batch', 'design', 'offtarget'],
        help='Command to execute'
    )
    parser.add_argument(
        '--guide',
        help='Guide RNA sequence'
    )
    parser.add_argument(
        '--target',
        help='Target DNA sequence'
    )
    parser.add_argument(
        '--guides',
        help='File containing guide sequences (one per line)'
    )
    parser.add_argument(
        '--k',
        type=float,
        default=K_STAR,
        help=f'Resolution exponent (default: {K_STAR})'
    )
    parser.add_argument(
        '--seed',
        type=int,
        help='Random seed for reproducibility'
    )
    parser.add_argument(
        '--output',
        help='Output JSON file'
    )
    
    args = parser.parse_args()
    
    # Initialize API
    api = GenomicDisruptionAPI(k=args.k, seed=args.seed)
    
    # Execute command
    if args.command == 'score':
        if not args.guide:
            print("Error: --guide required for score command")
            sys.exit(1)
        
        request = {'guide': args.guide}
        if args.target:
            request['target'] = args.target
        
        result = api.handle_score(request)
        print(json.dumps(result, indent=2))
        
        if args.output:
            with open(args.output, 'w') as f:
                json.dump(result, f, indent=2)
    
    elif args.command == 'design':
        if not args.target:
            print("Error: --target required for design command")
            sys.exit(1)
        
        request = {'target': args.target}
        result = api.handle_design(request)
        print(json.dumps(result, indent=2))
        
        if args.output:
            with open(args.output, 'w') as f:
                json.dump(result, f, indent=2)
    
    else:
        print(f"Command {args.command} not yet implemented in CLI")
        sys.exit(1)


if __name__ == '__main__':
    main()
