"""
CRISPR Guide Designer using Signal-Theoretic DNA Analysis

This module implements CRISPR guide RNA design using spectral disruption profiling
based on the wave-crispr-signal methodology.
"""

import numpy as np
from scipy.fft import fft
from scipy.stats import entropy
from collections import Counter
import re
from typing import List, Dict, Tuple, Optional, Union


class CRISPRGuideDesigner:
    """CRISPR guide RNA designer using signal-theoretic analysis."""
    
    # Base weights for wave function mapping
    WEIGHTS = {'A': 1 + 0j, 'T': -1 + 0j, 'C': 0 + 1j, 'G': 0 - 1j}
    
    def __init__(self, pam_pattern: str = "NGG", guide_length: int = 20):
        """
        Initialize CRISPR guide designer.
        
        Args:
            pam_pattern: PAM sequence pattern (default: NGG for Cas9)
            guide_length: Length of guide RNA (default: 20)
        """
        self.pam_pattern = pam_pattern.replace('N', '[ATCG]')
        self.guide_length = guide_length
        
    def build_waveform(self, seq: str, d: float = 0.34, zn_map: Optional[Dict[int, float]] = None) -> np.ndarray:
        """
        Build complex waveform from DNA sequence.
        
        Args:
            seq: DNA sequence
            d: Base spacing parameter (default: 0.34 nm)
            zn_map: Position-dependent scaling factors
            
        Returns:
            Complex waveform array
        """
        if zn_map is None:
            spacings = [d] * len(seq)
        else:
            spacings = [d * (1 + zn_map.get(i, 0)) for i in range(len(seq))]
        
        s = np.cumsum(spacings)
        wave = [self.WEIGHTS[base] * np.exp(2j * np.pi * s[i]) for i, base in enumerate(seq)]
        return np.array(wave)
    
    def compute_spectrum(self, waveform: np.ndarray) -> np.ndarray:
        """Compute frequency spectrum of waveform."""
        return np.abs(fft(waveform))
    
    def normalized_entropy(self, spectrum: np.ndarray) -> float:
        """Calculate normalized entropy of spectrum."""
        ps = spectrum / np.sum(spectrum)
        return entropy(ps, base=2)
    
    def count_sidelobes(self, spectrum: np.ndarray, threshold_ratio: float = 0.25) -> int:
        """Count spectral sidelobes above threshold."""
        peak = np.max(spectrum)
        return np.sum(spectrum > (threshold_ratio * peak))
    
    def calculate_on_target_score(self, guide_seq: str) -> float:
        """
        Calculate on-target efficiency score using spectral features.
        
        Args:
            guide_seq: Guide RNA sequence
            
        Returns:
            On-target efficiency score (0-1)
        """
        wave = self.build_waveform(guide_seq)
        spec = self.compute_spectrum(wave)
        
        # Spectral features for on-target prediction
        entropy_score = 1.0 - (self.normalized_entropy(spec) / 10.0)  # Lower entropy = higher efficiency
        sidelobe_score = 1.0 - (self.count_sidelobes(spec) / len(spec))  # Fewer sidelobes = higher efficiency
        
        # GC content factor (optimal ~40-60%)
        gc_content = (guide_seq.count('G') + guide_seq.count('C')) / len(guide_seq)
        gc_score = 1.0 - abs(gc_content - 0.5) * 2  # Penalize extreme GC content
        
        # Combined score
        return np.clip((entropy_score * 0.4 + sidelobe_score * 0.4 + gc_score * 0.2), 0, 1)
    
    def calculate_off_target_risk(self, guide_seq: str, target_seq: str, off_target_seq: str) -> float:
        """
        Calculate off-target risk using spectral signature comparison.
        
        Args:
            guide_seq: Guide RNA sequence  
            target_seq: On-target genomic sequence
            off_target_seq: Potential off-target sequence
            
        Returns:
            Off-target risk score (0-1)
        """
        target_wave = self.build_waveform(target_seq)
        off_target_wave = self.build_waveform(off_target_seq)
        
        target_spec = self.compute_spectrum(target_wave)
        off_target_spec = self.compute_spectrum(off_target_wave)
        
        # Spectral similarity using cross-correlation
        min_len = min(len(target_spec), len(off_target_spec))
        corr = np.corrcoef(target_spec[:min_len], off_target_spec[:min_len])[0, 1]
        
        # Sequence similarity (Hamming distance)
        min_seq_len = min(len(target_seq), len(off_target_seq))
        seq_similarity = sum(a == b for a, b in zip(target_seq[:min_seq_len], off_target_seq[:min_seq_len])) / min_seq_len
        
        # Combined risk score (high correlation + high sequence similarity = high risk)
        spectral_risk = max(0, corr) if not np.isnan(corr) else 0
        return np.clip((spectral_risk * 0.6 + seq_similarity * 0.4), 0, 1)
    
    def find_pam_sites(self, sequence: str) -> List[Tuple[int, str]]:
        """
        Find PAM sites in sequence.
        
        Args:
            sequence: DNA sequence to search
            
        Returns:
            List of (position, PAM_sequence) tuples
        """
        pam_sites = []
        pattern = re.compile(self.pam_pattern)
        
        for match in pattern.finditer(sequence.upper()):
            pam_sites.append((match.start(), match.group()))
            
        return pam_sites
    
    def design_guides(self, sequence: str, num_guides: int = 5) -> List[Dict]:
        """
        Design CRISPR guides for a target sequence.
        
        Args:
            sequence: Target DNA sequence
            num_guides: Number of top guides to return
            
        Returns:
            List of guide dictionaries with scores and positions
        """
        guides = []
        pam_sites = self.find_pam_sites(sequence)
        
        for pam_pos, pam_seq in pam_sites:
            # Extract guide sequence (upstream of PAM)
            guide_start = max(0, pam_pos - self.guide_length)
            guide_end = pam_pos
            
            if guide_end - guide_start >= self.guide_length:
                guide_seq = sequence[guide_start:guide_end].upper()
                
                # Calculate scores
                on_target_score = self.calculate_on_target_score(guide_seq)
                
                # Basic quality filters
                gc_content = (guide_seq.count('G') + guide_seq.count('C')) / len(guide_seq)
                has_poly_t = 'TTTT' in guide_seq  # Avoid poly-T stretches
                
                if 0.2 <= gc_content <= 0.8 and not has_poly_t:
                    guides.append({
                        'sequence': guide_seq,
                        'position': guide_start,
                        'pam_position': pam_pos,
                        'pam_sequence': pam_seq,
                        'on_target_score': on_target_score,
                        'gc_content': gc_content,
                        'length': len(guide_seq)
                    })
        
        # Sort by on-target score and return top guides
        guides.sort(key=lambda x: x['on_target_score'], reverse=True)
        return guides[:num_guides]
    
    def predict_repair_outcomes(self, guide_seq: str, target_seq: str) -> Dict:
        """
        Predict CRISPR repair outcomes using spectral stability analysis.
        
        Args:
            guide_seq: Guide RNA sequence
            target_seq: Target genomic sequence around cut site
            
        Returns:
            Dictionary with repair outcome predictions
        """
        # Find cut site (typically 3 bp upstream of PAM)
        cut_site = len(target_seq) // 2  # Assume cut site is in middle
        
        # Analyze sequence context around cut site
        context_seq = target_seq[max(0, cut_site-10):cut_site+10]
        wave = self.build_waveform(context_seq)
        spec = self.compute_spectrum(wave)
        
        entropy_val = self.normalized_entropy(spec)
        sidelobes = self.count_sidelobes(spec)
        
        # Predict repair pathway bias
        # Low entropy → NHEJ bias
        # High sidelobes → MMEJ bias  
        nhej_bias = 1.0 - (entropy_val / 10.0)
        mmej_bias = sidelobes / len(spec)
        hdr_efficiency = max(0, 1.0 - nhej_bias - mmej_bias)
        
        return {
            'nhej_probability': np.clip(nhej_bias, 0, 1),
            'mmej_probability': np.clip(mmej_bias, 0, 1),
            'hdr_efficiency': np.clip(hdr_efficiency, 0, 1),
            'spectral_entropy': entropy_val,
            'stability_score': 1.0 - entropy_val / 10.0
        }


def main():
    """Example usage of CRISPR guide designer."""
    # Example target sequence (PCSK9 exon 1)
    target_sequence = "ATGCTGCGGAGACCTGGAGAGAAAGCAGTGGCCGGGGCAGTGGGAGGAGGAGGAGCTGGAAGAGGAGAGAAAGGAGGAGCTGCAGGAGGAGAGGAGGAGGAGGGAGAGGAGGAGCTGGAGCTGAAGCTGGAGCTGGAGCTGGAGAGGAGAGAGGG"
    
    designer = CRISPRGuideDesigner()
    
    print("🧬 CRISPR Guide Designer - Signal-Theoretic Analysis")
    print("=" * 60)
    
    # Design guides
    guides = designer.design_guides(target_sequence, num_guides=5)
    
    print(f"Target sequence length: {len(target_sequence)} bp")
    print(f"Found {len(guides)} high-quality guide candidates:\n")
    
    for i, guide in enumerate(guides, 1):
        print(f"Guide #{i}:")
        print(f"  Sequence: {guide['sequence']}")
        print(f"  Position: {guide['position']}-{guide['position'] + guide['length']}")
        print(f"  PAM: {guide['pam_sequence']} (pos {guide['pam_position']})")
        print(f"  On-target score: {guide['on_target_score']:.3f}")
        print(f"  GC content: {guide['gc_content']:.1%}")
        
        # Predict repair outcomes
        context_start = max(0, guide['position'] - 20)
        context_end = min(len(target_sequence), guide['position'] + guide['length'] + 20)
        context_seq = target_sequence[context_start:context_end]
        
        repair = designer.predict_repair_outcomes(guide['sequence'], context_seq)
        print(f"  Repair outcomes - NHEJ: {repair['nhej_probability']:.2f}, MMEJ: {repair['mmej_probability']:.2f}, HDR: {repair['hdr_efficiency']:.2f}")
        print()


if __name__ == "__main__":
    main()