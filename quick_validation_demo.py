#!/usr/bin/env python3
"""
Quick Experimental Validation Demo: Biological vs Arbitrary DNA Encoding

This demonstrates the experimental framework for testing biological relevance
of spectral DNA encodings with reduced computational requirements for demo purposes.
"""

import numpy as np
import scipy.stats as stats
from scipy.fft import fft
from scipy.stats import entropy, pearsonr
import matplotlib.pyplot as plt
import pandas as pd
from collections import defaultdict
import random
from typing import Dict, List, Tuple, Optional

# --- Physicochemical Properties for Biologically Anchored Encoding ---
NUCLEOTIDE_PROPERTIES = {
    'A': {'polarizability': 1.65, 'hydrogen_bonds': 2, 'molecular_weight': 331.2, 'purine': 1},
    'T': {'polarizability': 1.52, 'hydrogen_bonds': 2, 'molecular_weight': 322.2, 'purine': 0},
    'C': {'polarizability': 1.44, 'hydrogen_bonds': 3, 'molecular_weight': 307.2, 'purine': 0},
    'G': {'polarizability': 1.72, 'hydrogen_bonds': 3, 'molecular_weight': 347.2, 'purine': 1}
}

class BiologicalAnchoredEncoder:
    """Encoder based on nucleotide physicochemical properties"""
    
    def __init__(self):
        self.encodings = {}
        for nucleotide, props in NUCLEOTIDE_PROPERTIES.items():
            # Real: normalized polarizability
            real_part = props['polarizability'] / 1.72
            # Imaginary: hydrogen bonding capacity  
            imaginary_part = props['hydrogen_bonds'] / 3.0
            self.encodings[nucleotide] = complex(real_part, imaginary_part)
    
    def encode_sequence(self, sequence: str) -> np.ndarray:
        return np.array([self.encodings[base] for base in sequence.upper()])

class ArbitraryEncoder:
    """Encoder with arbitrary/random mappings"""
    
    def __init__(self, seed: int = 42):
        random.seed(seed)
        np.random.seed(seed)
        
        self.encodings = {}
        nucleotides = ['A', 'T', 'C', 'G']
        for nuc in nucleotides:
            real_part = random.uniform(-1, 1)
            imaginary_part = random.uniform(-1, 1)
            self.encodings[nuc] = complex(real_part, imaginary_part)
    
    def encode_sequence(self, sequence: str) -> np.ndarray:
        return np.array([self.encodings[base] for base in sequence.upper()])

class QuickValidator:
    """Streamlined experimental validation"""
    
    def __init__(self):
        self.bio_encoder = BiologicalAnchoredEncoder()
        self.arb_encoder = ArbitraryEncoder()
    
    def generate_test_data(self, n_sequences: int = 200) -> Tuple[List[str], List[float]]:
        """Generate test CRISPR sequences and efficiency scores"""
        sequences = []
        efficiencies = []
        
        nucleotides = ['A', 'T', 'C', 'G']
        
        for _ in range(n_sequences):
            # Generate 20bp guide sequence
            guide = ''.join(random.choices(nucleotides, k=20))
            sequences.append(guide)
            
            # Simulate efficiency based on GC content and position effects
            gc_content = (guide.count('G') + guide.count('C')) / len(guide)
            efficiency = 0.5  # baseline
            
            # GC content effect
            if 0.4 <= gc_content <= 0.7:
                efficiency += 0.2
            else:
                efficiency -= 0.1
            
            # Add noise
            efficiency += random.gauss(0, 0.118)
            efficiency = max(0, min(1, efficiency))
            
            efficiencies.append(efficiency)
        
        return sequences, efficiencies
    
    def extract_spectral_features(self, sequence: str, encoder) -> Dict[str, float]:
        """Extract key spectral features"""
        try:
            encoded = encoder.encode_sequence(sequence)
            
            # Build waveform with position modulation
            positions = np.cumsum([0.34] * len(encoded))
            waveform = encoded * np.exp(2j * np.pi * positions)
            
            # Compute spectrum
            spectrum = np.abs(fft(waveform))
            
            return {
                'spectral_entropy': entropy(spectrum / np.sum(spectrum), base=2),
                'peak_magnitude': np.max(spectrum),
                'spectral_variance': np.var(spectrum),
                'spectral_centroid': np.sum(spectrum * np.arange(len(spectrum))) / np.sum(spectrum)
            }
        except:
            return {'spectral_entropy': 0, 'peak_magnitude': 0, 'spectral_variance': 0, 'spectral_centroid': 0}
    
    def run_validation(self) -> Dict:
        """Run quick validation experiment"""
        print("🧬 Quick Experimental Validation")
        print("=" * 50)
        
        # Generate test data
        print("📊 Generating test data...")
        sequences, efficiencies = self.generate_test_data(200)
        print(f"✅ Generated {len(sequences)} sequences")
        
        # Extract features for both encoders
        print("🔍 Extracting features...")
        bio_features = []
        arb_features = []
        
        for seq in sequences:
            bio_features.append(self.extract_spectral_features(seq, self.bio_encoder))
            arb_features.append(self.extract_spectral_features(seq, self.arb_encoder))
        
        bio_df = pd.DataFrame(bio_features)
        arb_df = pd.DataFrame(arb_features)
        
        # Test correlations
        print("📊 Testing correlations...")
        bio_correlations = {}
        arb_correlations = {}
        
        for feature in bio_df.columns:
            # Biological encoder correlations
            try:
                bio_r, bio_p = pearsonr(bio_df[feature], efficiencies)
                bio_correlations[feature] = (bio_r, bio_p)
            except:
                bio_correlations[feature] = (0.0, 1.0)
            
            # Arbitrary encoder correlations
            try:
                arb_r, arb_p = pearsonr(arb_df[feature], efficiencies)
                arb_correlations[feature] = (arb_r, arb_p)
            except:
                arb_correlations[feature] = (0.0, 1.0)
        
        # Analyze results
        bio_significant = sum(1 for r, p in bio_correlations.values() if abs(r) >= 0.3 and p < 0.05)
        arb_significant = sum(1 for r, p in arb_correlations.values() if abs(r) >= 0.3 and p < 0.05)
        
        results = {
            'biological_correlations': bio_correlations,
            'arbitrary_correlations': arb_correlations,
            'biological_significant': bio_significant,
            'arbitrary_significant': arb_significant,
            'hypothesis_supported': bio_significant > arb_significant
        }
        
        return results
    
    def print_results(self, results: Dict):
        """Print validation results"""
        print("\n" + "="*60)
        print("🧬 VALIDATION RESULTS")
        print("="*60)
        
        bio_corr = results['biological_correlations']
        arb_corr = results['arbitrary_correlations']
        
        print(f"\n🔬 BIOLOGICALLY ANCHORED ENCODER")
        print(f"{'Feature':<20} {'Correlation':<12} {'P-value':<10} {'Significant':<12}")
        print("-" * 60)
        
        for feature, (r, p) in bio_corr.items():
            sig = "✓" if abs(r) >= 0.3 and p < 0.05 else "✗"
            print(f"{feature:<20} {r:>11.4f} {p:>9.4f} {sig:>11}")
        
        print(f"\n🎲 ARBITRARY ENCODER")
        print(f"{'Feature':<20} {'Correlation':<12} {'P-value':<10} {'Significant':<12}")
        print("-" * 60)
        
        for feature, (r, p) in arb_corr.items():
            sig = "✓" if abs(r) >= 0.3 and p < 0.05 else "✗"
            print(f"{feature:<20} {r:>11.4f} {p:>9.4f} {sig:>11}")
        
        print(f"\n🎯 SUMMARY")
        print(f"Biological encoder significant correlations: {results['biological_significant']}")
        print(f"Arbitrary encoder significant correlations: {results['arbitrary_significant']}")
        
        if results['hypothesis_supported']:
            print("✅ HYPOTHESIS SUPPORTED: Biological encoding shows superior performance")
        else:
            print("❌ HYPOTHESIS NOT SUPPORTED: No clear advantage for biological encoding")
        
        # Create simple visualization
        self.create_visualization(results)
    
    def create_visualization(self, results: Dict):
        """Create results visualization"""
        bio_corr = results['biological_correlations']
        arb_corr = results['arbitrary_correlations']
        
        features = list(bio_corr.keys())
        bio_r_values = [bio_corr[f][0] for f in features]
        arb_r_values = [arb_corr[f][0] for f in features]
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # Correlation comparison
        x_pos = np.arange(len(features))
        width = 0.35
        
        ax1.bar(x_pos - width/2, bio_r_values, width, label='Biological', alpha=0.8, color='#2E8B57')
        ax1.bar(x_pos + width/2, arb_r_values, width, label='Arbitrary', alpha=0.8, color='#CD5C5C')
        ax1.axhline(y=0.3, color='red', linestyle='--', alpha=0.7, label='Significance threshold')
        ax1.axhline(y=-0.3, color='red', linestyle='--', alpha=0.7)
        ax1.set_ylabel('Pearson Correlation (r)')
        ax1.set_title('Correlation with CRISPR Efficiency')
        ax1.set_xticks(x_pos)
        ax1.set_xticklabels([f.replace('_', '\n') for f in features], rotation=45, ha='right')
        ax1.legend()
        ax1.grid(axis='y', alpha=0.3)
        
        # Summary comparison
        categories = ['Biological', 'Arbitrary']
        significant_counts = [results['biological_significant'], results['arbitrary_significant']]
        
        bars = ax2.bar(categories, significant_counts, color=['#2E8B57', '#CD5C5C'], alpha=0.8)
        ax2.set_ylabel('Significant Correlations')
        ax2.set_title('Performance Comparison')
        ax2.grid(axis='y', alpha=0.3)
        
        for bar, count in zip(bars, significant_counts):
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height + 0.05,
                    f'{count}', ha='center', va='bottom', fontweight='bold')
        
        plt.tight_layout()
        plt.savefig('/home/runner/work/wave-crispr-signal/wave-crispr-signal/quick_validation_results.png', 
                    dpi=150, bbox_inches='tight')
        print(f"\n📊 Visualization saved as 'quick_validation_results.png'")

def main():
    """Run quick validation demo"""
    validator = QuickValidator()
    results = validator.run_validation()
    validator.print_results(results)

if __name__ == "__main__":
    main()