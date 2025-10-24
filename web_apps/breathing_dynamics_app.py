#!/usr/bin/env python3
"""
Interactive Web Demo: DNA Breathing Dynamics for CRISPR Prediction

Demonstrates biophysically-grounded DNA breathing dynamics encoding for 
CRISPR activity prediction. Uses CZT/Goertzel for precise fractional-period
analysis at the helical period (10.5 bp).

Features:
- Real-time breathing spectrum visualization
- Power at 10.5 bp helical period
- Comparison with baseline models
- Interactive parameter tuning (temperature, Mg²⁺)
"""

from flask import Flask, render_template, request, jsonify
import numpy as np
import json
import sys
import os
from typing import Dict, List, Any
import logging
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import io
import base64

# Add parent directory for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from experiments.signal_theoretic_crispr.breathing_dynamics import (
    BreathingDynamicsEncoder,
    BreathingSpectralAnalyzer,
    ChirpZTransform,
    HELICAL_PERIOD_BP,
    BP_OPENING_LIFETIMES_MS
)

from experiments.signal_theoretic_crispr.ablation_tests import (
    RandomEncoder,
    AblationTester
)

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Create Flask app
app = Flask(__name__)


class BreathingDemonstrator:
    """Main class for breathing dynamics demonstration."""
    
    def __init__(self):
        """Initialize demonstrator."""
        self.default_temp = 37.0
        self.default_mg = 2.0
        self.logger = logging.getLogger(__name__)
    
    def validate_sequence(self, sequence: str) -> tuple[bool, str]:
        """
        Validate DNA sequence according to scientific gates.
        
        Args:
            sequence: DNA sequence string
            
        Returns:
            (is_valid, message)
        """
        if not sequence:
            return False, "Sequence cannot be empty"
        
        sequence = sequence.upper().strip()
        
        # Check for valid bases only (A, C, G, T, N)
        valid_bases = set('ACGTN')
        invalid_bases = set(sequence) - valid_bases
        
        if invalid_bases:
            return False, f"Invalid bases: {invalid_bases}. Only A/C/G/T/N allowed (human DNA)"
        
        # Length checks
        if len(sequence) < 15:
            return False, "Sequence must be at least 15 bases long"
        
        if len(sequence) > 200:
            return False, "Sequence too long (max 200 bases for web demo)"
        
        return True, "Valid"
    
    def analyze_sequence(self,
                        sequence: str,
                        temperature_c: float = 37.0,
                        mg_mm: float = 2.0,
                        use_czt: bool = True) -> Dict[str, Any]:
        """
        Perform complete breathing dynamics analysis.
        
        Args:
            sequence: DNA sequence
            temperature_c: Temperature in Celsius
            mg_mm: Mg²⁺ concentration in mM
            use_czt: Use CZT (True) or Goertzel (False)
            
        Returns:
            Analysis results dictionary
        """
        # Validate
        is_valid, msg = self.validate_sequence(sequence)
        if not is_valid:
            return {'error': msg}
        
        try:
            # Create analyzer
            analyzer = BreathingSpectralAnalyzer(
                temperature_c=temperature_c,
                mg_concentration_mm=mg_mm,
                use_czt=use_czt,
                seed=42
            )
            
            # Extract breathing features
            features = analyzer.extract_breathing_features(sequence, harmonics=3)
            
            # Get encoder info
            encoder_info = analyzer.encoder.get_encoding_info()
            
            # Create visualizations
            fig_spectrum = self.plot_spectrum(analyzer, sequence)
            fig_weights = self.plot_base_weights(analyzer.encoder)
            
            results = {
                'sequence': sequence,
                'length': len(sequence),
                'gc_content': features.get('breathing_gc_content', 0.0),
                'features': {
                    'czt_10.5bp_power': features.get('czt_period_10.5_total_power', 0.0),
                    'phase_coherence': features.get('breathing_phase_coherence', 0.0),
                    'amplitude_variance': features.get('breathing_amplitude_var', 0.0),
                    'fft_peak_power': features.get('fft_peak_power', 0.0),
                },
                'harmonics': {
                    'h1_power': features.get('czt_period_10.5_h1_power', 0.0),
                    'h2_power': features.get('czt_period_10.5_h2_power', 0.0),
                    'h3_power': features.get('czt_period_10.5_h3_power', 0.0),
                },
                'parameters': {
                    'temperature_c': temperature_c,
                    'mg_mm': mg_mm,
                    'helical_period_bp': HELICAL_PERIOD_BP,
                    'method': 'CZT' if use_czt else 'Goertzel'
                },
                'encoder_info': encoder_info,
                'plots': {
                    'spectrum': fig_spectrum,
                    'weights': fig_weights
                }
            }
            
            return results
            
        except Exception as e:
            # Log error but don't expose stack trace details
            self.logger.error(f"Analysis failed for sequence length {len(sequence)}")
            return {'error': 'Analysis failed. Please check your sequence and try again.'}
    
    def plot_spectrum(self, analyzer: BreathingSpectralAnalyzer, sequence: str) -> str:
        """
        Plot breathing spectrum.
        
        Returns:
            Base64-encoded PNG image
        """
        # Encode and analyze
        encoded = analyzer.encoder.encode_sequence(sequence, apply_helical_phase=True)
        encoded_dc = analyzer.remove_dc(encoded)
        encoded_windowed = analyzer.apply_window(encoded_dc, 'hamming')
        
        # Compute FFT
        from scipy.fft import fft, fftfreq
        spectrum = np.abs(fft(encoded_windowed))
        freqs = fftfreq(len(encoded_windowed), d=1.0)
        
        # Plot
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
        
        # Power spectrum
        ax1.plot(freqs[:len(freqs)//2], spectrum[:len(spectrum)//2], 'b-', linewidth=1)
        ax1.axvline(1/HELICAL_PERIOD_BP, color='r', linestyle='--', 
                   label=f'Helical period (1/{HELICAL_PERIOD_BP:.1f} bp⁻¹)')
        ax1.set_xlabel('Frequency (cycles/bp)')
        ax1.set_ylabel('Power')
        ax1.set_title('Breathing Dynamics Power Spectrum')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Waveform
        positions = np.arange(len(encoded_windowed))
        ax2.plot(positions, np.abs(encoded_windowed), 'g-', linewidth=1, label='Magnitude')
        ax2.plot(positions, np.angle(encoded_windowed), 'r-', linewidth=1, alpha=0.5, label='Phase')
        ax2.set_xlabel('Position (bp)')
        ax2.set_ylabel('Value')
        ax2.set_title('Encoded Waveform')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        # Convert to base64
        buf = io.BytesIO()
        plt.savefig(buf, format='png', dpi=100, bbox_inches='tight')
        buf.seek(0)
        img_base64 = base64.b64encode(buf.read()).decode('utf-8')
        plt.close(fig)
        
        return img_base64
    
    def plot_base_weights(self, encoder: BreathingDynamicsEncoder) -> str:
        """
        Plot base weights.
        
        Returns:
            Base64-encoded PNG image
        """
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
        
        bases = ['A', 'T', 'C', 'G']
        weights = [encoder.base_weights[b] for b in bases]
        
        # Real and imaginary components
        real_parts = [w.real for w in weights]
        imag_parts = [w.imag for w in weights]
        
        x = np.arange(len(bases))
        width = 0.35
        
        ax1.bar(x - width/2, real_parts, width, label='Real (Opening Rate)', color='skyblue')
        ax1.bar(x + width/2, imag_parts, width, label='Imag (Stability)', color='orange')
        ax1.set_xlabel('Base')
        ax1.set_ylabel('Weight Component')
        ax1.set_title('Biophysical Base Weights')
        ax1.set_xticks(x)
        ax1.set_xticklabels(bases)
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Complex plane
        for base, weight in zip(bases, weights):
            ax2.plot(weight.real, weight.imag, 'o', markersize=12, label=base)
            ax2.annotate(base, (weight.real, weight.imag), 
                        xytext=(5, 5), textcoords='offset points', fontsize=12)
        
        ax2.set_xlabel('Real (Opening Rate)')
        ax2.set_ylabel('Imaginary (Stability)')
        ax2.set_title('Complex Weight Encoding')
        ax2.axhline(0, color='k', linewidth=0.5, alpha=0.3)
        ax2.axvline(0, color='k', linewidth=0.5, alpha=0.3)
        ax2.grid(True, alpha=0.3)
        ax2.legend()
        
        plt.tight_layout()
        
        # Convert to base64
        buf = io.BytesIO()
        plt.savefig(buf, format='png', dpi=100, bbox_inches='tight')
        buf.seek(0)
        img_base64 = base64.b64encode(buf.read()).decode('utf-8')
        plt.close(fig)
        
        return img_base64


# Create global demonstrator
demonstrator = BreathingDemonstrator()


@app.route('/')
def index():
    """Render main page."""
    return render_template('breathing_demo.html')


@app.route('/analyze', methods=['POST'])
def analyze():
    """Analyze DNA sequence."""
    try:
        data = request.get_json()
        
        sequence = data.get('sequence', '').upper().strip()
        temperature = float(data.get('temperature', 37.0))
        mg_concentration = float(data.get('mg_concentration', 2.0))
        use_czt = data.get('use_czt', True)
        
        # Validate parameters
        if not (0 <= temperature <= 100):
            return jsonify({'error': 'Temperature must be between 0-100°C'}), 400
        
        if not (0.1 <= mg_concentration <= 100):
            return jsonify({'error': 'Mg²⁺ must be between 0.1-100 mM'}), 400
        
        # Analyze
        results = demonstrator.analyze_sequence(
            sequence,
            temperature_c=temperature,
            mg_mm=mg_concentration,
            use_czt=use_czt
        )
        
        if 'error' in results:
            return jsonify({'error': results['error']}), 400
        
        return jsonify(results)
        
    except Exception as e:
        # Log error but don't expose stack trace
        logger.error(f"Analysis request failed")
        return jsonify({'error': 'Analysis failed. Please try again.'}), 500


@app.route('/example/<example_type>')
def get_example(example_type):
    """Get example sequence."""
    examples = {
        'at_rich': 'AAATTTAAATTTAAATTTAT',
        'gc_rich': 'GGGCCCGGGCCCGGGCCCGC',
        'mixed': 'ATCGATCGATCGATCGATCG',
        'guide': 'GACGATCGATCGATCGATCGATCG'  # 24 bp typical guide
    }
    
    sequence = examples.get(example_type, examples['mixed'])
    return jsonify({'sequence': sequence})


@app.route('/info')
def get_info():
    """Get breathing dynamics information."""
    info = {
        'title': 'DNA Breathing Dynamics for CRISPR Prediction',
        'description': 'Biophysically-grounded encoding based on base-pair opening rates',
        'references': {
            'opening_lifetimes': 'PMC5393899',
            'helical_period': '~10.5 bp for B-DNA',
            'r_loop_formation': 'PNAS 1402597111'
        },
        'parameters': {
            'at_opening_ms': BP_OPENING_LIFETIMES_MS['AT'],
            'gc_opening_ms': BP_OPENING_LIFETIMES_MS['GC'],
            'helical_period_bp': HELICAL_PERIOD_BP
        },
        'features': [
            'CZT/Goertzel for fractional-period analysis',
            'Temperature and Mg²⁺ dependencies',
            'DC removal and windowing',
            'Phase-aware feature extraction'
        ]
    }
    
    return jsonify(info)


if __name__ == '__main__':
    # Check if running in production
    is_production = os.environ.get('FLASK_ENV') == 'production'
    
    if is_production:
        # Production settings
        app.run(host='0.0.0.0', port=int(os.environ.get('PORT', 5000)), debug=False)
    else:
        # Development settings - debug disabled for security
        print("\n" + "="*70)
        print("DNA BREATHING DYNAMICS DEMO")
        print("="*70)
        print("\nStarting web server on http://localhost:5000")
        print("\nExample sequences:")
        print("  AT-rich: AAATTTAAATTTAAATTTAT")
        print("  GC-rich: GGGCCCGGGCCCGGGCCCGC")
        print("  Mixed:   ATCGATCGATCGATCGATCG")
        print("\nPress Ctrl+C to stop")
        print("="*70 + "\n")
        
        # Run without debug mode for security
        app.run(host='0.0.0.0', port=5000, debug=False)
