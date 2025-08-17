"""
Static 3D Plot Generator for Wave-CRISPR-Signal Framework

This script generates static PNG versions of the 3D plots using matplotlib
for easier viewing and inclusion in documentation.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
from scipy.interpolate import griddata
from scipy.fft import fft
import warnings
warnings.filterwarnings('ignore')

try:
    from applications.crispr_guide_designer import CRISPRGuideDesigner
    from invariant_features import InvariantFeatureSet, ZetaUnfoldCalculator
except ImportError:
    import sys
    sys.path.append('.')
    from applications.crispr_guide_designer import CRISPRGuideDesigner
    from invariant_features import InvariantFeatureSet, ZetaUnfoldCalculator


class Static3DPlotGenerator:
    """Generate static PNG versions of 3D plots."""
    
    def __init__(self):
        """Initialize the static plot generator."""
        self.designer = CRISPRGuideDesigner()
        self.invariant_set = InvariantFeatureSet()
        
        # Set plotting style
        plt.style.use('seaborn-v0_8')
        sns.set_palette("viridis")
    
    def create_spectral_landscape_static(self, sequence: str, save_path: str = "3d_spectral_landscape.png"):
        """Create static 3D spectral landscape plot."""
        fig = plt.figure(figsize=(12, 9))
        ax = fig.add_subplot(111, projection='3d')
        
        # Generate spectral data
        window_size = 20
        positions = []
        frequencies = []
        magnitudes = []
        
        for i in range(0, len(sequence) - window_size + 1, 3):
            subseq = sequence[i:i + window_size]
            wave = self.designer.build_waveform(subseq)
            spectrum = np.abs(fft(wave))
            
            for f_idx, magnitude in enumerate(spectrum[:10]):  # First 10 frequencies
                positions.append(i + window_size//2)
                frequencies.append(f_idx)
                magnitudes.append(magnitude)
        
        # Create surface
        pos_unique = sorted(list(set(positions)))
        freq_unique = sorted(list(set(frequencies)))
        
        if len(pos_unique) > 1 and len(freq_unique) > 1:
            POS, FREQ = np.meshgrid(pos_unique, freq_unique)
            
            # Interpolate magnitudes
            points = np.column_stack((positions, frequencies))
            MAGNITUDES = griddata(points, magnitudes, (POS, FREQ), method='nearest', fill_value=0)
            
            surf = ax.plot_surface(POS, FREQ, MAGNITUDES, cmap='viridis', alpha=0.8)
            
            ax.set_xlabel('Sequence Position (bp)')
            ax.set_ylabel('Frequency Index')
            ax.set_zlabel('Spectral Magnitude')
            ax.set_title('3D Spectral Landscape Analysis')
            
            plt.colorbar(surf, shrink=0.5, aspect=5)
        
        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Static spectral landscape saved to {save_path}")
    
    def create_invariant_features_static(self, sequences: list, save_path: str = "3d_invariant_features.png"):
        """Create static 3D invariant features scatter plot."""
        fig = plt.figure(figsize=(12, 9))
        ax = fig.add_subplot(111, projection='3d')
        
        golden_proximities = []
        phase_bits = []
        delta_entropies = []
        
        for seq in sequences:
            features = self.invariant_set.calculate_complete_feature_set(seq)
            golden_proximities.append(features.get('delta_phi', 0))
            phase_bits.append(features.get('phase_bit', 0))
            delta_entropies.append(features.get('delta_phase_entropy', 0))
        
        # Color by phase bit
        colors = ['red' if pb == 0 else 'blue' for pb in phase_bits]
        
        scatter = ax.scatter(golden_proximities, phase_bits, delta_entropies, 
                           c=colors, s=100, alpha=0.7, edgecolors='black')
        
        ax.set_xlabel('Golden Proximity (δφ)')
        ax.set_ylabel('Phase Bit (π)')
        ax.set_zlabel('Phase Entropy Difference')
        ax.set_title('3D Invariant Feature Space')
        
        # Add legend
        red_patch = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10, label='Phase 0')
        blue_patch = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=10, label='Phase 1')
        ax.legend(handles=[red_patch, blue_patch])
        
        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Static invariant features plot saved to {save_path}")
    
    def create_golden_proximity_surface_static(self, save_path: str = "3d_golden_proximity.png"):
        """Create static 3D golden ratio proximity surface."""
        fig = plt.figure(figsize=(12, 9))
        ax = fig.add_subplot(111, projection='3d')
        
        # Generate parameter space
        a_range = np.linspace(0.3, 0.8, 25)
        b_range = np.linspace(8.0, 12.0, 25)
        c_fixed = 7.389
        
        A, B = np.meshgrid(a_range, b_range)
        Z_values = A * (B / c_fixed)
        
        # Golden ratio conjugate
        phi_conjugate = 0.618033988749895
        delta_phi_surface = np.abs(Z_values - phi_conjugate)
        
        surf = ax.plot_surface(A, B, delta_phi_surface, cmap='RdYlBu_r', alpha=0.8)
        
        # Add golden ratio plane
        ax.plot_surface(A, B, np.zeros_like(A), alpha=0.3, color='gold')
        
        ax.set_xlabel('Parameter A')
        ax.set_ylabel('Parameter B')
        ax.set_zlabel('Distance to φ⁻¹ (δφ)')
        ax.set_title('3D Golden Ratio Proximity Surface')
        
        plt.colorbar(surf, shrink=0.5, aspect=5)
        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Static golden proximity surface saved to {save_path}")
    
    def create_guide_landscape_static(self, sequence: str, save_path: str = "3d_guide_landscape.png"):
        """Create static 3D guide efficiency landscape."""
        fig = plt.figure(figsize=(12, 9))
        ax = fig.add_subplot(111, projection='3d')
        
        # Generate guides
        guides = self.designer.design_guides(sequence, num_guides=30, use_invariants=True)
        
        positions = [g['position'] for g in guides]
        gc_contents = [g['gc_content'] for g in guides]
        scores = [g['comprehensive_score'] if 'comprehensive_score' in g 
                 else g['on_target_score'] for g in guides]
        
        # Create scatter plot
        scatter = ax.scatter(positions, gc_contents, scores, c=scores, cmap='viridis', s=60, alpha=0.8)
        
        ax.set_xlabel('Sequence Position (bp)')
        ax.set_ylabel('GC Content')
        ax.set_zlabel('Comprehensive Score')
        ax.set_title('3D Guide Efficiency Landscape')
        
        plt.colorbar(scatter, shrink=0.5, aspect=5)
        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Static guide landscape saved to {save_path}")
    
    def create_phase_space_static(self, sequences: list, save_path: str = "3d_phase_space.png"):
        """Create static 3D phase space analysis."""
        fig = plt.figure(figsize=(12, 9))
        ax = fig.add_subplot(111, projection='3d')
        
        f_values = []
        z_values = []
        phase_bits = []
        iteration_numbers = []
        
        for seq_idx, seq in enumerate(sequences[:2]):  # Limit to 2 sequences
            a, b, c = 0.552, 9.061, 7.389
            calc = ZetaUnfoldCalculator(a, b, c)
            
            for iter_num in range(6):
                f_values.append(float(calc.F))
                z_values.append(float(calc.z))
                phase_bits.append(calc.get_phase_bit())
                iteration_numbers.append(iter_num + seq_idx * 6)
                
                calc = calc.unfold_next()
        
        # Color by phase bit
        colors = ['red' if pb == 0 else 'blue' for pb in phase_bits]
        
        # Plot trajectory
        ax.scatter(f_values, z_values, iteration_numbers, c=colors, s=80, alpha=0.8)
        
        # Connect points with lines
        for i in range(len(f_values)-1):
            if iteration_numbers[i+1] == iteration_numbers[i] + 1:
                ax.plot([f_values[i], f_values[i+1]], 
                       [z_values[i], z_values[i+1]], 
                       [iteration_numbers[i], iteration_numbers[i+1]], 
                       'gray', alpha=0.5)
        
        # Add phase boundaries
        iter_range = range(max(iteration_numbers) + 1)
        z_range = np.linspace(min(z_values), max(z_values), len(iter_range))
        
        ax.plot([0.096] * len(iter_range), z_range, iter_range, 'r-', linewidth=3, alpha=0.7, label='Phase 0 (F≈0.096)')
        ax.plot([0.517] * len(iter_range), z_range, iter_range, 'b-', linewidth=3, alpha=0.7, label='Phase 1 (F≈0.517)')
        
        ax.set_xlabel('F Value')
        ax.set_ylabel('Z Value')
        ax.set_zlabel('Iteration Number')
        ax.set_title('3D Phase Space Analysis: F Alternation Patterns')
        ax.legend()
        
        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Static phase space analysis saved to {save_path}")
    
    def create_z_framework_static(self, save_path: str = "3d_z_framework.png"):
        """Create static 3D Z framework parameter space."""
        fig = plt.figure(figsize=(12, 9))
        ax = fig.add_subplot(111, projection='3d')
        
        # Generate parameter combinations
        np.random.seed(42)  # For reproducibility
        a_values = np.random.uniform(0.3, 0.8, 100)
        b_values = np.random.uniform(8.0, 12.0, 100)
        c_values = np.random.uniform(6.0, 9.0, 100)
        
        z_values = a_values * (b_values / c_values)
        
        # Color by proximity to golden ratio
        phi_conjugate = 0.618033988749895
        golden_distances = np.abs(z_values - phi_conjugate)
        
        scatter = ax.scatter(a_values, b_values, c_values, c=golden_distances, 
                           cmap='RdYlBu_r', s=50, alpha=0.7)
        
        # Highlight optimal region
        optimal_mask = golden_distances < 0.1
        if np.any(optimal_mask):
            ax.scatter(a_values[optimal_mask], b_values[optimal_mask], c_values[optimal_mask], 
                      c='gold', s=100, alpha=1.0, marker='D', edgecolors='black', 
                      label='Optimal Region (δφ < 0.1)')
            ax.legend()
        
        ax.set_xlabel('Parameter A')
        ax.set_ylabel('Parameter B')
        ax.set_zlabel('Parameter C')
        ax.set_title('3D Z Framework Parameter Space')
        
        plt.colorbar(scatter, shrink=0.5, aspect=5, label='Distance to φ⁻¹')
        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Static Z framework plot saved to {save_path}")
    
    def create_curvature_disruption_static(self, sequence: str, save_path: str = "3d_curvature_disruption.png"):
        """Create static 3D curvature disruption analysis."""
        fig = plt.figure(figsize=(12, 9))
        ax = fig.add_subplot(111, projection='3d')
        
        positions = []
        disruption_scores = []
        pam_distances = []
        
        for pos in range(0, len(sequence) - 20, 3):
            subseq = sequence[pos:pos + 20]
            
            # Find PAM sites
            pam_positions = []
            for i in range(len(subseq) - 2):
                if subseq[i+1:i+3] == "GG":
                    pam_positions.append(i)
            
            if pam_positions:
                closest_pam = min(pam_positions, key=lambda x: abs(x - 10))
                pam_dist = abs(closest_pam - 10)
            else:
                pam_dist = 10
            
            # Calculate disruption score
            gc_content = (subseq.count('G') + subseq.count('C')) / len(subseq)
            complexity = len(set(subseq)) / 4.0
            disruption = abs(gc_content - 0.5) + (1 - complexity)
            
            positions.append(pos + 10)
            disruption_scores.append(disruption)
            pam_distances.append(pam_dist)
        
        # Create scatter plot
        scatter = ax.scatter(positions, pam_distances, disruption_scores, 
                           c=disruption_scores, cmap='Reds', s=60, alpha=0.8)
        
        ax.set_xlabel('Sequence Position (bp)')
        ax.set_ylabel('Distance to PAM Site (bp)')
        ax.set_zlabel('Disruption Score')
        ax.set_title('3D Curvature Disruption Analysis')
        
        plt.colorbar(scatter, shrink=0.5, aspect=5)
        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Static curvature disruption plot saved to {save_path}")
    
    def generate_all_static_plots(self, sequence: str, additional_sequences: list = None, 
                                output_dir: str = "./3d_plots/"):
        """Generate all static PNG versions of 3D plots."""
        import os
        
        if additional_sequences is None:
            additional_sequences = [
                "ATCGATCGATCGATCGATCGAAATTTGGGCCCAAACCCGGGAAATTTGGGCCCAAA",
                "GGCCGGCCGGCCGGCCGGCCAAATTTGGGCCCAAACCCGGGAAATTTGGGCCCAAA", 
                "AAATTTAAATTTAAATTTAAGGGGGGCCCCCCAAATTTGGGCCCAAACCCGGGAAA",
                "CCCCCCCCCCCCCCCCCCCCAAATTTGGGCCCAAACCCGGGAAATTTGGGCCCAAA"
            ]
        
        print("🎨 Generating static PNG versions of 3D plots...")
        print("=" * 60)
        
        # Generate all static plots
        self.create_spectral_landscape_static(sequence, os.path.join(output_dir, "static_3d_spectral_landscape.png"))
        self.create_invariant_features_static(additional_sequences, os.path.join(output_dir, "static_3d_invariant_features.png"))
        self.create_golden_proximity_surface_static(os.path.join(output_dir, "static_3d_golden_proximity.png"))
        self.create_phase_space_static(additional_sequences, os.path.join(output_dir, "static_3d_phase_space.png"))
        self.create_guide_landscape_static(sequence, os.path.join(output_dir, "static_3d_guide_landscape.png"))
        self.create_z_framework_static(os.path.join(output_dir, "static_3d_z_framework.png"))
        self.create_curvature_disruption_static(sequence, os.path.join(output_dir, "static_3d_curvature_disruption.png"))
        
        print("\n✅ All static 3D plots generated successfully!")
        print(f"📁 Saved to directory: {output_dir}")


def main():
    """Generate static PNG versions of all 3D plots."""
    target_sequence = "ATGCTGCGGAGACCTGGAGAGAAAGCAGTGGCCGGGGCAGTGGGAGGAGGAGGAGCTGGAAGAGGAGAGAAAGGAGGAGCTGCAGGAGGAGAGGAGGAGGAGGGAGAGGAGGAGCTGGAGCTGAAGCTGGAGCTGGAGCTGGAGAGGAGAGAGGG"
    
    generator = Static3DPlotGenerator()
    generator.generate_all_static_plots(target_sequence)
    

if __name__ == "__main__":
    main()