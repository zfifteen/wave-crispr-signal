"""
3D Visualizations for Wave-CRISPR-Signal Framework

This module creates comprehensive 3D visualizations showcasing the mathematical
and biological features of the wave-CRISPR-signal framework, including:

1. Spectral landscape analysis
2. Invariant feature space visualization  
3. Golden ratio proximity surfaces
4. Phase-space analysis
5. Z-framework parameter space
6. Guide efficiency landscapes
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from typing import List, Dict, Tuple, Optional
import seaborn as sns
from scipy.interpolate import griddata
from scipy.fft import fft
import warnings
warnings.filterwarnings('ignore')

try:
    from applications.crispr_guide_designer import CRISPRGuideDesigner
    from invariant_features import InvariantFeatureSet, ZetaUnfoldCalculator
    from applications.crispr_visualization import CRISPRVisualizer
except ImportError:
    import sys
    sys.path.append('.')
    from applications.crispr_guide_designer import CRISPRGuideDesigner
    from invariant_features import InvariantFeatureSet, ZetaUnfoldCalculator
    from applications.crispr_visualization import CRISPRVisualizer


class CRISPR3DVisualizer:
    """Advanced 3D visualization tools for CRISPR analysis."""
    
    def __init__(self):
        """Initialize the 3D visualizer with required components."""
        self.designer = CRISPRGuideDesigner()
        self.invariant_set = InvariantFeatureSet()
        self.visualizer = CRISPRVisualizer()
        
        # Set plotting style
        plt.style.use('seaborn-v0_8')
        sns.set_palette("viridis")
    
    def create_spectral_landscape_3d(self, sequence: str, 
                                   window_size: int = 20,
                                   save_path: Optional[str] = None) -> go.Figure:
        """
        Create 3D spectral landscape showing frequency content across sequence positions.
        
        Args:
            sequence: DNA sequence to analyze
            window_size: Sliding window size for spectral analysis
            save_path: Optional path to save the plot
            
        Returns:
            Plotly 3D surface figure
        """
        positions = []
        frequencies = []
        magnitudes = []
        
        # Sliding window spectral analysis
        for i in range(0, len(sequence) - window_size + 1, 2):
            subseq = sequence[i:i + window_size]
            wave = self.designer.build_waveform(subseq)
            spectrum = np.abs(fft(wave))
            
            for f_idx, magnitude in enumerate(spectrum[:len(spectrum)//2]):
                positions.append(i + window_size//2)
                frequencies.append(f_idx)
                magnitudes.append(magnitude)
        
        # Create 3D surface
        fig = go.Figure()
        
        # Convert to mesh grid format
        pos_unique = sorted(list(set(positions)))
        freq_unique = sorted(list(set(frequencies)))
        
        Z = np.zeros((len(freq_unique), len(pos_unique)))
        for p, f, m in zip(positions, frequencies, magnitudes):
            p_idx = pos_unique.index(p)
            f_idx = freq_unique.index(f)
            Z[f_idx, p_idx] = m
        
        fig.add_trace(go.Surface(
            x=pos_unique,
            y=freq_unique,
            z=Z,
            colorscale='Viridis',
            name='Spectral Magnitude'
        ))
        
        fig.update_layout(
            title='3D Spectral Landscape Analysis',
            scene=dict(
                xaxis_title='Sequence Position (bp)',
                yaxis_title='Frequency Index',
                zaxis_title='Spectral Magnitude',
                camera=dict(eye=dict(x=1.5, y=1.5, z=1.5))
            ),
            template='plotly_white'
        )
        
        if save_path:
            fig.write_html(save_path)
            print(f"3D spectral landscape saved to {save_path}")
        
        return fig
    
    def create_invariant_feature_space_3d(self, sequences: List[str], 
                                        save_path: Optional[str] = None) -> go.Figure:
        """
        Create 3D scatter plot of invariant feature space.
        
        Args:
            sequences: List of DNA sequences to analyze
            save_path: Optional path to save the plot
            
        Returns:
            Plotly 3D scatter figure
        """
        golden_proximities = []
        phase_bits = []
        delta_entropies = []
        sequence_labels = []
        
        for i, seq in enumerate(sequences):
            features = self.invariant_set.calculate_complete_feature_set(seq)
            golden_proximities.append(features.get('delta_phi', 0))
            phase_bits.append(features.get('phase_bit', 0))
            delta_entropies.append(features.get('delta_phase_entropy', 0))
            sequence_labels.append(f'Seq_{i+1}')
        
        fig = go.Figure()
        
        # Color by phase bit
        colors = ['red' if pb == 0 else 'blue' for pb in phase_bits]
        
        fig.add_trace(go.Scatter3d(
            x=golden_proximities,
            y=phase_bits,
            z=delta_entropies,
            mode='markers',
            marker=dict(
                size=8,
                color=colors,
                opacity=0.8,
                line=dict(width=2, color='black')
            ),
            text=sequence_labels,
            name='Invariant Features'
        ))
        
        fig.update_layout(
            title='3D Invariant Feature Space',
            scene=dict(
                xaxis_title='Golden Proximity (δφ)',
                yaxis_title='Phase Bit (π)',
                zaxis_title='Phase Entropy Difference',
                camera=dict(eye=dict(x=1.2, y=1.2, z=1.2))
            ),
            template='plotly_white'
        )
        
        if save_path:
            fig.write_html(save_path)
            print(f"3D invariant feature space saved to {save_path}")
        
        return fig
    
    def create_golden_ratio_proximity_surface(self, save_path: Optional[str] = None) -> go.Figure:
        """
        Create 3D surface showing golden ratio proximity landscape.
        
        Args:
            save_path: Optional path to save the plot
            
        Returns:
            Plotly 3D surface figure
        """
        # Generate parameter space for Z framework
        a_range = np.linspace(0.3, 0.8, 30)
        b_range = np.linspace(8.0, 12.0, 30)
        c_fixed = 7.389  # e^2
        
        A, B = np.meshgrid(a_range, b_range)
        
        # Calculate Z values and golden proximity
        Z_values = A * (B / c_fixed)
        
        # Golden ratio conjugate target (φ-1 ≈ 0.618)
        phi_conjugate = 0.618033988749895
        delta_phi_surface = np.abs(Z_values - phi_conjugate)
        
        fig = go.Figure()
        
        fig.add_trace(go.Surface(
            x=A,
            y=B,
            z=delta_phi_surface,
            colorscale='RdYlBu_r',
            name='Golden Proximity Distance'
        ))
        
        # Add golden ratio plane
        fig.add_trace(go.Surface(
            x=A,
            y=B,
            z=np.full_like(A, 0),
            colorscale=[[0, 'gold'], [1, 'gold']],
            opacity=0.3,
            name='Golden Ratio Target',
            showscale=False
        ))
        
        fig.update_layout(
            title='3D Golden Ratio Proximity Surface (δφ Landscape)',
            scene=dict(
                xaxis_title='Parameter A',
                yaxis_title='Parameter B', 
                zaxis_title='Distance to φ⁻¹ (δφ)',
                camera=dict(eye=dict(x=1.5, y=1.5, z=1.5))
            ),
            template='plotly_white'
        )
        
        if save_path:
            fig.write_html(save_path)
            print(f"3D golden ratio proximity surface saved to {save_path}")
        
        return fig
    
    def create_phase_space_analysis_3d(self, sequences: List[str],
                                     save_path: Optional[str] = None) -> go.Figure:
        """
        Create 3D phase space visualization showing F alternation patterns.
        
        Args:
            sequences: List of DNA sequences to analyze
            save_path: Optional path to save the plot
            
        Returns:
            Plotly 3D scatter figure
        """
        f_values = []
        z_values = []
        phase_bits = []
        iteration_numbers = []
        
        for seq_idx, seq in enumerate(sequences):
            features = self.invariant_set.calculate_complete_feature_set(seq)
            
            # Get initial Z framework values
            a, b, c = 0.552, 9.061, 7.389
            calc = ZetaUnfoldCalculator(a, b, c)
            
            # Perform several unfolding iterations
            for iter_num in range(8):
                f_values.append(float(calc.F))
                z_values.append(float(calc.z))
                phase_bits.append(calc.get_phase_bit())
                iteration_numbers.append(iter_num + seq_idx * 8)
                
                calc = calc.unfold_next()
        
        fig = go.Figure()
        
        # Color by phase bit
        colors = ['red' if pb == 0 else 'blue' for pb in phase_bits]
        
        fig.add_trace(go.Scatter3d(
            x=f_values,
            y=z_values,
            z=iteration_numbers,
            mode='markers+lines',
            marker=dict(
                size=6,
                color=colors,
                opacity=0.8
            ),
            line=dict(width=3, color='gray'),
            name='Phase Space Trajectory'
        ))
        
        # Add phase boundaries
        phase_0_f = 0.096
        phase_1_f = 0.517
        
        fig.add_trace(go.Scatter3d(
            x=[phase_0_f] * len(set(iteration_numbers)),
            y=np.linspace(min(z_values), max(z_values), len(set(iteration_numbers))),
            z=list(set(iteration_numbers)),
            mode='lines',
            line=dict(width=5, color='red'),
            name='Phase 0 Boundary (F≈0.096)',
            opacity=0.7
        ))
        
        fig.add_trace(go.Scatter3d(
            x=[phase_1_f] * len(set(iteration_numbers)),
            y=np.linspace(min(z_values), max(z_values), len(set(iteration_numbers))),
            z=list(set(iteration_numbers)),
            mode='lines',
            line=dict(width=5, color='blue'),
            name='Phase 1 Boundary (F≈0.517)',
            opacity=0.7
        ))
        
        fig.update_layout(
            title='3D Phase Space Analysis: F Alternation Patterns',
            scene=dict(
                xaxis_title='F Value',
                yaxis_title='Z Value',
                zaxis_title='Iteration Number',
                camera=dict(eye=dict(x=1.2, y=1.2, z=1.2))
            ),
            template='plotly_white'
        )
        
        if save_path:
            fig.write_html(save_path)
            print(f"3D phase space analysis saved to {save_path}")
        
        return fig
    
    def create_guide_efficiency_landscape_3d(self, sequence: str,
                                           save_path: Optional[str] = None) -> go.Figure:
        """
        Create 3D landscape of guide efficiency across position and GC content.
        
        Args:
            sequence: Target DNA sequence
            save_path: Optional path to save the plot
            
        Returns:
            Plotly 3D surface figure
        """
        guides = self.designer.design_guides(sequence, num_guides=50, use_invariants=True)
        
        positions = [g['position'] for g in guides]
        gc_contents = [g['gc_content'] for g in guides]
        scores = [g['comprehensive_score'] if 'comprehensive_score' in g 
                 else g['on_target_score'] for g in guides]
        
        # Create interpolated surface
        pos_range = np.linspace(min(positions), max(positions), 30)
        gc_range = np.linspace(min(gc_contents), max(gc_contents), 30)
        
        POS, GC = np.meshgrid(pos_range, gc_range)
        
        # Interpolate scores
        points = np.column_stack((positions, gc_contents))
        SCORES = griddata(points, scores, (POS, GC), method='cubic', fill_value=0)
        
        fig = go.Figure()
        
        # Add surface
        fig.add_trace(go.Surface(
            x=POS,
            y=GC,
            z=SCORES,
            colorscale='Viridis',
            name='Guide Efficiency'
        ))
        
        # Add actual guide points
        fig.add_trace(go.Scatter3d(
            x=positions,
            y=gc_contents,
            z=scores,
            mode='markers',
            marker=dict(
                size=5,
                color='red',
                opacity=0.8
            ),
            name='Actual Guides'
        ))
        
        fig.update_layout(
            title='3D Guide Efficiency Landscape',
            scene=dict(
                xaxis_title='Sequence Position (bp)',
                yaxis_title='GC Content',
                zaxis_title='Comprehensive Score',
                camera=dict(eye=dict(x=1.5, y=1.5, z=1.5))
            ),
            template='plotly_white'
        )
        
        if save_path:
            fig.write_html(save_path)
            print(f"3D guide efficiency landscape saved to {save_path}")
        
        return fig
    
    def create_z_framework_parameter_space_3d(self, save_path: Optional[str] = None) -> go.Figure:
        """
        Create 3D visualization of Z framework parameter space.
        
        Args:
            save_path: Optional path to save the plot
            
        Returns:
            Plotly 3D scatter figure
        """
        # Generate parameter combinations
        a_values = np.random.uniform(0.3, 0.8, 100)
        b_values = np.random.uniform(8.0, 12.0, 100)
        c_values = np.random.uniform(6.0, 9.0, 100)
        
        z_values = a_values * (b_values / c_values)
        
        # Color by proximity to golden ratio
        phi_conjugate = 0.618033988749895
        golden_distances = np.abs(z_values - phi_conjugate)
        
        fig = go.Figure()
        
        fig.add_trace(go.Scatter3d(
            x=a_values,
            y=b_values,
            z=c_values,
            mode='markers',
            marker=dict(
                size=5,
                color=golden_distances,
                colorscale='RdYlBu_r',
                opacity=0.8,
                colorbar=dict(title='Distance to φ⁻¹')
            ),
            text=[f'Z={z:.3f}, δφ={d:.3f}' for z, d in zip(z_values, golden_distances)],
            name='Parameter Space'
        ))
        
        # Add optimal region indicator
        optimal_mask = golden_distances < 0.1
        if np.any(optimal_mask):
            fig.add_trace(go.Scatter3d(
                x=a_values[optimal_mask],
                y=b_values[optimal_mask],
                z=c_values[optimal_mask],
                mode='markers',
                marker=dict(
                    size=8,
                    color='gold',
                    opacity=1.0,
                    symbol='diamond'
                ),
                name='Optimal Region (δφ < 0.1)'
            ))
        
        fig.update_layout(
            title='3D Z Framework Parameter Space',
            scene=dict(
                xaxis_title='Parameter A',
                yaxis_title='Parameter B',
                zaxis_title='Parameter C',
                camera=dict(eye=dict(x=1.2, y=1.2, z=1.2))
            ),
            template='plotly_white'
        )
        
        if save_path:
            fig.write_html(save_path)
            print(f"3D Z framework parameter space saved to {save_path}")
        
        return fig
    
    def create_curvature_disruption_3d(self, sequence: str, 
                                     save_path: Optional[str] = None) -> go.Figure:
        """
        Create 3D visualization of curvature disruption analysis.
        
        Args:
            sequence: DNA sequence to analyze
            save_path: Optional path to save the plot
            
        Returns:
            Plotly 3D surface figure
        """
        # Analyze curvature disruption across sequence
        positions = []
        disruption_scores = []
        pam_distances = []
        
        for pos in range(0, len(sequence) - 20, 2):
            subseq = sequence[pos:pos + 20]
            
            # Find PAM sites in the region
            pam_positions = []
            for i in range(len(subseq) - 2):
                if subseq[i+1:i+3] == "GG":  # NGG PAM
                    pam_positions.append(i)
            
            if pam_positions:
                closest_pam = min(pam_positions, key=lambda x: abs(x - 10))
                pam_dist = abs(closest_pam - 10)
            else:
                pam_dist = 10  # Max distance if no PAM
            
            # Calculate disruption score (simplified)
            gc_content = (subseq.count('G') + subseq.count('C')) / len(subseq)
            complexity = len(set(subseq)) / 4.0  # Normalized complexity
            disruption = abs(gc_content - 0.5) + (1 - complexity)
            
            positions.append(pos + 10)
            disruption_scores.append(disruption)
            pam_distances.append(pam_dist)
        
        # Create surface
        pos_range = np.linspace(min(positions), max(positions), 20)
        dist_range = np.linspace(0, max(pam_distances), 15)
        
        POS, DIST = np.meshgrid(pos_range, dist_range)
        
        # Interpolate disruption scores
        points = np.column_stack((positions, pam_distances))
        DISRUPTION = griddata(points, disruption_scores, (POS, DIST), 
                            method='linear', fill_value=np.mean(disruption_scores))
        
        fig = go.Figure()
        
        fig.add_trace(go.Surface(
            x=POS,
            y=DIST,
            z=DISRUPTION,
            colorscale='Reds',
            name='Curvature Disruption'
        ))
        
        # Add actual data points
        fig.add_trace(go.Scatter3d(
            x=positions,
            y=pam_distances,
            z=disruption_scores,
            mode='markers',
            marker=dict(
                size=4,
                color='blue',
                opacity=0.8
            ),
            name='Analysis Points'
        ))
        
        fig.update_layout(
            title='3D Curvature Disruption Analysis',
            scene=dict(
                xaxis_title='Sequence Position (bp)',
                yaxis_title='Distance to PAM Site (bp)',
                zaxis_title='Disruption Score',
                camera=dict(eye=dict(x=1.5, y=1.5, z=1.5))
            ),
            template='plotly_white'
        )
        
        if save_path:
            fig.write_html(save_path)
            print(f"3D curvature disruption analysis saved to {save_path}")
        
        return fig
    
    def generate_comprehensive_3d_plots(self, sequence: str, 
                                      additional_sequences: Optional[List[str]] = None,
                                      output_dir: str = "./3d_plots/") -> Dict[str, str]:
        """
        Generate all 3D plots for comprehensive analysis.
        
        Args:
            sequence: Primary DNA sequence for analysis
            additional_sequences: Additional sequences for comparative analysis
            output_dir: Directory to save plots
            
        Returns:
            Dictionary mapping plot names to file paths
        """
        import os
        os.makedirs(output_dir, exist_ok=True)
        
        if additional_sequences is None:
            additional_sequences = [
                "ATCGATCGATCGATCGATCGAAATTTGGGCCCAAACCCGGGAAATTTGGGCCCAAA",
                "GGCCGGCCGGCCGGCCGGCCAAATTTGGGCCCAAACCCGGGAAATTTGGGCCCAAA", 
                "AAATTTAAATTTAAATTTAAGGGGGGCCCCCCAAATTTGGGCCCAAACCCGGGAAA",
                "CCCCCCCCCCCCCCCCCCCCAAATTTGGGCCCAAACCCGGGAAATTTGGGCCCAAA"
            ]
        
        plot_files = {}
        
        print("🎨 Generating comprehensive 3D visualizations...")
        print("=" * 60)
        
        # 1. Spectral landscape
        print("📊 Creating 3D spectral landscape...")
        fig1 = self.create_spectral_landscape_3d(sequence)
        path1 = os.path.join(output_dir, "3d_spectral_landscape.html")
        fig1.write_html(path1)
        plot_files['spectral_landscape'] = path1
        
        # 2. Invariant feature space
        print("🔬 Creating 3D invariant feature space...")
        fig2 = self.create_invariant_feature_space_3d(additional_sequences)
        path2 = os.path.join(output_dir, "3d_invariant_features.html")
        fig2.write_html(path2)
        plot_files['invariant_features'] = path2
        
        # 3. Golden ratio proximity surface
        print("🌟 Creating 3D golden ratio proximity surface...")
        fig3 = self.create_golden_ratio_proximity_surface()
        path3 = os.path.join(output_dir, "3d_golden_proximity.html")
        fig3.write_html(path3)
        plot_files['golden_proximity'] = path3
        
        # 4. Phase space analysis
        print("🌀 Creating 3D phase space analysis...")
        fig4 = self.create_phase_space_analysis_3d(additional_sequences[:2])
        path4 = os.path.join(output_dir, "3d_phase_space.html")
        fig4.write_html(path4)
        plot_files['phase_space'] = path4
        
        # 5. Guide efficiency landscape
        print("🎯 Creating 3D guide efficiency landscape...")
        fig5 = self.create_guide_efficiency_landscape_3d(sequence)
        path5 = os.path.join(output_dir, "3d_guide_landscape.html")
        fig5.write_html(path5)
        plot_files['guide_landscape'] = path5
        
        # 6. Z framework parameter space
        print("📐 Creating 3D Z framework parameter space...")
        fig6 = self.create_z_framework_parameter_space_3d()
        path6 = os.path.join(output_dir, "3d_z_framework.html")
        fig6.write_html(path6)
        plot_files['z_framework'] = path6
        
        # 7. Curvature disruption analysis
        print("📈 Creating 3D curvature disruption analysis...")
        fig7 = self.create_curvature_disruption_3d(sequence)
        path7 = os.path.join(output_dir, "3d_curvature_disruption.html")
        fig7.write_html(path7)
        plot_files['curvature_disruption'] = path7
        
        print("\n✅ All 3D visualizations generated successfully!")
        print(f"📁 Saved to directory: {output_dir}")
        print("\n📋 Generated plots:")
        for name, path in plot_files.items():
            print(f"  • {name}: {path}")
        
        return plot_files


def main():
    """Generate comprehensive 3D plots for the wave-CRISPR-signal framework."""
    
    # Example target sequence
    target_sequence = "ATGCTGCGGAGACCTGGAGAGAAAGCAGTGGCCGGGGCAGTGGGAGGAGGAGGAGCTGGAAGAGGAGAGAAAGGAGGAGCTGCAGGAGGAGAGGAGGAGGAGGGAGAGGAGGAGCTGGAGCTGAAGCTGGAGCTGGAGCTGGAGAGGAGAGAGGG"
    
    # Create visualizer
    visualizer = CRISPR3DVisualizer()
    
    # Generate all plots
    plot_files = visualizer.generate_comprehensive_3d_plots(target_sequence)
    
    print(f"\n🎉 Complete! Generated {len(plot_files)} 3D visualizations.")
    print("🌐 Open the HTML files in a web browser to view interactive 3D plots.")
    

if __name__ == "__main__":
    main()