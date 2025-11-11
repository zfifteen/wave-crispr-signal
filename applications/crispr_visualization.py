"""
CRISPR Visualization Module

This module provides visualization tools for CRISPR guide design analysis
using signal-theoretic features and traditional metrics.
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from typing import List, Dict, Optional, Tuple

try:
    from .crispr_guide_designer import CRISPRGuideDesigner
except ImportError:
    # Handle relative import for direct execution
    import sys

    sys.path.append(".")
    from crispr_guide_designer import CRISPRGuideDesigner


class CRISPRVisualizer:
    """Visualization tools for CRISPR guide analysis."""

    def __init__(self, designer: Optional[CRISPRGuideDesigner] = None):
        """Initialize visualizer with optional designer instance."""
        self.designer = designer or CRISPRGuideDesigner()

        # Set plotting style
        plt.style.use("seaborn-v0_8")
        sns.set_palette("husl")

    def plot_sequence_spectrum(
        self,
        sequence: str,
        title: str = "DNA Sequence Spectrum",
        figsize: Tuple[int, int] = (12, 6),
    ) -> plt.Figure:
        """
        Plot frequency spectrum of DNA sequence.

        Args:
            sequence: DNA sequence
            title: Plot title
            figsize: Figure size

        Returns:
            Matplotlib figure
        """
        wave = self.designer.build_waveform(sequence)
        spec = self.designer.compute_spectrum(wave)

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

        # Waveform plot
        ax1.plot(np.real(wave), label="Real", alpha=0.8)
        ax1.plot(np.imag(wave), label="Imaginary", alpha=0.8)
        ax1.set_title("Complex Waveform")
        ax1.set_xlabel("Position")
        ax1.set_ylabel("Amplitude")
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # Spectrum plot
        ax2.plot(spec[: len(spec) // 2], linewidth=2)
        ax2.set_title("Frequency Spectrum")
        ax2.set_xlabel("Frequency Index")
        ax2.set_ylabel("Magnitude")
        ax2.grid(True, alpha=0.3)

        # Add spectral features
        entropy_val = self.designer.normalized_entropy(spec)
        sidelobes = self.designer.count_sidelobes(spec)
        ax2.text(
            0.02,
            0.98,
            f"Entropy: {entropy_val:.3f}\nSidelobes: {sidelobes}",
            transform=ax2.transAxes,
            verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
        )

        plt.suptitle(title)
        plt.tight_layout()
        return fig

    def plot_guide_scores(
        self,
        guides: List[Dict],
        title: str = "Guide Score Distribution",
        figsize: Tuple[int, int] = (10, 6),
    ) -> plt.Figure:
        """
        Plot distribution of guide scores.

        Args:
            guides: List of guide dictionaries
            title: Plot title
            figsize: Figure size

        Returns:
            Matplotlib figure
        """
        if not guides:
            raise ValueError("No guides provided for plotting")

        # Extract scores
        scores = [g["on_target_score"] for g in guides]
        gc_contents = [g["gc_content"] for g in guides]
        positions = [g["position"] for g in guides]

        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=figsize)

        # Score distribution
        ax1.hist(scores, bins=min(len(scores), 20), alpha=0.7, edgecolor="black")
        ax1.set_title("On-Target Score Distribution")
        ax1.set_xlabel("Score")
        ax1.set_ylabel("Frequency")
        ax1.axvline(
            np.mean(scores),
            color="red",
            linestyle="--",
            label=f"Mean: {np.mean(scores):.3f}",
        )
        ax1.legend()

        # GC content distribution
        ax2.hist(
            gc_contents,
            bins=min(len(gc_contents), 20),
            alpha=0.7,
            edgecolor="black",
            color="orange",
        )
        ax2.set_title("GC Content Distribution")
        ax2.set_xlabel("GC Content")
        ax2.set_ylabel("Frequency")
        ax2.axvline(
            np.mean(gc_contents),
            color="red",
            linestyle="--",
            label=f"Mean: {np.mean(gc_contents):.3f}",
        )
        ax2.legend()

        # Score vs GC content
        ax3.scatter(gc_contents, scores, alpha=0.7)
        ax3.set_xlabel("GC Content")
        ax3.set_ylabel("On-Target Score")
        ax3.set_title("Score vs GC Content")

        # Score vs position
        ax4.scatter(positions, scores, alpha=0.7, color="green")
        ax4.set_xlabel("Position")
        ax4.set_ylabel("On-Target Score")
        ax4.set_title("Score vs Position")

        plt.suptitle(title)
        plt.tight_layout()
        return fig

    def plot_guide_comparison(
        self,
        guides: List[Dict],
        max_guides: int = 10,
        figsize: Tuple[int, int] = (12, 8),
    ) -> plt.Figure:
        """
        Compare multiple guides across different metrics.

        Args:
            guides: List of guide dictionaries
            max_guides: Maximum number of guides to show
            figsize: Figure size

        Returns:
            Matplotlib figure
        """
        guides = guides[:max_guides]  # Limit number of guides

        # Prepare data
        guide_names = [f"Guide_{i+1}" for i in range(len(guides))]
        on_target_scores = [g["on_target_score"] for g in guides]
        gc_contents = [g["gc_content"] for g in guides]

        # Get repair outcomes if available
        nhej_probs = []
        mmej_probs = []
        hdr_effs = []

        for guide in guides:
            if "repair_outcomes" in guide:
                repair = guide["repair_outcomes"]
                nhej_probs.append(repair.get("nhej_probability", 0))
                mmej_probs.append(repair.get("mmej_probability", 0))
                hdr_effs.append(repair.get("hdr_efficiency", 0))

        # Create subplot layout
        if nhej_probs:  # Has repair data
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=figsize)
        else:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(figsize[0], figsize[1] // 2))

        # On-target scores
        bars1 = ax1.bar(guide_names, on_target_scores, alpha=0.8)
        ax1.set_title("On-Target Efficiency Scores")
        ax1.set_ylabel("Score")
        ax1.tick_params(axis="x", rotation=45)

        # Color bars by score
        for bar, score in zip(bars1, on_target_scores):
            if score > 0.7:
                bar.set_color("green")
            elif score > 0.5:
                bar.set_color("orange")
            else:
                bar.set_color("red")

        # GC content
        ax2.bar(guide_names, gc_contents, alpha=0.8, color="skyblue")
        ax2.set_title("GC Content")
        ax2.set_ylabel("GC Content")
        ax2.tick_params(axis="x", rotation=45)
        ax2.axhline(0.5, color="red", linestyle="--", alpha=0.7, label="Optimal")
        ax2.legend()

        # Repair outcomes if available
        if nhej_probs:
            # Stacked bar for repair pathways
            width = 0.6
            ax3.bar(guide_names, nhej_probs, width, label="NHEJ", alpha=0.8)
            ax3.bar(
                guide_names,
                mmej_probs,
                width,
                bottom=nhej_probs,
                label="MMEJ",
                alpha=0.8,
            )
            ax3.bar(
                guide_names,
                hdr_effs,
                width,
                bottom=[n + m for n, m in zip(nhej_probs, mmej_probs)],
                label="HDR",
                alpha=0.8,
            )
            ax3.set_title("Repair Pathway Predictions")
            ax3.set_ylabel("Probability")
            ax3.tick_params(axis="x", rotation=45)
            ax3.legend()

            # HDR efficiency specifically
            ax4.bar(guide_names, hdr_effs, alpha=0.8, color="purple")
            ax4.set_title("HDR Efficiency")
            ax4.set_ylabel("HDR Efficiency")
            ax4.tick_params(axis="x", rotation=45)

        plt.suptitle("CRISPR Guide Comparison")
        plt.tight_layout()
        return fig

    def plot_interactive_guide_map(
        self, sequence: str, guides: List[Dict], title: str = "CRISPR Guide Map"
    ) -> go.Figure:
        """
        Create interactive guide map using Plotly.

        Args:
            sequence: Target DNA sequence
            guides: List of guide dictionaries
            title: Plot title

        Returns:
            Plotly figure
        """
        fig = go.Figure()

        # Add sequence track
        positions = list(range(len(sequence)))

        # Convert sequence to numeric for plotting (A=1, T=2, C=3, G=4)
        seq_numeric = [
            {"A": 1, "T": 2, "C": 3, "G": 4}.get(base, 0) for base in sequence.upper()
        ]

        fig.add_trace(
            go.Scatter(
                x=positions,
                y=seq_numeric,
                mode="lines+markers",
                name="DNA Sequence",
                line=dict(width=2),
                marker=dict(size=3),
                hovertemplate="Position: %{x}<br>Base: %{text}<extra></extra>",
                text=list(sequence.upper()),
            )
        )

        # Add guide regions
        colors = px.colors.qualitative.Set1
        for i, guide in enumerate(guides[:10]):  # Limit to 10 guides for clarity
            start_pos = guide["position"]
            end_pos = start_pos + guide["length"]
            score = guide["on_target_score"]

            # Add guide region
            fig.add_trace(
                go.Scatter(
                    x=[start_pos, end_pos, end_pos, start_pos, start_pos],
                    y=[0.5, 0.5, 4.5, 4.5, 0.5],
                    fill="toself",
                    fillcolor=colors[i % len(colors)],
                    opacity=score,  # Opacity based on score
                    name=f"Guide {i+1} (Score: {score:.3f})",
                    line=dict(width=2),
                    hovertemplate=f'Guide {i+1}<br>Position: {start_pos}-{end_pos}<br>Score: {score:.3f}<br>Sequence: {guide["sequence"]}<extra></extra>',
                )
            )

            # Add PAM site marker
            pam_pos = guide["pam_position"]
            fig.add_trace(
                go.Scatter(
                    x=[pam_pos],
                    y=[5],
                    mode="markers",
                    marker=dict(
                        size=12,
                        symbol="diamond",
                        color=colors[i % len(colors)],
                        line=dict(width=2, color="black"),
                    ),
                    name=f"PAM {i+1}",
                    showlegend=False,
                    hovertemplate=f'PAM Site {i+1}<br>Position: {pam_pos}<br>Sequence: {guide["pam_sequence"]}<extra></extra>',
                )
            )

        fig.update_layout(
            title=title,
            xaxis_title="Position (bp)",
            yaxis_title="Base Type",
            yaxis=dict(
                tickmode="array", tickvals=[1, 2, 3, 4], ticktext=["A", "T", "C", "G"]
            ),
            hovermode="closest",
            template="plotly_white",
        )

        return fig

    def plot_spectral_heatmap(
        self, sequences: Dict[str, str], figsize: Tuple[int, int] = (12, 8)
    ) -> plt.Figure:
        """
        Create heatmap of spectral features across multiple sequences.

        Args:
            sequences: Dictionary mapping sequence names to sequences
            figsize: Figure size

        Returns:
            Matplotlib figure
        """
        # Calculate spectral features for each sequence
        features = []
        seq_names = []

        for name, seq in sequences.items():
            wave = self.designer.build_waveform(seq)
            spec = self.designer.compute_spectrum(wave)

            entropy_val = self.designer.normalized_entropy(spec)
            sidelobes = self.designer.count_sidelobes(spec)
            gc_content = (seq.count("G") + seq.count("C")) / len(seq)

            # Take first 50 spectral components
            spec_features = (
                spec[:50] if len(spec) >= 50 else np.pad(spec, (0, 50 - len(spec)))
            )

            features.append([entropy_val, sidelobes, gc_content] + list(spec_features))
            seq_names.append(name)

        # Create feature matrix
        feature_matrix = np.array(features)
        feature_names = ["Entropy", "Sidelobes", "GC_Content"] + [
            f"Freq_{i}" for i in range(50)
        ]

        # Create heatmap
        fig, ax = plt.subplots(figsize=figsize)

        # Normalize features for better visualization
        normalized_features = (
            feature_matrix - feature_matrix.mean(axis=0)
        ) / feature_matrix.std(axis=0)

        im = ax.imshow(normalized_features.T, cmap="RdBu_r", aspect="auto")

        # Set labels
        ax.set_xticks(range(len(seq_names)))
        ax.set_xticklabels(seq_names, rotation=45, ha="right")
        ax.set_yticks(range(len(feature_names)))
        ax.set_yticklabels(feature_names)

        # Add colorbar
        plt.colorbar(im, ax=ax, label="Normalized Feature Value")

        ax.set_title("Spectral Features Heatmap")
        ax.set_xlabel("Sequences")
        ax.set_ylabel("Features")

        plt.tight_layout()
        return fig

    def create_dashboard(
        self, sequence: str, guides: List[Dict], output_file: Optional[str] = None
    ) -> go.Figure:
        """
        Create comprehensive dashboard for CRISPR analysis.

        Args:
            sequence: Target DNA sequence
            guides: List of guide dictionaries
            output_file: Optional HTML output file

        Returns:
            Plotly figure
        """
        # Create subplots
        fig = make_subplots(
            rows=2,
            cols=2,
            subplot_titles=(
                "Guide Positions",
                "Score Distribution",
                "GC Content vs Score",
                "Sequence Spectrum",
            ),
            specs=[
                [{"secondary_y": False}, {"secondary_y": False}],
                [{"secondary_y": False}, {"secondary_y": False}],
            ],
        )

        # Guide positions
        positions = [g["position"] for g in guides]
        scores = [g["on_target_score"] for g in guides]

        fig.add_trace(
            go.Scatter(
                x=positions,
                y=scores,
                mode="markers+text",
                text=[f"G{i+1}" for i in range(len(guides))],
                textposition="top center",
                marker=dict(
                    size=10, colorscale="Viridis", color=scores, showscale=True
                ),
                name="Guides",
            ),
            row=1,
            col=1,
        )

        # Score distribution
        fig.add_trace(
            go.Histogram(x=scores, nbinsx=10, name="Score Distribution"), row=1, col=2
        )

        # GC content vs score
        gc_contents = [g["gc_content"] for g in guides]
        fig.add_trace(
            go.Scatter(
                x=gc_contents,
                y=scores,
                mode="markers",
                marker=dict(size=8),
                name="GC vs Score",
            ),
            row=2,
            col=1,
        )

        # Sequence spectrum
        wave = self.designer.build_waveform(sequence)
        spec = self.designer.compute_spectrum(wave)
        freq_range = list(range(min(50, len(spec))))

        fig.add_trace(
            go.Scatter(
                x=freq_range, y=spec[: len(freq_range)], mode="lines", name="Spectrum"
            ),
            row=2,
            col=2,
        )

        # Update layout
        fig.update_layout(
            title_text="CRISPR Guide Analysis Dashboard",
            title_x=0.5,
            showlegend=False,
            template="plotly_white",
        )

        # Update axes labels
        fig.update_xaxes(title_text="Position (bp)", row=1, col=1)
        fig.update_yaxes(title_text="Score", row=1, col=1)

        fig.update_xaxes(title_text="Score", row=1, col=2)
        fig.update_yaxes(title_text="Count", row=1, col=2)

        fig.update_xaxes(title_text="GC Content", row=2, col=1)
        fig.update_yaxes(title_text="Score", row=2, col=1)

        fig.update_xaxes(title_text="Frequency", row=2, col=2)
        fig.update_yaxes(title_text="Magnitude", row=2, col=2)

        # Save if requested
        if output_file:
            fig.write_html(output_file)
            print(f"Dashboard saved to {output_file}")

        return fig


def main():
    """Example usage of CRISPR visualization tools."""
    # Example sequence
    sequence = "ATGCTGCGGAGACCTGGAGAGAAAGCAGTGGCCGGGGCAGTGGGAGGAGGAGGAGCTGGAAGAGGAGAGAAAGGAGGAGCTGCAGGAGGAGAGGAGGAGGAGGGAGAGGAGGAGCTGGAGCTGAAGCTGGAGCTGGAGCTGGAGAGGAGAGAGGG"

    # Design guides
    designer = CRISPRGuideDesigner()
    guides = designer.design_guides(sequence, num_guides=5)

    # Add repair predictions
    for guide in guides:
        context_start = max(0, guide["position"] - 20)
        context_end = min(len(sequence), guide["position"] + guide["length"] + 20)
        context_seq = sequence[context_start:context_end]
        repair = designer.predict_repair_outcomes(guide["sequence"], context_seq)
        guide["repair_outcomes"] = repair

    # Create visualizer
    visualizer = CRISPRVisualizer(designer)

    print("🎨 Creating CRISPR visualizations...")

    # Plot sequence spectrum
    visualizer.plot_sequence_spectrum(sequence, "Target Sequence Spectrum")
    plt.savefig("sequence_spectrum.png", dpi=150, bbox_inches="tight")
    plt.show()

    # Plot guide scores
    visualizer.plot_guide_scores(guides, "Guide Analysis")
    plt.savefig("guide_scores.png", dpi=150, bbox_inches="tight")
    plt.show()

    # Plot guide comparison
    visualizer.plot_guide_comparison(guides, "Guide Comparison")
    plt.savefig("guide_comparison.png", dpi=150, bbox_inches="tight")
    plt.show()

    # Create interactive dashboard
    dashboard = visualizer.create_dashboard(sequence, guides, "crispr_dashboard.html")
    dashboard.show()

    print("✅ Visualizations complete!")
    print(
        "📁 Files saved: sequence_spectrum.png, guide_scores.png, guide_comparison.png, crispr_dashboard.html"
    )


if __name__ == "__main__":
    main()
