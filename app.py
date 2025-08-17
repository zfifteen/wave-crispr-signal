"""
Z Framework Empirical Proof Web Application

Flask web application demonstrating empirical validation of the Z Framework's
discrete domain form for DNA sequence analysis.

Features:
- Interactive DNA sequence analysis
- High-precision Z value calculations
- Convergence validation to golden ratio conjugate
- Falsification testing with random perturbations
- Cross-validation with CRISPR benchmarks
- Scientific documentation and methodology
"""

from flask import Flask, render_template, request, jsonify, flash, redirect, url_for
from flask_wtf import FlaskForm
from wtforms import TextAreaField, IntegerField, FloatField, BooleanField, SelectField
from wtforms.validators import DataRequired, Length, NumberRange
import json
import logging
from datetime import datetime
import os
import sys

# Add current directory to path for imports
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from z_framework import ZFrameworkCalculator, format_mpmath_for_json
from applications.crispr_guide_designer import CRISPRGuideDesigner
from applications.wave_crispr_metrics import WaveCRISPRMetrics

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Initialize Flask app
app = Flask(__name__)
app.config['SECRET_KEY'] = 'z-framework-empirical-proof-2025'
app.config['WTF_CSRF_ENABLED'] = False  # Disable CSRF for demo

# Initialize calculators
z_calc = ZFrameworkCalculator(precision_dps=50)
crispr_designer = CRISPRGuideDesigner()
crispr_metrics = WaveCRISPRMetrics()


class DNAAnalysisForm(FlaskForm):
    """Form for DNA sequence analysis input."""
    sequence = TextAreaField(
        'DNA Sequence',
        validators=[DataRequired(), Length(min=10, max=1000)],
        description='Enter DNA sequence (A, T, C, G only). Length: 10-1000 bases.',
        render_kw={"placeholder": "ATGCTGCGGAGACCTGGAGAGAAAGCAGTGG..."}
    )
    
    precision = IntegerField(
        'Calculation Precision (decimal places)',
        validators=[NumberRange(min=15, max=100)],
        default=50,
        description='Precision for mpmath calculations (15-100 decimal places)'
    )
    
    window_size = IntegerField(
        'Density Analysis Window Size',
        validators=[NumberRange(min=5, max=50)],
        default=10,
        description='Window size for density enhancement analysis'
    )
    
    num_perturbations = IntegerField(
        'Number of Falsification Tests',
        validators=[NumberRange(min=10, max=500)],
        default=100,
        description='Number of random perturbations for falsification testing'
    )
    
    perturbation_rate = FloatField(
        'Perturbation Rate',
        validators=[NumberRange(min=0.01, max=0.5)],
        default=0.1,
        description='Fraction of bases to randomly change (0.01-0.5)'
    )
    
    enable_crispr_validation = BooleanField(
        'Enable CRISPR Cross-Validation',
        default=True,
        description='Cross-validate results with CRISPR efficiency benchmarks'
    )


@app.route('/')
def index():
    """Main page with DNA analysis form."""
    form = DNAAnalysisForm()
    return render_template('index.html', form=form)


@app.route('/analyze', methods=['POST'])
def analyze_sequence():
    """Perform comprehensive Z Framework analysis on DNA sequence."""
    form = DNAAnalysisForm()
    
    if not form.validate_on_submit():
        flash('Please correct the form errors and try again.', 'error')
        return render_template('index.html', form=form)
    
    try:
        # Get form data
        sequence = form.sequence.data.upper().strip()
        precision = form.precision.data
        window_size = form.window_size.data
        num_perturbations = form.num_perturbations.data
        perturbation_rate = form.perturbation_rate.data
        enable_crispr = form.enable_crispr_validation.data
        
        # Validate DNA sequence
        valid_bases = set('ATCG')
        if not set(sequence).issubset(valid_bases):
            flash('Invalid DNA sequence. Only A, T, C, G bases are allowed.', 'error')
            return render_template('index.html', form=form)
        
        # Initialize calculator with requested precision
        calculator = ZFrameworkCalculator(precision_dps=precision)
        
        logger.info(f"Starting analysis for sequence of length {len(sequence)}")
        
        # Perform core Z Framework analysis
        z_results = calculator.calculate_z_values(sequence)
        
        # Calculate density enhancement
        density_results = calculator.calculate_density_enhancement(
            sequence, window_size=window_size
        )
        
        # Perform falsification testing
        falsification_results = calculator.perform_falsification_test(
            sequence, 
            num_perturbations=num_perturbations,
            perturbation_rate=perturbation_rate
        )
        
        # CRISPR cross-validation (if enabled)
        crispr_results = None
        if enable_crispr and len(sequence) >= 20:
            try:
                guides = crispr_designer.design_guides(sequence, num_guides=5)
                if guides:
                    # Analyze best guide with Z Framework
                    best_guide = guides[0]
                    guide_z_results = calculator.calculate_z_values(best_guide['sequence'])
                    
                    # Calculate comprehensive CRISPR metrics
                    comprehensive_score = crispr_metrics.calculate_comprehensive_score(
                        best_guide['sequence'], sequence
                    )
                    
                    crispr_results = {
                        'num_guides_found': len(guides),
                        'best_guide': best_guide,
                        'guide_z_analysis': format_mpmath_for_json(guide_z_results),
                        'comprehensive_score': format_mpmath_for_json(comprehensive_score)
                    }
                    
                    logger.info(f"CRISPR analysis completed: {len(guides)} guides found")
                else:
                    crispr_results = {'num_guides_found': 0, 'message': 'No CRISPR guides found'}
            except Exception as e:
                logger.warning(f"CRISPR analysis failed: {e}")
                crispr_results = {'error': str(e)}
        
        # Compile comprehensive results
        analysis_results = {
            'timestamp': datetime.now().isoformat(),
            'sequence_info': {
                'sequence': sequence,
                'length': len(sequence),
                'base_composition': {
                    base: sequence.count(base) for base in 'ATCG'
                },
                'gc_content': (sequence.count('G') + sequence.count('C')) / len(sequence)
            },
            'analysis_parameters': {
                'precision': precision,
                'window_size': window_size,
                'num_perturbations': num_perturbations,
                'perturbation_rate': perturbation_rate,
                'crispr_validation_enabled': enable_crispr
            },
            'z_framework_results': z_results,
            'density_enhancement': density_results,
            'falsification_testing': falsification_results,
            'crispr_cross_validation': crispr_results
        }
        
        # Convert mpmath values for JSON serialization
        z_results = format_mpmath_for_json(z_results)
        density_results = format_mpmath_for_json(density_results)
        falsification_results = format_mpmath_for_json(falsification_results)
        
        analysis_results = format_mpmath_for_json(analysis_results)
        
        logger.info("Analysis completed successfully")
        
        return render_template('results.html', results=analysis_results)
        
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        flash(f'Analysis failed: {str(e)}', 'error')
        return render_template('index.html', form=form)


@app.route('/api/analyze', methods=['POST'])
def api_analyze():
    """API endpoint for programmatic access to Z Framework analysis."""
    try:
        data = request.get_json()
        
        if not data or 'sequence' not in data:
            return jsonify({'error': 'DNA sequence is required'}), 400
        
        sequence = data['sequence'].upper().strip()
        precision = data.get('precision', 50)
        window_size = data.get('window_size', 10)
        num_perturbations = data.get('num_perturbations', 100)
        perturbation_rate = data.get('perturbation_rate', 0.1)
        
        # Validate inputs
        valid_bases = set('ATCG')
        if not set(sequence).issubset(valid_bases):
            return jsonify({'error': 'Invalid DNA sequence'}), 400
        
        if not (10 <= len(sequence) <= 1000):
            return jsonify({'error': 'Sequence length must be between 10 and 1000 bases'}), 400
        
        # Perform analysis
        calculator = ZFrameworkCalculator(precision_dps=precision)
        
        z_results = calculator.calculate_z_values(sequence)
        density_results = calculator.calculate_density_enhancement(sequence, window_size)
        falsification_results = calculator.perform_falsification_test(
            sequence, num_perturbations, perturbation_rate
        )
        
        results = {
            'sequence_length': len(sequence),
            'z_framework_results': z_results,
            'density_enhancement': density_results,
            'falsification_testing': falsification_results
        }
        
        # Convert for JSON
        results = format_mpmath_for_json(results)
        
        return jsonify(results)
        
    except Exception as e:
        logger.error(f"API analysis failed: {e}")
        return jsonify({'error': str(e)}), 500


@app.route('/documentation')
def documentation():
    """Scientific documentation and methodology page."""
    return render_template('documentation.html')


@app.route('/examples')
def examples():
    """Example DNA sequences and use cases."""
    examples = [
        {
            'name': 'PCSK9 Exon 1 (150 bp)',
            'sequence': 'ATGCTGCGGAGACCTGGAGAGAAAGCAGTGGCCGGGGCAGTGGGAGGAGGAGGAGCTGGAAGAGGAGAGAAAGGAGGAGCTGCAGGAGGAGAGGAGGAGGAGGGAGAGGAGGAGCTGGAGCTGAAGCTGGAGCTGGAGCTGGAGAGGAGAGAGGG',
            'description': 'PCSK9 gene exon 1 sequence for cholesterol regulation research'
        },
        {
            'name': 'High GC Content (60 bp)',
            'sequence': 'GGCCGGCCGGCCGGCCATGCGCGCGCGCATGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCC',
            'description': 'High GC content sequence for testing convergence properties'
        },
        {
            'name': 'Low Complexity (40 bp)', 
            'sequence': 'AAATTTAAATTTAAATTTAAATTTAAATTTAAATTTAAAA',
            'description': 'Low complexity sequence for falsification testing'
        },
        {
            'name': 'Balanced Composition (80 bp)',
            'sequence': 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG',
            'description': 'Balanced base composition for baseline measurements'
        }
    ]
    
    return render_template('examples.html', examples=examples)


@app.errorhandler(404)
def not_found(error):
    """404 error handler."""
    return render_template('error.html', 
                         error_code=404, 
                         error_message="Page not found"), 404


@app.errorhandler(500)
def internal_error(error):
    """500 error handler."""
    logger.error(f"Internal server error: {error}")
    return render_template('error.html',
                         error_code=500,
                         error_message="Internal server error"), 500


if __name__ == '__main__':
    # Create templates directory if it doesn't exist
    os.makedirs('templates', exist_ok=True)
    os.makedirs('static', exist_ok=True)
    
    logger.info("Starting Z Framework Empirical Proof Web Application")
    app.run(debug=True, host='0.0.0.0', port=8080)