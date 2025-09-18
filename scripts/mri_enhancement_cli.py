#!/usr/bin/env python3
"""
Command Line Interface for MRI Enhancement Tool

Provides CLI access to DICOM-based spine MRI analysis for differentiating
traumatic vs degenerative changes.

Metric rationale and appropriateness:
- The SpineMRIAnalyzer is designed to detect chronic sequelae of old traumatic injuries (e.g., healed compression fractures, post-traumatic syringomyelia) and distinguish them from degenerative spine disorders.
- Key metrics:
    * Vertebral height loss threshold (default 0.10 or 10%): Sensitive to subtle, stable residuals of old trauma (see Genant et al., https://pubs.rsna.org/doi/10.1148/radiol.1993.186.2.8421738).
    * Syrinx width ratio (default 0.70-0.90 or 70-90%): Targets significant post-traumatic syringomyelia, as supported by chronic spinal cord injury literature.
- These thresholds are partially appropriate for chronic injury detection: the height loss threshold is more sensitive than standard (20-25%), aiding subtle trauma detection but risking false positives; the syrinx ratio aligns with significant post-traumatic cavities.
- See detailed discussion in project documentation and:
    * https://pubmed.ncbi.nlm.nih.gov/8421738/ (Genant method)
    * https://pubmed.ncbi.nlm.nih.gov/15179352/ (PTS ratios)
"""

import argparse
import os
import sys
import logging
from pathlib import Path
from typing import Optional

# Add src to path for imports
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', '..'))

from src.applications.mri_enhancement import (
    MRIConfig, DICOMScanner, SpineMRIAnalyzer, 
    AnnotationGenerator, ReportGenerator
)


def setup_logging(log_level: str = "INFO", log_file: Optional[str] = None):
    """Setup logging configuration."""
    level = getattr(logging, log_level.upper(), logging.INFO)
    
    handlers = [logging.StreamHandler()]
    if log_file:
        handlers.append(logging.FileHandler(log_file))
        
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=handlers
    )


def main():
    """Main CLI entry point.

    Notes on metric appropriateness:
    - Height loss threshold (default 0.10):
        * Used to flag wedge/burst deformities from old axial trauma.
        * 10% is more sensitive than standard radiological criteria (Genant: 20-25% for mild VCFs), which may help detect subtle chronic trauma but can over-identify age-related or degenerative changes.
        * In chronic trauma, focal, stable loss is typical; in degeneration, loss is gradual, multilevel, and often <10-15% per level.
        * Literature: https://pubmed.ncbi.nlm.nih.gov/8421738/
    - Syrinx width ratio (default 0.70-0.90):
        * Measures syrinx diameter relative to cord width on axial T2.
        * 70-90% range is highly appropriate for chronic post-traumatic syringomyelia (PTS), as ratios >70% indicate significant, symptomatic cavities.
        * Chronic PTS: stable, high T2 signal, no enhancement; ratios in this range correlate with risk and are detectable in 12-22% of chronic SCI patients.
        * Literature: https://pubmed.ncbi.nlm.nih.gov/15179352/
    """
    parser = argparse.ArgumentParser(
        description="MRI Enhancement Tool - DICOM-based Spine Analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze single DICOM folder
  %(prog)s -i /path/to/dicom/folder -o /path/to/output
  
  # Analyze with custom configuration
  %(prog)s -i /path/to/dicom -o /path/to/output -c custom_config.json
  
  # Generate report only (skip image processing)
  %(prog)s -i /path/to/dicom -o /path/to/output --report-only
  
  # Specify spine region for optimized analysis
  %(prog)s -i /path/to/dicom -o /path/to/output --spine-region cervical
        """
    )
    
    # Input/Output arguments
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input folder path containing DICOM files (scanned recursively)'
    )
    
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output directory for results (annotations, report, etc.)'
    )
    
    # Configuration arguments
    parser.add_argument(
        '-c', '--config',
        help='Path to custom configuration JSON file'
    )
    
    parser.add_argument(
        '--spine-region',
        choices=['cervical', 'thoracic', 'lumbar', 'c-spine', 't-spine', 'l-spine'],
        help='Specific spine region for optimized analysis parameters'
    )
    
    # Processing options
    parser.add_argument(
        '--report-only',
        action='store_true',
        help='Generate report only (skip DICOM scanning and analysis)'
    )
    
    parser.add_argument(
        '--no-annotations',
        action='store_true',
        help='Skip annotation generation (analysis and report only)'
    )
    
    parser.add_argument(
        '--no-citations',
        action='store_true',
        help='Disable PubMed citation lookup'
    )
    
    # Detection thresholds
    parser.add_argument(
        '--height-loss-threshold',
        type=float,
        default=0.10,
        help=(
            'Height loss threshold for wedge fracture detection (default: 0.10). '
            'A lower threshold (10%) increases sensitivity for subtle, healed compression fractures from old trauma, '
            'but may over-detect mild degenerative changes. '
            'Standard radiological criteria (Genant) use 20-25% for mild VCFs. '
            'See: https://pubmed.ncbi.nlm.nih.gov/8421738/'
        )
    )
    
    parser.add_argument(
        '--syrinx-width-ratio',
        type=float,
        nargs=2,
        default=[0.70, 0.90],
        metavar=('MIN', 'MAX'),
        help=(
            'Syrinx width ratio range for detection (default: 0.70 0.90). '
            'Measures syrinx diameter relative to cord width on axial T2. '
            '70-90% targets significant post-traumatic syringomyelia (PTS) in chronic injury, '
            'as supported by SCI literature. '
            'See: https://pubmed.ncbi.nlm.nih.gov/15179352/'
        )
    )
    
    # HIPAA compliance
    parser.add_argument(
        '--no-anonymize',
        action='store_true',
        help='Disable patient data anonymization (REQUIRES: MRI_TOOL_ALLOW_PHI=AUTHORIZED_RESEARCH_USE environment variable for security)'
    )
    
    # Logging options
    parser.add_argument(
        '--log-level',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        default='INFO',
        help='Logging level (default: INFO)'
    )
    
    parser.add_argument(
        '--log-file',
        help='Path to log file (logs to console if not specified)'
    )
    
    # Performance options
    parser.add_argument(
        '--max-parallel-jobs',
        type=int,
        default=4,
        help='Maximum parallel processing jobs (default: 4)'
    )
    
    parser.add_argument(
        '--memory-limit',
        type=int,
        default=2048,
        help='Memory limit in MB (default: 2048)'
    )
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.log_level, args.log_file)
    logger = logging.getLogger(__name__)
    
    try:
        # Validate input/output paths
        if not os.path.exists(args.input):
            logger.error(f"Input path does not exist: {args.input}")
            return 1
            
        os.makedirs(args.output, exist_ok=True)
        
        # Load configuration
        if args.config:
            if not os.path.exists(args.config):
                logger.error(f"Configuration file not found: {args.config}")
                return 1
            config = MRIConfig.from_json(args.config)
        else:
            # Use default or spine-specific configuration
            if args.spine_region:
                from src.applications.mri_enhancement.config import get_spine_config
                config = get_spine_config(args.spine_region)
            else:
                config = MRIConfig()
                
        # Apply command line overrides
        # Height loss threshold (default 0.10) and syrinx width ratio (default 0.70-0.90) are set here.
        # These values are critical for distinguishing chronic traumatic findings from degenerative changes.
        config.detection_thresholds.height_loss_threshold = args.height_loss_threshold
        config.detection_thresholds.syrinx_width_ratio = tuple(args.syrinx_width_ratio)
        config.max_parallel_jobs = args.max_parallel_jobs
        config.memory_limit_mb = args.memory_limit
        
        # HIPAA SECURITY CONTROL: Require environment variable for anonymization bypass
        if args.no_anonymize:
            # Check for explicit authorization via environment variable
            auth_token = os.environ.get('MRI_TOOL_ALLOW_PHI')
            if auth_token != 'AUTHORIZED_RESEARCH_USE':
                logger.error("SECURITY: --no-anonymize requires MRI_TOOL_ALLOW_PHI=AUTHORIZED_RESEARCH_USE")
                logger.error("This flag disables patient data anonymization and exposes PHI.")
                logger.error("Only use in authorized research environments with IRB approval.")
                return 1
            else:
                logger.warning("SECURITY: Patient data anonymization DISABLED by environment authorization")
                logger.warning("Ensure compliance with HIPAA and institutional policies")
                
        config.anonymize_patient_data = not args.no_anonymize
        config.report_settings.pubmed_api_enabled = not args.no_citations
        
        # Validate configuration
        # Warn if thresholds are set outside literature-supported ranges.
        warnings = config.validate()
        for warning in warnings:
            logger.warning(f"Configuration warning: {warning}")
            
        logger.info("Starting MRI Enhancement Tool analysis")
        logger.info(f"Input: {args.input}")
        logger.info(f"Output: {args.output}")
        logger.info(f"Configuration: {args.config or 'default'}")
        
        if args.report_only:
            logger.info("Report-only mode: skipping DICOM scanning and analysis")
            # Would load pre-existing analysis results
            logger.error("Report-only mode not yet implemented")
            return 1
            
        # Step 1: Scan DICOM files
        logger.info("Step 1: Scanning DICOM files...")
        scanner = DICOMScanner(config)
        studies = scanner.scan_folder(args.input)
        
        if not studies:
            logger.error("No valid spine MRI studies found in input folder")
            return 1
            
        logger.info(f"Found {len(studies)} spine studies")
        
        # Step 2: Analyze each series
        logger.info("Step 2: Analyzing MRI series...")
        analyzer = SpineMRIAnalyzer(config)
        analysis_results = {}
        
        for study in studies:
            for series in study.get_spine_series():
                try:
                    logger.info(f"Analyzing series: {series.series_description}")
                    
                    # Load volume
                    volume, metadata = scanner.load_series_as_volume(series)
                    
                    # Analyze series
                    result = analyzer.analyze_series(volume, metadata, series)
                    analysis_results[series.series_uid] = result
                    
                    logger.info(f"Found {len(result.findings)} findings, "
                              f"overall assessment: {result.overall_assessment.value}")
                    
                except Exception as e:
                    logger.error(f"Error analyzing series {series.series_uid}: {e}")
                    continue
                    
        if not analysis_results:
            logger.error("No series could be analyzed successfully")
            return 1
            
        # Step 3: Generate annotations
        annotated_images = {}
        if not args.no_annotations:
            logger.info("Step 3: Generating annotated images...")
            annotation_generator = AnnotationGenerator(config)
            
            for study in studies:
                for series in study.get_spine_series():
                    if series.series_uid in analysis_results:
                        try:
                            # Use first file as representative
                            dicom_file = series.file_paths[0] if series.file_paths else None
                            if dicom_file:
                                annotations_dir = os.path.join(args.output, "annotations", series.series_uid)
                                generated_files = annotation_generator.generate_annotations(
                                    dicom_file, 
                                    analysis_results[series.series_uid],
                                    annotations_dir
                                )
                                annotated_images[series.series_uid] = generated_files
                                
                        except Exception as e:
                            logger.error(f"Error generating annotations for series {series.series_uid}: {e}")
                            annotated_images[series.series_uid] = []
                            
        # Step 4: Generate report
        logger.info("Step 4: Generating comprehensive report...")
        report_generator = ReportGenerator(config)
        
        report_path = os.path.join(args.output, "spine_mri_analysis_report.pdf")
        try:
            generated_report = report_generator.generate_report(
                studies, analysis_results, annotated_images, report_path
            )
            logger.info(f"Report generated: {generated_report}")
            
        except Exception as e:
            logger.error(f"Error generating report: {e}")
            return 1
            
        # Summary
        total_findings = sum(len(result.findings) for result in analysis_results.values())
        traumatic_findings = sum(
            sum(1 for finding in result.findings if finding.finding_type.name == 'TRAUMATIC')
            for result in analysis_results.values()
        )
        
        logger.info("Analysis completed successfully!")
        logger.info(f"Summary:")
        logger.info(f"  Studies processed: {len(studies)}")
        logger.info(f"  Series analyzed: {len(analysis_results)}")
        logger.info(f"  Total findings: {total_findings}")
        logger.info(f"  Traumatic findings: {traumatic_findings}")
        logger.info(f"  Report: {report_path}")
        
        if annotated_images:
            total_images = sum(len(images) for images in annotated_images.values())
            logger.info(f"  Annotated images: {total_images}")
            
        return 0
        
    except KeyboardInterrupt:
        logger.info("Analysis interrupted by user")
        return 1
        
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        if args.log_level == 'DEBUG':
            import traceback
            logger.debug(traceback.format_exc())
        return 1


if __name__ == '__main__':
    sys.exit(main())