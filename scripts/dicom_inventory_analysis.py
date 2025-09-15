#!/usr/bin/env python3
"""
DICOM Inventory Analysis Script for Patient A
===========================================

This script performs a comprehensive analysis of DICOM data for Patient A, 
examining metadata, technical parameters, and file organization across 
all modalities (CT, MRI, X-Ray).

IMPORTANT: This script analyzes ALL DICOM files in the dataset without 
sampling or skipping any files. Medical imaging data cannot be discarded 
or sampled as each file may contain unique technical parameters and 
metadata variations that are medically significant.

Usage:
    python scripts/dicom_inventory_analysis.py

Output: 
    - Complete inventory of all DICOM files
    - Technical parameters analysis by modality
    - Navigation guide for other LLMs
"""

import os
import sys
from pathlib import Path
from collections import defaultdict, Counter
import json
from datetime import datetime
import re

# Add project root to path for imports
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

import pydicom
import pandas as pd


class DICOMInventoryAnalyzer:
    """Comprehensive analyzer for DICOM data inventory."""
    
    def __init__(self, dicom_root: Path):
        self.dicom_root = dicom_root
        self.inventory = {
            'patient_id': 'Patient A',
            'analysis_timestamp': datetime.now().isoformat(),
            'modalities': {},
            'summary_stats': {},
            'technical_parameters': {},
            'file_organization': {},
            'navigation_guide': {}
        }
        
    def analyze_dicom_file(self, dcm_path: Path) -> dict:
        """Extract key metadata from a single DICOM file."""
        try:
            ds = pydicom.dcmread(str(dcm_path), stop_before_pixels=True)
            
            metadata = {
                'file_path': str(dcm_path.relative_to(self.dicom_root)),
                'file_size_mb': dcm_path.stat().st_size / (1024 * 1024),
                'sop_class_uid': str(getattr(ds, 'SOPClassUID', 'Unknown')),
                'modality': str(getattr(ds, 'Modality', 'Unknown')),
                'manufacturer': str(getattr(ds, 'Manufacturer', 'Unknown')),
                'manufacturer_model': str(getattr(ds, 'ManufacturerModelName', 'Unknown')),
                'study_date': str(getattr(ds, 'StudyDate', 'Unknown')),
                'series_date': str(getattr(ds, 'SeriesDate', 'Unknown')),
                'study_description': str(getattr(ds, 'StudyDescription', 'Unknown')),
                'series_description': str(getattr(ds, 'SeriesDescription', 'Unknown')),
                'body_part_examined': str(getattr(ds, 'BodyPartExamined', 'Unknown')),
                'image_type': getattr(ds, 'ImageType', []),
                'rows': getattr(ds, 'Rows', None),
                'columns': getattr(ds, 'Columns', None),
                'bits_allocated': getattr(ds, 'BitsAllocated', None),
                'bits_stored': getattr(ds, 'BitsStored', None),
                'pixel_representation': getattr(ds, 'PixelRepresentation', None),
                'photometric_interpretation': str(getattr(ds, 'PhotometricInterpretation', 'Unknown')),
                'pixel_spacing': getattr(ds, 'PixelSpacing', None),
                'slice_thickness': getattr(ds, 'SliceThickness', None),
                'spacing_between_slices': getattr(ds, 'SpacingBetweenSlices', None),
                'window_center': getattr(ds, 'WindowCenter', None),
                'window_width': getattr(ds, 'WindowWidth', None),
                'rescale_slope': getattr(ds, 'RescaleSlope', None),
                'rescale_intercept': getattr(ds, 'RescaleIntercept', None),
                'kvp': getattr(ds, 'KVP', None),
                'exposure_time': getattr(ds, 'ExposureTime', None),
                'tube_current': getattr(ds, 'XRayTubeCurrent', None),
                'slice_location': getattr(ds, 'SliceLocation', None),
                'image_position_patient': getattr(ds, 'ImagePositionPatient', None),
                'image_orientation_patient': getattr(ds, 'ImageOrientationPatient', None),
                'repetition_time': getattr(ds, 'RepetitionTime', None),
                'echo_time': getattr(ds, 'EchoTime', None),
                'flip_angle': getattr(ds, 'FlipAngle', None),
                'magnetic_field_strength': getattr(ds, 'MagneticFieldStrength', None),
                'sequence_name': str(getattr(ds, 'SequenceName', 'Unknown')),
                'protocol_name': str(getattr(ds, 'ProtocolName', 'Unknown')),
                'contrast_bolus_agent': str(getattr(ds, 'ContrastBolusAgent', 'None')),
            }
            return metadata
            
        except Exception as e:
            return {
                'file_path': str(dcm_path.relative_to(self.dicom_root)),
                'error': str(e),
                'file_size_mb': dcm_path.stat().st_size / (1024 * 1024)
            }
    
    def analyze_directory_structure(self):
        """Analyze the overall directory structure."""
        structure = {}
        
        for modality_dir in self.dicom_root.iterdir():
            if modality_dir.is_dir() and modality_dir.name != '__pycache__':
                modality_name = modality_dir.name
                structure[modality_name] = {
                    'path': str(modality_dir.relative_to(self.dicom_root)),
                    'series': {},
                    'dicom_files_count': 0,
                    'png_files_count': 0,
                    'metadata_files_count': 0,
                    'total_size_mb': 0
                }
                
                # Look for DICOM subdirectory
                dicom_subdir = modality_dir / 'DICOM'
                if dicom_subdir.exists():
                    for series_dir in dicom_subdir.iterdir():
                        if series_dir.is_dir():
                            series_name = series_dir.name
                            dcm_files = list(series_dir.glob('*.dcm'))
                            png_files = list(series_dir.glob('*.png'))
                            metadata_files = list(series_dir.glob('*_metadata.txt'))
                            
                            structure[modality_name]['series'][series_name] = {
                                'path': str(series_dir.relative_to(self.dicom_root)),
                                'dcm_files': len(dcm_files),
                                'png_files': len(png_files),
                                'metadata_files': len(metadata_files),
                                'example_dcm_file': str(dcm_files[0].name) if dcm_files else None,
                                'size_mb': sum(f.stat().st_size for f in series_dir.iterdir() if f.is_file()) / (1024 * 1024)
                            }
                            
                            structure[modality_name]['dicom_files_count'] += len(dcm_files)
                            structure[modality_name]['png_files_count'] += len(png_files)
                            structure[modality_name]['metadata_files_count'] += len(metadata_files)
                            structure[modality_name]['total_size_mb'] += structure[modality_name]['series'][series_name]['size_mb']
        
        return structure
    
    def extract_technical_parameters_by_modality(self, dicom_files_data):
        """Extract and categorize technical parameters by imaging modality."""
        modality_params = defaultdict(lambda: defaultdict(list))
        
        for file_data in dicom_files_data:
            if 'error' in file_data:
                continue
                
            modality = file_data.get('modality', 'Unknown')
            
            # Common parameters
            for param in ['rows', 'columns', 'bits_allocated', 'bits_stored', 
                         'pixel_spacing', 'slice_thickness']:
                value = file_data.get(param)
                if value is not None:
                    modality_params[modality][param].append(value)
            
            # CT-specific parameters
            if modality == 'CT':
                for param in ['kvp', 'tube_current', 'exposure_time', 'window_center', 'window_width']:
                    value = file_data.get(param)
                    if value is not None:
                        modality_params[modality][param].append(value)
            
            # MRI-specific parameters
            elif modality == 'MR':
                for param in ['repetition_time', 'echo_time', 'flip_angle', 
                             'magnetic_field_strength', 'sequence_name']:
                    value = file_data.get(param)
                    if value is not None:
                        modality_params[modality][param].append(value)
            
            # X-Ray/DX-specific parameters
            elif modality in ['DX', 'CR']:
                for param in ['kvp', 'exposure_time', 'window_center', 'window_width']:
                    value = file_data.get(param)
                    if value is not None:
                        modality_params[modality][param].append(value)
        
        # Calculate statistics for each parameter
        param_stats = {}
        for modality, params in modality_params.items():
            param_stats[modality] = {}
            for param, values in params.items():
                if values:
                    if isinstance(values[0], (int, float)):
                        param_stats[modality][param] = {
                            'min': min(values),
                            'max': max(values),
                            'unique_values': len(set(values)),
                            'total_count': len(values)
                        }
                    else:
                        param_stats[modality][param] = {
                            'unique_values': list(set(str(v) for v in values)),
                            'total_count': len(values)
                        }
        
        return param_stats
    
    def generate_navigation_guide(self, directory_structure, dicom_files_data):
        """Generate a comprehensive navigation guide for other LLMs."""
        guide = {
            'overview': {
                'patient': 'Patient A',
                'total_modalities': len(directory_structure),
                'modality_list': list(directory_structure.keys()),
                'total_series': sum(len(mod['series']) for mod in directory_structure.values()),
                'total_dicom_files': sum(mod['dicom_files_count'] for mod in directory_structure.values())
            },
            'modality_details': {},
            'anatomy_mapping': {},
            'finding_files_by': {
                'modality': {},
                'anatomy': {},
                'test_type': {},
                'parameters': {}
            },
            'file_examples': {}
        }
        
        # Extract anatomy information from file data
        anatomy_files = defaultdict(list)
        for file_data in dicom_files_data:
            if 'error' not in file_data:
                body_part = file_data.get('body_part_examined', 'Unknown')
                series_desc = file_data.get('series_description', '')
                study_desc = file_data.get('study_description', '')
                
                # Map anatomy
                anatomy_key = f"{body_part}_{series_desc}_{study_desc}".replace(' ', '_').upper()
                anatomy_files[anatomy_key].append(file_data['file_path'])
        
        # Organize by modality
        for modality_name, modality_data in directory_structure.items():
            guide['modality_details'][modality_name] = {
                'description': self._get_modality_description(modality_name),
                'series_count': len(modality_data['series']),
                'file_count': modality_data['dicom_files_count'],
                'series_list': list(modality_data['series'].keys()),
                'example_path': f"data/DICOM/{modality_name}/DICOM/",
                'typical_use_cases': self._get_modality_use_cases(modality_name)
            }
            
            # Add file finding guidance
            guide['finding_files_by']['modality'][modality_name] = {
                'base_path': f"data/DICOM/{modality_name}/DICOM/",
                'series_directories': list(modality_data['series'].keys()),
                'example_file': self._get_example_file(modality_data)
            }
        
        return guide
    
    def _get_modality_description(self, modality_name):
        """Get description for each modality."""
        descriptions = {
            'CT_HEAD_W_O_CONT': 'CT scan of head without contrast - computed tomography imaging for brain and skull evaluation',
            'MRI__CERVICAL_SPINE_W_O_CONT': 'MRI of cervical spine without contrast - magnetic resonance imaging of neck vertebrae',
            'MRI__THORACIC_SPINE_W_O_CONT': 'MRI of thoracic spine without contrast - magnetic resonance imaging of upper back vertebrae',
            'SPINE_THORACIC_2_VIEWS': 'X-ray of thoracic spine in 2 views - conventional radiography of upper back'
        }
        return descriptions.get(modality_name, f'Medical imaging study: {modality_name}')
    
    def _get_modality_use_cases(self, modality_name):
        """Get typical use cases for each modality."""
        use_cases = {
            'CT_HEAD_W_O_CONT': ['Brain injury assessment', 'Stroke evaluation', 'Headache investigation', 'Skull fracture detection'],
            'MRI__CERVICAL_SPINE_W_O_CONT': ['Neck pain evaluation', 'Disc herniation assessment', 'Spinal cord evaluation', 'Nerve compression'],
            'MRI__THORACIC_SPINE_W_O_CONT': ['Back pain investigation', 'Spinal deformity assessment', 'Disc disease evaluation', 'Tumor screening'],
            'SPINE_THORACIC_2_VIEWS': ['Spine alignment check', 'Fracture detection', 'Scoliosis screening', 'General spine evaluation']
        }
        return use_cases.get(modality_name, ['General medical imaging'])
    
    def _get_example_file(self, modality_data):
        """Get an example file from the modality."""
        for series_data in modality_data['series'].values():
            if series_data['example_dcm_file']:
                return series_data['example_dcm_file']
        return None
    
    def run_full_analysis(self):
        """Run the complete DICOM inventory analysis."""
        print("Starting DICOM Inventory Analysis for Patient A...")
        print(f"Analyzing directory: {self.dicom_root}")
        
        # 1. Analyze directory structure
        print("\n1. Analyzing directory structure...")
        directory_structure = self.analyze_directory_structure()
        self.inventory['file_organization'] = directory_structure
        
        # 2. Collect ALL DICOM files for comprehensive analysis
        print("\n2. Collecting ALL DICOM file metadata for comprehensive analysis...")
        dicom_files_data = []
        total_files = 0
        
        for modality_dir in self.dicom_root.iterdir():
            if modality_dir.is_dir() and modality_dir.name != '__pycache__':
                dicom_subdir = modality_dir / 'DICOM'
                if dicom_subdir.exists():
                    for series_dir in dicom_subdir.iterdir():
                        if series_dir.is_dir():
                            dcm_files = list(series_dir.glob('*.dcm'))
                            total_files += len(dcm_files)
                            
                            # Analyze ALL DICOM files - medical data cannot be sampled or discarded
                            for dcm_file in dcm_files:
                                file_data = self.analyze_dicom_file(dcm_file)
                                file_data['series_directory'] = series_dir.name
                                file_data['modality_directory'] = modality_dir.name
                                dicom_files_data.append(file_data)
                                
                                if len(dicom_files_data) % 100 == 0:
                                    print(f"  Processed {len(dicom_files_data)} files...")
        
        print(f"  Completed comprehensive analysis of ALL {len(dicom_files_data)} DICOM files (total: {total_files})")
        
        # 3. Extract technical parameters
        print("\n3. Extracting technical parameters by modality...")
        technical_params = self.extract_technical_parameters_by_modality(dicom_files_data)
        self.inventory['technical_parameters'] = technical_params
        
        # 4. Generate summary statistics
        print("\n4. Generating summary statistics...")
        self.inventory['summary_stats'] = self._generate_summary_stats(dicom_files_data, directory_structure)
        
        # 5. Create navigation guide
        print("\n5. Creating navigation guide...")
        navigation_guide = self.generate_navigation_guide(directory_structure, dicom_files_data)
        self.inventory['navigation_guide'] = navigation_guide
        
        print("\nAnalysis complete!")
        return self.inventory
    
    def _generate_summary_stats(self, dicom_files_data, directory_structure):
        """Generate overall summary statistics."""
        stats = {
            'total_modalities': len(directory_structure),
            'total_series': sum(len(mod['series']) for mod in directory_structure.values()),
            'total_dicom_files': sum(mod['dicom_files_count'] for mod in directory_structure.values()),
            'total_png_files': sum(mod['png_files_count'] for mod in directory_structure.values()),
            'total_metadata_files': sum(mod['metadata_files_count'] for mod in directory_structure.values()),
            'total_size_gb': sum(mod['total_size_mb'] for mod in directory_structure.values()) / 1024,
            'modalities_present': [],
            'manufacturers': Counter(),
            'study_dates': Counter(),
            'body_parts_examined': Counter()
        }
        
        for file_data in dicom_files_data:
            if 'error' not in file_data:
                modality = file_data.get('modality', 'Unknown')
                if modality not in stats['modalities_present']:
                    stats['modalities_present'].append(modality)
                
                stats['manufacturers'][file_data.get('manufacturer', 'Unknown')] += 1
                stats['study_dates'][file_data.get('study_date', 'Unknown')] += 1
                stats['body_parts_examined'][file_data.get('body_part_examined', 'Unknown')] += 1
        
        # Convert Counter objects to regular dicts for JSON serialization
        stats['manufacturers'] = dict(stats['manufacturers'])
        stats['study_dates'] = dict(stats['study_dates'])
        stats['body_parts_examined'] = dict(stats['body_parts_examined'])
        
        return stats


def main():
    """Main function to run the DICOM inventory analysis."""
    # Set up paths
    project_root = Path(__file__).parent.parent
    dicom_root = project_root / 'data' / 'DICOM'
    
    if not dicom_root.exists():
        print(f"Error: DICOM directory not found at {dicom_root}")
        return 1
    
    # Run analysis
    analyzer = DICOMInventoryAnalyzer(dicom_root)
    inventory = analyzer.run_full_analysis()
    
    # Save results
    output_file = project_root / 'DICOM_INVENTORY_ANALYSIS.json'
    with open(output_file, 'w') as f:
        json.dump(inventory, f, indent=2, default=str)
    
    print(f"\nResults saved to: {output_file}")
    
    # Generate summary report
    report_file = project_root / 'DICOM_ANALYSIS_REPORT.md'
    generate_markdown_report(inventory, report_file)
    print(f"Detailed report generated: {report_file}")
    
    return 0


def generate_markdown_report(inventory, output_file):
    """Generate a detailed markdown report from the inventory data."""
    with open(output_file, 'w') as f:
        f.write(f"""# DICOM Data Analysis Report for Patient A

**Analysis Timestamp:** {inventory['analysis_timestamp']}  
**Patient ID:** {inventory['patient_id']}

## Executive Summary

This comprehensive analysis examines the DICOM medical imaging data for Patient A, providing detailed technical specifications, file organization, and navigation guidance for automated processing and analysis.

### Key Statistics
- **Total Modalities:** {inventory['summary_stats']['total_modalities']}
- **Total Series:** {inventory['summary_stats']['total_series']}
- **Total DICOM Files:** {inventory['summary_stats']['total_dicom_files']:,}
- **Total PNG Files:** {inventory['summary_stats']['total_png_files']:,}
- **Total Metadata Files:** {inventory['summary_stats']['total_metadata_files']:,}
- **Total Dataset Size:** {inventory['summary_stats']['total_size_gb']:.2f} GB

### Imaging Modalities Present
{', '.join(inventory['summary_stats']['modalities_present'])}

## Directory Structure and Organization

### Overview
```
data/DICOM/
├── CT_HEAD_W_O_CONT/           # CT scan of head without contrast
├── MRI__CERVICAL_SPINE_W_O_CONT/   # MRI cervical spine without contrast  
├── MRI__THORACIC_SPINE_W_O_CONT/   # MRI thoracic spine without contrast
└── SPINE_THORACIC_2_VIEWS/          # X-ray thoracic spine 2 views
```

""")
        
        # Detailed modality breakdown
        for modality_name, modality_data in inventory['file_organization'].items():
            f.write(f"""### {modality_name}
**Description:** {inventory['navigation_guide']['modality_details'][modality_name]['description']}

- **Location:** `data/DICOM/{modality_name}/DICOM/`
- **Series Count:** {len(modality_data['series'])}
- **DICOM Files:** {modality_data['dicom_files_count']:,}
- **PNG Images:** {modality_data['png_files_count']:,}
- **Metadata Files:** {modality_data['metadata_files_count']:,}
- **Total Size:** {modality_data['total_size_mb']:.1f} MB

#### Series Breakdown:
""")
            for series_name, series_data in modality_data['series'].items():
                example_file = series_data['example_dcm_file'] or 'N/A'
                f.write(f"""- **{series_name}**
  - Path: `{series_data['path']}`
  - DICOM files: {series_data['dcm_files']}
  - PNG files: {series_data['png_files']}
  - Example file: `{example_file}`
  - Size: {series_data['size_mb']:.1f} MB

""")
        
        # Technical parameters section
        f.write(f"""## Technical Parameters by Modality

""")
        
        for modality, params in inventory['technical_parameters'].items():
            f.write(f"""### {modality} Parameters

""")
            for param_name, param_data in params.items():
                if isinstance(param_data, dict) and 'min' in param_data:
                    f.write(f"""- **{param_name}**: Range {param_data['min']} - {param_data['max']} ({param_data['unique_values']} unique values)
""")
                elif isinstance(param_data, dict) and 'unique_values' in param_data:
                    unique_vals = param_data['unique_values']
                    if len(unique_vals) <= 5:
                        f.write(f"""- **{param_name}**: {', '.join(map(str, unique_vals))}
""")
                    else:
                        f.write(f"""- **{param_name}**: {param_data['total_count']} values, {len(unique_vals)} unique
""")
            f.write("\n")
        
        # Navigation guide
        f.write(f"""## Navigation Guide for LLMs

### Finding Files by Modality

""")
        
        for modality_name, guide_data in inventory['navigation_guide']['finding_files_by']['modality'].items():
            f.write(f"""#### {modality_name}
- **Base Path:** `{guide_data['base_path']}`
- **Series Directories:** {', '.join(guide_data['series_directories'])}
- **Example File:** `{guide_data.get('example_file', 'N/A')}`
- **Use Cases:** {', '.join(inventory['navigation_guide']['modality_details'][modality_name]['typical_use_cases'])}

""")
        
        # Manufacturers and equipment
        f.write(f"""## Equipment and Manufacturers

""")
        for manufacturer, count in inventory['summary_stats']['manufacturers'].items():
            f.write(f"""- **{manufacturer}**: {count} files
""")
        
        # Study dates
        f.write(f"""

## Study Dates

""")
        for date, count in sorted(inventory['summary_stats']['study_dates'].items()):
            f.write(f"""- **{date}**: {count} files
""")
        
        # Body parts examined
        f.write(f"""

## Anatomical Regions Examined

""")
        for body_part, count in inventory['summary_stats']['body_parts_examined'].items():
            f.write(f"""- **{body_part}**: {count} files
""")
        
        f.write(f"""

## File Examples for Testing

### CT Head Example
```bash
# Example CT file location
data/DICOM/CT_HEAD_W_O_CONT/DICOM/SERIES_202/209585291.dcm
data/DICOM/CT_HEAD_W_O_CONT/DICOM/SERIES_202/209585291_metadata.txt
data/DICOM/CT_HEAD_W_O_CONT/DICOM/SERIES_202/209585291_s0000_view8.png
```

### MRI Cervical Spine Example  
```bash
# Example MRI cervical file location
data/DICOM/MRI__CERVICAL_SPINE_W_O_CONT/DICOM/SERIES_3/207815325.dcm
data/DICOM/MRI__CERVICAL_SPINE_W_O_CONT/DICOM/SERIES_3/207815325_metadata.txt
data/DICOM/MRI__CERVICAL_SPINE_W_O_CONT/DICOM/SERIES_3/207815325_s0000_view8.png
```

### X-Ray Thoracic Spine Example
```bash
# Example X-ray file location  
data/DICOM/SPINE_THORACIC_2_VIEWS/DICOM/SERIES_357/204792524.dcm
data/DICOM/SPINE_THORACIC_2_VIEWS/DICOM/SERIES_357/204792524_metadata.txt
data/DICOM/SPINE_THORACIC_2_VIEWS/DICOM/SERIES_357/204792524_s0000_view8.png
```

## Quick Reference Commands

### Count Files by Type
```bash
# Count DICOM files
find data/DICOM -name "*.dcm" | wc -l

# Count PNG images  
find data/DICOM -name "*.png" | wc -l

# Count metadata files
find data/DICOM -name "*_metadata.txt" | wc -l
```

### Find Files by Modality
```bash
# CT files
find data/DICOM/CT_HEAD_W_O_CONT -name "*.dcm"

# MRI cervical files
find data/DICOM/MRI__CERVICAL_SPINE_W_O_CONT -name "*.dcm"

# MRI thoracic files  
find data/DICOM/MRI__THORACIC_SPINE_W_O_CONT -name "*.dcm"

# X-ray files
find data/DICOM/SPINE_THORACIC_2_VIEWS -name "*.dcm"
```

---
*Generated by DICOM Inventory Analysis Script - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}*
""")


if __name__ == '__main__':
    sys.exit(main())