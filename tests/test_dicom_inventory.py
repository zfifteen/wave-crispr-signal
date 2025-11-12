#!/usr/bin/env python3
"""
Test script to validate DICOM inventory analysis results
"""

import sys
from pathlib import Path
import json

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

def test_dicom_analysis():
    """Test that our DICOM analysis results are accurate."""
    
    # Load the analysis results
    analysis_file = project_root / 'DICOM_INVENTORY_ANALYSIS.json'
    if not analysis_file.exists():
        print("❌ Analysis file not found. Run the analysis first.")
        return False
    
    with open(analysis_file, 'r') as f:
        inventory = json.load(f)
    
    print("Testing DICOM Inventory Analysis Results...")
    
    # Test 1: Check file counts
    dicom_root = project_root / 'data' / 'DICOM'
    actual_dcm_files = len(list(dicom_root.rglob('*.dcm')))
    actual_png_files = len(list(dicom_root.rglob('*.png')))
    actual_metadata_files = len(list(dicom_root.rglob('*_metadata.txt')))
    
    expected_dcm = inventory['summary_stats']['total_dicom_files']
    expected_png = inventory['summary_stats']['total_png_files'] 
    expected_metadata = inventory['summary_stats']['total_metadata_files']
    
    print(f"\n📊 File Count Validation:")
    print(f"  DICOM files: {actual_dcm_files} (expected: {expected_dcm}) {'✅' if actual_dcm_files == expected_dcm else '❌'}")
    print(f"  PNG files: {actual_png_files} (expected: {expected_png}) {'✅' if actual_png_files == expected_png else '❌'}")
    print(f"  Metadata files: {actual_metadata_files} (expected: {expected_metadata}) {'✅' if actual_metadata_files == expected_metadata else '❌'}")
    
    # Test 2: Check that modalities are correctly identified
    expected_modalities = set(['CT', 'MR', 'DX'])
    actual_modalities = set(inventory['summary_stats']['modalities_present'])
    
    print(f"\n🔬 Modality Validation:")
    print(f"  Expected modalities: {expected_modalities}")
    print(f"  Found modalities: {actual_modalities}")
    print(f"  Modalities match: {'✅' if expected_modalities == actual_modalities else '❌'}")
    
    # Test 3: Check example files exist
    # Dynamically extract example DICOM file paths from the inventory
    example_files = []
    # Try to get one example file from up to 3 different series
    if 'series' in inventory:
        for series in list(inventory['series'].values())[:3]:
            if 'files' in series and series['files']:
                # Prepend the relative path to the file
                example_files.append(series['files'][0])
    if not example_files:
        print("⚠️  No example files found in inventory for validation.")
    
    print(f"\n📁 Example Files Validation:")
    all_examples_exist = True
    for example_file in example_files:
        file_path = project_root / example_file
        exists = file_path.exists()
        print(f"  {example_file}: {'✅' if exists else '❌'}")
        if not exists:
            all_examples_exist = False
    
    # Test 4: Check directory structure
    expected_dirs = ['CT_HEAD_W_O_CONT', 'MRI__CERVICAL_SPINE_W_O_CONT', 
                     'MRI__THORACIC_SPINE_W_O_CONT', 'SPINE_THORACIC_2_VIEWS']
    
    print(f"\n📂 Directory Structure Validation:")
    all_dirs_exist = True
    for dir_name in expected_dirs:
        dir_path = dicom_root / dir_name
        exists = dir_path.exists() and dir_path.is_dir()
        print(f"  {dir_name}: {'✅' if exists else '❌'}")
        if not exists:
            all_dirs_exist = False
    
    # Test 5: Check manufacturers
    manufacturers = inventory['summary_stats']['manufacturers']
    expected_manufacturers = {'Philips', 'Siemens Healthineers', 'GE Healthcare'}
    actual_manufacturers = set(manufacturers.keys())
    
    print(f"\n🏭 Manufacturer Validation:")
    print(f"  Expected: {expected_manufacturers}")
    print(f"  Found: {actual_manufacturers}")
    print(f"  Contains expected manufacturers: {'✅' if expected_manufacturers.issubset(actual_manufacturers) else '❌'}")
    
    # Overall assessment
    all_tests_passed = (
        actual_dcm_files == expected_dcm and
        actual_png_files == expected_png and
        actual_metadata_files == expected_metadata and
        expected_modalities == actual_modalities and
        all_examples_exist and
        all_dirs_exist and
        expected_manufacturers.issubset(actual_manufacturers)
    )
    
    print(f"\n🎯 Overall Test Result: {'✅ ALL TESTS PASSED' if all_tests_passed else '❌ SOME TESTS FAILED'}")

    # Additional statistics
    print(f"\n📈 Summary Statistics:")
    print(f"  Total DICOM dataset size: {inventory['summary_stats']['total_size_gb']:.2f} GB")
    print(f"  Number of series: {inventory['summary_stats']['total_series']}")
    print(f"  Patient ID: {inventory['patient_id']}")

    assert all_tests_passed, "DICOM analysis validation failed"

if __name__ == '__main__':
    try:
        test_dicom_analysis()
        print("✅ All tests passed")
        sys.exit(0)
    except AssertionError as e:
        print(f"❌ Test failed: {e}")
        sys.exit(1)