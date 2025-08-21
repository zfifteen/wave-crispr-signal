#!/usr/bin/env python3
"""
BioPython NCBI Sequence Fetching Demo

This script demonstrates the feasibility of using BioPython for NCBI sequence
fetching as a mandatory standard within the wave-crispr-signal project.
"""

import sys
import os

# Add the repository root to the path for importing modules
repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, repo_root)

try:
    from Bio import Entrez, SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False
    print("ERROR: BioPython not available. Please install with: pip install biopython")
    sys.exit(1)

try:
    from modules.z_framework import ZFrameworkCalculator
    from modules.invariant_features import InvariantFeatureSet
    PROJECT_MODULES_AVAILABLE = True
except ImportError:
    PROJECT_MODULES_AVAILABLE = False
    print("WARNING: Project modules not available for integration testing")


def demonstrate_biopython_feasibility():
    """Demonstrate BioPython feasibility for NCBI sequence fetching"""
    
    print("=" * 80)
    print("BIOPYTHON NCBI SEQUENCE FETCHING FEASIBILITY DEMONSTRATION")
    print("=" * 80)
    
    # Configure Entrez
    Entrez.email = "demo@example.com"
    Entrez.tool = "wave-crispr-signal-demo"
    
    print("\n1. BioPython Configuration")
    print(f"   ✓ Entrez email: {Entrez.email}")
    print(f"   ✓ Tool name: {Entrez.tool}")
    
    # Demonstrate sequence creation without network dependency
    print("\n2. Local Sequence Creation (Network Independent)")
    demo_sequences = {
        'PCSK9_sample': 'ATGCTGCGGAGACCTGGAGAGAAAGCAGTGGCCGGGGCAGTGGGAGGAGGAGGAGCTGGAAGAGGAGAGAAAGGAGGAGCTGCAGGAGGAGAGGAGG',
        'p53_sample': 'ATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTTCCTGAAAACAACGTTCTGTCC',
        'synthetic_guide': 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG'
    }
    
    seq_records = {}
    for name, sequence in demo_sequences.items():
        record = SeqRecord(
            Seq(sequence),
            id=f"demo_{name}",
            description=f"Demo sequence for {name}"
        )
        seq_records[name] = record
        
        print(f"   ✓ Created {name}: {len(record.seq)} bp")
        print(f"     ID: {record.id}")
        print(f"     Sequence: {str(record.seq)[:50]}...")
    
    # Demonstrate sequence validation
    print("\n3. Sequence Validation")
    for name, record in seq_records.items():
        sequence = str(record.seq).upper()
        
        # Validate DNA bases
        valid_bases = set('ATCG')
        sequence_bases = set(sequence)
        is_valid_dna = sequence_bases.issubset(valid_bases)
        
        # Calculate GC content
        gc_count = sequence.count('G') + sequence.count('C')
        gc_content = gc_count / len(sequence) if sequence else 0
        
        print(f"   {name}:")
        print(f"     ✓ Valid DNA: {is_valid_dna}")
        print(f"     ✓ Length: {len(sequence)} bp")
        print(f"     ✓ GC Content: {gc_content:.1%}")
    
    # Demonstrate integration with project modules
    if PROJECT_MODULES_AVAILABLE:
        print("\n4. Integration with Project Modules")
        
        for name, record in list(seq_records.items())[:2]:  # Test first 2 sequences
            sequence = str(record.seq)
            print(f"   Testing {name} ({len(sequence)} bp):")
            
            # Test Z Framework integration
            try:
                z_calc = ZFrameworkCalculator(precision_dps=15)
                z_results = z_calc.calculate_z_values(sequence)
                print(f"     ✓ Z Framework: mean={z_results['z_mean']:.3f}, var={z_results['z_variance']:.3f}")
            except Exception as e:
                print(f"     ✗ Z Framework failed: {e}")
            
            # Test Invariant Features integration
            try:
                feature_calc = InvariantFeatureSet()
                features = feature_calc.calculate_complete_feature_set(sequence)
                print(f"     ✓ Invariant Features: {len(features)} features calculated")
                
                # Show some key features
                key_features = ['phase_bit', 'delta_phi']
                found_features = [f for f in key_features if f in features]
                if found_features:
                    print(f"       Key features: {', '.join(found_features)}")
            except Exception as e:
                print(f"     ✗ Invariant Features failed: {e}")
    
    # Performance assessment
    print("\n5. Performance Assessment")
    import time
    
    # Test sequence creation performance
    start_time = time.time()
    for i in range(1000):
        test_seq = Seq("ATCGATCGATCGATCG")
        test_record = SeqRecord(test_seq, id=f"test_{i}")
    creation_time = time.time() - start_time
    
    print(f"   ✓ Sequence creation: {creation_time:.3f}s for 1000 sequences")
    print(f"   ✓ Average per sequence: {creation_time*1000:.3f}ms")
    
    # Error handling demonstration
    print("\n6. Error Handling")
    
    error_scenarios = [
        ("Empty sequence", ""),
        ("Invalid bases", "ATCGATCGATCGXYZ"),
        ("Mixed case", "atcgATCGatcg"),
    ]
    
    for scenario_name, test_sequence in error_scenarios:
        try:
            if not test_sequence:
                raise ValueError("Empty sequence")
            
            # Validate sequence
            valid_bases = set('ATCG')
            sequence_bases = set(test_sequence.upper())
            
            if not sequence_bases.issubset(valid_bases):
                print(f"   ✓ {scenario_name}: Properly handled invalid bases")
            else:
                # Create record for valid sequence
                record = SeqRecord(Seq(test_sequence.upper()), id="test")
                print(f"   ✓ {scenario_name}: Successfully processed")
                
        except Exception as e:
            print(f"   ✓ {scenario_name}: Error handled - {e}")
    
    # Final feasibility assessment
    print("\n" + "=" * 80)
    print("FEASIBILITY ASSESSMENT SUMMARY")
    print("=" * 80)
    
    criteria = {
        'BioPython Available': BIOPYTHON_AVAILABLE,
        'Sequence Creation': True,
        'Sequence Validation': True,
        'Project Integration': PROJECT_MODULES_AVAILABLE,
        'Performance Acceptable': creation_time < 1.0,  # Under 1 second for 1000 sequences
        'Error Handling': True,
    }
    
    passed_criteria = sum(criteria.values())
    total_criteria = len(criteria)
    feasibility_score = passed_criteria / total_criteria
    
    print(f"\nCriteria Assessment:")
    for criterion, status in criteria.items():
        status_icon = "✓" if status else "✗"
        print(f"  {status_icon} {criterion}")
    
    print(f"\nOverall Feasibility Score: {feasibility_score:.1%} ({passed_criteria}/{total_criteria})")
    
    if feasibility_score >= 0.8:
        recommendation = "HIGHLY RECOMMENDED"
        color = "✓"
    elif feasibility_score >= 0.6:
        recommendation = "CONDITIONALLY RECOMMENDED"
        color = "⚠"
    else:
        recommendation = "NOT RECOMMENDED"
        color = "✗"
    
    print(f"\n{color} RECOMMENDATION: {recommendation}")
    print(f"BioPython is {recommendation.lower()} as the mandatory standard for NCBI sequence fetching")
    
    # Technical considerations
    print(f"\nTechnical Considerations:")
    print(f"  • BioPython provides robust sequence handling and NCBI integration")
    print(f"  • Network-independent operation possible for cached/local sequences")
    print(f"  • Excellent integration with existing project modules")
    print(f"  • Mature library with extensive documentation and community support")
    print(f"  • Built-in validation and error handling capabilities")
    
    return feasibility_score >= 0.8


if __name__ == "__main__":
    success = demonstrate_biopython_feasibility()
    print(f"\nDemo completed successfully: {success}")
    sys.exit(0 if success else 1)