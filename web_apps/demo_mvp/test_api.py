#!/usr/bin/env python3
"""
Simple test for demo MVP API
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from fastapi.testclient import TestClient
from web_apps.demo_mvp.app import app

client = TestClient(app)

def test_health():
    """Test health endpoint."""
    response = client.get("/health")
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "healthy"
    print("✓ Health check passed")

def test_info():
    """Test info endpoint."""
    response = client.get("/info")
    assert response.status_code == 200
    data = response.json()
    assert "helical_period" in data
    assert data["helical_period"] == 10.5
    print("✓ Info endpoint passed")

def test_analyze():
    """Test analyze endpoint with valid sequences."""
    payload = {
        "guide": "GACGAUCGAUCGAUCGAUCG",  # RNA with U
        "target": "ATGCGATCGATCGATCGATCGCTAGCTAGCTA",  # DNA with T
        "rate_ratio": 20.0,
        "weight_mode": "rate_ratio"
    }
    
    response = client.post("/analyze", json=payload)
    assert response.status_code == 200
    
    data = response.json()
    assert data["success"] == True
    assert "peak" in data
    assert "p10" in data["peak"]
    assert "p5_25" in data["peak"]
    assert "p3_5" in data["peak"]
    assert "rotation" in data
    assert "breathing_features" in data
    
    print("✓ Analyze endpoint passed")
    print(f"  10.5 bp peak: {data['peak']['p10']:.4f}")
    print(f"  5.25 bp peak: {data['peak']['p5_25']:.4f}")
    print(f"  3.5 bp peak: {data['peak']['p3_5']:.4f}")
    print(f"  GC content: {data['breathing_features']['gc_content']:.2%}")

def test_invalid_guide():
    """Test that invalid guide sequences are rejected."""
    payload = {
        "guide": "GACGATCGATCGATCGATCG",  # DNA (T) not RNA - should fail
        "target": "ATGCGATCGATCGATCGATCGCTAGCTAGCTA",
        "rate_ratio": 20.0
    }
    
    response = client.post("/analyze", json=payload)
    # Should fail validation because guide should be RNA with U
    assert response.status_code == 422  # Validation error
    print("✓ Invalid guide rejection passed")

def test_invalid_target():
    """Test that invalid target sequences are rejected."""
    payload = {
        "guide": "GACGAUCGAUCGAUCGAUCG",  # Valid RNA
        "target": "ATGCGXYZGATCGATCGATCGCTAGCTAGCTA",  # Invalid bases
        "rate_ratio": 20.0
    }
    
    response = client.post("/analyze", json=payload)
    assert response.status_code == 422  # Validation error
    print("✓ Invalid target rejection passed")

if __name__ == "__main__":
    print("Testing Demo MVP API...")
    print()
    
    test_health()
    test_info()
    test_analyze()
    test_invalid_guide()
    test_invalid_target()
    
    print()
    print("✓ All tests passed!")
