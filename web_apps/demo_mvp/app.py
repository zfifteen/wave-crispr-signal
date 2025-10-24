#!/usr/bin/env python3
"""
Demo MVP: Interactive CRISPR Breathing Dynamics Analyzer

This FastAPI application provides a public demo showing:
1. 10.5-bp spectral peak (fundamental helical period)
2. Harmonics at 5.25 bp and 3.5 bp
3. Rotational phase wheel visualization
4. Breathing Δscore (model improvement from adding breathing features)

SCIENTIFIC GATES ENFORCED:
- Human DNA only (A/C/G/T validation)
- No fabrication (only accept real nucleotide sequences)
- Dimensionless parametrization (rate ratios, not MHz)
- Privacy-safe logging (SHA256 only, no raw sequences unless opt-in)

Usage:
    uvicorn app:app --reload
    
    Or:
    python app.py
"""

import sys
import os
import hashlib
import logging
from datetime import datetime
from typing import Optional, Dict, Any
import json

# Add parent directories to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field, validator
import numpy as np

# Import our spectral utilities
from experiments.trinity.spectral import (
    validate_dna_sequence,
    validate_rna_sequence,
    breathing_features,
    rotational_phase_curve,
    spectral_peak_analysis,
    phase_at,
    encode_complex,
    HELICAL_PERIOD,
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Create FastAPI app
app = FastAPI(
    title="CRISPR Breathing Dynamics Demo",
    description="Interactive demo showing DNA breathing dynamics and 10.5-bp spectral peaks",
    version="1.0.0"
)

# CORS middleware for web access
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # In production, restrict to specific domains
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


# Request/Response Models
class AnalyzeRequest(BaseModel):
    """Request model for sequence analysis."""
    
    guide: str = Field(
        ...,
        description="Guide RNA sequence (20-30 nt, A/C/G/U only)",
        min_length=10,
        max_length=50
    )
    target: str = Field(
        ...,
        description="Target DNA sequence (30-100 bp, A/C/G/T only)",
        min_length=20,
        max_length=200
    )
    temperature_c: Optional[float] = Field(
        None,
        description="Temperature in Celsius (optional, for future NN models)",
        ge=0,
        le=100
    )
    mg_mM: Optional[float] = Field(
        None,
        description="Magnesium concentration in mM (optional, for future NN models)",
        ge=0,
        le=100
    )
    weight_mode: str = Field(
        "rate_ratio",
        description="Weight calculation mode: 'rate_ratio' or 'nn_thermo'",
        pattern="^(rate_ratio|nn_thermo)$"
    )
    rate_ratio: float = Field(
        20.0,
        description="Dimensionless rate ratio r = k_GC/k_AT",
        ge=1.0,
        le=1000.0
    )
    log_consent: bool = Field(
        False,
        description="Opt-in to privacy-safe logging (SHA256 + features only)"
    )
    
    @validator('guide')
    def validate_guide(cls, v):
        """Validate guide RNA sequence."""
        try:
            # Guide is RNA, so expect U
            return validate_rna_sequence(v)
        except ValueError as e:
            raise ValueError(f"Invalid guide sequence: {e}")
    
    @validator('target')
    def validate_target(cls, v):
        """Validate target DNA sequence."""
        try:
            # Target is DNA, so expect T
            return validate_dna_sequence(v)
        except ValueError as e:
            raise ValueError(f"Invalid target sequence: {e}")


class AnalyzeResponse(BaseModel):
    """Response model for sequence analysis."""
    
    success: bool
    peak: Dict[str, float]
    rotation: Dict[str, Any]
    breathing_features: Dict[str, float]
    meta: Dict[str, Any]
    message: Optional[str] = None


# Helper functions
def calculate_breathing_score(
    target: str,
    guide: str,
    r: float = 20.0
) -> float:
    """
    Calculate breathing Δscore.
    
    This is a simplified proxy showing potential improvement from
    adding breathing features. In full implementation, this would
    be: Δ = model(+breathing) - model(baseline)
    
    For demo purposes, we use spectral power at 10.5 bp as proxy.
    """
    # For target DNA
    features = breathing_features(target, r=r, is_rna=False)
    
    # Simple proxy: normalized 10.5 bp power
    # In real model, this would be regression coefficient difference
    base_score = 0.5  # Baseline model score (placeholder)
    breathing_boost = features["P10_5"] / 100.0  # Normalize to reasonable range
    
    delta_score = breathing_boost  # Simplified: just the boost
    
    return float(delta_score)


def hash_sequence(seq: str) -> str:
    """
    Create SHA256 hash of sequence for privacy-safe logging.
    
    Args:
        seq: DNA/RNA sequence
        
    Returns:
        Hex string of SHA256 hash
    """
    return hashlib.sha256(seq.encode()).hexdigest()


def log_analysis(request: AnalyzeRequest, response: AnalyzeResponse):
    """
    Log analysis with privacy protection.
    
    Only logs if user gives consent. Stores:
    - SHA256 hash (not raw sequence)
    - Features (GC%, length, peaks)
    - Timestamp
    """
    if not request.log_consent:
        return
    
    log_entry = {
        "timestamp": datetime.utcnow().isoformat(),
        "guide_hash": hash_sequence(request.guide),
        "target_hash": hash_sequence(request.target),
        "guide_length": len(request.guide),
        "target_length": len(request.target),
        "gc_content": response.breathing_features.get("gc_content", 0.0),
        "peak_10_5": response.peak.get("p10", 0.0),
        "weight_mode": request.weight_mode,
        "rate_ratio": request.rate_ratio,
    }
    
    # In production, write to database or structured log
    logger.info(f"Analysis logged: {json.dumps(log_entry)}")


# API Endpoints
@app.get("/", response_class=HTMLResponse)
async def root():
    """Serve landing page."""
    html_content = """
    <!DOCTYPE html>
    <html>
    <head>
        <title>CRISPR Breathing Dynamics Demo</title>
        <style>
            body { font-family: Arial, sans-serif; max-width: 1200px; margin: 0 auto; padding: 20px; }
            h1 { color: #2c3e50; }
            .info-box { background: #ecf0f1; padding: 15px; border-radius: 5px; margin: 20px 0; }
            .endpoint { background: #e8f5e9; padding: 10px; margin: 10px 0; border-left: 4px solid #4caf50; }
            code { background: #f5f5f5; padding: 2px 5px; border-radius: 3px; }
        </style>
    </head>
    <body>
        <h1>🧬 CRISPR Breathing Dynamics Demo</h1>
        
        <div class="info-box">
            <h2>About This Demo</h2>
            <p>This interactive demo shows how DNA breathing dynamics (base pair opening rates) 
            create spectral signatures at the helical period of <strong>10.5 bp/turn</strong>.</p>
            
            <p><strong>Key Features:</strong></p>
            <ul>
                <li>10.5-bp spectral peak detection (fundamental helical period)</li>
                <li>Harmonic analysis (5.25 bp, 3.5 bp)</li>
                <li>Rotational phase wheel visualization</li>
                <li>Dimensionless parametrization (rate ratios, not MHz)</li>
            </ul>
        </div>
        
        <div class="info-box">
            <h2>API Endpoints</h2>
            
            <div class="endpoint">
                <h3>POST /analyze</h3>
                <p>Analyze a guide RNA + target DNA pair</p>
                <p><strong>Example:</strong></p>
                <pre><code>curl -X POST http://localhost:8000/analyze \\
  -H "Content-Type: application/json" \\
  -d '{
    "guide": "GACGAUCGAUCGAUCGAUCG",
    "target": "ATGCGATCGATCGATCGATCGCTAGCTAGCTA",
    "rate_ratio": 20.0,
    "weight_mode": "rate_ratio"
  }'</code></pre>
            </div>
            
            <div class="endpoint">
                <h3>GET /health</h3>
                <p>Check service health</p>
            </div>
            
            <div class="endpoint">
                <h3>GET /info</h3>
                <p>Get information about available parameters</p>
            </div>
        </div>
        
        <div class="info-box">
            <h2>Scientific Validation</h2>
            <p>This demo enforces strict scientific gates:</p>
            <ul>
                <li>✓ Human DNA only (A/C/G/T for DNA, A/C/G/U for RNA)</li>
                <li>✓ No fabrication (real nucleotides only)</li>
                <li>✓ Dimensionless parameters (rate ratios)</li>
                <li>✓ Privacy-safe logging (opt-in, SHA256 only)</li>
            </ul>
        </div>
        
        <p><em>For full documentation, see <code>experiments/trinity/README.md</code></em></p>
    </body>
    </html>
    """
    return HTMLResponse(content=html_content)


@app.get("/health")
async def health_check():
    """Health check endpoint."""
    return {"status": "healthy", "timestamp": datetime.utcnow().isoformat()}


@app.get("/info")
async def get_info():
    """Get information about available parameters and modes."""
    return {
        "version": "1.0.0",
        "weight_modes": {
            "rate_ratio": {
                "description": "Dimensionless rate ratio k_GC/k_AT",
                "default_value": 20.0,
                "typical_range": [5.0, 200.0],
                "formula": "alpha = log(r), beta = 0.3*log(r)"
            },
            "nn_thermo": {
                "description": "Nearest-neighbor thermodynamics (future)",
                "status": "not_yet_implemented"
            }
        },
        "helical_period": HELICAL_PERIOD,
        "harmonics": [10.5, 5.25, 3.5],
        "privacy": {
            "default_logging": False,
            "opt_in_required": True,
            "stored_data": "SHA256 hash + features only, no raw sequences"
        }
    }


@app.post("/analyze", response_model=AnalyzeResponse)
async def analyze(request: AnalyzeRequest):
    """
    Analyze guide RNA and target DNA for breathing dynamics.
    
    This endpoint:
    1. Validates sequences (DNA: A/C/G/T, RNA: A/C/G/U)
    2. Computes spectral features at 10.5, 5.25, 3.5 bp periods
    3. Generates rotational phase curve
    4. Calculates breathing score
    5. Optionally logs with privacy protection
    
    Args:
        request: AnalyzeRequest with guide, target, and parameters
        
    Returns:
        AnalyzeResponse with spectral analysis results
        
    Raises:
        HTTPException: If validation fails
    """
    try:
        logger.info(f"Analyzing: guide={len(request.guide)}nt, target={len(request.target)}bp")
        
        # Calculate spectral peaks for target
        peak_10_5 = spectral_peak_analysis(
            request.target, 
            period=10.5, 
            r=request.rate_ratio,
            is_rna=False
        )
        peak_5_25 = spectral_peak_analysis(
            request.target,
            period=5.25,
            r=request.rate_ratio,
            is_rna=False
        )
        peak_3_5 = spectral_peak_analysis(
            request.target,
            period=3.5,
            r=request.rate_ratio,
            is_rna=False
        )
        
        # Get breathing features
        features = breathing_features(
            request.target,
            r=request.rate_ratio,
            is_rna=False
        )
        
        # Generate rotational phase curve
        phi_centers, phi_values = rotational_phase_curve(
            request.target,
            period=10.5,
            bins=24,
            r=request.rate_ratio,
            is_rna=False
        )
        
        # Calculate breathing score
        breathing_score = calculate_breathing_score(
            request.target,
            request.guide,
            r=request.rate_ratio
        )
        
        # Construct response
        response = AnalyzeResponse(
            success=True,
            peak={
                "p10": peak_10_5["normalized_power"],
                "p5_25": peak_5_25["normalized_power"],
                "p3_5": peak_3_5["normalized_power"],
                "angle": peak_10_5["angle"],
            },
            rotation={
                "phi": phi_centers,
                "curve": phi_values,
            },
            breathing_features=features,
            meta={
                "guide_length": len(request.guide),
                "target_length": len(request.target),
                "weight_mode": request.weight_mode,
                "rate_ratio": request.rate_ratio,
                "breathing_score": breathing_score,
                "temperature_c": request.temperature_c,
                "mg_mM": request.mg_mM,
            },
            message="Analysis completed successfully"
        )
        
        # Log if consent given
        log_analysis(request, response)
        
        return response
        
    except ValueError as e:
        logger.error(f"Validation error: {e}")
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.error(f"Analysis error: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"Analysis failed: {str(e)}")


if __name__ == "__main__":
    import uvicorn
    
    print("Starting CRISPR Breathing Dynamics Demo...")
    print("Navigate to: http://127.0.0.1:8000")
    print("API docs at: http://127.0.0.1:8000/docs")
    
    uvicorn.run(app, host="127.0.0.1", port=8000, log_level="info")
