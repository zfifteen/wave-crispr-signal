# Models

This directory contains pretrained weights, model definitions, and exported model artifacts.

## Contents

This directory is intended for:

- **Model Architectures** - Definitions of neural network or ML model structures
- **Pretrained Weights** - Trained model parameters (not committed to git)
- **ONNX Exports** - Exported models in ONNX format for deployment
- **Model Metadata** - Configuration files describing model specifications

## File Organization

- **architectures/** - Model architecture definitions (Python files)
- **weights/** - Model weights and checkpoints (excluded from git)
- **exports/** - Exported model files for deployment
- **configs/** - Model-specific configuration files

## .gitignore Policy

Large model files should be excluded from git:

```gitignore
# Model weights and checkpoints
*.pth
*.pkl
*.h5
*.weights
checkpoints/

# Large exports
*.onnx
*.pb
```

Only model definitions, configurations, and small reference files should be committed.

## Usage

```python
from models.architectures import CRISPREfficiencyModel
from models.configs import model_config

# Load model architecture
model = CRISPREfficiencyModel(model_config)

# Load pretrained weights (if available)
# model.load_weights('models/weights/crispr_model.pth')
```

## Note

This directory structure is prepared for future model development. Currently, the project focuses on signal processing and mathematical analysis rather than machine learning models.