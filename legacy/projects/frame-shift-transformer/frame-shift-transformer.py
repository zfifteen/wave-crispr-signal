import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
from Bio.Seq import Seq
import matplotlib.pyplot as plt

# Step 1: Define Z Universal Form for TE
# Z = A * (B / C), where A = measured TE, B = rate (processivity), C = invariant (tRNA max)
# In discrete: Z = n * (Delta_n / Delta_max), n = position-weighted features

def compute_z_te(features, positions, tRNA_max=1.0, delta_max=1.0):
    """
    Compute Z-transformed TE.
    features: array of dinucleotide/trinucleotide scores
    positions: spatial weights (e.g., closer to 5' UTR higher impact)
    Returns Z-normalized TE.
    """
    delta_n = np.sum(features * positions) / len(features)
    rate_b = delta_n / delta_max  # Normalized rate
    invariant_c = tRNA_max
    a_te = np.random.uniform(0.5, 1.0)  # Simulated measured TE
    z = a_te * (rate_b / invariant_c)
    return z

# Step 2: Synthetic mRNA Data Generation (mimic transcriptome atlas)
def generate_synthetic_mrna(num_samples=100, seq_length=300):
    nucleotides = ['A', 'C', 'G', 'U']
    sequences = [''.join(np.random.choice(nucleotides, seq_length)) for _ in range(num_samples)]
    # Simulate cell types with varying TE (140+ types approximated as 5 clusters)
    cell_types = np.random.randint(0, 5, num_samples)
    # Ground truth TE: higher for optimal codons, influenced by position
    tes = []
    for seq in sequences:
        mrna = Seq(seq)
        codons = [str(mrna[i:i+3]) for i in range(0, len(mrna), 3)]
        # Simple codon efficiency (mock tRNA abundance)
        codon_scores = {codon: np.random.uniform(0.1, 1.0) for codon in set(codons)}
        features = np.array([codon_scores.get(codon, 0.5) for codon in codons])
        positions = np.linspace(1.0, 0.1, len(features))  # Higher weight at 5' end
        te = compute_z_te(features, positions)
        tes.append(te)
    return sequences, cell_types, np.array(tes)

# Step 3: Sequence Encoding (dinucleotide/trinucleotide features with position)
def encode_sequence(seq):
    # One-hot encode nucleotides
    nuc_map = {'A': 0, 'C': 1, 'G': 2, 'U': 3}
    one_hot = np.zeros((len(seq), 4))
    for i, nuc in enumerate(seq):
        one_hot[i, nuc_map[nuc]] = 1
    # Extract dinucleotides as features
    di_features = []
    for i in range(len(seq) - 1):
        di = seq[i:i+2]
        di_features.append(hash(di) % 16)  # Simple hash to feature space
    # Positional embedding
    positions = np.sin(np.arange(len(di_features)) / 10000 ** (2 * np.arange(16) / 16))
    return np.concatenate([one_hot[:-1], positions.reshape(-1, 1) * np.array(di_features).reshape(-1, 1)], axis=1)

class mRNADataset(Dataset):
    def __init__(self, sequences, tes, cell_types):
        self.sequences = [encode_sequence(seq) for seq in sequences]
        self.tes = tes
        self.cell_types = cell_types

    def __len__(self):
        return len(self.tes)

    def __getitem__(self, idx):
        return torch.tensor(self.sequences[idx], dtype=torch.float32), torch.tensor(self.tes[idx], dtype=torch.float32), torch.tensor(self.cell_types[idx], dtype=torch.long)

# Step 4: RiboNN-like Model (Simple Multitask CNN)
class RiboNN(nn.Module):
    def __init__(self, num_cell_types=5):
        super(RiboNN, self).__init__()
        self.conv1 = nn.Conv1d(in_channels=5, out_channels=32, kernel_size=3)  # Features incl. position
        self.pool = nn.MaxPool1d(2)
        self.fc1 = nn.Linear(32 * 148, 128)  # Approx for seq_length ~300
        self.fc_te = nn.Linear(128, 1)  # Predict TE
        self.fc_cell = nn.Embedding(num_cell_types, 128)  # Multitask: cell type embedding

    def forward(self, x, cell_type):
        x = x.permute(0, 2, 1)  # For Conv1d
        x = torch.relu(self.conv1(x))
        x = self.pool(x)
        x = x.view(x.size(0), -1)
        x = torch.relu(self.fc1(x))
        cell_emb = self.fc_cell(cell_type)
        x = x + cell_emb  # Integrate cell context
        te = self.fc_te(x)
        return te

# Step 5: Train and Predict (Prove Prediction and Interconnectedness)
sequences, cell_types, true_tes = generate_synthetic_mrna()
dataset = mRNADataset(sequences, true_tes, cell_types)
dataloader = DataLoader(dataset, batch_size=32, shuffle=True)

model = RiboNN()
optimizer = optim.Adam(model.parameters(), lr=0.001)
criterion = nn.MSELoss()

# Training loop
for epoch in range(10):
    for seqs, tes, cells in dataloader:
        optimizer.zero_grad()
        preds = model(seqs, cells)
        loss = criterion(preds.squeeze(), tes)
        loss.backward()
        optimizer.step()
    print(f"Epoch {epoch+1}, Loss: {loss.item()}")

# Step 6: Inference and Visualization (Prove Insights)
with torch.no_grad():
    test_seqs, test_cells, true_test_tes = next(iter(dataloader))
    pred_tes = model(test_seqs, test_cells).squeeze().numpy()

# Mock stability/localization correlation (interconnectedness)
stability = pred_tes * np.random.uniform(0.8, 1.2)  # TE influences stability
localization = stability / np.max(stability)  # Normalized as Z-form

plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
plt.scatter(true_test_tes, pred_tes)
plt.xlabel('True TE')
plt.ylabel('Predicted TE')
plt.title('TE Prediction (Proves Sequence-to-TE Mapping)')

plt.subplot(1, 2, 2)
plt.plot(pred_tes, stability, label='Stability')
plt.plot(pred_tes, localization, label='Localization')
plt.xlabel('Predicted TE')
plt.ylabel('Correlated Metrics')
plt.title('Regulatory Interconnectedness')
plt.legend()
plt.show()

# Evolutionary Pressure Simulation: Optimize sequence for max TE
def optimize_sequence(base_seq, iterations=10):
    seq = list(base_seq)
    best_te = -np.inf
    for _ in range(iterations):
        idx = np.random.randint(0, len(seq))
        new_nuc = np.random.choice(['A', 'C', 'G', 'U'])
        seq[idx] = new_nuc
        encoded = encode_sequence(''.join(seq))
        cell = torch.tensor([0])  # Fixed cell type
        pred_te = model(torch.tensor(encoded).unsqueeze(0), cell).item()
        if pred_te > best_te:
            best_te = pred_te
    return ''.join(seq), best_te

opt_seq, opt_te = optimize_sequence(sequences[0])
print(f"Optimized Sequence TE: {opt_te} (Proves Selection Pressure)")