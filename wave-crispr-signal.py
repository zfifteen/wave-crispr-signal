import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft
from scipy.stats import entropy
from collections import Counter

# --- Base Weights (Wave Function Mapping) ---
weights = {'A': 1 + 0j, 'T': -1 + 0j, 'C': 0 + 1j, 'G': 0 - 1j}

# --- Mock PCSK9 Exon 1 Sequence (150 bp) ---
sequence = "ATGCTGCGGAGACCTGGAGAGAAAGCAGTGGCCGGGGCAGTGGGAGGAGGAGGAGCTGGAAGAGGAGAGAAAGGAGGAGCTGCAGGAGGAGAGGAGGAGGAGGGAGAGGAGGAGCTGGAGCTGAAGCTGGAGCTGGAGCTGGAGAGGAGAGAGGG"

def build_waveform(seq, d=0.34, zn_map=None):
    if zn_map is None:
        spacings = [d]*len(seq)
    else:
        spacings = [d*(1+zn_map.get(i, 0)) for i in range(len(seq))]
    s = np.cumsum(spacings)
    wave = [weights[base]*np.exp(2j * np.pi * s[i]) for i, base in enumerate(seq)]
    return np.array(wave)

def compute_spectrum(waveform):
    return np.abs(fft(waveform))

def normalized_entropy(spectrum):
    ps = spectrum / np.sum(spectrum)
    return entropy(ps, base=2)

def count_sidelobes(spectrum, threshold_ratio=0.25):
    peak = np.max(spectrum)
    return np.sum(spectrum > (threshold_ratio * peak))

def mutate_and_analyze(seq, pos, new_base):
    if seq[pos] == new_base:
        return None
    mutated = list(seq)
    mutated[pos] = new_base
    zn = pos / len(seq)
    zn_map = {pos: zn}

    base_wave = build_waveform(seq)
    mut_wave = build_waveform(mutated, zn_map=zn_map)

    base_spec = compute_spectrum(base_wave)
    mut_spec = compute_spectrum(mut_wave)

    f1_index = 10
    delta_f1 = 100 * (mut_spec[f1_index] - base_spec[f1_index]) / base_spec[f1_index]

    side_lobe_delta = count_sidelobes(mut_spec) - count_sidelobes(base_spec)
    entropy_jump = normalized_entropy(mut_spec) - normalized_entropy(base_spec)

    hotspot_score = zn * abs(delta_f1) + side_lobe_delta + entropy_jump

    return {
        "position": pos,
        "wt": seq[pos],
        "edit": new_base,
        "zn": zn,
        "delta_f1 (%)": delta_f1,
        "side_lobe_Δ": side_lobe_delta,
        "entropy_Δ": entropy_jump,
        "hotspot_score": hotspot_score
    }

# --- Run Analysis ---
results = []
for pos in range(0, len(sequence), 15):
    for b in "ATCG":
        r = mutate_and_analyze(sequence, pos, b)
        if r:
            results.append(r)

# --- Display Top Edits ---
results.sort(key=lambda r: -r["hotspot_score"])
for r in results[:6]:
    print(f"n={r['position']:>3}  {r['wt']}→{r['edit']}  "
          f"Δf₁={r['delta_f1 (%)']:+.1f}%  "
          f"SideLobesΔ={r['side_lobe_Δ']}  "
          f"Score={r['hotspot_score']:.2f}")

# --- Optional: Plot baseline spectrum ---
base_wave = build_waveform(sequence)
spec = compute_spectrum(base_wave)
plt.plot(spec)
plt.title("Baseline Vibrational Spectrum (PCSK9 Exon 1)")
plt.xlabel("Frequency Index")
plt.ylabel("Magnitude")
plt.grid(True)
plt.show()
