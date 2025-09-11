# experiments/smoke_mri_zbin.py
# Minimal smoke test: handles mixed slice shapes by grouping and classifies each group.
# Usage: python experiments/smoke_mri_zbin.py /path/to/DICOM
#
# Dependencies: pydicom, numpy
#   pip install pydicom numpy
#
# Output: For each shape group, prints a single-line verdict:
#   [shape (H, W)] <degenerative-like|old-trauma-like|ambiguous> | segments=<n> | widest=<rows> rows | coverage=<0..1> | focality=<ratio> | slices=<n>/<total>

import sys, os, glob
import numpy as np

try:
    import pydicom
except ImportError:
    print("Please: pip install pydicom numpy")
    sys.exit(1)

def iter_slices(dicom_dir):
    files = sorted(glob.glob(os.path.join(dicom_dir, "**", "*.dcm"), recursive=True))
    if not files:
        raise RuntimeError(f"No .dcm files found under {dicom_dir}")

    def sort_key(fp):
        try:
            ds = pydicom.dcmread(fp, stop_before_pixels=True, force=True)
            return int(getattr(ds, "InstanceNumber", 0))
        except Exception:
            return 0
    files.sort(key=sort_key)

    for fp in files:
        try:
            ds = pydicom.dcmread(fp, force=True)
            arr = ds.pixel_array.astype(np.float32)
            slope = float(getattr(ds, "RescaleSlope", 1.0))
            inter = float(getattr(ds, "RescaleIntercept", 0.0))
            arr = arr * slope + inter
            yield arr
        except Exception:
            # skip unreadable
            continue

def group_by_shape(dicom_dir):
    groups = {}
    total = 0
    for arr in iter_slices(dicom_dir):
        total += 1
        groups.setdefault(arr.shape, []).append(arr)
    return groups, total

def minmax01(img):
    # robust normalize to [0,1] using 1–99th percentiles; tolerate constant arrays
    img = np.asarray(img, dtype=np.float32)
    if not np.isfinite(img).any():
        return np.zeros_like(img, dtype=np.float32)
    lo, hi = np.nanpercentile(img, [1, 99])
    if not np.isfinite(lo): lo = np.nanmin(img)
    if not np.isfinite(hi): hi = np.nanmax(img)
    if hi <= lo:
        return np.zeros_like(img, dtype=np.float32)
    img = np.clip(img, lo, hi)
    return (img - lo) / (hi - hi + 1e-9 if hi == lo else (hi - lo))

def theta_prime(x, k=0.3):
    # θ′(x,k) = φ * frac(x/φ)^k ; for x∈[0,1], frac(x/φ)=x/φ
    phi = (1 + 5 ** 0.5) / 2.0
    x = np.clip(x, 0.0, 1.0).astype(np.float32)
    # Avoid domain issues for very small x
    base = x / phi
    base = np.clip(base, 0.0, 1.0)
    return phi * (base ** k)

def smooth1d(v, w=7):
    if w <= 1:
        return v
    k = np.ones(w, dtype=np.float32) / float(w)
    return np.convolve(v, k, mode="same")

def count_segments(v, thr=0.6, min_len_frac=0.02):
    n = len(v)
    min_len = max(1, int(n * min_len_frac))
    runs, cur = [], 0
    for val in v:
        if val >= thr:
            cur += 1
        else:
            if cur >= min_len:
                runs.append(cur)
            cur = 0
    if cur >= min_len:
        runs.append(cur)
    return runs

def classify_profile(profile):
    # Returns: label, nseg, maxw, foc_ratio, coverage
    v = smooth1d(profile, w=7)
    vmax, vmin = float(np.max(v)), float(np.min(v))
    if not np.isfinite(vmax) or vmax - vmin <= 0:
        return "ambiguous", 0, 0, 1.0, 0.0
    v = (v - vmin) / (vmax - vmin + 1e-9)
    runs = count_segments(v, thr=0.6, min_len_frac=0.02)
    nseg = len(runs)
    length = len(v)
    maxw = max(runs) if runs else 0
    coverage = (sum(runs) / float(length)) if length else 0.0
    vmean = float(np.mean(v)) if np.isfinite(v).all() else 1.0
    foc_ratio = float(np.max(v) / (vmean + 1e-9))

    # Simple explicit rules
    if nseg >= 3 and coverage >= 0.30:
        label = "degenerative-like"
    elif nseg <= 2 and foc_ratio >= 1.55 and (maxw / float(length)) <= 0.20:
        label = "old-trauma-like"
    else:
        label = "ambiguous"
    return label, nseg, maxw, foc_ratio, coverage

def analyze_stack(stack):
    # Take mid slice, build craniocaudal darkness profile, classify
    mid = stack[len(stack)//2]
    img = minmax01(mid)
    h, w = img.shape
    c0, c1 = int(w*0.35), int(w*0.65)
    if c1 <= c0:
        c0, c1 = 0, w
    strip = img[:, c0:c1]
    # Darkness = 1 - intensity (low T2 → dark)
    darkness = 1.0 - np.mean(strip, axis=1)
    # θ′ transform and re-normalize
    tp = theta_prime(darkness, k=0.3)
    tmax, tmin = float(np.max(tp)), float(np.min(tp))
    if tmax - tmin <= 0:
        tp = np.zeros_like(tp, dtype=np.float32)
    else:
        tp = (tp - tmin) / (tmax - tmin + 1e-9)
    return classify_profile(tp)

def main():
    if len(sys.argv) < 2:
        print("Usage: python experiments/smoke_mri_zbin.py /path/to/DICOM")
        sys.exit(2)

    dicom_dir = sys.argv[1]
    groups, total = group_by_shape(dicom_dir)

    shapes = list(groups.keys())
    shapes_sorted = sorted(shapes, key=lambda s: len(groups[s]), reverse=True)
    counts = [(s, len(groups[s])) for s in shapes_sorted]
    if len(shapes_sorted) > 1:
        print(f"Warning: mixed slice shapes {shapes_sorted}, processing each group separately ({sum(c for _, c in counts)}/{total} slices)")

    for shape in shapes_sorted:
        slices = groups[shape]
        try:
            stack = np.stack(slices, axis=0)
        except Exception:
            print(f"[shape {shape}] ambiguous | segments=0 | widest=0 rows | coverage=0.00 | focality=1.00 | slices={len(slices)}/{total}")
            continue

        label, nseg, maxw, foc, cov = analyze_stack(stack)
        print(f"[shape {shape}] {label} | segments={nseg} | widest={maxw} rows | coverage={cov:.2f} | focality={foc:.2f} | slices={len(slices)}/{total}")

if __name__ == "__main__":
    main()
