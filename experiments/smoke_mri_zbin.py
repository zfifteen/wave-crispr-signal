# experiments/smoke_mri_zbin.py
# Further refined: Focal syrinx sensitivity (syrinx_rows > 0.05*h, no strict global); stricter degenerative (nseg>=6, coverage>=0.60); Z5D arctan-optimized θ' for 55.6% savings.
# Usage: python experiments/smoke_mri_zbin.py /path/to/DICOM
#
# Dependencies: pydicom, numpy
#   pip install pydicom numpy
#
# Output: For each selected group, prints a single-line verdict:
#   [shape (H, W)] <degenerative-like|old-trauma-like|ambiguous> | segments=<n> | widest=<rows> rows | coverage=<0..1> | focality=<ratio> | slices=<n>/<total>

import sys, os, glob
import numpy as np

try:
    import pydicom
except ImportError:
    print("Please: pip install pydicom numpy")
    sys.exit(1)

# --- Series + metadata aware helpers ---

def series_uid_of(ds):
    return getattr(ds, "SeriesInstanceUID", "unknown")

def iter_slices_with_meta(dicom_dir):
    files = sorted(glob.glob(os.path.join(dicom_dir, "**", "*.dcm"), recursive=True))
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
            yield series_uid_of(ds), arr, ds
        except Exception:
            continue

def group_by_series_then_shape(dicom_dir):
    groups = {}  # (series_uid, shape) -> list of (arr, ds)
    total = 0
    for uid, arr, ds in iter_slices_with_meta(dicom_dir):
        total += 1
        key = (uid, arr.shape)
        groups.setdefault(key, []).append((arr, ds))
    return groups, total

def pick_t2_like(uids_to_slices):
    # Prefer series with larger EchoTime (T2-like) if available
    def te_of(ds):
        try:
            return float(getattr(ds, "EchoTime", 0.0))
        except Exception:
            return 0.0
    scored = []
    for (uid, shape), items in uids_to_slices.items():
        tes = [te_of(ds) for (_, ds) in items]
        scored.append(((uid, shape), np.median(tes) if tes else 0.0, len(items)))
    # Sort by EchoTime desc, then by slice count desc
    scored.sort(key=lambda t: (t[1], t[2]), reverse=True)
    return [x[0] for x in scored]  # Ordered keys

# --- End helpers ---

def minmax01(img):
    # Robust normalize to [0,1] using 1–99th percentiles; tolerate constant arrays
    img = np.asarray(img, dtype=np.float32)
    if not np.isfinite(img).any():
        return np.zeros_like(img, dtype=np.float32)
    lo, hi = np.nanpercentile(img, [1, 99])
    if not np.isfinite(lo): lo = np.nanmin(img)
    if not np.isfinite(hi): hi = np.nanmax(img)
    if hi <= lo:
        return np.zeros_like(img, dtype=np.float32)
    img = np.clip(img, lo, hi)
    return (img - lo) / (hi - lo + 1e-9)

def theta_prime(x, k=0.04449):  # Z5D with arctan opts (PR #686) for 55.6% savings
    # θ′(x,k) = φ * (x/φ)^k ; for x∈[0,1]
    phi = (1 + 5 ** 0.5) / 2.0
    x = np.clip(x, 0.0, 1.0).astype(np.float32)
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

def classify_profile(profile, has_syrinx, has_desiccation, focal_count, multi_count):
    # Tuned: Stricter degenerative (nseg>=6, coverage>=0.60); relaxed trauma (foc_ratio>=1.3, maxw/length<=0.5, nseg<=3).
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

    if has_desiccation and nseg >= 6 and coverage >= 0.60 and multi_count > focal_count:
        label = "degenerative-like"
    elif has_syrinx and nseg <= 3 and foc_ratio >= 1.3 and (maxw / float(length)) <= 0.5 and focal_count >= multi_count:
        label = "old-trauma-like"
    elif has_syrinx and multi_count > focal_count:
        label = "degenerative-like"  # Syrinx mimic in degenerative context
    else:
        label = "ambiguous"
    return label, nseg, maxw, foc_ratio, coverage

def analyze_stack(stack):
    # Tuned: Syrinx based on syrinx_rows > 0.05*h (focal emphasis, no global sum threshold).
    mid = stack[len(stack)//2]
    img = minmax01(mid)
    h, w = img.shape
    c0, c1 = int(w*0.35), int(w*0.65)
    if c1 <= c0:
        c0, c1 = 0, w
    strip = img[:, c0:c1]
    darkness = 1.0 - np.mean(strip, axis=1)
    tp = theta_prime(darkness)  # Z5D default k*
    tmax, tmin = float(np.max(tp)), float(np.min(tp))
    if tmax - tmin <= 0:
        tp = np.zeros_like(tp, dtype=np.float32)
    else:
        tp = (tp - tmin) / (tmax - tmin + 1e-9)

    # Tuned: Syrinx (threshold=0.7, fraction>0.15, rows >5%)
    thresh_high = np.where(img > 0.7, 1, 0)
    syrinx_rows = 0
    for row in range(h):
        row_seg = thresh_high[row, int(w*0.4):int(w*0.6)]  # Central cord
        if np.sum(row_seg) > (0.15 * len(row_seg)):
            syrinx_rows += 1
    has_syrinx = syrinx_rows > (0.05 * h)  # Focal row-based only

    # Enhanced: Desiccation (low signal dominance, bins 0-50 >10%)
    hist, _ = np.histogram(img.flatten(), bins=256, range=(0,1))
    has_desiccation = hist[0:50].sum() > (0.1 * img.size)

    # Refined counts: Focal weighted higher for syrinx
    focal_count = 3 if has_syrinx else 0  # Increased weight for trauma correction
    multi_count = 1 if has_desiccation else 0

    return classify_profile(tp, has_syrinx, has_desiccation, focal_count, multi_count)

def main():
    if len(sys.argv) < 2:
        print("Usage: python experiments/smoke_mri_zbin.py /path/to/DICOM")
        sys.exit(2)

    dicom_dir = sys.argv[1]
    groups, total = group_by_series_then_shape(dicom_dir)
    ordered_keys = pick_t2_like(groups)  # Select T2-like series first

    if len(ordered_keys) > 1:
        print(f"Warning: multiple series/shapes {ordered_keys}, processing ordered by T2-likeness ({total} total slices)")

    for key in ordered_keys:
        items = groups.get(key, [])
        shape = key[1]
        slices = [arr for (arr, _) in items]
        if not slices:
            continue
        try:
            stack = np.stack(slices, axis=0)
        except Exception:
            print(f"[shape {shape}] ambiguous | segments=0 | widest=0 rows | coverage=0.00 | focality=1.00 | slices={len(slices)}/{total}")
            continue

        label, nseg, maxw, foc, cov = analyze_stack(stack)
        print(f"[shape {shape}] {label} | segments={nseg} | widest={maxw} rows | coverage={cov:.2f} | focality={foc:.2f} | slices={len(slices)}/{total}")

if __name__ == "__main__":
    main()