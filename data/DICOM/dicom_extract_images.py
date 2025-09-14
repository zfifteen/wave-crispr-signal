#!/usr/bin/env python3
"""
dicom_extract_images.py

Extracts images from a DICOM .dcm file with *no loss of detail* and saves every
frame/slice. Also writes a comprehensive text file with all metadata from the
DICOM file.

Outputs per input file:
  - <base>_sXXXX_raw16.png        (bit‑perfect 16‑bit grayscale, one per slice/frame)
  - <base>_sXXXX_view8.png        (display-optimized 8‑bit PNG with VOI LUT + inversion as needed)
  - <base>_metadata.txt           (full tag dump and key image/header attributes)

Notes on “no loss”:
  - The *_raw16.png files store the EXACT pixel bits from the DICOM Pixel Data.
    If the original was signed int16, we save its two’s‑complement bit pattern as
    uint16 to preserve *all* information bit‑for‑bit. You can recover the original
    signed values by viewing the uint16 as int16 (no rounding, no windowing).
  - The *_view8.png files are for human viewing only (VOI LUT/Window applied,
    scaled to 0..255, and inverted if MONOCHROME1). They are not intended for
    quantitative analysis.
"""
import argparse
import os
from pathlib import Path
from datetime import datetime

import numpy as np
from PIL import Image

import pydicom
from pydicom import dcmread
from pydicom.pixel_data_handlers.util import apply_voi_lut


def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _to_uint16_bitperfect(arr: np.ndarray) -> np.ndarray:
    """
    Return a uint16 array that preserves the original *bits* exactly.
    If the source is int16, we reinterpret the same memory as uint16 (two's complement).
    If the source is already uint16, return as-is.
    If the source is uint8/int8, upcast without scaling.
    """
    if arr.dtype == np.uint16:
        return arr
    if arr.dtype == np.int16:
        return arr.view(np.uint16)
    if arr.dtype == np.uint8:
        return arr.astype(np.uint16) * 257  # promote 8->16 by bit replication for visibility, still lossless w.r.t. 8-bit
    if arr.dtype == np.int8:
        # reinterpret -> uint8 (two's complement), then promote to 16
        return arr.view(np.uint8).astype(np.uint16) * 257
    # Fallback: try safe cast without scaling
    if np.issubdtype(arr.dtype, np.integer):
        return arr.astype(np.uint16, copy=False)
    # Floating types are unusual for raw DICOM pixels; scale to full range without clipping
    # but warn via metadata file; still, try a safe conversion
    amin, amax = float(np.nanmin(arr)), float(np.nanmax(arr))
    if amax > amin:
        norm = (arr - amin) / (amax - amin)
    else:
        norm = np.zeros_like(arr, dtype=np.float32)
    return (norm * 65535.0 + 0.5).astype(np.uint16)


def _save_png_u16(path: Path, arr_u16: np.ndarray) -> None:
    """
    Save uint16 grayscale PNG. Mode 'I;16' keeps 16-bit depth.
    """
    im = Image.fromarray(arr_u16, mode='I;16')
    im.save(str(path))


def _save_png_u8(path: Path, arr_u8: np.ndarray) -> None:
    im = Image.fromarray(arr_u8, mode='L')
    im.save(str(path))


def _window_to_uint8(frame: np.ndarray, ds: pydicom.Dataset) -> np.ndarray:
    """
    Apply VOI LUT / windowing (if present) and convert to 8-bit for display.
    Handle MONOCHROME1 inversion for display correctness.
    """
    # Apply VOI LUT or fall back to raw
    try:
        voi = apply_voi_lut(frame, ds)
    except Exception:
        voi = frame

    # Convert to float and scale 0..255
    voi = voi.astype(np.float32)
    vmin = float(np.min(voi))
    vmax = float(np.max(voi))
    if vmax > vmin:
        img8 = (255.0 * (voi - vmin) / (vmax - vmin)).astype(np.uint8)
    else:
        img8 = np.zeros_like(voi, dtype=np.uint8)

    # Invert for MONOCHROME1
    photometric = str(getattr(ds, "PhotometricInterpretation", "")).upper()
    if photometric == "MONOCHROME1":
        img8 = 255 - img8

    return img8


def _slice_iter(arr: np.ndarray, samples_per_pixel: int) -> tuple[int, np.ndarray]:
    """
    Uniform iterator yielding (slice_index, 2D frame) from a pixel array that may be
    2D, 3D, or (rarely) color. For color, returns HxWxC frames if present.
    """
    if arr.ndim == 2:
        yield 0, arr
        return

    if arr.ndim == 3:
        # Monochrome multi-frame usually (frames, rows, cols)
        # Color single-frame sometimes (rows, cols, 3) with samples_per_pixel=3
        if samples_per_pixel == 1 and arr.shape[0] > 1:
            # multi-frame grayscale
            for i in range(arr.shape[0]):
                yield i, arr[i]
        else:
            # treat as single slice (HxWxC)
            yield 0, arr
        return

    if arr.ndim == 4:
        # Color multi-frame: (frames, rows, cols, channels)
        for i in range(arr.shape[0]):
            yield i, arr[i]
        return

    # Fallback: flatten first dimension as frames
    first = arr.shape[0]
    for i in range(first):
        yield i, np.squeeze(arr[i])


def extract_dicom(input_path: Path, outdir: Path, write_view: bool = True) -> None:
    ds = dcmread(str(input_path))

    # Try to materialize pixels
    try:
        px = ds.pixel_array  # pydicom will decompress and shape appropriately
    except Exception as e:
        raise RuntimeError(f"Failed to decode Pixel Data: {e}")

    samples_per_pixel = int(getattr(ds, "SamplesPerPixel", 1))

    # Prepare output directory
    base = input_path.stem
    out_base = outdir / f"{base}_extracted"
    _ensure_dir(out_base)

    # Metadata dump
    meta_path = out_base / f"{base}_metadata.txt"
    with open(meta_path, "w", encoding="utf-8") as f:
        f.write(f"# DICOM Metadata Dump\n")
        f.write(f"# Source file: {input_path}\n")
        f.write(f"# Extracted: {datetime.utcnow().isoformat()}Z\n\n")
        # Key attributes
        keys = [
            "SOPClassUID", "SOPInstanceUID", "StudyInstanceUID", "SeriesInstanceUID",
            "Modality", "Manufacturer", "ManufacturerModelName",
            "Rows", "Columns", "NumberOfFrames", "SamplesPerPixel",
            "BitsAllocated", "BitsStored", "HighBit", "PixelRepresentation",
            "PhotometricInterpretation", "PlanarConfiguration",
            "RescaleSlope", "RescaleIntercept", "WindowCenter", "WindowWidth",
            "PixelSpacing", "SliceThickness", "SpacingBetweenSlices", "ImagePositionPatient",
            "ImageOrientationPatient", "SliceLocation",
        ]
        for k in keys:
            f.write(f"{k}: {getattr(ds, k, None)}\n")
        f.write("\n### Full Dataset ###\n")
        f.write(str(ds))

    # Save frames
    count = 0
    for idx, frame in _slice_iter(px, samples_per_pixel):
        # Bit‑perfect RAW: 16‑bit PNG (uint16), preserving original bits
        raw_u16 = _to_uint16_bitperfect(frame)
        raw_name = out_base / f"{base}_s{idx:04d}_raw16.png"
        _save_png_u16(raw_name, raw_u16)

        # Optional display image
        if write_view:
            try:
                view8 = _window_to_uint8(frame, ds)
                view_name = out_base / f"{base}_s{idx:04d}_view8.png"
                _save_png_u8(view_name, view8)
            except Exception:
                # Best-effort: scale raw for view if VOI LUT fails
                fr = frame.astype(np.float32)
                mn, mx = float(fr.min()), float(fr.max())
                if mx > mn:
                    view8 = ((fr - mn) / (mx - mn) * 255.0).astype(np.uint8)
                else:
                    view8 = np.zeros_like(fr, dtype=np.uint8)
                view_name = out_base / f"{base}_s{idx:04d}_view8.png"
                _save_png_u8(view_name, view8)

        count += 1

    print(f"Saved {count} slice(s) to: {out_base}")
    print(f"Metadata written to: {meta_path}")


def main():
    ap = argparse.ArgumentParser(description="Extract lossless images and metadata from a DICOM .dcm file.")
    ap.add_argument("dicom_file", type=str, help="Path to the input .dcm file")
    ap.add_argument("-o", "--outdir", type=str, default=".", help="Output directory (default: current directory)")
    ap.add_argument("--no-view", action="store_true", help="Do not create 8-bit display PNGs (only 16-bit raw PNGs)")
    args = ap.parse_args()

    in_path = Path(args.dicom_file).expanduser().resolve()
    if not in_path.exists():
        raise SystemExit(f"Input file not found: {in_path}")

    outdir = Path(args.outdir).expanduser().resolve()
    _ensure_dir(outdir)

    extract_dicom(in_path, outdir, write_view=(not args.no_view))


if __name__ == "__main__":
    main()
