#!/usr/bin/env python3
"""
dicom_extract_images_recursive_same_dir.py

Recursively extracts images and full metadata from DICOM files with *no loss of detail*,
writing outputs **in the same directory as each source file** (no subfolders). It detects
DICOM files by signature (not just by extension), so files *without* .dcm are handled too.

Per-file outputs:
  - <base>_sXXXX_raw16.png   (bit-perfect 16-bit grayscale, one per slice/frame)
  - <base>_sXXXX_view8.png   (display-optimized 8-bit PNG with VOI LUT + inversion if MONOCHROME1)
  - <base>_metadata.txt      (full tag dump + key attributes)
"""
import argparse
from pathlib import Path
from datetime import datetime

import numpy as np
from PIL import Image

import pydicom
from pydicom import dcmread
from pydicom.pixel_data_handlers.util import apply_voi_lut

# Prefer official helper if available
try:
    from pydicom.misc import is_dicom as _pydicom_is_dicom  # type: ignore
except Exception:
    _pydicom_is_dicom = None


def is_dicom_file(path: Path) -> bool:
    """
    Robust DICOM detection:
      1) Try pydicom.misc.is_dicom if present.
      2) Fallback to Part 10 signature check: 128-byte preamble + 'DICM'.
      3) As last resort, try a very fast header read via dcmread(stop_before_pixels=True).
    """
    try:
        if _pydicom_is_dicom is not None:
            return bool(_pydicom_is_dicom(str(path)))
    except Exception:
        pass

    try:
        with path.open("rb") as f:
            header = f.read(132)
            if len(header) >= 132 and header[128:132] == b"DICM":
                return True
    except Exception:
        return False

    # Last resort: attempt lightweight read (can catch non-Part10 DICOMs)
    try:
        dcmread(str(path), stop_before_pixels=True, force=True)
        return True
    except Exception:
        return False


def _to_uint16_bitperfect(arr: np.ndarray) -> np.ndarray:
    """
    Return a uint16 array that preserves the original *bits* exactly.
    If the source is int16, reinterpret the same memory as uint16 (two's complement).
    If uint8/int8, upcast without scaling (replicate bits for visibility 8->16).
    """
    if arr.dtype == np.uint16:
        return arr
    if arr.dtype == np.int16:
        return arr.view(np.uint16)
    if arr.dtype == np.uint8:
        return arr.astype(np.uint16) * 257
    if arr.dtype == np.int8:
        return arr.view(np.uint8).astype(np.uint16) * 257
    if np.issubdtype(arr.dtype, np.integer):
        return arr.astype(np.uint16, copy=False)
    # Floating fallback (rare). Normalize without clipping.
    amin = float(np.nanmin(arr))
    amax = float(np.nanmax(arr))
    if amax > amin:
        norm = (arr - amin) / (amax - amin)
    else:
        norm = np.zeros_like(arr, dtype=np.float32)
    return (norm * 65535.0 + 0.5).astype(np.uint16)


def _save_png_u16(path: Path, arr_u16: np.ndarray) -> None:
    im = Image.fromarray(arr_u16, mode='I;16')
    im.save(str(path))


def _save_png_u8(path: Path, arr_u8: np.ndarray) -> None:
    im = Image.fromarray(arr_u8, mode='L')
    im.save(str(path))


def _window_to_uint8(frame, ds: pydicom.Dataset) -> np.ndarray:
    """
    Apply VOI LUT / windowing (if present) and convert to 8-bit for display.
    Handle MONOCHROME1 inversion for display correctness.
    """
    try:
        voi = apply_voi_lut(frame, ds)
    except Exception:
        voi = frame

    voi = np.asarray(voi, dtype=np.float32)
    vmin = float(np.min(voi))
    vmax = float(np.max(voi))
    if vmax > vmin:
        img8 = (255.0 * (voi - vmin) / (vmax - vmin)).astype(np.uint8)
    else:
        img8 = np.zeros_like(voi, dtype=np.uint8)

    photometric = str(getattr(ds, "PhotometricInterpretation", "")).upper()
    if photometric == "MONOCHROME1":
        img8 = 255 - img8

    return img8


def _slice_iter(arr: np.ndarray, samples_per_pixel: int):
    """
    Yield (slice_index, 2D frame) for arrays that may be 2D, 3D, or 4D.
    """
    if arr.ndim == 2:
        yield 0, arr
        return
    if arr.ndim == 3:
        if samples_per_pixel == 1 and arr.shape[0] > 1:
            for i in range(arr.shape[0]):
                yield i, arr[i]
        else:
            yield 0, arr
        return
    if arr.ndim == 4:
        for i in range(arr.shape[0]):
            yield i, arr[i]
        return
    first = arr.shape[0]
    for i in range(first):
        yield i, np.squeeze(arr[i])


def extract_dicom_same_dir(input_path: Path, write_view: bool = True) -> None:
    ds = dcmread(str(input_path))

    try:
        px = ds.pixel_array
    except Exception as e:
        raise RuntimeError(f"Failed to decode Pixel Data for {input_path}: {e}")

    samples_per_pixel = int(getattr(ds, "SamplesPerPixel", 1))

    base = input_path.stem
    out_dir = input_path.parent  # SAME directory as the .dcm

    # Metadata
    meta_path = out_dir / f"{base}_metadata.txt"
    with open(meta_path, "w", encoding="utf-8") as f:
        f.write(f"# DICOM Metadata Dump\n")
        f.write(f"# Source file: {input_path}\n")
        f.write(f"# Extracted: {datetime.utcnow().isoformat()}Z\n\n")
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

    count = 0
    for idx, frame in _slice_iter(px, samples_per_pixel):
        raw_u16 = _to_uint16_bitperfect(frame)
        raw_name = out_dir / f"{base}_s{idx:04d}_raw16.png"
        _save_png_u16(raw_name, raw_u16)

        if write_view:
            try:
                view8 = _window_to_uint8(frame, ds)
            except Exception:
                fr = np.asarray(frame, dtype=np.float32)
                mn, mx = float(fr.min()), float(fr.max())
                if mx > mn:
                    view8 = ((fr - mn) / (mx - mn) * 255.0).astype(np.uint8)
                else:
                    view8 = np.zeros_like(fr, dtype=np.uint8)
            view_name = out_dir / f"{base}_s{idx:04d}_view8.png"
            _save_png_u8(view_name, view8)

        count += 1

    print(f"[OK] {input_path} | slices: {count} | outputs in: {out_dir}")


def _iter_candidate_files(root: Path):
    """
    Recursively yield all regular files under root. We do NOT restrict to *.dcm;
    we detect DICOMs by signature to catch files without extensions.
    """
    for p in root.rglob("*"):
        if p.is_file():
            yield p


def main():
    ap = argparse.ArgumentParser(description="Recursively extract lossless images + metadata in the SAME directory as each DICOM. Accepts file or directory.")
    ap.add_argument("input_path", type=str, help="Path to a DICOM file or a directory to scan recursively")
    ap.add_argument("--no-view", action="store_true", help="Do not create 8-bit display PNGs (only 16-bit raw PNGs)")
    args = ap.parse_args()

    in_path = Path(args.input_path).expanduser().resolve()

    if in_path.is_file():
        if not is_dicom_file(in_path):
            raise SystemExit(f"Not a DICOM file (by signature): {in_path}")
        extract_dicom_same_dir(in_path, write_view=(not args.no_view))
        return

    if in_path.is_dir():
        any_found = False
        for f in _iter_candidate_files(in_path):
            if not is_dicom_file(f):
                continue
            any_found = True
            try:
                extract_dicom_same_dir(f, write_view=(not args.no_view))
            except Exception as e:
                print(f"[WARN] Skipping {f}: {e}")
        if not any_found:
            print(f"[INFO] No DICOM files found under: {in_path}")
        return

    raise SystemExit(f"Path not found: {in_path}")


if __name__ == "__main__":
    main()
