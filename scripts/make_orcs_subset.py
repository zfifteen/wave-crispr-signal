#!/usr/bin/env python3
# scripts/make_orcs_subset.py
# Build a small, CI-friendly subset of the BioGRID ORCS 'screens' table using only local data.
# - Input (first existing): data/BIOGRID-ORCS-ALL-homo-sapiens-1.1.17.screens
#                           data/BIOGRID-ORCS-ALL-homo_sapiens-LATEST.screens
# - Output: data/ORCS_subset.screens
# - Filters (defaults; can be overridden by optional CLI flags): cell lines, screen types, throughput
# - No network, no changes to data/*. Read-only.

from __future__ import annotations
import os, sys, csv, argparse
from typing import List

DEFAULT_INPUTS = [
    "data/BIOGRID-ORCS-ALL-homo-sapiens-1.1.17.screens",
    "data/BIOGRID-ORCS-ALL-homo_sapiens-LATEST.screens",
]

DEFAULT_CELL_LINES = ["HCT 116", "DLD-1", "HeLa", "hTERT-RPE1"]
DEFAULT_SCREEN_TYPES = ["Negative Selection", "Positive Selection"]
DEFAULT_THROUGHPUT = "High Throughput"

OUT_PATH = "data/ORCS_subset.screens"

def find_input() -> str:
    for p in DEFAULT_INPUTS:
        if os.path.exists(p):
            return p
    sys.stderr.write(
        "[error] No ORCS screens file found. Expected one of:\n  - "
        + "\n  - ".join(DEFAULT_INPUTS)
        + "\n"
    )
    sys.exit(2)

def parse_args():
    ap = argparse.ArgumentParser(description="Make a filtered subset of BioGRID ORCS screens.")
    ap.add_argument("--cell-lines", nargs="*", default=DEFAULT_CELL_LINES,
                    help="Cell lines to keep (exact match on CELL_LINE column).")
    ap.add_argument("--screen-types", nargs="*", default=DEFAULT_SCREEN_TYPES,
                    help="Screen types to keep (exact match on SCREEN_TYPE).")
    ap.add_argument("--throughput", default=DEFAULT_THROUGHPUT,
                    help="Substring required in THROUGHPUT column (default: 'High Throughput').")
    ap.add_argument("--in", dest="in_path", default=None,
                    help="Override input path. If not given, auto-detects.")
    ap.add_argument("--out", dest="out_path", default=OUT_PATH,
                    help="Output subset path (default: data/ORCS_subset.screens).")
    return ap.parse_args()

def main():
    args = parse_args()
    in_path = args.in_path or find_input()
    out_path = args.out_path

    # Read file: preserve comment lines (#...), detect header, then TSV rows
    comments: List[str] = []
    rows: List[List[str]] = []
    header: List[str] | None = None
    total = 0

    with open(in_path, "r", encoding="utf-8", newline="") as f:
        for line in f:
            if line.startswith("#"):
                # Special case: header line starts with #SCREEN_ID
                if line.startswith("#SCREEN_ID"):
                    comments.append(line.rstrip("\n"))
                    if header is None:
                        # Remove the # prefix and use as header
                        header = line.rstrip("\n")[1:].split("\t")
                        continue
                else:
                    comments.append(line.rstrip("\n"))
                    continue
            if header is None:
                header = line.rstrip("\n").split("\t")
                continue
            total += 1
            rows.append(line.rstrip("\n").split("\t"))

    if header is None:
        sys.stderr.write("[error] Could not find header row in ORCS screens file.\n")
        sys.exit(3)

    # Build column index map (robust to column order)
    idx = {name: i for i, name in enumerate(header)}

    required_cols = ["CELL_LINE", "THROUGHPUT", "SCREEN_TYPE"]
    for c in required_cols:
        if c not in idx:
            sys.stderr.write(f"[error] Required column '{c}' not found in header.\n")
            sys.exit(4)

    keep_rows: List[List[str]] = []
    for r in rows:
        cell_line = r[idx["CELL_LINE"]]
        throughput = r[idx["THROUGHPUT"]]
        screen_type = r[idx["SCREEN_TYPE"]]

        if args.cell_lines and cell_line not in args.cell_lines:
            continue
        if args.screen_types and screen_type not in args.screen_types:
            continue
        if args.throughput and (args.throughput not in throughput):
            continue
        keep_rows.append(r)

    # Write subset: preserve comments and header
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", encoding="utf-8", newline="") as out:
        for c in comments:
            out.write(c + "\n")
        out.write("\t".join(header) + "\n")
        for r in keep_rows:
            out.write("\t".join(r) + "\n")

    print(f"rows_in={total} rows_out={len(keep_rows)} input={in_path} output={out_path}")

if __name__ == "__main__":
    main()