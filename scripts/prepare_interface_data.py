#!/usr/bin/env python3
"""Rebuild the public UI bundle from the checked-in data."""

from __future__ import annotations

import sys
import time
from pathlib import Path

from tqdm import tqdm

from global_aggregator import export_ui_data_js, update_global_files

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    print("\n" + "=" * 70)
    print("ESM-2 Interface Data Bundle Regeneration")
    print("=" * 70 + "\n")
    
    # Step 1: Build global CSV files
    print("Step 1/2: Building global CSV files from alignment data...")
    print("  • Reading proteins metadata")
    print("  • Reading alignment files (MMseqs2, Needleman-Wunsch, TM-align, Foldseek)")
    print("  • Reading embedding similarity data")
    print()
    
    with tqdm(total=100, desc="Building global files", unit="%", ncols=70) as pbar:
        # Simulate progress tracking
        pbar.update(10)
        result = update_global_files(ROOT)
        pbar.update(90)
    
    print(f"✓ Updated proteins CSV: {result['proteins_csv']}")
    print(f"✓ Updated pairs CSV: {result['pairs_csv']}")
    print()
    
    # Step 2: Export UI data bundle
    print("Step 2/2: Generating UI data bundle (JavaScript)...")
    print("  • Loading CSV data")
    print("  • Coercing numeric values")
    print("  • Creating JSON payload")
    print()
    
    with tqdm(total=100, desc="Exporting UI bundle", unit="%", ncols=70) as pbar:
        pbar.update(20)
        out_js = export_ui_data_js(ROOT)
        pbar.update(80)
    
    print(f"✓ Wrote UI bundle: {out_js}")
    print()
    print("=" * 70)
    print("✓ Bundle regeneration complete!")
    print("=" * 70 + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())