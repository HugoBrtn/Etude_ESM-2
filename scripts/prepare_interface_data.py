#!/usr/bin/env python3
"""Rebuild the public UI bundle from the checked-in data."""

from __future__ import annotations

from pathlib import Path

from global_aggregator import export_ui_data_js, update_global_files

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    result = update_global_files(ROOT)
    out_js = export_ui_data_js(ROOT)
    print(f"[INFO] Updated proteins CSV: {result['proteins_csv']}")
    print(f"[INFO] Updated pairs CSV: {result['pairs_csv']}")
    print(f"[INFO] Wrote UI bundle: {out_js}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())