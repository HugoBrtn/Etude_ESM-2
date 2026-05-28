#!/usr/bin/env python3
"""Build UI data bundle from GLOBAL CSV files."""

from __future__ import annotations

from pathlib import Path

from global_aggregator import export_ui_data_js

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    out_js = export_ui_data_js(ROOT)
    print(f"[INFO] Wrote {out_js}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
