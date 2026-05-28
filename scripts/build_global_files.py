#!/usr/bin/env python3
"""CLI wrapper to rebuild global aggregated files."""

from __future__ import annotations

from pathlib import Path

from global_aggregator import update_global_files

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    result = update_global_files(ROOT)
    print(f"[INFO] Updated proteins CSV: {result['proteins_csv']}")
    print(f"[INFO] Updated pairs CSV: {result['pairs_csv']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
