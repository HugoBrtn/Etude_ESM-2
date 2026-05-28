#!/usr/bin/env python3
"""Launch the Flask UI for the public dataset."""

from __future__ import annotations

from pathlib import Path
import runpy

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    runpy.run_path(str(ROOT / "ui" / "app.py"), run_name="__main__")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())