#!/usr/bin/env python3
"""Cleanup legacy per-pair artifacts to save disk space.

Safe defaults:
- Keep alignment.txt.gz files used by the UI.
- Only remove embedding_similarity per-pair folders if a summary exists.
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"


def _remove_path(path: Path, dry_run: bool) -> None:
    if not path.exists():
        return
    if dry_run:
        print(f"[DRY] remove {path}")
        return
    if path.is_dir():
        shutil.rmtree(path, ignore_errors=True)
    else:
        path.unlink(missing_ok=True)


def cleanup_alignment_dir(root: Path, dry_run: bool) -> None:
    if not root.exists():
        return
    for pair_dir in sorted([p for p in root.iterdir() if p.is_dir()]):
        _remove_path(pair_dir / "tmp", dry_run)
        _remove_path(pair_dir / "alignment.tsv", dry_run)
        _remove_path(pair_dir / "metadata.json", dry_run)
        _remove_path(pair_dir / "alignment.txt", dry_run)
        if not (pair_dir / "alignment.txt.gz").exists():
            # Drop empty pair folders that no longer keep alignment text
            if dry_run:
                print(f"[DRY] remove {pair_dir}")
            else:
                try:
                    pair_dir.rmdir()
                except OSError:
                    pass


def cleanup_embedding_similarity(root: Path, dry_run: bool) -> None:
    if not root.exists():
        return
    summary_path = root / "pairwise_summary.csv"
    if not summary_path.exists():
        print("[WARN] No embedding summary found; keeping per-pair folders")
        return
    for pair_dir in sorted([p for p in root.iterdir() if p.is_dir()]):
        _remove_path(pair_dir, dry_run)


def main() -> int:
    parser = argparse.ArgumentParser(description="Cleanup legacy alignment artifacts to save disk space")
    parser.add_argument("--dry-run", action="store_true", help="Preview deletions without modifying files")
    args = parser.parse_args()

    cleanup_alignment_dir(DATA_DIR / "alignment_mmseq2", args.dry_run)
    cleanup_alignment_dir(DATA_DIR / "alignment_needleman_wunsh", args.dry_run)
    cleanup_alignment_dir(DATA_DIR / "alignment_structure_tmscore", args.dry_run)
    cleanup_embedding_similarity(DATA_DIR / "embedding_similarity", args.dry_run)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
