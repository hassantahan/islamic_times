"""Compatibility shim for mapper workflows.

This module keeps the historical `mapper.py` entrypoint while delegating to the
new package implementation under `islamic_times.mapper`.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Support running from repository root after migrating to src-layout.
ROOT_DIR = Path(__file__).resolve().parent
SRC_DIR = ROOT_DIR / "src"
if SRC_DIR.exists():
    sys.path.insert(0, str(SRC_DIR))

from islamic_times.mapper import _MAPPER_IMPORT_ERROR, main, require_mapper_dependencies
from islamic_times.mapper.cli import main as _cli_main

__all__ = ["_MAPPER_IMPORT_ERROR", "require_mapper_dependencies", "main"]


if __name__ == "__main__":
    raise SystemExit(_cli_main())

