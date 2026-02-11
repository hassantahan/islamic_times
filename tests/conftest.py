from __future__ import annotations

import sys
from datetime import datetime, timezone
from pathlib import Path

import pytest

ROOT_DIR = Path(__file__).resolve().parents[1]
SRC_DIR = ROOT_DIR / "src"
if SRC_DIR.exists():
    sys.path.insert(0, str(SRC_DIR))


@pytest.fixture()
def fixed_datetime_utc() -> datetime:
    return datetime(2025, 6, 1, 12, 0, 0, tzinfo=timezone.utc)


@pytest.fixture()
def toronto_observer_kwargs(fixed_datetime_utc: datetime) -> dict[str, float | datetime | bool | str]:
    return {
        "latitude": 43.651070,
        "longitude": -79.347015,
        "elevation": 10.0,
        "temperature": 15.0,
        "pressure": 101.325,
        "date": fixed_datetime_utc,
        "method": "ISNA",
        "find_local_tz": False,
    }
