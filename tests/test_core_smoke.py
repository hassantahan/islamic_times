from __future__ import annotations

from datetime import datetime, timezone

from islamic_times.islamic_times import ITLocation
from islamic_times.it_dataclasses import ObserverInfo, PrayerTimes, Visibilities


def test_public_api_smoke() -> None:
    location = ITLocation(
        latitude=43.651070,
        longitude=-79.347015,
        elevation=10.0,
        temperature=15.0,
        pressure=101.325,
        date=datetime(2025, 6, 1, 12, 0, 0, tzinfo=timezone.utc),
        method="ISNA",
        find_local_tz=False,
    )

    observer = location.observer()
    prayers = location.prayer_times()

    assert isinstance(observer, ObserverInfo)
    assert isinstance(prayers, PrayerTimes)
    assert prayers.fajr.name.lower() == "fajr"
    assert prayers.sunrise.time is not None


def test_visibilities_smoke() -> None:
    location = ITLocation(
        latitude=43.651070,
        longitude=-79.347015,
        elevation=10.0,
        date=datetime(2025, 6, 1, 12, 0, 0, tzinfo=timezone.utc),
        find_local_tz=False,
    )
    visibilities = location.visibilities(days=1, criterion=1)

    assert isinstance(visibilities, Visibilities)
    assert len(visibilities.dates) == 1
    assert len(visibilities.q_values) == 1
    assert len(visibilities.classifications) == 1


def test_visibilities_smoke_shaukat() -> None:
    location = ITLocation(
        latitude=43.651070,
        longitude=-79.347015,
        elevation=10.0,
        date=datetime(2025, 6, 1, 12, 0, 0, tzinfo=timezone.utc),
        find_local_tz=False,
    )
    visibilities = location.visibilities(days=1, criterion=2)

    assert isinstance(visibilities, Visibilities)
    assert visibilities.criterion == "Shaukat"
    assert len(visibilities.dates) == 1
    assert len(visibilities.q_values) == 1
    assert len(visibilities.classifications) == 1
