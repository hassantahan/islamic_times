from __future__ import annotations

from datetime import datetime, timezone

import pytest

from islamic_times.islamic_times import ITLocation
from islamic_times.it_dataclasses import ObserverInfo


def test_init_rejects_non_numeric_latitude(toronto_observer_kwargs: dict[str, object]) -> None:
    kwargs = dict(toronto_observer_kwargs)
    kwargs["latitude"] = "43.6"
    with pytest.raises(TypeError, match="latitude"):
        ITLocation(**kwargs)  # type: ignore[arg-type]


def test_init_rejects_out_of_range_latitude(toronto_observer_kwargs: dict[str, object]) -> None:
    kwargs = dict(toronto_observer_kwargs)
    kwargs["latitude"] = 95.0
    with pytest.raises(ValueError, match="Latitudes must be between -90° and 90°"):
        ITLocation(**kwargs)


def test_init_rejects_out_of_range_longitude(toronto_observer_kwargs: dict[str, object]) -> None:
    kwargs = dict(toronto_observer_kwargs)
    kwargs["longitude"] = 190.0
    with pytest.raises(ValueError, match="Longitudes must be between -180° and 180°"):
        ITLocation(**kwargs)


def test_init_rejects_non_datetime_date(toronto_observer_kwargs: dict[str, object]) -> None:
    kwargs = dict(toronto_observer_kwargs)
    kwargs["date"] = "2025-06-01"
    with pytest.raises(TypeError, match="must be of type `datetime`"):
        ITLocation(**kwargs)  # type: ignore[arg-type]


def test_init_rejects_invalid_find_local_tz_type(toronto_observer_kwargs: dict[str, object]) -> None:
    kwargs = dict(toronto_observer_kwargs)
    kwargs["find_local_tz"] = "false"
    with pytest.raises(ValueError, match="must be of type `bool`"):
        ITLocation(**kwargs)  # type: ignore[arg-type]


def test_init_rejects_invalid_numeric_find_local_tz(toronto_observer_kwargs: dict[str, object]) -> None:
    kwargs = dict(toronto_observer_kwargs)
    kwargs["find_local_tz"] = 2
    with pytest.raises(TypeError, match="must be of type `bool`"):
        ITLocation(**kwargs)  # type: ignore[arg-type]


def test_init_rejects_unknown_prayer_method(toronto_observer_kwargs: dict[str, object]) -> None:
    kwargs = dict(toronto_observer_kwargs)
    kwargs["method"] = "NOT_A_METHOD"
    with pytest.raises(ValueError, match="Invalid prayer method"):
        ITLocation(**kwargs)


def test_init_rejects_invalid_asr_type(toronto_observer_kwargs: dict[str, object]) -> None:
    kwargs = dict(toronto_observer_kwargs)
    kwargs["asr_type"] = 3
    with pytest.raises(ValueError, match="'asr_type' must be either 0 or 1"):
        ITLocation(**kwargs)


def test_init_manual_mode_sets_modified_flags(toronto_observer_kwargs: dict[str, object]) -> None:
    kwargs = dict(toronto_observer_kwargs)
    kwargs["auto_calculate"] = False
    location = ITLocation(**kwargs)

    assert location.auto_calculate is False
    assert location.datetime_modified is True
    assert location.prayers_modified is True


def test_init_auto_mode_produces_observer_info(toronto_observer_kwargs: dict[str, object]) -> None:
    location = ITLocation(**toronto_observer_kwargs)
    observer = location.observer()
    assert isinstance(observer, ObserverInfo)
    assert observer.latitude.decimal == pytest.approx(43.651070)
    assert observer.longitude.decimal == pytest.approx(-79.347015)


def test_get_timezone_uses_utc_for_naive_datetime() -> None:
    location = ITLocation(date=datetime(2025, 6, 1, 12, 0, 0), find_local_tz=False)
    assert location.dates_times().utc_offset == pytest.approx(0.0)


def test_init_find_local_tz_tracks_named_timezone() -> None:
    location = ITLocation(
        latitude=43.651070,
        longitude=-79.347015,
        date=datetime(2025, 1, 15, 12, 0, 0, tzinfo=timezone.utc),
        find_local_tz=True,
    )
    dt_info = location.dates_times()

    assert dt_info.timezone is not None
    assert "/" in dt_info.timezone
    assert dt_info.utc_offset == pytest.approx(-5.0)
    assert dt_info.date.hour == 7


def test_init_find_local_tz_rejects_nonexistent_local_time() -> None:
    with pytest.raises(ValueError, match="does not exist"):
        ITLocation(
            latitude=43.651070,
            longitude=-79.347015,
            date=datetime(2025, 3, 9, 2, 30, 0),
            find_local_tz=True,
        )
