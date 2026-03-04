from __future__ import annotations

from datetime import datetime, timedelta, timezone

import pytest

from islamic_times.islamic_times import ITLocation
from islamic_times.it_dataclasses import DateTimeInfo, PrayerTimes


def test_manual_mode_blocks_access_until_recompute(toronto_observer_kwargs: dict[str, object]) -> None:
    kwargs = dict(toronto_observer_kwargs)
    kwargs["auto_calculate"] = False
    location = ITLocation(**kwargs)

    with pytest.raises(ValueError):
        location.dates_times()
    with pytest.raises(ValueError):
        location.prayer_times()
    with pytest.raises(ValueError):
        location.sun()
    with pytest.raises(ValueError):
        location.moon()

    location.calculate_astro()
    assert isinstance(location.dates_times(), DateTimeInfo)

    location.calculate_prayer_times()
    assert isinstance(location.prayer_times(), PrayerTimes)


def test_manual_mode_calculate_prayer_times_requires_astro_first(toronto_observer_kwargs: dict[str, object]) -> None:
    kwargs = dict(toronto_observer_kwargs)
    kwargs["auto_calculate"] = False
    location = ITLocation(**kwargs)

    with pytest.raises(ValueError, match="prayer times cannot be calculated"):
        location.calculate_prayer_times()


def test_update_time_requires_datetime_argument(toronto_observer_kwargs: dict[str, object]) -> None:
    location = ITLocation(**toronto_observer_kwargs)
    with pytest.raises(TypeError, match="must be of type `datetime`"):
        location.update_time("2025-06-01")  # type: ignore[arg-type]


def test_update_time_preserves_existing_timezone_for_naive_input(toronto_observer_kwargs: dict[str, object]) -> None:
    location = ITLocation(**toronto_observer_kwargs)

    location.update_time(datetime(2025, 6, 2, 6, 0, 0))
    updated_date = location.dates_times().date
    assert updated_date.tzinfo is not None
    assert updated_date.utcoffset() == timedelta(0)


def test_update_time_with_local_tz_refreshes_dst_offset() -> None:
    location = ITLocation(
        latitude=43.651070,
        longitude=-79.347015,
        date=datetime(2025, 1, 15, 12, 0, 0, tzinfo=timezone.utc),
        find_local_tz=True,
    )
    winter_offset = location.dates_times().utc_offset

    location.update_time(datetime(2025, 7, 15, 12, 0, 0, tzinfo=timezone.utc))
    summer_offset = location.dates_times().utc_offset

    assert winter_offset == pytest.approx(-5.0)
    assert summer_offset == pytest.approx(-4.0)


def test_update_time_with_local_tz_rejects_ambiguous_naive_wall_time() -> None:
    location = ITLocation(
        latitude=43.651070,
        longitude=-79.347015,
        date=datetime(2025, 1, 15, 12, 0, 0, tzinfo=timezone.utc),
        find_local_tz=True,
    )

    with pytest.raises(ValueError, match="ambiguous"):
        location.update_time(datetime(2025, 11, 2, 1, 30, 0))


def test_manual_mode_update_marks_state_as_stale(toronto_observer_kwargs: dict[str, object]) -> None:
    kwargs = dict(toronto_observer_kwargs)
    kwargs["auto_calculate"] = False
    location = ITLocation(**kwargs)
    location.calculate_astro()
    location.calculate_prayer_times()
    _ = location.prayer_times()

    location.update_time(datetime(2025, 6, 3, 12, 0, 0, tzinfo=timezone.utc))

    with pytest.raises(ValueError):
        location.prayer_times()


def test_manual_mode_setters_keep_prayer_state_dirty(toronto_observer_kwargs: dict[str, object]) -> None:
    kwargs = dict(toronto_observer_kwargs)
    kwargs["auto_calculate"] = False
    location = ITLocation(**kwargs)
    assert location.prayers_modified is True

    location.set_asr_type(1)
    assert location.prayers_modified is True

    location.set_midnight_type(1)
    assert location.prayers_modified is True

    location.set_extreme_latitude_rule("ANGLEBASED")
    assert location.prayers_modified is True
