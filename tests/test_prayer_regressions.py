from datetime import datetime, timedelta, timezone
from math import atan, degrees
from zoneinfo import ZoneInfo

import pytest
import pytz

from islamic_times.islamic_times import ITLocation
from islamic_times.it_dataclasses import Angle
from islamic_times.prayer_times import asr_time
from islamic_times.time_equations import localize_or_convert_datetime


@pytest.mark.parametrize("asr_type, shadow_factor", [(0, 1), (1, 2)])
def test_asr_at_zero_noon_shadow(asr_type: int, shadow_factor: int) -> None:
    noon = datetime(2026, 3, 20, 12, tzinfo=timezone.utc)
    expected = noon + timedelta(hours=degrees(atan(shadow_factor)) / 15)
    assert asr_time(noon, Angle(0), Angle(0), asr_type) == expected


@pytest.mark.parametrize("month, day", [(3, 7), (10, 31)])
@pytest.mark.parametrize("tz", [ZoneInfo("America/Toronto"), pytz.timezone("America/Toronto")])
@pytest.mark.parametrize("method, endpoint", [("ISNA", "sunrise"), ("JAFARI", "fajr")])
def test_midnight_is_elapsed_midpoint_on_dst_nights(month: int, day: int, tz, method: str, endpoint: str) -> None:
    wall_noon = datetime(2026, month, day, 12)
    dates = [localize_or_convert_datetime(wall_noon + timedelta(days=i), tz) for i in range(2)]
    first, second = [
        ITLocation(latitude=43.653661, longitude=-79.382924, date=date, method=method).prayer_times()
        for date in dates
    ]
    sunset = first.sunset.when
    target = getattr(second, endpoint).when
    midnight = first.midnight.when
    assert sunset is not None and target is not None and midnight is not None
    assert sunset.utcoffset() != target.utcoffset()
    start, end, actual = [event.astimezone(timezone.utc) for event in (sunset, target, midnight)]
    assert (actual - (start + (end - start) / 2)).total_seconds() == pytest.approx(0, abs=0.1)


@pytest.mark.parametrize("latitude, longitude", [(69.65, 18.96), (78.22, 15.65)])
def test_makkah_polar_isha_uses_adjusted_maghrib(latitude: float, longitude: float) -> None:
    with pytest.warns(UserWarning, match="Extreme latitude"):
        observer = ITLocation(
            latitude=latitude, longitude=longitude,
            date=datetime(2026, 6, 21, 12, tzinfo=ZoneInfo("Europe/Oslo")), method="MAKKAH",
        )
    result = observer.prayer_times()
    assert result.maghrib.when is not None and result.isha.when is not None
    assert result.isha.when - result.maghrib.when == timedelta(minutes=90)
    assert result.isha.status == "estimated"
    assert result.isha.fallback_rule == "ANGLEBASED"
    assert result.zuhr.status == "calculated"
    assert result.isha.to_dict()["status"] == "estimated"

    with pytest.warns(UserWarning, match="Extreme latitude"):
        observer.set_extreme_latitude_rule("NONE")
    result = observer.prayer_times()
    assert result.isha.when is None
    assert result.isha.status == "unavailable"
    assert not result.extreme_latitude_applied
    assert result.zuhr.status == "calculated"


def test_makkah_ramadan_interval_survives_asr_override() -> None:
    observer = ITLocation(
        latitude=21.42, longitude=39.82,
        date=datetime(2026, 2, 25, 12, tzinfo=ZoneInfo("Asia/Riyadh")), method="MAKKAH",
    )
    observer.set_asr_type(1)
    result = observer.prayer_times()
    assert result.isha.when is not None and result.maghrib.when is not None
    assert result.isha.when - result.maghrib.when == timedelta(minutes=120)


def test_missing_reference_solar_event_does_not_fail_the_day() -> None:
    observer = ITLocation(
        latitude=78, longitude=18.96, date=datetime(2026, 4, 21, 12, tzinfo=timezone.utc),
        method="MAKKAH", auto_calculate=False,
    )
    with pytest.warns(UserWarning, match="Extreme latitude"):
        result = observer.prayer_schedule(observer.observer_dateinfo.date.date(),
                                          observer.observer_dateinfo.date.date())[0].prayers
    assert result.zuhr.status == "calculated"
    assert result.extreme_latitude_reason is not None
