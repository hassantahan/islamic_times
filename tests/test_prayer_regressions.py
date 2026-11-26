from dataclasses import replace
from datetime import datetime, timedelta, timezone
from math import atan, degrees
from zoneinfo import ZoneInfo

import pytest
import pytz

from islamic_times import sun_equations
from islamic_times.islamic_times import ITLocation
from islamic_times.it_dataclasses import Angle, AsrConvention
from islamic_times.prayer_profiles import MosqueProfile
from islamic_times.prayer_times import asr_time, get_prayer_method
from islamic_times.time_equations import localize_or_convert_datetime, observer_datetime_info


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


@pytest.mark.parametrize(
    "method, fajr_angle, isha_angle",
    [("fcna canada", 13, 13), ("KUWAIT", 18, 17.5)],
)
def test_published_angle_presets_follow_twilight_and_sunset(
    method: str, fajr_angle: float, isha_angle: float,
) -> None:
    observer = ITLocation(
        latitude=43.65, longitude=-79.38,
        date=datetime(2026, 11, 25, 12, tzinfo=ZoneInfo("America/Toronto")), method=method,
    )
    result = observer.prayer_times()
    assert result.method == get_prayer_method(method)
    assert result.method.asr_type == AsrConvention.STANDARD
    assert result.maghrib.when == result.sunset.when
    for prayer, angle in ((result.fajr, fajr_angle), (result.isha, isha_angle)):
        assert prayer.when is not None
        sun = sun_equations.sunpos(observer_datetime_info(prayer.when), observer.observer_info)
        assert sun.true_altitude.decimal == pytest.approx(-angle, abs=0.01)


@pytest.mark.parametrize("day", [datetime(2026, 2, 25, 12), datetime(2026, 7, 25, 12)])
def test_qatar_keeps_fixed_isha_interval_in_and_out_of_ramadan(day: datetime) -> None:
    observer = ITLocation(
        latitude=25.29, longitude=51.53, date=day.replace(tzinfo=ZoneInfo("Asia/Qatar")), method="QATAR",
    )
    result = observer.prayer_times()
    assert result.maghrib.when is not None and result.isha.when is not None
    assert result.maghrib.when == result.sunset.when
    assert result.isha.when - result.maghrib.when == timedelta(minutes=90)


def test_qatar_polar_maghrib_and_isha_follow_the_estimated_sunset() -> None:
    with pytest.warns(UserWarning, match="Extreme latitude"):
        observer = ITLocation(
            latitude=69.65, longitude=18.96,
            date=datetime(2026, 6, 21, 12, tzinfo=ZoneInfo("Europe/Oslo")), method="QATAR",
        )
    result = observer.prayer_times()
    assert result.maghrib.when is not None and result.isha.when is not None
    assert result.maghrib.when == result.sunset.when
    assert result.isha.when - result.maghrib.when == timedelta(minutes=90)
    assert result.sunset.fallback_rule == result.maghrib.fallback_rule == result.isha.fallback_rule == "MIDDLENIGHT"


def test_sunset_offset_survives_profile_serialization_and_old_profiles_load() -> None:
    profile = MosqueProfile.from_preset("Qatar method", "QATAR")
    assert MosqueProfile.from_dict(profile.to_dict()) == profile
    previous = MosqueProfile.from_preset("Previous profile", "MAKKAH")
    payload = previous.to_dict()
    del payload["method"]["maghrib_offset_minutes"]
    assert MosqueProfile.from_dict(payload) == previous


def test_explicit_sunset_offset_is_applied_before_fixed_isha() -> None:
    observer = ITLocation(
        latitude=25.29, longitude=51.53, date=datetime(2026, 7, 25, 12, tzinfo=ZoneInfo("Asia/Qatar")),
        method="QATAR",
    )
    method = replace(get_prayer_method("QATAR"), maghrib_offset_minutes=1)
    result = observer.prayer_times(profile=MosqueProfile("Sunset plus one minute", method))
    assert result.sunset.when is not None and result.maghrib.when is not None and result.isha.when is not None
    assert result.maghrib.when - result.sunset.when == timedelta(minutes=1)
    assert result.isha.when - result.maghrib.when == timedelta(minutes=90)


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
