from __future__ import annotations

import math
from datetime import datetime, timedelta, timezone

from islamic_times.it_dataclasses import (
    Angle,
    DateTimeInfo,
    Distance,
    DistanceUnit,
    DistanceUnits,
    IslamicDateInfo,
    MeccaInfo,
    MoonInfo,
    ObserverInfo,
    Prayer,
    PrayerMethod,
    PrayerTimes,
    RightAscension,
    SunInfo,
    Visibilities,
)


def test_angle_properties_and_string_representation() -> None:
    angle = Angle(-12.5)
    assert angle.dms == (-12, 30, 0.0)
    assert math.isclose(angle.radians, math.radians(-12.5), rel_tol=0.0, abs_tol=1e-12)
    assert "°" in str(angle)


def test_right_ascension_conversions() -> None:
    ra = RightAscension(5.5)
    assert ra.hms[0] == 5
    assert isinstance(ra.decimal_degrees, Angle)
    assert isinstance(ra.radians, float)
    assert "h" in ra.hms_str()


def test_distance_unit_conversion_and_formatting() -> None:
    custom_unit = DistanceUnit("unit", "u", 2.0)
    assert math.isclose(custom_unit.convert_to(10, DistanceUnits.METRE), 20.0, abs_tol=1e-12)

    distance = Distance(1, DistanceUnits.KILOMETRE)
    assert math.isclose(distance.in_unit(DistanceUnits.METRE), 1000.0, abs_tol=1e-12)
    converted = distance.to(DistanceUnits.MILE)
    assert isinstance(converted, Distance)
    assert converted.unit == DistanceUnits.MILE


def test_islamic_date_helpers() -> None:
    islamic = IslamicDateInfo(1447, 9, 1)
    assert islamic.hijri_month_name == "Ramaḍān"
    assert islamic.hijri_day_of_week_name("Friday") == "al-Jumuʿah"
    assert "Ramaḍān" in islamic.full_date("Friday")


def test_observerinfo_defaults_are_physical() -> None:
    observer = ObserverInfo(latitude=Angle(0), longitude=Angle(0), elevation=Distance(0))
    assert math.isclose(observer.pressure, 101.325, abs_tol=1e-12)
    assert math.isclose(observer.temperature, 10.0, abs_tol=1e-12)


def test_datetimeinfo_properties() -> None:
    tz = timezone(timedelta(hours=-5))
    dt = datetime(2025, 1, 1, 12, 0, 0, tzinfo=tz)
    info = DateTimeInfo(date=dt, jd=2460000.5, deltaT=69.0, hijri=IslamicDateInfo(1446, 7, 1))

    assert info.gregorian_date.endswith("2025")
    assert info.timezone is not None
    assert math.isclose(info.utc_offset, -5.0, abs_tol=1e-12)
    assert math.isclose(info.jde, 2460000.5 + 69.0 / 86400.0, rel_tol=0.0, abs_tol=1e-12)
    assert info.format_utc_offset() == "UTC-05:00"
    assert "Time & Date" in str(info)


def test_prayer_and_prayertimes_string_output() -> None:
    method = PrayerMethod(name="Test", fajr_angle=Angle(18), isha_angle=Angle(17))
    now = datetime(2025, 6, 1, 12, 0, 0)

    prayers = [
        Prayer("Fajr", now, method),
        Prayer("Sunrise", now, method),
        Prayer("Ẓuhr", now, method),
        Prayer("ʿAṣr", now, method),
        Prayer("Sunset", now, method),
        Prayer("Maghrib", now, method),
        Prayer("ʿIshāʾ", now, method),
        Prayer("Midnight", now, method),
    ]
    prayer_times = PrayerTimes(
        method=method,
        fajr=prayers[0],
        sunrise=prayers[1],
        zuhr=prayers[2],
        asr=prayers[3],
        sunset=prayers[4],
        maghrib=prayers[5],
        isha=prayers[6],
        midnight=prayers[7],
    )

    assert prayers[0].time_str.endswith("01-06-2025")
    assert "Prayer Times at Observer Timezone" in str(prayer_times)
    assert "Midnight" in str(prayer_times)


def test_mecca_sun_moon_and_visibility_strings() -> None:
    now = datetime(2025, 6, 1, 12, 0, 0)
    mecca = MeccaInfo(distance=Distance(1000, DistanceUnits.KILOMETRE), angle=Angle(45), cardinal="NE")
    assert "Mecca" in str(mecca)

    sun = SunInfo(
        sunrise=now,
        sun_transit=now,
        sunset=now,
        apparent_altitude=Angle(10),
        true_azimuth=Angle(180),
        geocentric_distance=Distance(1, DistanceUnits.AU),
        apparent_declination=Angle(20),
        apparent_right_ascension=RightAscension(6),
        greenwich_hour_angle=Angle(30),
        local_hour_angle=Angle(40),
    )
    assert "The Sun" in str(sun)

    moon = MoonInfo(
        moonrise=now,
        moon_transit=now,
        moonset="Moonset does not exist.",
        illumination=0.5,
        apparent_altitude=Angle(12),
        true_azimuth=Angle(200),
        geocentric_distance=Distance(384400, DistanceUnits.KILOMETRE),
        parallax=Angle(1),
        topocentric_declination=Angle(5),
        topocentric_right_ascension=RightAscension(7),
        greenwich_hour_angle=Angle(20),
        local_hour_angle=Angle(25),
    )
    assert "The Moon" in str(moon)

    vis = Visibilities(
        criterion="Yallop",
        dates=(now, now + timedelta(days=1)),
        q_values=(0.1234, -0.5678),
        classifications=("A: Easily visible.", "F: Not visible."),
    )
    rendered = str(vis)
    assert "Visibility of New Moon Crescent" in rendered
    assert "Criterion" in rendered


def test_visibilities_str_with_zero_q_should_include_context_lines() -> None:
    vis = Visibilities(
        criterion="Yallop",
        dates=(datetime(2025, 6, 1, 12, 0, 0, tzinfo=timezone.utc),),
        q_values=(0.0,),
        classifications=("A: Easily visible.",),
    )
    rendered = str(vis)
    assert "Visibility of New Moon Crescent" in rendered
    assert "Criterion" in rendered
