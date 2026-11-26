from csv import DictReader
from datetime import date, datetime, timedelta, timezone
from io import StringIO
from math import cos, sin
from zoneinfo import ZoneInfo

import pytest

from islamic_times import sun_equations as se
from islamic_times import time_equations as te
from islamic_times.islamic_times import ITLocation
from islamic_times.it_dataclasses import Angle, AsrConvention
from islamic_times.prayer_calendar import prayer_schedule_to_csv, prayer_schedule_to_ics
from islamic_times.prayer_fadilah import _western_redness_delay
from islamic_times.prayer_profiles import MosqueProfile


@pytest.mark.parametrize(
    ("lat", "lon", "local_date"),
    [
        (43.65, -79.38, datetime(2026, 3, 8, 12, tzinfo=ZoneInfo("America/Toronto"))),
        (21.42, 39.83, datetime(2026, 5, 28, 12, tzinfo=ZoneInfo("Asia/Riyadh"))),
        (-37.81, 144.96, datetime(2026, 11, 24, 12, tzinfo=ZoneInfo("Australia/Melbourne"))),
    ],
)
def test_sistani_windows_follow_the_observers_afternoon_shadow(
    lat: float, lon: float, local_date: datetime,
) -> None:
    observer = ITLocation(latitude=lat, longitude=lon, date=local_date)
    windows = observer.sistani_daylight_fadilah()
    zuhr, asr = windows["zuhr"], windows["asr"]
    sun = observer.sun()
    assert isinstance(sun.sunrise, datetime) and isinstance(sun.sunset, datetime)
    assert zuhr.start == te.time_midpoint(sun.sunrise, sun.sunset)
    assert zuhr.preferred_by == asr.start
    assert zuhr.end == asr.preferred_by
    assert zuhr.start < asr.start < zuhr.end < asr.end < sun.sunset

    def shadow_ratio(instant: datetime) -> float:
        altitude = se.sunpos(te.observer_datetime_info(instant), observer.observer_info).true_altitude.radians
        return cos(altitude) / sin(altitude)

    baseline = shadow_ratio(zuhr.start)
    for instant, parts in ((asr.start, 2), (zuhr.end, 4), (asr.end, 6)):
        assert instant.tzinfo is local_date.tzinfo
        assert shadow_ratio(instant) - baseline == pytest.approx(parts / 7, abs=0.0006)

    observer.set_prayer_method("ISNA")
    assert observer.sistani_daylight_fadilah() == windows


def test_polar_day_does_not_use_a_fallback_for_fadilah() -> None:
    local_date = datetime(2026, 6, 21, 12, tzinfo=ZoneInfo("Europe/Oslo"))
    with pytest.warns(UserWarning, match="Extreme latitude"):
        observer = ITLocation(latitude=69.65, longitude=18.96, date=local_date)
    assert observer.sistani_daylight_fadilah() == {}
    assert observer.sistani_evening_fadilah() == {}
    assert observer.prayer_times().asr.when is None
    assert observer.prayer_times().asr.fadilah is None


@pytest.mark.parametrize("method", ["JAFARI", "TEHRAN"])
def test_shia_presets_schedule_preferred_asr_and_allow_explicit_overrides(method: str) -> None:
    local_date = datetime(2026, 11, 25, 12, tzinfo=ZoneInfo("America/Toronto"))
    observer = ITLocation(latitude=43.65, longitude=-79.38, date=local_date, method=method)
    window = observer.sistani_daylight_fadilah()["asr"]
    asr = observer.prayer_times().asr
    assert asr.method.asr_type == AsrConvention.SISTANI_PREFERRED
    assert asr.when == window.start
    assert asr.fadilah == window
    assert asr.to_dict()["fadilah"]["end"] == window.end.isoformat()

    profile = MosqueProfile.from_preset("Local", method)
    restored = MosqueProfile.from_dict(profile.to_dict())
    day = local_date.date()
    schedule = observer.prayer_schedule(day, day, profile=restored)
    assert schedule[0].prayers.asr.fadilah == window
    assert "preferred ʿAṣr window" in prayer_schedule_to_ics(schedule).replace("\r\n ", "")

    observer.set_prayer_method(method, asr_type=AsrConvention.STANDARD)
    standard = observer.prayer_times().asr
    observer.set_asr_type(AsrConvention.HANAFI)
    hanafi = observer.prayer_times().asr
    assert window.end < standard.when < hanafi.when
    assert standard.fadilah is None
    assert hanafi.fadilah is None


@pytest.mark.parametrize(
    ("local_day", "latitude", "minutes"),
    [
        (date(2026, 12, 21), 55, 79.4),
        (date(2026, 6, 21), -55, 79.4),
        (date(2024, 6, 21), -55, 79.4),
        (date(2026, 3, 22), 55, 54.84),
        (date(2026, 6, 22), 55, 81.44),
        (date(2026, 1, 12), 0, 62),
    ],
)
def test_redness_model_tracks_latitude_and_the_local_season(local_day, latitude, minutes) -> None:
    delay = _western_redness_delay(local_day, latitude)
    assert delay.total_seconds() == pytest.approx(minutes * 60, abs=0.5)


def test_evening_fadilah_uses_elapsed_night_across_dst() -> None:
    zone = ZoneInfo("America/Toronto")
    local_date = datetime(2026, 3, 7, 12, tzinfo=zone)
    observer = ITLocation(latitude=43.65, longitude=-79.38, date=local_date, method="JAFARI")
    tomorrow = ITLocation(latitude=43.65, longitude=-79.38, date=local_date + timedelta(days=1), method="JAFARI")
    windows = observer.sistani_evening_fadilah()
    maghrib, isha = windows["maghrib"], windows["isha"]
    sunset, next_dawn = observer.prayer_times().sunset.when, tomorrow.prayer_times().fajr.when

    assert sunset.utcoffset() != next_dawn.utcoffset()
    assert maghrib.start == observer.prayer_times().maghrib.when
    assert maghrib.end == isha.start == observer.prayer_times().isha.when
    assert maghrib.estimated and isha.estimated and isha.preferred_by is None
    assert "preferred_by" not in isha.to_dict()
    elapsed = next_dawn.astimezone(timezone.utc) - sunset.astimezone(timezone.utc)
    expected_end = sunset.astimezone(timezone.utc) + elapsed / 3
    assert isha.end.astimezone(timezone.utc) == expected_end


@pytest.mark.parametrize("method", ["JAFARI", "TEHRAN"])
def test_shia_presets_schedule_estimated_isha_and_attach_window(method: str) -> None:
    local_date = datetime(2026, 11, 26, 12, tzinfo=ZoneInfo("America/Toronto"))
    observer = ITLocation(latitude=43.65, longitude=-79.38, date=local_date, method=method)
    window = observer.sistani_evening_fadilah()["isha"]
    isha = observer.prayer_times().isha
    assert isha.when == window.start
    assert isha.fadilah == window
    assert isha.to_dict()["fadilah"] == window.to_dict()
    assert isha.status == "estimated"

    profile = MosqueProfile.from_preset("Local", method, offset_minutes={"isha": 5})
    serialized = profile.to_dict()
    restored = MosqueProfile.from_dict(serialized)
    day = local_date.date()
    schedule = observer.prayer_schedule(day, day, profile=restored)
    scheduled = schedule[0].prayers.isha
    assert scheduled.when == window.start + timedelta(minutes=5)
    assert scheduled.fadilah == window
    row = next(DictReader(StringIO(prayer_schedule_to_csv(schedule))))
    assert row["isha"] == scheduled.when.isoformat(timespec="seconds")
    assert "isha" in row["estimated"].split()
    ics = prayer_schedule_to_ics(schedule).replace("\r\n ", "")
    assert f"DTSTART:{scheduled.when.astimezone(timezone.utc):%Y%m%dT%H%M%SZ}" in ics
    assert "estimated preferred ʿIshāʾ window" in ics

    observer.set_custom_prayer_angles(isha_angle=14)
    angle_result = observer.prayer_times().isha
    assert not angle_result.method.use_isha_fadilah
    assert angle_result.fadilah is None
    assert angle_result.status == "calculated"
    assert window.start < angle_result.when
    assert observer.sistani_evening_fadilah()["isha"] == window

    # A profile saved before this setting represented an explicit 14° timetable.
    del serialized["method"]["use_isha_fadilah"]
    legacy = MosqueProfile.from_dict(serialized)
    assert observer.prayer_times(profile=legacy).isha.when == angle_result.when + timedelta(minutes=5)


def test_evening_maghrib_is_estimated_when_dawn_angle_is_unavailable() -> None:
    local_date = datetime(1988, 6, 12, 12, tzinfo=ZoneInfo("Europe/London"))
    with pytest.warns(UserWarning, match="Extreme latitude"):
        observer = ITLocation(latitude=53.75, longitude=-2.48333, date=local_date, method="JAFARI")
    times = observer.prayer_times()
    assert times.maghrib.fallback_rule
    assert times.isha.status == "unavailable"
    assert times.isha.fallback_rule is None
    windows = observer.sistani_evening_fadilah()
    assert set(windows) == {"maghrib"}
    maghrib = windows["maghrib"]
    assert maghrib.estimated
    assert times.sunset.when < maghrib.start < maghrib.end
    assert maghrib.end - times.sunset.when == timedelta(minutes=78, seconds=16)


def test_evening_fadilah_omits_locations_outside_the_model() -> None:
    local_date = datetime(2026, 12, 21, 12, tzinfo=ZoneInfo("Europe/Oslo"))
    observer = ITLocation(latitude=59.91, longitude=10.75, date=local_date, method="JAFARI")
    assert observer.prayer_times().maghrib.when is not None
    assert observer.sistani_evening_fadilah() == {}
    assert observer.prayer_times().isha.when is None


@pytest.mark.parametrize("method", ["JAFARI", "TEHRAN"])
def test_solar_descent_propagates_to_both_windows_and_mosque_schedules(method: str) -> None:
    local = datetime(2021, 3, 11, 12, tzinfo=ZoneInfo("Asia/Jakarta"))
    observer = ITLocation(latitude=-6.303291, longitude=106.872089, date=local, method=method)
    seasonal = observer.prayer_times().isha.when
    observer.set_western_redness_model("solar_descent")
    times = observer.prayer_times()
    windows = observer.sistani_evening_fadilah()
    assert windows["maghrib"].end == windows["isha"].start == times.isha.when
    assert times.isha.when < seasonal
    assert times.isha.status == "estimated"
    delay = (times.isha.when - times.sunset.when).total_seconds() / 60
    assert delay == pytest.approx(50.53, abs=0.05)

    payload = MosqueProfile("Local", observer.method).to_dict()
    profile = MosqueProfile.from_dict(payload)
    observer.set_western_redness_model("seasonal")
    assert observer.prayer_times().isha.when == seasonal
    schedule = observer.prayer_schedule(local.date(), local.date(), profile=profile)
    assert schedule[0].prayers.isha.fadilah == windows["isha"]
    assert schedule[0].prayers.isha.when == times.isha.when
    del payload["method"]["western_redness_model"]
    assert observer.prayer_times(profile=MosqueProfile.from_dict(payload)).isha.when == seasonal


def test_solar_descent_on_a_shallow_summer_night() -> None:
    local = datetime(1988, 6, 12, 12, tzinfo=ZoneInfo("Europe/London"))
    with pytest.warns(UserWarning, match="Extreme latitude"):
        observer = ITLocation(latitude=53.75, longitude=-2.48333, date=local, method="JAFARI")
    seasonal = observer.sistani_evening_fadilah()
    with pytest.warns(UserWarning, match="Extreme latitude"):
        observer.set_western_redness_model("solar_descent")
    # The Sun never reaches the model's reference altitude that night.
    assert observer.sistani_evening_fadilah() == seasonal
    assert set(seasonal) == {"maghrib"}
    assert observer.prayer_times().isha.when is None


def test_redness_model_preserves_an_angle_based_timetable() -> None:
    local = datetime(2021, 3, 11, 12, tzinfo=ZoneInfo("Asia/Jakarta"))
    observer = ITLocation(latitude=-6.303291, longitude=106.872089, date=local, method="ISNA")
    seasonal = observer.sistani_evening_fadilah()
    isha = observer.prayer_times().isha.when
    observer.set_western_redness_model("solar_descent")
    assert observer.sistani_evening_fadilah()["isha"].start != seasonal["isha"].start
    assert observer.prayer_times().isha.when == isha
    with pytest.raises(ValueError, match="western_redness_model"):
        observer.set_western_redness_model("unknown")
    assert observer.method.western_redness_model == "solar_descent"


def test_experimental_blend_uses_the_local_sunset_and_restores_from_profile() -> None:
    local = datetime(2026, 11, 26, 12, tzinfo=ZoneInfo("America/Toronto"))
    observer = ITLocation(latitude=43.65, longitude=-79.38, date=local, method="JAFARI")
    seasonal = observer.sistani_evening_fadilah()["maghrib"].end
    fixed = se.find_proper_suntime(observer.observer_dateinfo, observer.observer_info, "set", Angle(14))
    observer.set_western_redness_model("blend")
    expected = seasonal + (fixed - seasonal) * 0.37
    windows = observer.sistani_evening_fadilah()
    assert expected == windows["maghrib"].end == windows["isha"].start
    profile = MosqueProfile.from_dict(MosqueProfile("Local", observer.method).to_dict())
    observer.set_western_redness_model("seasonal")
    assert observer.prayer_times(profile=profile).isha.when == expected


@pytest.mark.parametrize("method", ["JAFARI", "TEHRAN"])
def test_fajr_estimates_the_end_and_supports_a_local_crossing(method) -> None:
    local = datetime(2026, 3, 8, 12, tzinfo=ZoneInfo("America/Toronto"))
    observer = ITLocation(latitude=43.65, longitude=-79.38, date=local, method=method)
    dawn = observer.prayer_times().fajr
    civil_dawn = se.find_proper_suntime(observer.observer_dateinfo, observer.observer_info, "rise", Angle(6))
    assert dawn.fadilah.start == dawn.when
    assert dawn.fadilah.end == civil_dawn < observer.prayer_times().sunrise.when
    assert dawn.fadilah.estimated and dawn.status == "calculated"
    assert dawn.fadilah.to_dict()["end_sign"] == "whole_horizon_lit"
    observer.set_fajr_fadilah_end_angle(10)
    result = observer.prayer_times()
    target = se.find_proper_suntime(observer.observer_dateinfo, observer.observer_info, "rise", Angle(10))
    assert result.fajr.when == dawn.when < target < result.sunrise.when
    assert result.fajr.fadilah.end == target and result.fajr.fadilah.estimated
    assert result.fajr.status == "calculated"
    saved = MosqueProfile("Local", observer.method, offset_minutes={"fajr": 3}).to_dict()
    shifted = observer.prayer_times(profile=MosqueProfile.from_dict(saved)).fajr
    assert shifted.when == dawn.when + timedelta(minutes=3)
    assert shifted.fadilah == result.fajr.fadilah
    del saved["method"]["fajr_fadilah_end_degrees"]
    assert observer.prayer_times(profile=MosqueProfile.from_dict(saved)).fajr.fadilah == dawn.fadilah
    with pytest.raises(ValueError, match="smaller than the Fajr angle"):
        observer.set_fajr_fadilah_end_angle(observer.method.fajr_angle.decimal)
    observer.set_fajr_fadilah_end_angle(None)
    assert observer.prayer_times().fajr.fadilah == dawn.fadilah
    observer.set_prayer_method("ISNA")
    assert observer.prayer_times().fajr.fadilah is None


def test_changing_fajr_angle_updates_the_opening_only() -> None:
    local = datetime(2026, 11, 26, 12, tzinfo=ZoneInfo("America/Toronto"))
    observer = ITLocation(latitude=43.65, longitude=-79.38, date=local, method="JAFARI")
    original = observer.prayer_times().fajr.fadilah
    observer.set_custom_prayer_angles(fajr_angle=18)
    changed = observer.prayer_times().fajr.fadilah
    assert changed.start < original.start < original.end == changed.end
    observer.set_custom_prayer_angles(fajr_angle=5)
    assert observer.prayer_times().fajr.when > original.end
    assert observer.prayer_times().fajr.fadilah is None


def test_absent_local_dawn_and_redness_crossing() -> None:
    local = datetime(1988, 6, 12, 12, tzinfo=ZoneInfo("Europe/London"))
    with pytest.warns(UserWarning, match="Extreme latitude"):
        observer = ITLocation(latitude=53.75, longitude=-2.48, date=local, method="JAFARI")
        seasonal = observer.sistani_evening_fadilah()
        observer.set_western_redness_model("blend")
    assert observer.prayer_times().fajr.fadilah is None
    assert observer.sistani_evening_fadilah() == seasonal


def test_blend_uses_the_crossing_after_sunset_when_it_falls_tomorrow() -> None:
    local = datetime(1988, 6, 1, 12, tzinfo=ZoneInfo("Europe/London"))
    with pytest.warns(UserWarning, match="Extreme latitude"):
        observer = ITLocation(latitude=53.75, longitude=-2.48, date=local, method="JAFARI")
        seasonal = observer.sistani_evening_fadilah()["maghrib"].end
        observer.set_western_redness_model("blend")
    crossing = se.find_proper_suntime(
        te.observer_datetime_info(local + timedelta(days=1)), observer.observer_info, "set", Angle(14),
    )
    assert crossing.date() == local.date() + timedelta(days=1)
    assert observer.sistani_evening_fadilah()["maghrib"].end == seasonal + (crossing - seasonal) * 0.37
