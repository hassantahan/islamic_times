from datetime import datetime
from math import cos, sin
from zoneinfo import ZoneInfo

import pytest

from islamic_times import sun_equations as se
from islamic_times import time_equations as te
from islamic_times.islamic_times import ITLocation
from islamic_times.it_dataclasses import AsrConvention
from islamic_times.prayer_calendar import prayer_schedule_to_ics
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
