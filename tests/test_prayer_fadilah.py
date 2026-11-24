from datetime import datetime
from math import cos, sin
from zoneinfo import ZoneInfo

import pytest

from islamic_times import sun_equations as se
from islamic_times import time_equations as te
from islamic_times.islamic_times import ITLocation


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
