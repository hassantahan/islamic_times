from datetime import datetime, time
from zoneinfo import ZoneInfo

import pytest

from islamic_times.islamic_times import ITLocation


@pytest.mark.parametrize(
    "day, latitude, longitude, zone, reference",
    [
        ("2025-04-28", 51.5074, -0.1278, "Europe/London", ("05:38", "20:19", "05:30", "22:01")),
        ("2026-11-24", -36.85, 174.76, "Pacific/Auckland", ("05:58", "20:18", "20:04", "04:55")),
    ],
)
def test_rise_and_set_against_usno(
    day: str, latitude: float, longitude: float, zone: str, reference: tuple[str, ...],
) -> None:
    # USNO Complete Sun and Moon Data for One Day, at sea level and rounded to the minute.
    # https://aa.usno.navy.mil/data/api#rstt
    local_noon = datetime.fromisoformat(day).replace(hour=12, tzinfo=ZoneInfo(zone))
    location = ITLocation(latitude=latitude, longitude=longitude, date=local_noon)
    sun, moon = location.sun(), location.moon()
    for actual, expected in zip((sun.sunrise, sun.sunset, moon.moonrise, moon.moonset), reference):
        assert isinstance(actual, datetime)
        target = datetime.combine(local_noon.date(), time.fromisoformat(expected), local_noon.tzinfo)
        assert abs((actual - target).total_seconds()) <= 75
