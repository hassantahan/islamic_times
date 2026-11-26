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


@pytest.mark.parametrize(
    "day, latitude, longitude, zone, rise, setting",
    [
        (19, -5.0775, 114.6014, "Asia/Jakarta", "05:12:34", "17:18:07"),
        (20, -5.0775, 114.6014, "Asia/Jakarta", "05:12:04", "17:17:54"),
        (21, -5.0775, 114.6014, "Asia/Jakarta", "05:11:35", "17:17:41"),
        (20, 1.872, -157.4278, "Pacific/Kiritimati", "06:19:46", "18:26:31"),
        (20, -13.83, -171.75, "UTC", "17:17:27", "05:22:52"),
    ],
)
def test_solar_events_belong_to_requested_civil_day(day, latitude, longitude, zone, rise, setting):
    # Independent PyEphem 4.2.1, solar centre at -50 arcminutes, zero refraction.
    noon = datetime(2020, 9, day, 12, tzinfo=ZoneInfo(zone))
    observer = ITLocation(latitude=latitude, longitude=longitude, date=noon, method="JAFARI")
    sun = observer.sun()
    for actual, clock in ((sun.sunrise, rise), (sun.sunset, setting)):
        expected = datetime.combine(noon.date(), time.fromisoformat(clock), noon.tzinfo)
        assert abs((actual - expected).total_seconds()) < 10
    if zone == "Asia/Jakarta":
        evening = observer.sistani_evening_fadilah()
        assert evening["maghrib"].end == evening["isha"].start == observer.prayer_times().isha.when
        assert evening["isha"].start.date() == noon.date()


def test_sunset_near_midnight_on_a_clock_change_day():
    # Express equatorial Pacific events on the Toronto clock to put sunset
    # in the hour gained or lost. PyEphem gives two sunsets on the long day.
    zone = ZoneInfo("America/Toronto")
    spring = ITLocation(latitude=0, longitude=-155, date=datetime(2020, 3, 8, 12, tzinfo=zone))
    assert not isinstance(spring.sun().sunset, datetime)
    autumn = ITLocation(latitude=0, longitude=-155, date=datetime(2020, 11, 1, 12, tzinfo=zone))
    first = datetime(2020, 11, 1, 0, 7, tzinfo=zone)
    assert abs(autumn.sun().sunset.timestamp() - first.timestamp()) < 10
    for observer in (spring, autumn):
        forecast = observer.crescent_forecast(observer.observer_dateinfo.date, days=1)
        expected = observer.sun().sunset
        assert forecast.evenings[0].sunset == (expected if isinstance(expected, datetime) else None)
