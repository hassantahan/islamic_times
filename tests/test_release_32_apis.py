from datetime import date, datetime, timedelta, timezone
from zoneinfo import ZoneInfo

import pytest

from islamic_times import astro_core, moon_equations
from islamic_times.islamic_times import ITLocation
from islamic_times.prayer_times import prayer_method_catalog
from islamic_times.time_equations import civil_hijri_to_gregorian, gregorian_to_civil_hijri


def test_civil_hijri_round_trips_one_cycle_and_supported_limit() -> None:
    epoch = date(622, 7, 19)
    assert civil_hijri_to_gregorian(1, 1, 1) == epoch
    assert civil_hijri_to_gregorian(1446, 1, 1) == date(2024, 7, 8)
    for value in [*(epoch + timedelta(days=i) for i in range(10631)), date.max]:
        hijri = gregorian_to_civil_hijri(value)
        assert civil_hijri_to_gregorian(hijri.hijri_year, hijri.hijri_month, hijri.hijri_day) == value
    assert civil_hijri_to_gregorian(2, 12, 30) == epoch + timedelta(days=708)
    with pytest.raises(ValueError, match="day"):
        civil_hijri_to_gregorian(1, 12, 30)
    with pytest.raises(ValueError, match="epoch"):
        gregorian_to_civil_hijri(epoch - timedelta(days=1))


def test_prayer_schedule_reuses_the_daily_solar_calculation(monkeypatch: pytest.MonkeyPatch) -> None:
    zone = ZoneInfo("America/Toronto")
    location_kwargs = {"latitude": 43.65, "longitude": -79.38, "method": "ISNA"}
    start = date(2026, 3, 7)
    dates = (start, start + timedelta(days=1))
    expected = [
        ITLocation(**location_kwargs, date=datetime(day.year, day.month, day.day, 12, tzinfo=zone)).prayer_times()
        for day in dates
    ]
    observer = ITLocation(**location_kwargs, date=datetime(2026, 1, 1), auto_calculate=False)

    def fail_if_moon_is_computed(*args: object) -> None:
        pytest.fail("A prayer-only schedule must not calculate lunar positions.")

    monkeypatch.setattr(moon_equations, "moonpos", fail_if_moon_is_computed)
    schedule = observer.prayer_schedule(*dates, time_zone="America/Toronto")
    assert [entry.prayers for entry in schedule] == expected
    assert tuple(entry.date for entry in schedule) == dates
    assert schedule[1].to_dict()["date"] == "2026-03-08"
    with pytest.raises(ValueError):
        observer.prayer_times()
    assert any(method.keys and method.keys[0] == "MAKKAH" for method in prayer_method_catalog())


def test_upcoming_lunar_phases_are_chronological_and_use_event_offsets() -> None:
    zone = ZoneInfo("America/Toronto")
    start = datetime(2026, 11, 24, 12, tzinfo=zone)
    observer = ITLocation(date=start, auto_calculate=False)
    phases = observer.upcoming_moon_phases(5)
    assert [name for name, _ in phases] == ["Last Quarter", "New Moon", "First Quarter", "Full Moon", "Last Quarter"]
    instants = [event.astimezone(timezone.utc) for _, event in phases]
    assert start < instants[0] < instants[1] < instants[2] < instants[3] < instants[4]
    # USNO Dates of Primary Phases for 2026, rounded to the minute.
    assert abs((instants[0] - datetime(2026, 12, 1, 6, 8, tzinfo=timezone.utc)).total_seconds()) < 60
    assert observer.upcoming_moon_phases(1, after=phases[0][1]) == phases[1:2]
    march = observer.upcoming_moon_phases(1, after=datetime(2026, 3, 7, 12, tzinfo=zone))[0]
    assert march[1].utcoffset() == timedelta(hours=-4)


def test_detailed_crescent_forecast_agrees_with_the_same_conjunction_in_legacy_api() -> None:
    local = datetime(2025, 4, 27, 12, tzinfo=ZoneInfo("America/Toronto"))
    observer = ITLocation(latitude=43.65, longitude=-79.38, date=local)
    conjunction = astro_core.next_phases_of_moon_utc(local)[0].replace(tzinfo=timezone.utc)
    legacy = observer.visibilities()
    forecast = observer.crescent_forecast(conjunction)
    assert forecast.conjunction == conjunction
    assert forecast.criterion == legacy.criterion
    assert tuple(entry.date for entry in forecast.evenings) == tuple(local.date() + timedelta(days=i) for i in range(3))
    assert tuple(entry.score for entry in forecast.evenings) == pytest.approx(legacy.q_values)
    for entry in forecast.evenings:
        assert entry.best_time is not None and entry.moonset is not None and entry.sunset is not None
        assert entry.sunset < entry.best_time < entry.moonset
        assert entry.lag_minutes > 0
    assert forecast.to_dict()["evenings"][0]["status"] == "evaluated"


def test_crescent_without_sunset_has_no_numeric_score() -> None:
    zone = ZoneInfo("Europe/Oslo")
    observer = ITLocation(latitude=78.22, longitude=15.65, date=datetime(2026, 6, 15, tzinfo=zone),
                          auto_calculate=False)
    conjunction = datetime(2026, 6, 15, 3, tzinfo=timezone.utc)
    for entry in observer.crescent_forecast(conjunction).evenings:
        assert entry.status in ("no_sunset", "no_sunset_or_moonset")
        assert entry.sunset is None and entry.best_time is None and entry.score is None
