from __future__ import annotations

import csv
import json
from dataclasses import replace
from datetime import date, datetime, timedelta, timezone
from io import StringIO
from zoneinfo import ZoneInfo

import pytest

from islamic_times.islamic_times import ITLocation
from islamic_times.it_dataclasses import PRAYER_KEYS, PrayerDay
from islamic_times.prayer_calendar import prayer_schedule_to_csv, prayer_schedule_to_ics
from islamic_times.prayer_profiles import MosqueProfile
from islamic_times.prayer_times import get_prayer_method


@pytest.fixture
def toronto() -> ITLocation:
    return ITLocation(
        latitude=43.65, longitude=-79.38,
        date=datetime(2026, 3, 8, 12, tzinfo=ZoneInfo("America/Toronto")),
        method="ISNA",
    )


def _ical_events(calendar: str) -> list[dict[str, str]]:
    """Read unfolded VEVENT properties for assertions."""
    unfolded = calendar.replace("\r\n ", "")
    return [
        dict(line.split(":", 1) for line in event.split("\r\n", 1)[1].split("END:VEVENT")[0].splitlines())
        for event in unfolded.split("BEGIN:VEVENT")[1:]
    ]


def test_profile_roundtrip_preserves_settings_and_snapshots_offsets() -> None:
    method = replace(get_prayer_method(" mecca "), asr_type=1, extreme_lats="ONESEVENTH")
    offsets = {"fajr": -3, "maghrib": 4}
    profile = MosqueProfile("Al-Nūr", method, offset_minutes=offsets)
    offsets["maghrib"] = 0

    payload = json.loads(json.dumps(profile.to_dict()))
    restored = MosqueProfile.from_dict(payload)
    assert restored == profile
    assert restored.method.ramadan_isha_offset_minutes == 120
    assert restored.offset_minutes["maghrib"] == 4
    payload["offset_minutes"]["maghrib"] = 10
    assert restored.offset_minutes["maghrib"] == 4


@pytest.mark.parametrize("invalid", [{"sunset": 3}, {"isha": True}, {"fajr": 2.5}, {"asr": 10**20}])
def test_profile_rejects_invalid_offsets(invalid: dict[str, int | float]) -> None:
    with pytest.raises((TypeError, ValueError)):
        MosqueProfile.from_preset("Local", "ISNA", offset_minutes=invalid)


def test_profile_rejects_invalid_serialized_method() -> None:
    payload = MosqueProfile.from_preset("Local", "ISNA").to_dict()
    payload["method"]["fajr_angle_degrees"] = float("nan")
    with pytest.raises(ValueError, match="fajr_angle"):
        MosqueProfile.from_dict(payload)

    payload["schema_version"] = "nonexistent"
    with pytest.raises(ValueError, match="schema version"):
        MosqueProfile.from_dict(payload)


def test_profile_offsets_apply_after_calculation_without_changing_the_observer(toronto: ITLocation) -> None:
    original = toronto.prayer_times()
    fixed_method = MosqueProfile.from_preset("Makkah base", "MAKKAH")
    profile = MosqueProfile.from_preset("Local", "MAKKAH", offset_minutes={"maghrib": 5, "asr": -4})
    base = toronto.prayer_times(profile=fixed_method)
    adjusted = toronto.prayer_times(profile=profile)

    assert adjusted.method == profile.method
    assert adjusted.maghrib.when - base.maghrib.when == timedelta(minutes=5)
    assert adjusted.asr.when - base.asr.when == timedelta(minutes=-4)
    assert adjusted.isha == base.isha
    assert adjusted.midnight == base.midnight
    assert adjusted.sunset == base.sunset
    assert adjusted.to_dict()["mosque_profile_name"] == "Local"
    assert adjusted.to_dict()["prayers"]["asr"]["offset_minutes"] == -4
    assert "mosque_profile_name" not in original.to_dict()
    assert "offset_minutes" not in original.asr.to_dict()
    assert toronto.prayer_times() is original
    civil_day = toronto.dates_times().date.date()
    assert toronto.prayer_schedule(civil_day, civil_day, profile=profile)[0].prayers == adjusted


def test_daily_profile_in_manual_mode_requires_only_current_astronomy() -> None:
    observer = ITLocation(date=datetime(2026, 3, 8, tzinfo=timezone.utc), auto_calculate=False)
    profile = MosqueProfile.from_preset("Local", "ISNA")
    with pytest.raises(ValueError, match="calculate_astro"):
        observer.prayer_times(profile=profile)
    observer.calculate_astro()
    assert observer.prayer_times(profile=profile).mosque_profile_name == "Local"
    with pytest.raises(ValueError, match="calculate_prayer_times"):
        observer.prayer_times()


def test_csv_preserves_dst_offset_and_rollover(toronto: ITLocation) -> None:
    start, end = date(2026, 3, 7), date(2026, 3, 8)
    profile = MosqueProfile.from_preset("Local", "ISNA", offset_minutes={"fajr": -420})
    original = toronto.prayer_schedule(start, end)
    adjusted = toronto.prayer_schedule(start, end, profile=profile)

    csv_text = prayer_schedule_to_csv(adjusted)
    first, second = list(csv.DictReader(StringIO(csv_text)))
    assert (first["date"], second["date"]) == (start.isoformat(), end.isoformat())
    assert (first["isha"][-6:], second["isha"][-6:]) == ("-05:00", "-04:00")
    assert second["hijri_date"] == "1447-09-19"
    assert second["fajr"] == adjusted[1].prayers.fajr.when.isoformat(timespec="seconds")
    # The adjusted Fajr crosses both the change in offset and civil midnight.
    shifted = datetime.fromisoformat(second["fajr"])
    baseline = original[1].prayers.fajr.when
    assert shifted.date() == start and shifted.utcoffset() == timedelta(hours=-5)
    assert shifted.astimezone(timezone.utc) - baseline.replace(microsecond=0).astimezone(timezone.utc) == (
        timedelta(minutes=-420)
    )
    assert csv_text.endswith("\r\n")


def test_ics_preserves_instants_and_stable_ids_across_revisions(toronto: ITLocation) -> None:
    start, end = date(2026, 3, 7), date(2026, 3, 8)
    schedule = toronto.prayer_schedule(start, end)
    stamp = datetime(2026, 2, 1, tzinfo=timezone.utc)
    ics = prayer_schedule_to_ics(schedule, calendar_id="toronto-prayers", generated_at=stamp)
    events = _ical_events(ics)
    assert len(events) == 10
    assert len({event["UID"] for event in events}) == 10
    assert all(event["DTSTAMP"] == "20260201T000000Z" for event in events)
    expected = [
        getattr(day.prayers, key).when.astimezone(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        for day in schedule for key in PRAYER_KEYS
    ]
    assert [event["DTSTART"] for event in events] == expected

    profile = MosqueProfile.from_preset("Local", "ISNA", offset_minutes={"isha": 3})
    revised = toronto.prayer_schedule(start, end, profile=profile)
    revision = _ical_events(prayer_schedule_to_ics(revised, calendar_id="toronto-prayers"))
    for event, updated in zip(events, revision, strict=True):
        assert event["UID"] == updated["UID"]
        if updated["SUMMARY"] == "ʿIshāʾ":
            assert event["DTSTART"] != updated["DTSTART"]
            assert "Mosque offset: +3 minutes." in updated["DESCRIPTION"]
        else:
            assert event["DTSTART"] == updated["DTSTART"]


def test_ics_escapes_and_folds_unicode_calendar_text(toronto: ITLocation) -> None:
    day = date(2026, 3, 8)
    title = "Ṣalāh" * 20 + ",;\\\r\nBEGIN:VEVENT"
    schedule = toronto.prayer_schedule(day, day)
    ics = prayer_schedule_to_ics(schedule, calendar_name=title)
    assert all(len(line.encode("utf-8")) <= 75 for line in ics.split("\r\n"))
    assert any(line.startswith(" ") for line in ics.split("\r\n"))
    assert ics.count("\r\nBEGIN:VEVENT\r\n") == 5
    assert "Ṣalāh" * 20 + "\\,\\;\\\\\\nBEGIN:VEVENT" in ics.replace("\r\n ", "")


def test_exports_handle_estimated_and_unavailable_prayers() -> None:
    civil_day = date(2026, 6, 21)
    observer = ITLocation(
        latitude=78.22, longitude=15.65, date=datetime(2026, 6, 21, 12, tzinfo=ZoneInfo("Europe/Oslo")),
        auto_calculate=False,
    )
    profile = MosqueProfile.from_preset("Local", "MAKKAH", offset_minutes={"isha": 5})
    with pytest.warns(UserWarning, match="Extreme latitude"):
        estimated = observer.prayer_schedule(civil_day, civil_day, profile=profile)[0]
    assert "isha" in next(csv.DictReader(StringIO(prayer_schedule_to_csv([estimated]))))["estimated"].split()
    isha = next(event for event in _ical_events(prayer_schedule_to_ics([estimated])) if event["SUMMARY"] == "ʿIshāʾ")
    assert "Estimated using ANGLEBASED." in isha["DESCRIPTION"]

    direct = replace(profile, method=replace(profile.method, extreme_lats="NONE"))
    with pytest.warns(UserWarning, match="Extreme latitude"):
        unavailable_prayers = observer.prayer_schedule(civil_day, civil_day, profile=direct)[0]
    assert unavailable_prayers.prayers.isha.offset_minutes == 0
    assert next(csv.DictReader(StringIO(prayer_schedule_to_csv([unavailable_prayers]))))["isha"] == ""
    events = _ical_events(prayer_schedule_to_ics([unavailable_prayers]))
    assert all(event["SUMMARY"] != "ʿIshāʾ" for event in events)
    assert any(event["SUMMARY"] == "Ẓuhr" for event in events)


def test_exports_reject_naive_times(toronto: ITLocation) -> None:
    day = date(2026, 3, 8)
    schedule = toronto.prayer_schedule(day, day)[0]
    prayer = replace(schedule.prayers.fajr, time=datetime(2026, 3, 8, 6))
    naive = PrayerDay(day, schedule.hijri, replace(schedule.prayers, fajr=prayer))
    for export in (prayer_schedule_to_csv, prayer_schedule_to_ics):
        with pytest.raises(ValueError, match="timezone-aware"):
            export([naive])
