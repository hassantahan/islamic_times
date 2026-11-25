"""Export daily prayer schedules as CSV or iCalendar."""

import csv
import json
from collections.abc import Iterable
from datetime import datetime, timezone
from hashlib import sha256
from io import StringIO

from islamic_times.it_dataclasses import PRAYER_EVENT_KEYS, PRAYER_KEYS, Prayer, PrayerDay


def _aware(value: datetime) -> datetime:
    """Require a timestamp with a UTC offset so exports have unambiguous instants."""
    if not isinstance(value, datetime) or value.utcoffset() is None:
        raise ValueError("Calendar exports require timezone-aware timestamps.")
    return value


def prayer_schedule_to_csv(schedule: Iterable[PrayerDay]) -> str:
    """Return a CSV timetable with one row per observer-local civil date.

    Parameters
    ----------
    schedule : Iterable[PrayerDay]
        The days supplied by ``ITLocation.prayer_schedule()``.

    Returns
    -------
    str
        Gregorian and Hijri dates, local ISO timestamps with UTC offsets,
        and space-separated event keys in the estimated column. Empty time
        fields represent unavailable events.
    """
    output = StringIO(newline="")
    writer = csv.writer(output, lineterminator="\r\n")
    writer.writerow(("date", "hijri_date", *PRAYER_EVENT_KEYS, "estimated"))
    for day in schedule:
        hijri = day.hijri
        hijri_date = f"{hijri.hijri_year:04d}-{hijri.hijri_month:02d}-{hijri.hijri_day:02d}"
        times: list[str] = []
        estimated: list[str] = []
        for key in PRAYER_EVENT_KEYS:
            prayer: Prayer = getattr(day.prayers, key)
            when = prayer.when
            times.append(_aware(when).isoformat(timespec="seconds") if when is not None else "")
            if prayer.status == "estimated":
                estimated.append(key)
        writer.writerow((day.date.isoformat(), hijri_date, *times, " ".join(estimated)))
    return output.getvalue()


def _calendar_text(value: str) -> str:
    """Escape iCalendar TEXT values, including untrusted names and descriptions."""
    if any((ord(char) < 32 and char not in "\t\r\n") or ord(char) == 127 for char in value):
        raise ValueError("Calendar text contains an unsupported control character.")
    value = value.replace("\r\n", "\n").replace("\r", "\n")
    return value.replace("\\", "\\\\").replace("\n", "\\n").replace(";", "\\;").replace(",", "\\,")


def _fold_line(line: str) -> str:
    """Fold at 75 UTF-8 octets, including the continuation's leading space."""
    parts: list[str] = []
    start, used = 0, 0
    for index, char in enumerate(line):
        length = len(char.encode("utf-8"))
        if used + length > 75:
            parts.append(line[start:index])
            start, used = index, 1
        used += length
    parts.append(line[start:])
    return "\r\n ".join(parts)


def _utc_timestamp(value: datetime) -> str:
    """Format an aware instant as a UTC iCalendar DATE-TIME."""
    utc = _aware(value).astimezone(timezone.utc)
    return f"{utc.year:04d}{utc:%m%dT%H%M%SZ}"


def _event_lines(
    prayer: Prayer, start: str, stamp: str, uid: str, profile_name: str | None,
) -> list[str]:
    """Serialize one prayer at its adjusted time and attach provenance when relevant."""
    lines = [
        "BEGIN:VEVENT",
        f"UID:{uid}",
        f"DTSTAMP:{stamp}",
        f"DTSTART:{start}",
        f"SUMMARY:{_calendar_text(prayer.name)}",
        "TRANSP:TRANSPARENT",
    ]
    notes = []
    if profile_name is not None:
        notes.append(f"Mosque profile: {profile_name}")
    if prayer.fadilah is not None:
        notes.append("Timed from the start of Sayyid Sistani’s preferred ʿAṣr window.")
    if prayer.fallback_rule is not None:
        notes.append(f"Estimated using {prayer.fallback_rule}.")
    if prayer.offset_minutes:
        notes.append(f"Mosque offset: {prayer.offset_minutes:+d} minutes.")
    if notes:
        lines.append("DESCRIPTION:" + _calendar_text("\n".join(notes)))
    lines.append("END:VEVENT")
    return lines


def prayer_schedule_to_ics(
    schedule: Iterable[PrayerDay],
    *,
    calendar_name: str = "Prayer Times",
    calendar_id: str | None = None,
    generated_at: datetime | None = None,
) -> str:
    """Return an iCalendar document containing the available daily prayers.

    Parameters
    ----------
    schedule : Iterable[PrayerDay]
        The days supplied by ``ITLocation.prayer_schedule()``.
    calendar_name : str, optional
        The calendar title.
    calendar_id : str, optional
        A unique identifier for a timetable, such as a mosque and city slug.
        Reuse it to keep event IDs stable across revised exports. By default
        IDs are stable across exports of the same prayer times and titles.
    generated_at : datetime, optional
        An aware publication timestamp. Defaults to the current UTC time.

    Notes
    -----
    Offsets and daylight saving changes are preserved by exporting each
    event's instant in UTC. Encode the returned text as UTF-8 and write it
    without newline translation.

    Raises
    ------
    ValueError
        If there are no available prayers or a timestamp is naive.
    """
    if not isinstance(calendar_name, str) or not calendar_name.strip():
        raise ValueError("'calendar_name' must be a non-empty string.")
    if calendar_id is not None and (not isinstance(calendar_id, str) or not calendar_id.strip()):
        raise ValueError("'calendar_id' must be a non-empty string.")
    stamp = _utc_timestamp(datetime.now(timezone.utc) if generated_at is None else generated_at)
    lines = [
        "BEGIN:VCALENDAR",
        "VERSION:2.0",
        "PRODID:-//Islamic Times//Prayer Calendar//EN",
        "CALSCALE:GREGORIAN",
        f"NAME:{_calendar_text(calendar_name)}",
    ]
    events: list[tuple[str, list[str]]] = []
    for day in schedule:
        for key in PRAYER_KEYS:
            prayer: Prayer = getattr(day.prayers, key)
            when = prayer.when
            if when is not None:
                start = _utc_timestamp(when)
                profile_name = day.prayers.mosque_profile_name
                identity: tuple[str | None, ...] = (day.date.isoformat(), key, calendar_id)
                if calendar_id is None:
                    identity += (calendar_name, profile_name, start)
                digest = sha256(json.dumps(identity).encode("utf-8")).hexdigest()[:32]
                uid = f"{digest}@islamictimes.org"
                events.append((start, _event_lines(prayer, start, stamp, uid, profile_name)))
    if not events:
        raise ValueError("There are no available prayer times to export.")
    for _, event in sorted(events, key=lambda item: item[0]):
        lines.extend(event)
    lines.append("END:VCALENDAR")
    return "\r\n".join(_fold_line(line) for line in lines) + "\r\n"
