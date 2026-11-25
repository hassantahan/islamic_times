# Calculation APIs

## Prayer schedules

```python
from datetime import date, datetime

from islamic_times.islamic_times import ITLocation
from islamic_times.prayer_times import prayer_method_catalog

location = ITLocation(
    latitude=43.65,
    longitude=-79.38,
    date=datetime(2026, 3, 1),
    method="ISNA",
    auto_calculate=False,
)
month = location.prayer_schedule(
    date(2026, 3, 1),
    date(2026, 3, 31),
    time_zone="America/Toronto",
)
print(month[0].prayers.fajr.when)
print(month[0].prayers.fajr.status)
print([method.keys[0] for method in prayer_method_catalog()])
```

The start and end dates are inclusive. Each day uses the selected timezone's
offset on that date and shares the solar calculation used by `prayer_times()`.
The schedule uses local noon as the reference instant and leaves the
`ITLocation` date and cached results intact. Use `auto_calculate=False` to
avoid computing lunar positions when creating a schedule-only observer.

`Prayer.when` is a datetime or `None`. `Prayer.status` is `calculated`,
`estimated`, or `unavailable`. Estimated events carry the name of their
`fallback_rule`. The existing `Prayer.time` field and serialized time payload
still support existing consumers. Selecting `NONE` as a high-latitude rule
preserves available direct calculations and leaves missing events unavailable.

The Makkah preset exposes its ordinary 90-minute and Ramadan 120-minute
ʿIshāʾ intervals as method metadata. They are applied to a resolved Maghrib.
An explicit custom ʿIshāʾ angle switches the method to angle-based calculation.

## Mosque profiles and calendar export

```python
import json
from datetime import date, datetime
from pathlib import Path

from islamic_times.islamic_times import ITLocation
from islamic_times.prayer_calendar import prayer_schedule_to_csv, prayer_schedule_to_ics
from islamic_times.prayer_profiles import MosqueProfile

profile = MosqueProfile.from_preset(
    "Local mosque",
    "ISNA",
    offset_minutes={"fajr": -2, "zuhr": 5, "maghrib": 3},
)
profile_file = Path("mosque-profile.json")
profile_file.write_text(json.dumps(profile.to_dict(), ensure_ascii=False, indent=2), encoding="utf-8")
profile = MosqueProfile.from_dict(json.loads(profile_file.read_text(encoding="utf-8")))

location = ITLocation(
    latitude=43.65,
    longitude=-79.38,
    date=datetime(2026, 3, 1),
    auto_calculate=False,
)
schedule = location.prayer_schedule(
    date(2026, 3, 1),
    date(2026, 3, 31),
    time_zone="America/Toronto",
    profile=profile,
)
Path("prayer-times.csv").write_bytes(prayer_schedule_to_csv(schedule).encode("utf-8"))
calendar = prayer_schedule_to_ics(schedule, calendar_name=profile.name, calendar_id="local-mosque-toronto")
Path("prayer-times.ics").write_bytes(calendar.encode("utf-8"))
```

Profiles can also wrap a custom `PrayerMethod`, such as the
`method` returned by `location.prayer_times()` after configuring an observer.
Use the same profile across locations and date ranges. For a single day,
`location.prayer_times(profile=profile)` uses the observer's current astronomy.
When `auto_calculate=False`, call `location.calculate_astro()` first for that
single-day calculation.

Offsets are signed whole minutes for `fajr`, `zuhr`, `asr`, `maghrib`, and
`isha`. Each offset is applied independently after the base timetable is
calculated. A time can cross midnight or a timezone transition. Each adjusted
prayer exposes its `offset_minutes` alongside its calculation status.

CSV includes Gregorian and Hijri dates, local ISO timestamps with UTC offsets,
and the keys of estimated events. An empty time field indicates that an event
is unavailable. The row's date identifies the original timetable day and the
timestamp includes the adjusted event's actual date.

iCalendar includes each available daily prayer at its actual UTC instant.
Estimated prayers carry the fallback rule in their descriptions. Use the same
`calendar_id` when exporting revisions of a timetable to preserve event IDs.
The exported text uses CRLF line endings for compatibility with
[RFC 5545](https://www.rfc-editor.org/rfc/rfc5545.html).

## Faḍīlah windows for Ẓuhr and ʿAṣr

`JAFARI` and `TEHRAN` use the start of Sayyid Sistani’s preferred ʿAṣr
window for the timetable entry. Read `observer.prayer_times().asr.fadilah`
for the whole window. Under Sistani’s ruling, ʿAṣr can be offered as soon as
Ẓuhr has been prayed. [Prayer timing, ruling 717](https://www.sistani.org/english/book/48/2209/).

| `asr_type` | Timetable convention | Source |
| --- | --- | --- |
| `0` | One object height added to the noon shadow, the Standard convention | [Mālikī, Shāfiʿī, and Ḥanbalī rule](https://islamweb.net/ar/fatwa/print.php?id=117832) |
| `1` | Two object heights added to the noon shadow, the Ḥanafī convention | [Al-Marghīnānī, al-Hidāyah](https://parallelquran.com/fiqh/al-hidaya/9/) |
| `2` | Start of Sistani’s preferred window | [Sistani, rulings 966–967](https://www.sistani.org/persian/book/26575/6324/) |

Omit `asr_type` to use the preset. The selector also accepts the corresponding
`AsrConvention` enum from `islamic_times.it_dataclasses`. CSV and iCalendar use
the selected entry, with a mosque offset if configured. Serialized mosque profiles
retain their selected convention. The window bounds always describe the
astronomical estimate before offsets.

`observer.sistani_daylight_fadilah()` returns the available preferred windows
for the observer’s current date. Each result includes `start`, `preferred_by`
and `end`, with timezone-aware datetimes. The accessor works with any prayer preset.

Legal noon is the elapsed-time midpoint between sunrise and sunset. The
calculation uses the observer’s geometric solar position to estimate the
shadow at legal noon and finds when the afternoon shadow grows by 2/7, 4/7
and 6/7 of the height of a vertical object. Ẓuhr runs from noon to the 4/7
crossing, preferably by 2/7. ʿAṣr runs from the 2/7 to the 6/7 crossing,
preferably by 4/7.

The calculation covers ordinary conditions with a level horizon. A window
is omitted if the underlying sunrise, sunset or shadow crossing is unavailable.
The rules are from Sayyid Sistani’s office.
[Legal noon, rulings 953–954](https://www.sistani.org/persian/book/26575/6319/).
[Preferred windows, rulings 966–967](https://www.sistani.org/persian/book/26575/6324/).

## Civil Hijri conversion

```python
from datetime import date

from islamic_times.time_equations import civil_hijri_to_gregorian, gregorian_to_civil_hijri

hijri = gregorian_to_civil_hijri(date(2024, 7, 8))
original = civil_hijri_to_gregorian(hijri.hijri_year, hijri.hijri_month, hijri.hijri_day)
```

Both directions use the tabular 30-year type II calendar and the Friday epoch
of 19 July 622 in the proleptic Gregorian calendar. Days begin at civil midnight.
The converter calculates tabular dates. Observed month starts depend on the
calendar or community used. Dates before the epoch and invalid month days
raise `ValueError`. The legacy tuple converter `gregorian_to_hijri()` remains
available with its existing behavior. Existing prayer and mapper date labels,
including the Makkah Ramadan rule, continue to use that legacy conversion.

## Lunar phases and crescent forecasts

```python
from datetime import datetime, timezone

from islamic_times.islamic_times import ITLocation

observer = ITLocation(
    latitude=43.65,
    longitude=-79.38,
    date=datetime(2026, 11, 24, tzinfo=timezone.utc),
    find_local_tz=True,
    auto_calculate=False,
)
phases = observer.upcoming_moon_phases()
conjunction = next(when for name, when in phases if name == "New Moon")
forecast = observer.crescent_forecast(conjunction, days=3, criterion=1)
print(forecast.to_dict())
```

`upcoming_moon_phases()` returns the requested number of phases strictly after
the reference instant, with each timestamp localized at its own offset.
`moonphases()` retains its legacy nearest-phase selection and naive timestamps.

`crescent_forecast()` begins on the conjunction's local civil date. Each evening
reports available sunset and moonset times, the model's best evaluation time,
the signed lag from sunset to moonset, a status, and a score when evaluated.
Pass the same UTC conjunction to `islamic-times-map generate --conjunction ...`.
The global mapper retains its grid and configured fixed offset, while the local
forecast uses the observer's timezone for each evening. The visibility criteria
predict viewing conditions using astronomical geometry.
