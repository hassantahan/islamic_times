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

### Additional prayer presets

| Method key | Fajr | ʿIshāʾ | High-latitude default | Published rules |
| --- | --- | --- | --- | --- |
| `FCNA_CANADA` | 13° | 13° | Angle-based | [FCNA twilight recommendation](https://fiqhcouncil.org/the-suggested-calculation-method-for-fajr-and-isha/) and [prayer boundaries](https://fiqhcouncil.org/fifteen-or-eighteen-degrees-calculating-prayer-fasting-times-in-islam/) |
| `KUWAIT` | 18° | 17.5° | Half the night | [Batoul Apps Adhan calculation parameters](https://github.com/batoulapps/adhan-js/blob/d2993d56fac7d8a06abba4837347bd3fc9cbd188/src/CalculationMethod.ts#L59-L67) |
| `QATAR` | 18° | 90 minutes after Maghrib year-round | Half the night | [Batoul Apps Adhan calculation parameters](https://github.com/batoulapps/adhan-js/blob/d2993d56fac7d8a06abba4837347bd3fc9cbd188/src/CalculationMethod.ts#L59-L67) |

These presets use solar transit for Ẓuhr, the one-shadow convention for ʿAṣr,
and calculated sunset for Maghrib. The library reports the midpoint between
sunset and the next sunrise as Islamic midnight. It uses its own solar solver,
timestamp precision and polar fallback. Kuwait and Qatar follow the published
Adhan variants, including its [ordinary high-latitude night fraction](https://github.com/batoulapps/adhan-js/blob/d2993d56fac7d8a06abba4837347bd3fc9cbd188/src/CalculationParameters.ts#L13-L17).
Qatar's fixed ʿIshāʾ interval stays in force when Maghrib is estimated. Local
mosque timetable adjustments can be applied with a `MosqueProfile`.

`PrayerMethod.maghrib_offset_minutes=0` uses calculated sunset when
`maghrib_angle` is zero. Set it to `1` for sunset plus one minute. Omitted or
null offsets in saved profiles retain the legacy one-minute adjustment.
With a fixed ʿIshāʾ interval, that interval starts at the resulting Maghrib.

Najaf’s Astronomical Research Center publishes 18° for Fajr and 4.5° for
Maghrib. Its calculator publishes shared prayer landmarks and omits a
separate ʿIshāʾ timestamp. Sayyid Sistani prescribes the evening prayers in
sequence. A complete five-entry Najaf preset requires a separately sourced
convention for the ʿIshāʾ entry. The two published angles can be applied to a
custom profile. [Najaf calculator](https://english.nojumi.org/prayertimes) ·
[Sayyid Sistani’s evening rulings](https://www.sistani.org/english/book/48/2212/).

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

## Fajr faḍīlah

With Jaʿfarī and Tehran, `observer.prayer_times().fajr.fadilah` opens at
the preset's local dawn. Its `end_sign` is `"whole_horizon_lit"`. Ruling 965
describes dawn light reaching the sky and the entire horizon, apparently before
red appears in the east. Sunrise marks the end of the prayer's valid time.
[Preferred Fajr, ruling 965](https://www.sistani.org/persian/book/26575/6323/).

The default estimates the end at civil dawn, when the Sun's geometric centre
rises to 6° below the horizon. At civil twilight the horizon is generally clear
in good conditions. This is a proxy for the observed sign and has no calibration
against timings of the whole horizon brightening. Fajr observational studies
such as OpenFajr record the first light spreading across the eastern horizon.
Changing the Fajr angle changes the opening only.
[U.S. Naval Observatory twilight definitions](https://aa.usno.navy.mil/faq/RST_defs).
[OpenFajr observational study](https://openfajr.org/docs/research_paper.pdf).

Colour all-sky frames at Kolonica Observatory, Slovakia, show diffuse blue
dawn light around the visible skyline by 02:16 UTC on 12 and 14 August 2026
and 04:51 UTC on 19 November 2026 (about 27, 30 and 20 minutes before 6°).
These selected frames have exposures of 2.1, 2.2 and 2.8 seconds. The camera
adjusts exposure, and trees and structures obscure parts of the lowest sky.
The frames show how early the camera can record the spread of dawn light.
They cannot time the disappearance of darkness to the unaided eye.
[12 August frame](https://astronomy.science.upjs.sk/allsky/archive/20260811-20260812/20260812_0218.jpg).
[14 August frame](https://astronomy.science.upjs.sk/allsky/archive/20260813-20260814/20260814_0219.jpg).
[19 November frame](https://astronomy.science.upjs.sk/allsky/archive/20261118-20261119/20261119_0451.jpg).

Aboushelib and Morcos recorded sky brightness every minute on four August 2022
mornings at Kottamia Observatory, Egypt. Digitizing the published curves gives
a median overnight baseline from 01:00 to 01:30 UTC for each direction.
The west was the last of the three sampled directions to reach both comparison
levels. It reached 10 times its overnight brightness 12–16 minutes before the
calculated 6° crossing, and 100 times its overnight brightness 2–4 minutes before.

| Date | West, 10× | West, 100× | Calculated 6° |
| --- | ---: | ---: | ---: |
| 2022-08-02 | 02:32 | 02:43 | 02:46 |
| 2022-08-03 | 02:31 | 02:43 | 02:47 |
| 2022-08-04 | 02:33 | 02:44 | 02:48 |
| 2022-08-07 | 02:38 | 02:48 | 02:50 |

Times are UTC and rounded to the nearest minute. The two levels compare the
shape of the recorded brightening with civil dawn. Establishing the end sign
would require human observations around the full horizon. The study measured
only east, west and zenith, with artificial skyglow in the western baseline.
[Aboushelib and Morcos, figures 13–16 and table 5](https://doi.org/10.1038/s41598-026-74791-y).

An observational check used the darkest sampled horizon directions in
Koomen et al.'s clear-sky twilight measurements. At a solar depression of 6°,
those directions were about 0.84 cd/m² at Sacramento Peak and 0.97 cd/m² in
Maryland, compared with 0.039 and 0.034 cd/m² at 9° depression. The tables
describe gradual brightening and assume symmetry between sunrise and sunset.
They give brightness by direction and solar altitude. A timing error for the
proxy would require an independent observation of the end sign.
[Koomen et al., tables I and II](https://doi.org/10.1364/JOSA.42.000353).

The Oldenburg Eye2Sky archive also covers the whole sky. On 20 June 2022,
15 December 2022 and 20 March 2023, the first morning frames in the OLUOL
archive came about 42, 35 and 29 minutes after civil dawn, respectively.
Those dates lack frames bracketing the estimate.
[Eye2Sky public dataset](https://zenodo.org/records/12804613).

SkyCam-O at the Liverpool Telescope on La Palma recorded frames either side
of civil dawn (07:01 UTC on 12 March 2026 and 07:00 UTC on 13 March). Its
automatic exposure and monochrome output make the frames useful for sky
condition checks. The archived frames supply no naked-eye endpoint or calibrated
luminance series.
[SkyCam instrument description](https://telescope.livjm.ac.uk/TelInst/Inst/SkyCam/).
[12 March frames](https://telescope.livjm.ac.uk/data/archive/webfiles/Skycam/2026/20260311/index_o.shtml).
[13 March frames](https://telescope.livjm.ac.uk/data/archive/webfiles/Skycam/2026/20260312/index_o.shtml).

A community with a locally measured solar depression for that sign can call
`observer.set_fajr_fadilah_end_angle(degrees)`. The library finds the morning
crossing between local Fajr and sunrise and marks that end as estimated.
Pass `None` to restore the default. This setting travels with a mosque profile.
An unavailable closing crossing leaves `end=None` and the observed sign.
The library omits the preferred window if the computed end is at or before
the configured Fajr opening, if Fajr needs a substitute, or if local sunrise
is unavailable.

## Evening faḍīlah estimates

`observer.sistani_evening_fadilah()` returns available `maghrib` and `isha`
windows for a non-traveller. Each window has timezone-aware bounds and
`estimated=True`.

By default the library uses Adhan’s red-twilight formula for its Moonsighting Committee
method. It estimates elapsed minutes after local sunset from latitude and time of year.
The library applies it through 55° north or south. This estimates the fading
of western redness and does not model local clouds, haze or horizon obstructions.
[Moonsighting Committee observations](https://moonsighting.com/faq_pt.html).
[Adhan’s red-twilight formula](https://github.com/batoulapps/adhan-js/blob/d2993d56fac7d8a06abba4837347bd3fc9cbd188/src/Astronomical.ts#L453-L496).

A comparison covering 25 evenings in Indonesia and Québec held out each observing
campaign in turn. The seasonal estimate missed the reported time or photographic
interval by 11.6 minutes on average, giving each campaign equal weight. A fitted
blend with the 14° crossing reduced this to 8.2 minutes. On separate observations
the blend was closer at Barombong but about 12 minutes late at St Joseph, where
the seasonal estimate fell within the observed interval. These results support
keeping the boundary labelled as an estimate. The field records, their sources
and the fitting procedure are in the [calibration example](../examples/calibrate_western_redness.py).

### Experimental solar-descent model

Call `observer.set_western_redness_model("solar_descent")` to try a geometry-based
estimate for the same evening boundary. Mosque profiles retain this choice.
Use `"seasonal"` to return to the default.

The model measures how far the Sun descends in the first 30 minutes after
local sunset. If it drops by 7.5°, the fitted delay is 50 minutes. At half
that rate the unweighted estimate is 100 minutes. The model blends this
estimate with the seasonal one using the Sun’s vertical speed near 13.3°
below the horizon. A steep path gives the descent estimate more weight.
If the Sun never gets that far below the horizon, the seasonal estimate
gets full weight. The reference angle controls the blend. It does not
define when visible redness disappears.

The fit used 11 visual or photographic evenings across five campaigns, with
equal weight per campaign. Refitting with each campaign held out gave a mean
absolute error of 8.5 minutes, compared with 9.4 for seasonal and 7.6 for a
fixed 14° crossing. Error is the distance to the reported instant or interval.
Two coastal locations held out of fitting gave opposite results. At Maron
the model was 7.6 minutes early and at Benteng Portugis it was 3.3 minutes
past the observed interval. Neither estimator accounts for cloud, haze or
horizon visibility on an individual evening.
[Coastal field observations, pp. 77 and 82–83](https://eprints.walisongo.ac.id/id/eprint/1863/4/092111122_Bab3.pdf).

### Western redness at new observing sites

Four previously unused sites provided 47 dated, naked-eye observations of the
last visible western red. The US logs cover July to September 2002. At
Yağlıpınar near Ankara, the team recorded red and residual horizontal light
separately. Its 31 evening red measurements span summer, autumn and winter.
The comparison uses local clock times with historical daylight-saving offsets.
[US observers’ field logs](https://birka.nur.nu/prayertimes/prayertimes-references/quoted/omar-afzal_Fajr%20Isha%20Observations.doc).
[Ankara field tables 6.3 and 6.4a, pp. 49–51](https://isamveri.org/pdfdrg/D184688/2013/2012_OZDEMIRS.pdf#page=19).

The table gives mean absolute error in minutes. The experimental blend was
fitted on 13 earlier evenings at seven other sites. It places the prediction
37% of the way from the seasonal time to the 14° crossing. The fitting rule
minimizes the largest increase in a training site's error relative to its
better baseline. If the crossing is absent, the blend uses seasonal alone.
Select it with `observer.set_western_redness_model("blend")`. Mosque profiles
retain this setting. The default is `"seasonal"`.

| Held-out site | Evenings | Fixed 14° | Seasonal | Solar descent | Experimental blend |
| --- | ---: | ---: | ---: | ---: | ---: |
| Tucson | 6 | 3.3 | 4.7 | 4.4 | 4.1 |
| Ithaca | 6 | 27.1 | 12.4 | 19.5 | 17.9 |
| Chicago area | 4 | 22.0 | 7.6 | 12.4 | 12.5 |
| Yağlıpınar | 31 | 3.3 | 8.3 | 5.8 | 5.3 |
| Each site weighted equally | 47 | 14.0 | 8.2 | 10.5 | 9.9 |
| Each evening weighted equally | 47 | 8.0 | 8.3 | 8.0 | 7.3 |

The 14° crossing had the smallest error in each of the three observed seasons
at Yağlıpınar. In a second check, we fitted each experimental model on two
regions and evaluated the third. The blend beat 14° in North America
(10.7 vs 16.0 minutes), while 14° was closer in Indonesia (5.85 vs 6.25)
and Türkiye (3.3 vs 5.4). At Ithaca the seasonal estimate was 12 minutes late
on average. At Yağlıpınar it was eight minutes early. A fitted depression
angle that varied with solar descent missed by 19.1, 9.9 and 9.6 minutes
in those respective regions.
Some Indonesian results are distances to photographic intervals and represent
lower bounds on the timing error. The study sites were held out of our fits.
The original training sites for the published seasonal formula are incompletely documented.

Two new Indonesian beaches provided six enhanced-image intervals across January
and March. Seasonal averaged 2.6 minutes to those intervals and 14° averaged
5.7. Blackburn's monthly observations provide a separate high-latitude check.
In June 1988 observers reported 81–84 minutes after sunset, and the Sun reached
14° below the horizon on only one of the 30 evenings. Blackburn contributed
to the committee's observational work. The research examples keep these
monthly ranges separate from nightly scores.
[Indonesian field study, tables 4.4–4.6](https://eprints.walisongo.ac.id/id/eprint/22652/1/Skripsi_1902046035_Ibanez_Sofadella_Agil_Aswindana_Lengkap.pdf#page=97).
[Blackburn observations, table 4](https://astronomycenter.net/ref/yaqub.pdf#page=123).
[Reproduce the comparisons](../examples/validate_western_redness.py).

### Prayer windows and timetable entries

The selected prayer method supplies Maghrib and the next day's Fajr angle.
Maghrib's preferred window ends at the estimated twilight boundary. ʿIshāʾ's
preferred window begins there and ends one-third of the elapsed time from
sunset to the next dawn. Jaʿfarī and Tehran use the opening of this preferred
window as the ʿIshāʾ timetable entry and include the bounds in
`prayer_times().isha.fadilah`. The timetable and calendar exports mark this
entry as estimated.
[Sistani's Maghrib ruling](https://www.sistani.org/persian/book/26575/6325/).
[Sistani's ʿIshāʾ ruling](https://www.sistani.org/persian/book/26575/6326/).

The calculation requires local sunset and Maghrib. It omits the ʿIshāʾ
window if the dawn crossing is unavailable or if the estimated boundary
falls at or after the end of the first third of the night.
The preferred ʿIshāʾ timetable entry is unavailable when there is no window.
To use a fixed angle, call `set_custom_prayer_angles(isha_angle=14)` or pass
the desired angle. New mosque profiles made from either preset use the
preferred window. Profiles serialized before the `use_isha_fadilah` setting
was added retain their stored angle or interval.

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
raise `ValueError`. The tuple converter `gregorian_to_hijri()` uses the same
tabular dates from the civil Hijri epoch onward. The prayer timetable and map
pipeline use this converter for date labels and the Makkah Ramaḍān interval.

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
