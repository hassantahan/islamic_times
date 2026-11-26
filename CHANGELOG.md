# Changelog

## 3.2.0

- Fixed sunrise, sunset and twilight crossings sometimes returning the following
  civil day's event, which could also suppress the evening prayer window.
- Corrected the Ḥanafī ʿAṣr shadow factor and calculated Islamic midnight from
  elapsed time across timezone offset changes.
- Fixed the Makkah preset at polar latitudes and exposed which prayer times were
  estimated by a high-latitude rule.
- Validate custom prayer angles before applying any of the requested changes.
- Added date-range prayer schedules, a prayer method catalog, two-way tabular
  Hijri conversion, chronological lunar phases, and explicit conjunction forecasts.
- Unified Gregorian-to-Hijri conversion with the tabular civil calendar. This
  fixes day-zero dates and can change Hijri labels and the dates of the Makkah
  Ramaḍān interval.
- Added reusable mosque profiles with JSON round-tripping and per-prayer minute
  adjustments, plus CSV and iCalendar prayer schedule exports.
- Added estimated Ẓuhr and ʿAṣr faḍīlah windows based on Sayyid Sistani’s
  ordinary daylight shadow rules.
- Added Maghrib and ʿIshāʾ faḍīlah estimates for non-travellers using a
  seasonal model for western twilight, with experimental solar-descent and
  seasonal/14° blend options.
- Added the Fajr preferred window for Jaʿfarī and Tehran, with civil dawn as
  the default closing estimate and an optional local solar angle.
- Jaʿfarī and Tehran use the opening of the estimated preferred ʿIshāʾ window
  as the default timetable entry and attach the window to the result.
- Jaʿfarī and Tehran use the opening of the preferred ʿAṣr window as the
  default timetable entry and expose the full window on the prayer result.
- Added the FCNA Canada preset and the published Kuwait and Qatar Adhan
  variants, with Maghrib at calculated sunset.
- Reused adaptive map samples, reduced working memory for native visibility
  batches, and assembled multiprocessing results as chunks complete.
- Added the `islamic-times-map` command, an explicit conjunction option, and a
  choice of dense or adaptive category computation.
- Updated the wheel builder and added checks for Python 3.10 through 3.13 wheel coverage.
- Linux wheels target x86-64. The source distribution remains available for
  32-bit Linux.
