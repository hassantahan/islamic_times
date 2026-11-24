# Changelog

## 3.2.0

- Corrected the Ḥanafī ʿAṣr shadow factor and calculated Islamic midnight from
  elapsed time across timezone offset changes.
- Fixed the Makkah preset at polar latitudes and exposed which prayer times were
  estimated by a high-latitude rule.
- Validate custom prayer angles before applying any of the requested changes.
- Added date-range prayer schedules, a prayer method catalog, two-way tabular
  Hijri conversion, chronological lunar phases, and explicit conjunction forecasts.
- Added reusable mosque profiles with JSON round-tripping and per-prayer minute
  adjustments, plus CSV and iCalendar prayer schedule exports.
- Reused adaptive map samples, reduced working memory for native visibility
  batches, and assembled multiprocessing results as chunks complete.
- Added the `islamic-times-map` command, an explicit conjunction option, and a
  choice of dense or adaptive category computation.
- Updated the wheel builder and added checks for Python 3.10 through 3.13 wheel coverage.
