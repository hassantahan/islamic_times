"""Predict new observing sites using only previously studied sites for fitting.

Run with the library installed: python examples/validate_western_redness.py
Times are local. Errors are in minutes, with equal weight per location. The
primary endpoint is the last visible red. Enhanced images and published
monthly ranges are separate checks with different measurement precision.
"""

from calendar import monthrange
from collections import Counter, defaultdict
from datetime import date, datetime, time, timedelta
from functools import lru_cache
from statistics import mean
from zoneinfo import ZoneInfo

from calibrate_western_redness import COASTAL_VISUAL, EXTERNAL, FIELD, SITES, fit_solar_descent, records
from islamic_times.islamic_times import ITLocation
from islamic_times.it_dataclasses import Angle
from islamic_times.prayer_fadilah import _western_redness_delay
from islamic_times.sun_equations import find_proper_suntime


LOCATIONS = {
    **SITES,
    "Tucson": (32.4204, -110.9645, "America/Phoenix"),  # Observer's coordinates.
    "Ithaca": (42.44, -76.50, "America/New_York"),  # Hilltop, approximate.
    "Chicago": (41.88, -87.98, "America/Chicago"),  # Islamic Foundation, Villa Park area.
    "Yağlıpınar": (39.64, 32.808, "Europe/Istanbul"),  # Türksat hill, near Ankara.
    "Trisik": (-7.9744, 110.1933, "Asia/Jakarta"),
    "Parangtritis": (-8.0244, 110.33, "Asia/Jakarta"),  # Beach; thesis prints 91 arcminutes.
}
ELEVATIONS = {"Tucson": 842}

# Jim Stamm's log, Omar Afzal's Ithaca log, and the Islamic Foundation report
# forwarded by Syed Quadri. Sunset is retained for a sensitivity check. Use
# the final disappearance after the red returns on Sep 12 and 17 in Tucson.
# https://birka.nur.nu/prayertimes/prayertimes-references/quoted/omar-afzal_Fajr%20Isha%20Observations.doc
VISUAL = {
    "Tucson": [
        ("2002-09-12", "18:29", "19:25"), ("2002-09-13", "18:27", "19:32"),
        ("2002-09-14", "18:26", "19:34"), ("2002-09-15", "18:25", "19:34"),
        ("2002-09-16", "18:25", "19:32"), ("2002-09-17", "18:23", "19:32"),
    ],
    "Ithaca": [
        ("2002-09-05", "19:35", "20:21"), ("2002-09-12", "19:23", "20:14"),
        ("2002-09-14", "19:20", "20:01"), ("2002-09-15", "19:18", "19:53"),
        ("2002-09-17", "19:13", "19:59"), ("2002-09-18", "19:10", "19:58"),
    ],
    "Chicago": [
        ("2002-07-29", "20:13", "21:06"), ("2002-08-07", "20:03", "21:06"),
        ("2002-09-03", "19:22", "20:12"), ("2002-09-05", "19:18", "20:20"),
    ],
}

# Özdemir et al., tables 6.3 and 6.4a (printed pp. 49–51). Mean of the
# observers' clocks for the end of redness. In 6.4a the red column is on the
# right; horizontal light and camera measurements use different endpoints.
# The journal summary also lists Gerede, but supplies no separate nightly clocks.
# https://isamveri.org/pdfdrg/D184688/2013/2012_OZDEMIRS.pdf#page=19
ANKARA = {"Yağlıpınar": [
    ("2011-07-31", "21:14:56"), ("2011-08-10", "21:13:17"), ("2011-09-26", "19:54:00"),
    ("2011-09-27", "19:46:30"), ("2011-09-29", "19:40:45"), ("2011-10-03", "19:39:40"),
    ("2011-10-04", "19:36:20"), ("2011-10-05", "19:29:30"), ("2011-10-06", "19:33:00"),
    ("2011-10-19", "19:19:00"), ("2011-10-20", "19:17:24"), ("2011-10-27", "19:03:26"),
    ("2011-11-01", "17:54:24"), ("2011-11-22", "17:39:45"), ("2011-11-23", "17:33:12"),
    ("2011-11-25", "17:42:00"), ("2011-11-28", "17:41:30"), ("2011-11-29", "17:40:00"),
    ("2011-11-30", "17:41:36"), ("2011-12-01", "17:42:48"), ("2011-12-05", "17:42:24"),
    ("2011-12-15", "17:35:00"), ("2012-01-03", "17:52:00"), ("2012-06-11", "21:40"),
    ("2012-06-12", "21:37"), ("2012-06-13", "21:38"), ("2012-06-18", "21:48"),
    ("2012-06-19", "21:52"), ("2012-06-20", "21:46"), ("2012-06-21", "21:47"),
    ("2012-07-10", "21:43"),
]}

# Aswindana, ch. III tables 3.1–3.6 and ch. IV pp. 88–91, tables 4.4–4.6.
# Enhanced long exposures. Use conservative brackets from frames at the
# thresholds explicitly described as still red / already dark. Round outward
# across discrepancies in the summaries. Clouds obscured three transitions.
# https://eprints.walisongo.ac.id/id/eprint/22652/1/Skripsi_1902046035_Ibanez_Sofadella_Agil_Aswindana_Lengkap.pdf
PHOTOS = {
    "Trisik": [("2023-01-19", "19:14", "19:22"), ("2023-01-23", "19:14", "19:22"),
               ("2023-03-24", "18:31", "18:57")],
    "Parangtritis": [("2023-01-21", "19:14", "19:22"), ("2023-01-24", "18:48", "19:16"),
                     ("2023-03-28", "18:27", "18:57")],
}

# Rifa'i, table 4.1 p. 88. Same beach as the earlier Jepara development site.
# On May 26 the text says 18:39 and the table says 18:36, so retain both bounds.
# https://eprints.walisongo.ac.id/id/eprint/1360/5/082111063_Bab4.pdf
REPEAT_SITE = {"Jepara": [
    ("2012-04-16", "18:39"), ("2012-04-18", "18:38"), ("2012-04-21", "18:38"),
    ("2012-04-22", "18:38"), ("2012-04-23", "18:39"), ("2012-04-24", "18:37"),
    ("2012-04-25", "18:38"), ("2012-04-26", "18:38"), ("2012-05-25", "18:36"),
    ("2012-05-26", "18:36", "18:39"),
]}

# Blackburn table 4, Shafaqe Ahmar only. One row spans March/April. Keep
# reported periods intact. No times were observed in the other rows.
# https://astronomycenter.net/ref/yaqub.pdf#page=123
BLACKBURN = [
    (1987, 9, 9, 55, 60), (1987, 10, 10, 58, 68), (1987, 11, 11, 68, 68),
    (1987, 12, 12, 53, 80), (1988, 1, 1, 50, 66), (1988, 3, 4, 56, 56),
    (1988, 5, 5, 57, 75), (1988, 6, 6, 81, 84), (1988, 8, 8, 67, 68),
]


@lru_cache(maxsize=None)
def observer(site, day):
    lat, lon, tz = LOCATIONS[site]
    noon = datetime.combine(day, time(12), ZoneInfo(tz))
    return ITLocation(latitude=lat, longitude=lon, elevation=ELEVATIONS.get(site, 0),
                      date=noon, method="JAFARI")


@lru_cache(maxsize=None)
def crossing(site, day, degrees):
    """Find the descending crossing belonging to this sunset, including after midnight."""
    loc = observer(site, day)
    sunset = loc.sun().sunset
    for delta in (0, 1):
        context = observer(site, day + timedelta(days=delta)).observer_dateinfo
        try:
            event = find_proper_suntime(context, loc.observer_info, "set", Angle(degrees))
        except ArithmeticError:
            continue
        delay = (event.timestamp() - sunset.timestamp()) / 60
        if 0 < delay < 12 * 60:
            return delay
    return None


def adjusted_crossing(row, base, slope):
    """Let effective depression change with the measured post-sunset descent rate."""
    degrees = base + slope * (1 - row.geometry.descent_ratio)
    degrees = max(8, min(18, round(degrees * 4) / 4))
    delay = crossing(row.site, row.day, degrees)
    return row.seasonal if delay is None else delay


def fit_angle(rows, *, varying):
    """Fit to development locations only. Prefer a constant, then 14°, on ties."""
    counts = Counter(row.site for row in rows)
    slopes = range(-6, 7) if varying else (0,)
    candidates = [(step / 4, slope) for step in range(44, 69) for slope in slopes]

    def loss(params):
        total = sum(abs(row.error(adjusted_crossing(row, *params))) / counts[row.site] for row in rows)
        base, slope = params
        return round(total, 9), abs(slope), abs(base - 14)

    return min(candidates, key=loss)


def blended_delay(row, weight):
    fixed = crossing(row.site, row.day, 14)
    if fixed is None:
        return row.seasonal
    return row.seasonal + weight * (fixed - row.seasonal)


def fit_compromise(rows):
    """Minimize the largest site's excess MAE over its better baseline, using training sites only."""
    by_site = defaultdict(list)
    for row in rows:
        by_site[row.site].append(row)

    def site_error(site_rows, weight):
        return mean(abs(row.error(blended_delay(row, weight))) for row in site_rows)

    baselines = {site: min(site_error(group, 0), site_error(group, 1)) for site, group in by_site.items()}

    def worst_excess(weight):
        return max(site_error(group, weight) - baselines[site] for site, group in by_site.items())

    return min((i / 100 for i in range(101)), key=worst_excess)


def report(rows, constant, varying, weight, descent=(50, 1), *, night_average=False):
    scores = []
    all_nights = []
    for site in dict.fromkeys(row.site for row in rows):
        signed = []
        for row in (r for r in rows if r.site == site):
            baseline = crossing(row.site, row.day, 14)
            assert baseline is not None  # Same nights and events for every model scored here.
            signed.append([row.error(value) for value in (
                baseline, row.seasonal, row.geometry.delay_minutes(row.seasonal, *descent),
                adjusted_crossing(row, *constant), adjusted_crossing(row, *varying), blended_delay(row, weight),
            )])
        errors = [mean(map(abs, column)) for column in zip(*signed)]
        scores.append(errors)
        all_nights.extend(signed)
        print(f"{site:14} {len(signed):2}  " + "  ".join(f"{v:6.2f}" for v in errors))
    print("Equal-site mean    " + "  ".join(f"{mean(col):6.2f}" for col in zip(*scores)))
    if night_average:
        print("Equal-night mean   " + "  ".join(f"{mean(map(abs, col)):6.2f}" for col in zip(*all_nights)))


def region_holdout(rows):
    """Refit on other regions. Keep repeat nights, sites and campaigns in the held-out region together."""
    regions = {
        "North America": {"Quebec", "St Joseph", *VISUAL},
        "Indonesia": {"Jepara", "Situbondo", "Barombong", *COASTAL_VISUAL, *PHOTOS},
        "Türkiye": set(ANKARA),
    }
    assert set.union(*regions.values()) == {row.site for row in rows}
    for region, sites in regions.items():
        training = [row for row in rows if row.site not in sites]
        held_out = [row for row in rows if row.site in sites]
        constant = fit_angle(training, varying=False)
        varying = fit_angle(training, varying=True)
        descent = fit_solar_descent(training)
        weight = fit_compromise(training)
        print(f"\nHold out {region}; angle fits {constant}, {varying}; descent {descent}; blend weight {weight}")
        report(held_out, constant, varying, weight, descent)


def blackburn_ranges():
    print("\nBlackburn: observed, seasonal and 14° monthly ranges (minutes), and 14° coverage.")
    for year, start, end, low, high in BLACKBURN:
        days = [date(year, month, day) for month in range(start, end + 1)
                for day in range(1, monthrange(year, month)[1] + 1)]
        seasonal = [_western_redness_delay(day, LOCATIONS["Blackburn"][0]).total_seconds() / 60 for day in days]
        fixed = [v for day in days if (v := crossing("Blackburn", day, 14)) is not None]
        angle_range = f"{min(fixed):.0f}–{max(fixed):.0f}" if fixed else "unreachable"
        print(f"{year}-{start:02}..{end:02}   {low}–{high}   {min(seasonal):.0f}–{max(seasonal):.0f}"
              f"   {angle_range}   {len(fixed)}/{len(days)}")


def main():
    training = records({**{site: FIELD[site] for site in ("Jepara", "Situbondo", "Quebec")},
                        **{site: EXTERNAL[site] for site in ("Barombong", "St Joseph")}, **COASTAL_VISUAL})
    constant = fit_angle(training, varying=False)
    varying = fit_angle(training, varying=True)
    weight = fit_compromise(training)
    print(f"{len(training)} development evenings; fitted (base degrees, slope): {constant}, {varying}; blend {weight}")
    print("Mean absolute error, minutes: 14°, seasonal, solar descent, fitted angle, varying angle, compromise.")
    clocks = {site: [(day, end) for day, _, end in nights] for site, nights in VISUAL.items()}
    visual = records(clocks, sites=LOCATIONS, elevations=ELEVATIONS)
    turkey = records(ANKARA, sites=LOCATIONS)
    print("\nNew sites: reported clock of final visible red")
    report(visual + turkey, constant, varying, weight, night_average=True)
    print("\nYağlıpınar by meteorological season (same site and observers)")
    for season, months in (("Summer", (6, 7, 8)), ("Autumn", (9, 10, 11)), ("Winter", (12, 1, 2))):
        print(season)
        report([row for row in turkey if row.day.month in months], constant, varying, weight)
    delays = {site: [(day, (datetime.fromisoformat(day + 'T' + end)
                            - datetime.fromisoformat(day + 'T' + start)).total_seconds() / 60)
                      for day, start, end in nights] for site, nights in VISUAL.items()}
    print("\nU.S. observations, using the field observer's sunset as the origin")
    report(records(delays, sites=LOCATIONS, elevations=ELEVATIONS), constant, varying, weight)
    photos = records(PHOTOS, sites=LOCATIONS)
    print("\nNew sites: enhanced photographic brackets (one observing team)")
    report(photos, constant, varying, weight)
    repeated = records(REPEAT_SITE)
    print("\nRepeat site: separate Jepara campaign in April–May")
    report(repeated, constant, varying, weight)
    region_holdout(training + visual + turkey + photos + repeated)
    blackburn_ranges()


if __name__ == "__main__":
    main()
