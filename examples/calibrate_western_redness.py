"""Evaluate western-redness estimates by holding out entire observing campaigns.

Run with islamic_times 3.2.0 installed: python examples/calibrate_western_redness.py
All fit and error calculations are in minutes after sunset. A prediction inside
a photographic bracket has zero error to the bracket. This is a lower bound
on the error to the unknown instant of disappearance.
"""

from collections import Counter
from dataclasses import dataclass, replace
from datetime import date, datetime, time
from statistics import mean
from zoneinfo import ZoneInfo

from islamic_times.islamic_times import ITLocation
from islamic_times.it_dataclasses import Angle
from islamic_times.prayer_fadilah import _TwilightGeometry
from islamic_times.sun_equations import find_proper_suntime

# A clock is local civil time. Numeric values are delays from the authors' sunset.
# Each campaign gets one vote, regardless of the number of observation nights.
SITES = {
    "Jepara": (-6.612128, 110.650919, "Asia/Jakarta"),
    "Masalembu": (-5.0775, 114.601389, "Asia/Jakarta"),
    "Tuban": (-7.059722, 111.688056, "Asia/Jakarta"),
    "Situbondo": (-7.6053, 114.0401, "Asia/Jakarta"),
    "Quebec": (45.5, -74.0, "America/Toronto"),  # Approximate Montreal-area geometry
    "St Joseph": (42.10, -86.48, "America/Detroit"),  # City centre
    "Blackburn": (53.75, -2.48, "Europe/London"),
    "Barombong": (-5.219, 119.384, "Asia/Makassar"),  # Approximate beach coordinates
    "East Jakarta": (-6.303291, 106.872089, "Asia/Jakarta"),
    "Maron": (-6.954314, 110.35945, "Asia/Jakarta"),
    "Benteng Portugis": (-6.406472, 110.917917, "Asia/Jakarta"),
}

FIELD = {
    # Last red / first dark, chapter III. Omit the fifth night: chapters
    # III and IV disagree on whether the photographs were taken Oct 10 or 11.
    # https://eprints.walisongo.ac.id/id/eprint/15773/
    "Jepara": [
        ("2020-10-05", "18:17", "18:34"), ("2020-10-06", "18:13", "18:33"),
        ("2020-10-07", "18:14", "18:37"), ("2020-10-08", "18:24", "18:36"),
    ],
    # Table 2, broadband Sky Quality Meter (SQM) red-disappearance estimates.
    # https://eprints.walisongo.ac.id/id/eprint/15779/
    "Masalembu": [
        ("2020-09-20", "18:17"), ("2020-09-21", "18:18"), ("2020-10-03", "18:17"),
        ("2020-10-11", "18:17"), ("2020-10-13", "18:17"), ("2020-10-14", "18:18"),
        ("2020-10-16", "18:18"),
    ],
    # Table 4.2, p. 89. Onset of fading inferred from the SQM curve.
    # https://eprints.walisongo.ac.id/id/eprint/29206/
    "Tuban": [
        ("2024-06-19", "18:20:24"), ("2024-06-20", "18:20:12"), ("2024-09-20", "18:20:47"),
        ("2024-09-21", "18:20:02"), ("2024-09-22", "18:19:34"), ("2024-12-11", "18:20:27"),
        ("2024-12-19", "18:29:19"), ("2024-12-20", "18:29:29"), ("2024-12-21", "18:29:48"),
    ],
    # Authors report the disappearance of red at approximately 18:06, white at 18:23.
    # https://www.aimmah.ac.id/blog/2026/09/19/observasi-syafaq-di-pantai-tangsi-situbondo-ungkap-perubahan-cahaya-senja-dan-kaitannya-dengan-waktu-shalat/
    "Situbondo": [("2026-08-30", "18:06")],
    # Field sheet pp. 20-21. June 16 arrival was late. June 20 has no usable result.
    # Keep nearby Quebec locations and the common observing team in the same fold.
    # https://birka.nur.nu/prayertimes/prayertimes-references/quoted/sidq_GhauriU.pdf#page=21
    "Quebec": [("2004-06-05", 96), ("2004-06-10", 95), ("2004-06-19", 100), ("2004-06-23", 100)],
}

# Reserved in the original offset/blend comparison. Blackburn
# was among the original committee's sites. Its June range is a consistency
# check, not two independently timed nightly observations.
EXTERNAL = {
    # Tables 4.1-4.2, July 18. At 18:58 both faint red and white persisted.
    # At 19:01 the observer recorded darkness.
    # https://journal.uin-alauddin.ac.id/index.php/hisabuna/article/download/23014/14600
    "Barombong": [("2021-07-18", "18:58", "19:01")],
    # https://moonsighting.com/how-we.html (Tayebi's May 19 field account)
    "St Joseph": [("2009-05-19", "22:11", "22:15")],
    # https://astronomycenter.net/ref/yaqub.pdf (Table 4 and June notes)
    "Blackburn": [("1988-06-12", 81, 84), ("1988-06-13", 81, 84)],
}

# Septianto et al., Table 3. Drone-image brightness cutoffs, not visual sightings.
# Use March 16 from the dated table (the prose conclusion says March 18).
# The paper transposes the latitude/longitude labels, but identifies the campus.
# https://ejournal.unib.ac.id/kumparan_fisika/article/download/18302/9407
JAKARTA = {"East Jakarta": [
    ("2021-03-11", "18:59:19"), ("2021-03-12", "18:57:34"), ("2021-03-13", "18:51:30"),
    ("2021-03-16", "18:38:18"), ("2021-03-25", "18:45:17"), ("2021-04-01", "18:51:38"),
    ("2021-04-02", "18:46:15"),
]}

# A separate observing campaign. Maron gives a direct disappearance time.
# At Benteng the sunset was clouded; red was clear at 18:12 and the sky dark by 18:19.
# Marina and Mrican report darkness or arrival after the red was gone.
# https://eprints.walisongo.ac.id/id/eprint/1863/4/092111122_Bab3.pdf (pp. 77, 82-83)
COASTAL_VISUAL = {"Maron": [("2013-05-30", "18:31")],
                  "Benteng Portugis": [("2013-05-01", "18:12", "18:19")]}


@dataclass(frozen=True)
class Record:
    site: str
    day: date
    low: float
    high: float
    seasonal: float
    fixed: float | None
    geometry: _TwilightGeometry

    def error(self, prediction):
        """Signed distance to the reported time or interval (late is positive)."""
        return prediction - min(max(prediction, self.low), self.high)


def records(campaigns, *, sites=SITES, elevations=None):
    rows = []
    for site, evenings in campaigns.items():
        lat, lon, zone_name = sites[site]
        zone = ZoneInfo(zone_name)
        for stamp, *observed in evenings:
            day = date.fromisoformat(stamp)
            noon = datetime.combine(day, time(12), zone)
            elevation = elevations.get(site, 0) if elevations else 0
            loc = ITLocation(latitude=lat, longitude=lon, elevation=elevation, date=noon, method="JAFARI")
            sunset = loc.sun().sunset
            assert sunset.date() == day
            boundary = loc.sistani_evening_fadilah()["maghrib"].end
            seasonal = (boundary - sunset).total_seconds() / 60
            try:
                deep = find_proper_suntime(loc.observer_dateinfo, loc.observer_info, "set", Angle(14))
                fixed = (deep - sunset).total_seconds() / 60
                assert fixed > 0  # All available 14-degree events in these campaigns precede midnight.
            except ArithmeticError:
                fixed = None
            bounds = []
            for label in observed:
                if isinstance(label, str):
                    clock = datetime.combine(day, time.fromisoformat(label), zone)
                    bounds.append((clock - sunset).total_seconds() / 60)
                else:
                    bounds.append(float(label))
            geometry = _TwilightGeometry.at_sunset(sunset, loc.observer_info)
            rows.append(Record(site, day, bounds[0], bounds[-1], seasonal, fixed, geometry))
    return rows


def feature(row, model):
    if model == "offset":
        return 1.0
    # The blend falls back to the seasonal estimate on nights with no 14-degree crossing.
    return 0.0 if row.fixed is None else row.fixed - row.seasonal


def fit(training, model):
    """Minimize mean campaign absolute interval error with one free parameter."""
    counts = Counter(row.site for row in training)
    # Absolute interval loss is piecewise linear. Its minimum lies at a bound or
    # a breakpoint. In a tie use the smallest departure from the original model.
    candidates = {0.0}
    if model == "blend":
        candidates.add(1.0)
    for row in training:
        x = feature(row, model)
        if x:
            for bound in (row.low, row.high):
                value = (bound - row.seasonal) / x
                if model == "offset" or 0 <= value <= 1:
                    candidates.add(value)

    def score(value):
        total = sum(abs(row.error(row.seasonal + value * feature(row, model))) / counts[row.site]
                    for row in training)
        return round(total, 10), abs(value), value

    return min(candidates, key=score)


def errors(training, evaluation):
    offset = fit(training, "offset")
    blend = fit(training, "blend")
    results = []
    for row in evaluation:
        predictions = [row.seasonal, row.seasonal + offset,
                       row.seasonal + blend * feature(row, "blend"), row.fixed]
        results.append([None if p is None else row.error(p) for p in predictions])
    return results, offset, blend


def summarize(values):
    # Refuse to improve a metric by quietly dropping a missing event.
    return [mean(abs(e) for e in column) if None not in column else None for column in zip(*values)]


def display(label, values):
    print(f"{label:23}", "  ".join(f"{v:7.2f}" if v is not None else "    n/a" for v in values))


def holdout(rows):
    site_errors = []
    for site in dict.fromkeys(r.site for r in rows):
        test = [r for r in rows if r.site == site]
        training = [r for r in rows if r.site != site]
        values, offset, blend = errors(training, test)
        scores = summarize(values)
        display(f"{site} (n={len(test)})", scores)
        print(f"  training offset {offset:+.2f} min, blend {blend:.2f}")
        site_errors.append(scores)
    display("Site-balanced mean", summarize(site_errors))


def fit_solar_descent(rows):
    """Fit T / rate**p by equal-site absolute distance to the reported intervals."""
    site_counts = Counter(row.site for row in rows)
    best_fit = None
    best_key = None
    for step in range(101):
        exponent = step / 100
        factors = [row.geometry.descent_ratio ** -exponent for row in rows]
        # For fixed p the loss is piecewise linear in T. All slope changes
        # occur at a reported boundary divided by its night's multiplier.
        scales = sorted({bound / factor for row, factor in zip(rows, factors) for bound in (row.low, row.high)})
        losses = [
            round(sum(abs(row.error(scale * factor)) / site_counts[row.site]
                      for row, factor in zip(rows, factors)), 10)
            for scale in scales
        ]
        minimum = min(losses)
        optimum = [scale for scale, loss in zip(scales, losses) if loss == minimum]
        # Resolve flat optima at the centre of the minimizing scale interval.
        # Prefer the middle of the exponent range in case of an exact tie.
        key = (minimum, abs(exponent - 0.5), exponent)
        if best_key is None or key < best_key:
            best_key = key
            best_fit = ((optimum[0] + optimum[-1]) / 2, exponent)
    scale, exponent = best_fit
    return round(scale), exponent  # Clock observations do not support a sub-minute reference scale.


def compare_descent(training, rows):
    scale, exponent = fit_solar_descent(training)
    blend = fit(training, "blend")
    reference_interval = fit_constant_delay(training)
    differences = []
    for row in rows:
        estimates = [row.seasonal, row.fixed, row.seasonal + blend * feature(row, "blend"),
                     reference_interval, scale / row.geometry.descent_ratio**exponent,
                     row.geometry.delay_minutes(row.seasonal, scale, exponent)]
        differences.append([None if value is None else row.error(value) for value in estimates])
    return summarize(differences), scale, exponent


def fit_constant_delay(training):
    """Same training target and site weights, with no solar geometry."""
    # The existing offset solver fits a constant when its baseline is zero.
    return fit([replace(row, seasonal=0.0) for row in training], "offset")


def solar_descent_experiment():
    # Refits use only campaigns that actually tracked visible red. The Barombong
    # and Michigan accounts were examined in the earlier comparison. They can
    # participate in cross-validation but are not fresh holdouts for model design.
    visual = records({**{k: FIELD[k] for k in ("Jepara", "Situbondo", "Quebec")},
                      **{k: EXTERNAL[k] for k in ("Barombong", "St Joseph")}})
    print("\nColumns: seasonal, fixed 14 degrees, refitted blend, constant, raw descent, weighted descent.")
    held_out = []
    for site in dict.fromkeys(r.site for r in visual):
        scores, scale, exponent = compare_descent([r for r in visual if r.site != site],
                                                  [r for r in visual if r.site == site])
        display(site, scores)
        print(f"  fit without {site}: T={scale:.3f} min, p={exponent:.2f}")
        held_out.append(scores)
    display("Site-balanced mean", summarize(held_out))

    indonesia = {"Jepara", "Situbondo", "Barombong"}
    for label, train_sites in (("Train in Indonesia", indonesia), ("Train in North America", {"Quebec", "St Joseph"})):
        training = [r for r in visual if r.site in train_sites]
        print(label)
        for site in sorted({r.site for r in visual} - train_sites):
            display(site, compare_descent(training, [r for r in visual if r.site == site])[0])

    scale, exponent = fit_solar_descent(visual)
    print(f"\nFit to all five visual campaigns: T={scale} min, p={exponent:.2f}")
    print("Transfer checks: meter/camera brightness, monthly range (Blackburn), and 2013 coastal observations.")
    for dataset in ({"Masalembu": FIELD["Masalembu"]}, {"Tuban": FIELD["Tuban"]}, JAKARTA,
                    {"Blackburn": EXTERNAL["Blackburn"]}, *({s: v} for s, v in COASTAL_VISUAL.items())):
        rows = records(dataset)
        display(rows[0].site, compare_descent(visual, rows)[0])
        for row in rows:
            delay = row.geometry.delay_minutes(row.seasonal, scale, exponent)
            print(f"  {row.day}: weighted estimate {delay:.2f}, "
                  f"observation [{row.low:.2f},{row.high:.2f}] min after sunset")


def main():
    rows = records(FIELD)
    print("Minutes of absolute error. Columns: current, fitted offset, fitted blend, fixed 14 degrees.")
    print("Hold out each observing campaign")
    holdout(rows)
    print("\nSensitivity: visual and photographic evidence only (remove both SQM campaigns, then refit)")
    holdout([r for r in rows if r.site not in {"Tuban", "Masalembu"}])
    print("\nRegional transfer: train in Quebec, hold out Indonesia together")
    for site in sorted(FIELD.keys() - {"Quebec"}):
        values, _, _ = errors([r for r in rows if r.site == "Quebec"], [r for r in rows if r.site == site])
        display(site, summarize(values))
    print("\nExternal observations: fit once to all five campaigns")
    offset, blend = fit(rows, "offset"), fit(rows, "blend")
    print(f"Full-fit offset {offset:+.2f} min, blend {blend:.2f}")
    for row in records(EXTERNAL):
        values, _, _ = errors(rows, [row])
        display(f"{row.site} {row.day}", values[0])  # Signed to show early versus late.
        print(f"  observed [{row.low:.2f}, {row.high:.2f}] min after sunset; current {row.seasonal:.2f}")


if __name__ == "__main__":
    main()
    solar_descent_experiment()
