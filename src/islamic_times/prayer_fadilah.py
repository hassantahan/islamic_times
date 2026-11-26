"""Preferred prayer windows from Sayyid Sistani's published rulings."""

from calendar import isleap
from dataclasses import dataclass
from datetime import date, datetime, timedelta, timezone
from math import cos, inf, isfinite, radians, sin, sqrt, tan

from islamic_times import sun_equations as se
from islamic_times import time_equations as te
from islamic_times.it_dataclasses import Angle, DateTimeInfo, FadilahWindow, ObserverInfo, SunInfo


_CIVIL_DAWN_DEGREES = 6

# Red twilight offsets from Adhan's Moonsighting Committee method.
# Each coefficient scales with absolute latitude at the given seasonal
# anchor, counted in days from the local winter solstice.
# https://github.com/batoulapps/adhan-js/blob/d2993d56fac7d8a06abba4837347bd3fc9cbd188/src/Astronomical.ts#L453-L496
_RED_TWILIGHT_ANCHORS = (
    (0, 17.4), (91, -7.16), (137, 5.12), (183, 19.44),
    (229, 5.12), (275, -7.16), (366, 17.4),
)
_RED_TWILIGHT_MAX_LATITUDE = 55


def _western_redness_delay(local_day: date, latitude: float) -> timedelta:
    """Estimate elapsed time from sunset to the fading of western redness."""
    leap_day = int(isleap(local_day.year))
    winter_offset = 10 if latitude >= 0 else -(172 + leap_day)
    winter_day = (local_day.timetuple().tm_yday + winter_offset) % (365 + leap_day)
    for (left, before), (right, after) in zip(_RED_TWILIGHT_ANCHORS, _RED_TWILIGHT_ANCHORS[1:]):
        if winter_day < right:
            fraction = (winter_day - left) / (right - left)
            coefficient = before + fraction * (after - before)
            return timedelta(seconds=round(60 * (62 + coefficient * abs(latitude) / 55)))
    raise AssertionError("The seasonal anchors must cover the entire year.")


@dataclass(frozen=True)
class _TwilightGeometry:
    """Solar descent at sunset, with latitude and declination in radians."""

    descent_ratio: float
    latitude: float
    declination: float

    @classmethod
    def at_sunset(cls, sunset: datetime, observer: ObserverInfo) -> "_TwilightGeometry":
        """Average vertical speed over 30 elapsed minutes, relative to 15°/hour."""
        start = sunset.astimezone(timezone.utc)
        end = start + timedelta(minutes=30)
        setting = se.sunpos(te.observer_datetime_info(start), observer)
        later = se.sunpos(te.observer_datetime_info(end), observer)
        ratio = (setting.true_altitude.decimal - later.true_altitude.decimal) / 7.5
        return cls(ratio, observer.latitude.radians, setting.apparent_declination.radians)

    def delay_minutes(self, seasonal: float, scale: float = 50.0, exponent: float = 1.0) -> float:
        """Blend the fitted descent delay with the seasonal estimate on shallow paths."""
        # At equatorial reference speed, the training delay corresponds to this altitude.
        altitude = radians(-(5 / 6 + 0.25 * scale))
        # Differentiate sin(h) = sin(phi) sin(delta) + cos(phi) cos(delta) cos(H).
        # The normalized vertical speed tends to zero at grazing incidence.
        radial = cos(self.latitude) * cos(self.declination)
        along = sin(altitude) - sin(self.latitude) * sin(self.declination)
        speed_squared = radial**2 - along**2
        if speed_squared <= 0 or self.descent_ratio <= 0:
            return seasonal
        weight = min(1.0, sqrt(speed_squared) / cos(altitude))
        descent = scale / self.descent_ratio**exponent
        return seasonal + weight * (descent - seasonal)


def _shadow_ratio(instant: datetime, observer: ObserverInfo) -> float:
    """Calculate the geometric shadow of a vertical object per unit height."""
    altitude = se.sunpos(te.observer_datetime_info(instant), observer).true_altitude.radians
    return 1 / tan(altitude) if altitude > 0 else inf


def _shadow_crossing(
    start: datetime, end: datetime, observer: ObserverInfo, target: float,
) -> datetime:
    """Find the afternoon instant the shadow reaches a bracketed ratio."""
    while end - start > timedelta(seconds=1):
        middle = start + (end - start) / 2
        if _shadow_ratio(middle, observer) < target:
            start = middle
        else:
            end = middle
    return end


def calculate_sistani_daylight_fadilah(observer: ObserverInfo, sun: SunInfo) -> dict[str, FadilahWindow]:
    """Estimate Ẓuhr and ʿAṣr faḍīlah windows under ordinary conditions.

    Parameters
    ----------
    observer : ObserverInfo
        Coordinates and elevation used for the solar position.
    sun : SunInfo
        Actual sunrise and sunset for the observer's civil day.

    Returns
    -------
    dict[str, FadilahWindow]
        Available windows keyed by ``zuhr`` and ``asr`` in the observer's
        timezone. A window is omitted if its boundaries cannot be found.

    Notes
    -----
    Ruling 953 defines legal noon as halfway between sunrise and sunset.
    Rulings 966–967 define the ordinary windows in terms of 2/7, 4/7, and
    6/7 of the object's height added to its shadow at that instant. See
    https://www.sistani.org/persian/book/26575/6319/ and
    https://www.sistani.org/persian/book/26575/6324/.

    The estimate uses geometric solar altitude and a level horizon. It
    does not use high-latitude timetable substitutions or infer intense heat.
    """
    sunrise, sunset = sun.sunrise, sun.sunset
    if not isinstance(sunrise, datetime) or not isinstance(sunset, datetime):
        return {}
    if sunrise.utcoffset() is None or sunset.utcoffset() is None:
        return {}
    start, end = sunrise.astimezone(timezone.utc), sunset.astimezone(timezone.utc)
    if not timedelta(0) < end - start < timedelta(days=1):
        return {}

    noon = start + (end - start) / 2
    baseline = _shadow_ratio(noon, observer)
    if not isfinite(baseline):
        return {}
    sunset_shadow = _shadow_ratio(end, observer)
    zone = sunrise.tzinfo
    crossings = {}
    for parts in (2, 4, 6):
        target = baseline + parts / 7
        if target <= sunset_shadow:
            crossings[parts] = _shadow_crossing(noon, end, observer, target).astimezone(zone)

    windows = {}
    if 2 in crossings and 4 in crossings:
        windows["zuhr"] = FadilahWindow(noon.astimezone(zone), crossings[2], crossings[4])
    if 4 in crossings and 6 in crossings:
        windows["asr"] = FadilahWindow(crossings[2], crossings[4], crossings[6])
    return windows


def calculate_sistani_fajr_fadilah(
    context: DateTimeInfo, observer: ObserverInfo, dawn: datetime, sunrise: datetime,
    end_degrees: float | None = None,
) -> FadilahWindow | None:
    """Open at local dawn. Ruling 965 ends the window when dawn lights the entire horizon.

    Civil dawn is a proxy for the end sign. Its 6° solar crossing is independent
    of the prayer method's Fajr angle. A local calibration can replace it.
    An unavailable crossing leaves the observable sign for the end.
    If the crossing precedes the configured opening, there is no window.
    """
    if dawn.utcoffset() is None or sunrise.utcoffset() is None:
        return None
    start, rise = dawn.astimezone(timezone.utc), sunrise.astimezone(timezone.utc)
    if not timedelta(0) < rise - start < timedelta(days=1):
        return None
    degrees = _CIVIL_DAWN_DEGREES if end_degrees is None else end_degrees
    end = None
    try:
        crossing = se.find_proper_suntime(context, observer, "rise", Angle(degrees))
    except ArithmeticError:
        pass
    else:
        crossing_utc = crossing.astimezone(timezone.utc)
        if crossing_utc <= start:
            return None
        if crossing_utc < rise:
            end = crossing
    return FadilahWindow(dawn, None, end, estimated=end is not None, end_sign="whole_horizon_lit")


def calculate_sistani_evening_fadilah(
    local_day: date, latitude: float, sunset: datetime, maghrib: datetime, next_dawn: datetime | None,
    *, geometry: _TwilightGeometry | None = None, blend_crossing: datetime | None = None,
) -> dict[str, FadilahWindow]:
    """Estimate the ordinary non-traveller windows using local solar events.

    Rulings 968–969 use the disappearance of western redness as the shared
    boundary. The seasonal estimate uses sunset, latitude and civil date.
    Supplying sunset geometry selects the experimental solar-descent estimate.
    A supplied 14° crossing selects the experimental 37% blend. A missing
    crossing leaves the seasonal estimate intact. All three estimates
    cover latitudes through 55° north or south.

    The end of the ʿIshāʾ window is one-third of the elapsed night from sunset
    to next dawn. When the dawn crossing is unavailable, only the Maghrib
    window can be estimated.
    """
    if abs(latitude) > _RED_TWILIGHT_MAX_LATITUDE:
        return {}

    sunset_utc = sunset.astimezone(timezone.utc)
    delay = _western_redness_delay(local_day, latitude)
    if geometry is not None:
        delay = timedelta(seconds=round(60 * geometry.delay_minutes(delay.total_seconds() / 60)))
    elif blend_crossing is not None:
        fixed_delay = blend_crossing.astimezone(timezone.utc) - sunset_utc
        delay += (fixed_delay - delay) * 0.37
    redness_end = sunset_utc + delay
    if not sunset_utc <= maghrib.astimezone(timezone.utc) < redness_end:
        return {}

    first_third = None
    if next_dawn is not None:
        night = next_dawn.astimezone(timezone.utc) - sunset_utc
        if timedelta(0) < night < timedelta(days=1):
            if redness_end > sunset_utc + night / 2:
                return {}
            first_third = sunset_utc + night / 3

    local_boundary = redness_end.astimezone(sunset.tzinfo)
    windows = {"maghrib": FadilahWindow(maghrib, None, local_boundary, estimated=True)}
    if first_third is not None and redness_end < first_third:
        windows["isha"] = FadilahWindow(
            local_boundary, None, first_third.astimezone(sunset.tzinfo), estimated=True,
        )
    return windows
