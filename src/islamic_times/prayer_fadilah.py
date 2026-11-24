"""Preferred daytime prayer windows from Sayyid Sistani's published rulings."""

from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from math import inf, isfinite, tan

from islamic_times import sun_equations as se
from islamic_times import time_equations as te
from islamic_times.it_dataclasses import ObserverInfo, SunInfo


@dataclass(frozen=True, slots=True)
class FadilahWindow:
    """A preferred period and its earlier recommended limit, in local time."""

    start: datetime
    preferred_by: datetime
    end: datetime


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
