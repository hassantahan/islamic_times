"""
Astronomical time equations and calendar conversion helpers.

This module provides shared temporal primitives used by the solar/lunar
calculation layers, including:

- Gregorian <-> Julian Day conversion,
- Hijri civil-date conversion helpers,
- Delta T approximations,
- timezone offset lookup from coordinates,
- utility helpers for midpoint and offset formatting.
"""

import math
import atexit
from typing import Tuple, Dict, Optional, Any
from warnings import warn
from datetime import datetime, time
from islamic_times.it_dataclasses import Angle

# Constants used in astronomical calculations.
J2000: float                = 2451545.0
HIJRI_JD: float             = 1948439.389675
HIJRI_STANDARD_LNG: float   = 39.71658
JULIAN_CENTURY: float       = 36525.0
JULIAN_MILLENNIUM: float    = JULIAN_CENTURY * 10
ASTRONOMICAL_UNIT: float    = 149597870.7
TROPICAL_YEAR: float        = 365.24219878
EARTH_RADIUS_KM: float      = 6378.14
MECCA_LAT: float            = 21.420164986
MECCA_LONG: float           = 39.822330044

ISLAMIC_MONTHS: Dict[int, str] = {
        1: 'Muḥarram',
        2: 'Ṣaffar',
        3: 'Rabīʿ al-Awwal',
        4: 'Rabīʿ al-Thānī',
        5: 'Jumādā al-Ūlā',
        6: 'Jumādā al-Thāniyah',
        7: 'Rajab',
        8: 'Shaʿbān',
        9: 'Ramaḍān',
        10: 'Shawwāl',
        11: 'Dhū al-Qaʿdah',
        12: 'Dhū al-Ḥijjah',
    }

ISLAMIC_DAYS = {
        'Sunday': 'al-Aḥad',
        'Monday': 'al-Ithnayn',
        'Tuesday': 'al-Thulāthāʾ',
        'Wednesday': 'al-Arbiʿāʾ',
        'Thursday': 'al-Khamīs',
        'Friday': 'al-Jumuʿah',
        'Saturday': 'al-Sabt',
    }

_TZ_FINDER: Optional[Any] = None
_TZ_NAME_CACHE: Dict[Tuple[float, float], str] = {}
_UTC_OFFSET_CACHE: Dict[Tuple[str, int], float] = {}
_TZ_COORD_ROUND_DIGITS = 4


def _get_timezone_finder() -> Any:
    """Return a lazily-initialized TimezoneFinder singleton."""
    global _TZ_FINDER
    if _TZ_FINDER is None:
        from timezonefinder import TimezoneFinder
        _TZ_FINDER = TimezoneFinder()
    return _TZ_FINDER


def _cleanup_timezone_finder() -> None:
    """Release timezonefinder file resources before interpreter teardown."""
    global _TZ_FINDER
    tf = _TZ_FINDER
    _TZ_FINDER = None
    if tf is None:
        return

    for attr_name in ("boundaries", "holes"):
        polygon_array = getattr(tf, attr_name, None)
        coord_accessor = getattr(polygon_array, "coordinates", None)
        cleanup = getattr(coord_accessor, "cleanup", None)
        if callable(cleanup):
            try:
                cleanup()
            except Exception:
                pass


atexit.register(_cleanup_timezone_finder)

def fraction_of_day(date: datetime) -> float:
    """Return elapsed fraction of the civil day for a datetime.

    Parameters
    ----------
    date : datetime
        Datetime whose elapsed day fraction is required.

    Returns
    -------
    float
        Fraction in the range [0, 1).
    """
    return (date - datetime.combine(date.date(), time(0, tzinfo=date.tzinfo))).total_seconds() / (3600 * 24)

# Taken from pg. 88 of AA 
# (DOES NOT TAKE JDE)
def greenwich_mean_sidereal_time(julian_day: float) -> Angle:
    """Compute Greenwich Mean Sidereal Time (GMST) for a Julian Day.

    Parameters
    ----------
    julian_day : float
        Julian Day.

    Returns
    -------
    Angle
        GMST angle normalized to [0, 360) degrees.
    """
    warn('This particular function will no longer be supported in python. The proper function is in the C extension as "islamic_times.astro_core.greenwich_mean_sidereal_time()".', DeprecationWarning)

    from islamic_times._legacy_py_impl import time_equations as legacy_te
    return legacy_te.greenwich_mean_sidereal_time(julian_day)

# Look to Jean Meeus' "Astronomical Algorithms"
def gregorian_to_jd(date: datetime, zone: float = 0) -> float:
    """Convert a Gregorian datetime to Julian Day.

    Parameters
    ----------
    date : datetime
        Gregorian datetime.
    zone : float, default 0
        Timezone offset in hours applied as ``jd - zone / 24``.

    Returns
    -------
    float
        Julian Day value.
    """
    warn('This particular function will no longer be supported in python. The proper function is in the C extension as "islamic_times.astro_core.gregorian_to_jd()".', DeprecationWarning)

    from islamic_times._legacy_py_impl import time_equations as legacy_te
    return legacy_te.gregorian_to_jd(date, zone)

# Look to Jean Meeus' "Astronomical Algorithms" pg. 
def jd_to_gregorian(jd: float, adjust_for_tz_diff: float = 0) -> datetime:
    """Convert Julian Day to Gregorian datetime.

    Parameters
    ----------
    jd : float
        Julian Day.
    adjust_for_tz_diff : float, default 0
        Hour offset adjustment applied before conversion.

    Returns
    -------
    datetime
        Gregorian datetime. Returns ``datetime.min`` when ``jd == -1``.
    """
    warn('This particular function will no longer be supported in python. The proper function is in the C extension as "islamic_times.astro_core.jd_to_gregorian()".', DeprecationWarning)

    from islamic_times._legacy_py_impl import time_equations as legacy_te
    return legacy_te.jd_to_gregorian(jd, adjust_for_tz_diff)

# Look to Jean Meeus' "Astronomical Algorithms".
# Limitation: this civil-date conversion assumes midnight day boundaries.
# Religious sunset-based day boundaries are handled at higher application layers.
def gregorian_to_hijri(year: int, month: int, day: float, zone: int = 0) -> Tuple[int, int, int]:
    """Convert Gregorian civil date to Hijri civil date tuple.

    Parameters
    ----------
    year : int
        Gregorian year.
    month : int
        Gregorian month.
    day : float
        Gregorian day value.
    zone : int, default 0
        Retained for compatibility. Currently unused.

    Returns
    -------
    tuple[int, int, int]
        Hijri year, month, and day.
    """
    x = year
    m = month

    if m < 3:
        x -= 1
        m += 12
    
    alpha = int(x / 100)
    beta = 2 - alpha + int(alpha / 4)
 
    day += beta
    if day == 0:
        m -= 1
        if m in [1, 3, 5, 7, 8, 10, 12]:
            day = 31
        elif m == 2:
            day = 28
        else:
            day = 30

    b = int(365.25 * x) + int(30.6001 * (m + 1)) + day + 1722519 + beta
    c = int((b - 122.1) / 365.25)
    d = int(365.25 * c)
    e = int((b - d) / 30.6001)

    if x % 4 == 0:
        w = 1
    else:
        w = 2
    
    n = int(275 * m / 9) - w * int((m + 9) / 12) + day - 30
    A = x - 623
    B = int(A / 4)
    C = A % 4
    C_1 = 365.2501 * C
    C_2 = int(C_1)

    if C_1 - C_2 > 0.5:
        C_2 += 1

    D_prime = 1461 * B + 170 + C_2
    q = int(D_prime / 10631)
    r = D_prime % 10631
    j = int(r / 354)
    k = r % 354
    o = int((11 * j + 14) / 30)
    h = 30 * q + j + 1
    jj = k - o + n - 1

    if jj > 354:
        cl = h % 30
        dl = (11 * cl + 3) % 30
        if dl < 19:
            jj -= 354
            h += 1
        elif dl > 18:
            jj -= 355
            h += 1
    elif jj == 0:
        jj = 355
        h -= 1
    
    s = int((jj - 1) / 29.5)
    if jj != 355:
        i_month = 1 + s
        i_day = int(jj - 29.5 * s)
    else:
        i_month = 12
        i_day = 30
    i_year = h

    return (i_year, i_month, i_day)

def get_islamic_month(month: int) -> str:
    """Return Islamic month name for month number.

    Parameters
    ----------
    month : int
        Islamic month number in the range 1-12.

    Returns
    -------
    str
        Month name, or ``"Invalid month number"`` for invalid input.
    """

    return ISLAMIC_MONTHS.get(month, "Invalid month number")

def get_islamic_day(day: str) -> str:
    """Return Arabic Islamic weekday name for an English weekday key.

    Parameters
    ----------
    day : str
        English weekday name (for example, ``"Friday"``).

    Returns
    -------
    str
        Islamic weekday name, or ``"Invalid day"`` for invalid input.
    """

    return ISLAMIC_DAYS.get(day, "Invalid day")

# Based off of https://eclipse.gsfc.nasa.gov/LEcat5/deltatpoly.html
def delta_t_approx(year: int, month: int) -> float:
    """Approximate Delta T (TT - UT) in seconds.

    Parameters
    ----------
    year : int
        Gregorian year.
    month : int
        Gregorian month.

    Returns
    -------
    float
        Approximate Delta T in seconds.

    References
    ----------
    https://eclipse.gsfc.nasa.gov/LEcat5/deltatpoly.html
    """
    warn('This particular function will no longer be supported in python. The proper function is in the C extension as "islamic_times.astro_core.delta_t_approx()".', DeprecationWarning)
    
    from islamic_times._legacy_py_impl import time_equations as legacy_te
    return legacy_te.delta_t_approx(year, month)

# The most computationally expensive function
# Finds the UTC offset given a date and coordinates
def find_utc_offset(lat: float, long: float, day: datetime) -> Tuple[str, float]:
    """Determine timezone name and UTC offset for a coordinate/date.

    Parameters
    ----------
    lat : float
        Latitude in decimal degrees.
    long : float
        Longitude in decimal degrees.
    day : datetime
        Date used for offset resolution (DST-aware).

    Returns
    -------
    tuple[str, float]
        IANA timezone name and UTC offset in decimal hours.

    Notes
    -----
    This function is practically expensive because it performs polygon-based
    timezone lookup plus timezone-localization work. Call it sparingly and cache
    results when many calculations reuse the same location/date.
    """
    from pytz import timezone as pytz_timezone

    if not isinstance(day, datetime):
        raise TypeError(f"'day' must be of type `datetime`, but got `{type(day).__name__}`.")

    day_date = day.date()
    coord_key = (round(lat, _TZ_COORD_ROUND_DIGITS), round(long, _TZ_COORD_ROUND_DIGITS))
    timezone_str = _TZ_NAME_CACHE.get(coord_key)

    if timezone_str is None:
        tf = _get_timezone_finder()
        timezone_str = tf.certain_timezone_at(lat=lat, lng=long)
        if timezone_str is None:
            timezone_str = tf.timezone_at(lat=lat, lng=long)
        if not timezone_str:
            raise ValueError(
                f"Unable to determine timezone for latitude={lat}, longitude={long}, date={day_date.isoformat()}."
            )
        _TZ_NAME_CACHE[coord_key] = timezone_str

    offset_key = (timezone_str, day_date.toordinal())
    cached_offset = _UTC_OFFSET_CACHE.get(offset_key)
    if cached_offset is not None:
        return timezone_str, cached_offset

    timezone = pytz_timezone(timezone_str)
    localized_datetime = timezone.localize(datetime.combine(day_date, datetime.min.time()))
    utc_offset = localized_datetime.utcoffset()
    if utc_offset is None:
        raise ValueError(
            f"Unable to resolve UTC offset for timezone='{timezone_str}' on date={day_date.isoformat()}."
        )

    offset_hours = utc_offset.total_seconds() / 3600.0
    _UTC_OFFSET_CACHE[offset_key] = offset_hours
    return timezone_str, offset_hours

# Finds the middle time between two datetimes. Used to find islamic midnight (usually either between sunset & sunrise, or sunrise & fajr).
def time_midpoint(datetime1: datetime, datetime2: datetime) -> datetime:
    """Compute midpoint datetime between two datetime values.

    Parameters
    ----------
    datetime1 : datetime
        First datetime.
    datetime2 : datetime
        Second datetime.

    Returns
    -------
    datetime
        Midpoint between the two inputs.

    Raises
    ------
    TypeError
        If either input is not a ``datetime`` instance.
    ValueError
        If either value represents a missing event sentinel.
    """
    if type(datetime1) is not datetime or type(datetime2) is not datetime:
        raise TypeError
    
    if datetime1 == math.inf or datetime2 == math.inf:
        raise ValueError

    # Calculate the difference
    difference = datetime2 - datetime1

    # Calculate the midpoint
    midpoint = datetime1 + difference / 2

    return midpoint

# Simple way to make the offset look nice
def format_utc_offset(utc_offset: float) -> str:
    """Format UTC offset hours as ``UTC+HH:MM``.

    Parameters
    ----------
    utc_offset : float
        UTC offset in decimal hours.

    Returns
    -------
    str
        Formatted offset string.
    """
    hours = int(utc_offset)
    minutes = int((utc_offset - hours) * 60)

    # Since minutes are absolute, take the absolute value
    minutes = abs(minutes)

    # Use the built-in `format` function to convert to the required format
    offset_str = "UTC{:+03d}:{:02d}".format(hours, minutes)
    return offset_str
