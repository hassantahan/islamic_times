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

import atexit
from math import floor
from typing import Tuple, Dict, Optional, Any
from warnings import warn
from datetime import date, datetime, time, timedelta, timezone
from islamic_times.it_dataclasses import Angle, DateTimeInfo, IslamicDateInfo

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
_TZ_OBJECT_CACHE: Dict[str, Any] = {}
_PYTZ_TIMEZONE_FACTORY: Optional[Any] = None
_TZ_COORD_ROUND_DIGITS = 4


def _get_timezone_finder() -> Any:
    """Return a lazily-initialized TimezoneFinder singleton."""
    global _TZ_FINDER
    if _TZ_FINDER is None:
        from timezonefinder import TimezoneFinder
        _TZ_FINDER = TimezoneFinder()
    return _TZ_FINDER


def _get_pytz_timezone_factory() -> Any:
    """Return cached `pytz.timezone` callable."""
    global _PYTZ_TIMEZONE_FACTORY
    if _PYTZ_TIMEZONE_FACTORY is None:
        from pytz import timezone as pytz_timezone
        _PYTZ_TIMEZONE_FACTORY = pytz_timezone
    return _PYTZ_TIMEZONE_FACTORY


def _get_timezone_object(timezone_str: str) -> Any:
    """Return cached timezone object for an IANA timezone name."""
    timezone = _TZ_OBJECT_CACHE.get(timezone_str)
    if timezone is not None:
        return timezone

    timezone_factory = _get_pytz_timezone_factory()
    timezone = timezone_factory(timezone_str)
    _TZ_OBJECT_CACHE[timezone_str] = timezone
    return timezone


def _cleanup_timezone_finder() -> None:
    """Release timezonefinder file resources before interpreter teardown."""
    global _TZ_FINDER
    tf = _TZ_FINDER
    _TZ_FINDER = None
    _TZ_OBJECT_CACHE.clear()
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


def find_timezone_name(lat: float, long: float, day: datetime) -> str:
    """Determine IANA timezone name for a coordinate/date context."""
    if not isinstance(day, datetime):
        raise TypeError(f"'day' must be of type `datetime`, but got `{type(day).__name__}`.")

    day_date = day.date()
    coord_key = (round(lat, _TZ_COORD_ROUND_DIGITS), round(long, _TZ_COORD_ROUND_DIGITS))
    timezone_str = _TZ_NAME_CACHE.get(coord_key)
    if timezone_str is not None:
        return timezone_str

    tf = _get_timezone_finder()
    timezone_str = tf.certain_timezone_at(lat=lat, lng=long)
    if timezone_str is None:
        timezone_str = tf.timezone_at(lat=lat, lng=long)
    if not timezone_str:
        raise ValueError(
            f"Unable to determine timezone for latitude={lat}, longitude={long}, date={day_date.isoformat()}."
        )

    _TZ_NAME_CACHE[coord_key] = timezone_str
    return timezone_str


def localize_or_convert_datetime(day: datetime, timezone_obj: Any) -> datetime:
    """Return datetime normalized into ``timezone_obj`` with DST-safe semantics.

    For naive datetimes, this performs strict localization (`is_dst=None`) when
    supported by the timezone object so ambiguous and nonexistent wall times are
    surfaced as explicit errors.
    """
    if not isinstance(day, datetime):
        raise TypeError(f"'day' must be of type `datetime`, but got `{type(day).__name__}`.")

    if day.tzinfo is None:
        localize = getattr(timezone_obj, "localize", None)
        if callable(localize):
            from pytz.exceptions import AmbiguousTimeError, NonExistentTimeError

            try:
                return localize(day, is_dst=None)
            except NonExistentTimeError as exc:
                raise ValueError(
                    "The provided local time does not exist because of a DST transition. "
                    f"datetime={day.isoformat()}, timezone={timezone_obj}."
                ) from exc
            except AmbiguousTimeError as exc:
                raise ValueError(
                    "The provided local time is ambiguous because of a DST transition. "
                    f"datetime={day.isoformat()}, timezone={timezone_obj}."
                ) from exc
        return day.replace(tzinfo=timezone_obj)

    converted = day.astimezone(timezone_obj)
    normalize = getattr(timezone_obj, "normalize", None)
    if callable(normalize):
        return normalize(converted)
    return converted


def find_timezone(lat: float, long: float, day: datetime) -> Any:
    """Resolve and return timezone object for a coordinate/date context."""
    timezone_str = find_timezone_name(lat, long, day)
    return _get_timezone_object(timezone_str)


def observer_datetime_info(value: datetime) -> DateTimeInfo:
    """Build the Julian-date, time-scale, and Hijri context for an aware instant."""
    from islamic_times import astro_core

    offset = value.utcoffset()
    if offset is None:
        raise ValueError("Input datetime must include timezone information.")
    jd = astro_core.gregorian_to_jd(value, offset.total_seconds() / 3600)
    tt_utc, ut1_utc, delta_t = resolve_time_scales(jd)
    return DateTimeInfo(
        date=value,
        jd=jd,
        deltaT=delta_t,
        tt_minus_utc=tt_utc,
        ut1_minus_utc=ut1_utc,
        hijri=IslamicDateInfo(*gregorian_to_hijri(value.year, value.month, value.day)),
    )


def _civil_midnight_offset(midnight: datetime, context: datetime) -> float:
    """Effective UTC offset at the first instant of a civil date, including clock changes."""
    offset = context.utcoffset()
    if offset is None:
        return 0.0
    nominal = midnight.replace(tzinfo=timezone.utc)
    guess = nominal - offset
    zone = context.tzinfo
    if (guess.astimezone(zone).replace(tzinfo=None) == midnight
            and (guess - timedelta(seconds=1)).astimezone(zone).date() < midnight.date()):
        return offset.total_seconds() / 3600

    # Searching in UTC also handles zones that skip midnight or repeat it.
    # ZoneInfo and pytz both convert unambiguously in this direction.
    lo, hi = -86400, 86400
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if (nominal + timedelta(seconds=mid)).astimezone(zone).date() < midnight.date():
            lo = mid
        else:
            hi = mid
    return -hi / 3600


def civil_day_offsets(context: datetime) -> tuple[float, float]:
    """Return the UTC offsets delimiting the observer's local calendar day."""
    if context.tzinfo is None or isinstance(context.tzinfo, timezone):
        offset = context.utcoffset() or timedelta(0)
        return (offset.total_seconds() / 3600,) * 2
    midnight = datetime.combine(context.date(), datetime.min.time())
    return _civil_midnight_offset(midnight, context), _civil_midnight_offset(midnight + timedelta(days=1), context)


def local_event_datetime(value: datetime, context: DateTimeInfo) -> datetime:
    """Convert a native wall time from its solver offset to the observer timezone."""
    offset = context.date.utcoffset()
    if offset is None:
        return value
    instant = value.replace(tzinfo=timezone(offset)).astimezone(timezone.utc)
    return instant.astimezone(context.date.tzinfo)


def add_elapsed_time(value: datetime, duration: timedelta) -> datetime:
    """Add elapsed time, including across timezone offset changes."""
    if value.utcoffset() is None:
        return value + duration
    return (value.astimezone(timezone.utc) + duration).astimezone(value.tzinfo)


# Tabular Islamic calendar, Friday epoch and 30-year type II leap cycle.
_CIVIL_HIJRI_EPOCH = date(622, 7, 19)


def _hijri_year_start(year: int) -> int:
    """Return the number of days from the civil Hijri epoch to year start."""
    return 354 * (year - 1) + (3 + 11 * year) // 30


def _hijri_month_start(month: int) -> int:
    """Return the day index of the first day of a month within a Hijri year."""
    return 29 * (month - 1) + month // 2


def gregorian_to_civil_hijri(value: date) -> IslamicDateInfo:
    """Convert a Gregorian civil date to the tabular Islamic calendar.

    Uses a midnight day boundary and the type II leap cycle from 1 Muḥarram
    1 AH (19 July 622 in the proleptic Gregorian calendar).
    """
    if not isinstance(value, date) or isinstance(value, datetime):
        raise TypeError("'value' must be a Gregorian date.")
    days = (value - _CIVIL_HIJRI_EPOCH).days
    if days < 0:
        raise ValueError("'value' must be on or after the civil Hijri epoch.")
    year = (30 * days + 10646) // 10631
    day_of_year = days - _hijri_year_start(year)
    month = 1
    while month < 12 and day_of_year >= _hijri_month_start(month + 1):
        month += 1
    return IslamicDateInfo(year, month, day_of_year - _hijri_month_start(month) + 1)


def civil_hijri_to_gregorian(year: int, month: int, day: int) -> date:
    """Convert a tabular Hijri date to a proleptic Gregorian civil date."""
    if any(type(part) is not int for part in (year, month, day)):
        raise TypeError("Hijri year, month, and day must be integers.")
    if not 1 <= year <= gregorian_to_civil_hijri(date.max).hijri_year or not 1 <= month <= 12:
        raise ValueError("Hijri year or month is out of range.")
    year_start = _hijri_year_start(year)
    month_start = _hijri_month_start(month)
    if month == 12:
        month_length = _hijri_year_start(year + 1) - year_start - month_start
    else:
        month_length = _hijri_month_start(month + 1) - month_start
    if not 1 <= day <= month_length:
        raise ValueError("Hijri day is out of range for the selected month.")
    ordinal = _CIVIL_HIJRI_EPOCH.toordinal() + year_start + month_start + day - 1
    return date.fromordinal(ordinal)

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

def gregorian_to_hijri(year: int, month: int, day: float, zone: int = 0) -> Tuple[int, int, int]:
    """Convert Gregorian civil date to Hijri civil date tuple.

    Parameters
    ----------
    year : int
        Gregorian year.
    month : int
        Gregorian month.
    day : float
        Gregorian day value. The fractional part is ignored for civil dates.
    zone : int, default 0
        Retained for compatibility. Currently unused.

    Returns
    -------
    tuple[int, int, int]
        Hijri year, month, and day.
    """
    try:
        civil_date = date(year, month, 1) + timedelta(days=floor(day) - 1)
    except (TypeError, ValueError, OverflowError):
        civil_date = None
    if civil_date is not None and civil_date >= _CIVIL_HIJRI_EPOCH:
        hijri = gregorian_to_civil_hijri(civil_date)
        return hijri.hijri_year, hijri.hijri_month, hijri.hijri_day

    # Retain the legacy Meeus arithmetic for inputs outside the civil calendar.
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
    
    import islamic_times.astro_core as fast_astro
    return fast_astro.delta_t_approx(year, month)


def resolve_time_scales(julian_day_utc: float) -> tuple[float, float, float]:
    """Resolve TT-UTC, UT1-UTC, and Delta T for a UTC Julian Day.

    Parameters
    ----------
    julian_day_utc : float
        UTC-based Julian Day.

    Returns
    -------
    tuple[float, float, float]
        ``(tt_minus_utc_seconds, ut1_minus_utc_seconds, delta_t_seconds)``.
    """
    import islamic_times.astro_core as fast_astro
    return fast_astro.resolve_time_scales(julian_day_utc)

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
    if not isinstance(day, datetime):
        raise TypeError(f"'day' must be of type `datetime`, but got `{type(day).__name__}`.")

    day_date = day.date()
    timezone_str = find_timezone_name(lat, long, day)
    offset_key = (timezone_str, day_date.toordinal())
    cached_offset = _UTC_OFFSET_CACHE.get(offset_key)
    if cached_offset is not None:
        return timezone_str, cached_offset

    timezone_obj = _get_timezone_object(timezone_str)
    localized_datetime = localize_or_convert_datetime(datetime.combine(day_date, datetime.min.time()), timezone_obj)
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
        If either input is not a datetime or only one input is timezone-aware.
    """
    if type(datetime1) is not datetime or type(datetime2) is not datetime:
        raise TypeError
    
    if datetime1.utcoffset() is None and datetime2.utcoffset() is None:
        return datetime1 + (datetime2 - datetime1) / 2
    if datetime1.utcoffset() is None or datetime2.utcoffset() is None:
        raise TypeError("Both midpoint endpoints must be naive or both timezone-aware.")

    start_utc = datetime1.astimezone(timezone.utc)
    end_utc = datetime2.astimezone(timezone.utc)
    return (start_utc + (end_utc - start_utc) / 2).astimezone(datetime1.tzinfo)

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
