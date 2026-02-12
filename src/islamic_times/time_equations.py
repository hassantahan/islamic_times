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
import pytz
from typing import Tuple, Dict
from warnings import warn
from datetime import datetime, time
from timezonefinder import TimezoneFinder
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

    t = (julian_day - J2000) / JULIAN_CENTURY
    t2 = t ** 2
    t3 = t ** 3

    theta_zero = 280.46061837 + 360.98564736629 * (julian_day - J2000) + \
                 0.000387933 * t2 - t3 / 38710000

    return Angle(theta_zero % 360)

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

    y = date.year
    m = date.month
    day = date.day + fraction_of_day(date)

    if date.month <= 2:
        y -= 1
        m += 12
    
    a = int(math.floor(y / 100))
    b = 2 - a + int(math.floor(a / 4))

    jd = int(math.floor(365.25 * (y + 4716))) + int(math.floor(30.6001 * (m + 1))) + day + b - 1524.5 - zone / 24
    return jd

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

    if jd == -1:
        return datetime.min

    jd = jd + 0.5 - adjust_for_tz_diff / 24
    z = int(jd)
    f = jd - z

    if z < 2299161:
        a = z
    else:
        alpha = int(math.floor((z - 1867216.25) / 36524.25))
        a = z + 1 + alpha - int(math.floor(alpha / 4))

    b = a + 1524
    c = int(math.floor((b - 122.1) / 365.25))
    d = int(math.floor(365.25 * c))
    e = int(math.floor((b - d) / 30.6001))

    day = b - d - int(30.6001 * e) + f
    if e < 14:
        month = e - 1
    else:
        month = e - 13
    if month > 2:
        year = c - 4716
    else:
        year = c - 4715

    # calculate hours, minutes, and seconds
    f_day = day - int(math.floor(day))
    hour = f_day * 24
    minute = (hour - int(math.floor(hour))) * 60
    second = (minute - int(math.floor(minute))) * 60
    microsec = (second - int(math.floor(second))) * 1000000

    return datetime(int(year), int(month), int(day), int(hour), int(minute), int(second), int(microsec))

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
    
    y = year + (month - 0.5) / 12

    if year < -500:
        u = (year - 1820) / 100
        return -20 + 32 * u**2
    elif year < 500:
        u = y / 100
        return 10583.6 - 1014.41 * u + 33.78311 * u**2 - 5.952053 * u**3 - 0.1798452 * u**4 + 0.022174192 * u**5 + 0.0090316521 * u**6
    elif year < 1600:
        u = (y - 1000) / 100
        return 1574.2 - 556.01 * u + 71.23472 * u**2 + 0.319781 * u**3 - 0.8503463 * u**4 - 0.005050998 * u**5 + 0.0083572073 * u**6
    elif year < 1700:
        t = y - 1600
        return 120 - 0.9808 * t - 0.01532 * t**2 + t**3 / 7129
    elif year < 1800:
        t = y - 1700
        return 8.83 + 0.1603 * t - 0.0059285 * t**2 + 0.00013336 * t**3 - t**4 / 1174000
    elif year < 1860:
        t = y - 1800
        return 13.72 - 0.332447 * t + 0.0068612 * t**2 + 0.0041116 * t**3 - 0.00037436 * t**4 + 0.0000121272 * t**5 - 0.0000001699 * t**6 + 0.000000000875 * t**7
    elif year < 1900:
        t = y - 1860
        return 7.62 + 0.5737 * t - 0.251754 * t**2 + 0.01680668 * t**3 -0.0004473624 * t**4 + t**5 / 233174
    elif year < 1920:
        t = y - 1900
        return -2.79 + 1.494119 * t - 0.0598939 * t**2 + 0.0061966 * t**3 - 0.000197 * t**4
    elif year < 1941:
        t = y - 1920
        return 21.20 + 0.84493 * t - 0.076100 * t**2 + 0.0020936 * t**3
    elif year < 1961:
        t = y - 1950
        return 29.07 + 0.407 * t - t**2 / 233 + t**3 / 2547
    elif year < 1986:
        t = y - 1975
        return 45.45 + 1.067 * t - t**2 / 260 - t**3 / 718
    elif year < 2005:
        t = y - 2000
        return 63.86 + 0.3345 * t - 0.060374 * t**2 + 0.0017275 * t**3 + 0.000651814 * t**4 + 0.00002373599 * t**5
    elif year < 2050:
        t = y - 2000
        return 62.92 + 0.32217 * t + 0.005589 * t**2
    elif year < 2150:
        return -20 + 32 * ((y - 1820) / 100)**2 - 0.5628 * (2150 - y)
    else:
        u = (year - 1820) / 100
        return -20 + 32 * u**2

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
    # Create a TimezoneFinder object
    tf = TimezoneFinder()

    # Get the timezone name for the coordinates
    timezone_str = tf.certain_timezone_at(lat=lat, lng=long)

    # Get the timezone object for the given name
    timezone = pytz.timezone(timezone_str) # type: ignore

    # Localize the given date with the timezone
    localized_datetime = timezone.localize(datetime.combine(day, datetime.min.time()))

    # Get the UTC offset for the given timezone and date
    utc_offset = localized_datetime.utcoffset()

    # Calculate the UTC offset in hours
    return (timezone_str, utc_offset.total_seconds() / 3600) # type: ignore

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
