"""Deprecated pure-Python time-equation implementations."""

from __future__ import annotations

import math
from datetime import datetime, time

from islamic_times.it_dataclasses import Angle

J2000 = 2451545.0
JULIAN_CENTURY = 36525.0


def fraction_of_day(date: datetime) -> float:
    return (date - datetime.combine(date.date(), time(0, tzinfo=date.tzinfo))).total_seconds() / (3600 * 24)


def greenwich_mean_sidereal_time(julian_day: float) -> Angle:
    t = (julian_day - J2000) / JULIAN_CENTURY
    t2 = t ** 2
    t3 = t ** 3

    theta_zero = 280.46061837 + 360.98564736629 * (julian_day - J2000) + 0.000387933 * t2 - t3 / 38710000
    return Angle(theta_zero % 360)


def gregorian_to_jd(date: datetime, zone: float = 0) -> float:
    y = date.year
    m = date.month
    day = date.day + fraction_of_day(date)

    if date.month <= 2:
        y -= 1
        m += 12

    a = int(math.floor(y / 100))
    b = 2 - a + int(math.floor(a / 4))
    return int(math.floor(365.25 * (y + 4716))) + int(math.floor(30.6001 * (m + 1))) + day + b - 1524.5 - zone / 24


def jd_to_gregorian(jd: float, adjust_for_tz_diff: float = 0) -> datetime:
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
    month = e - 1 if e < 14 else e - 13
    year = c - 4716 if month > 2 else c - 4715

    f_day = day - int(math.floor(day))
    hour = f_day * 24
    minute = (hour - int(math.floor(hour))) * 60
    second = (minute - int(math.floor(minute))) * 60
    microsec = (second - int(math.floor(second))) * 1000000

    return datetime(int(year), int(month), int(day), int(hour), int(minute), int(second), int(microsec))


def delta_t_approx(year: int, month: int) -> float:
    y = year + (month - 0.5) / 12

    if year < -500:
        u = (year - 1820) / 100
        return -20 + 32 * u**2
    if year < 500:
        u = y / 100
        return 10583.6 - 1014.41 * u + 33.78311 * u**2 - 5.952053 * u**3 - 0.1798452 * u**4 + 0.022174192 * u**5 + 0.0090316521 * u**6
    if year < 1600:
        u = (y - 1000) / 100
        return 1574.2 - 556.01 * u + 71.23472 * u**2 + 0.319781 * u**3 - 0.8503463 * u**4 - 0.005050998 * u**5 + 0.0083572073 * u**6
    if year < 1700:
        t = y - 1600
        return 120 - 0.9808 * t - 0.01532 * t**2 + t**3 / 7129
    if year < 1800:
        t = y - 1700
        return 8.83 + 0.1603 * t - 0.0059285 * t**2 + 0.00013336 * t**3 - t**4 / 1174000
    if year < 1860:
        t = y - 1800
        return 13.72 - 0.332447 * t + 0.0068612 * t**2 + 0.0041116 * t**3 - 0.00037436 * t**4 + 0.0000121272 * t**5 - 0.0000001699 * t**6 + 0.000000000875 * t**7
    if year < 1900:
        t = y - 1860
        return 7.62 + 0.5737 * t - 0.251754 * t**2 + 0.01680668 * t**3 - 0.0004473624 * t**4 + t**5 / 233174
    if year < 1920:
        t = y - 1900
        return -2.79 + 1.494119 * t - 0.0598939 * t**2 + 0.0061966 * t**3 - 0.000197 * t**4
    if year < 1941:
        t = y - 1920
        return 21.20 + 0.84493 * t - 0.076100 * t**2 + 0.0020936 * t**3
    if year < 1961:
        t = y - 1950
        return 29.07 + 0.407 * t - t**2 / 233 + t**3 / 2547
    if year < 1986:
        t = y - 1975
        return 45.45 + 1.067 * t - t**2 / 260 - t**3 / 718
    if year < 2005:
        t = y - 2000
        return 63.86 + 0.3345 * t - 0.060374 * t**2 + 0.0017275 * t**3 + 0.000651814 * t**4 + 0.00002373599 * t**5
    if year < 2050:
        t = y - 2000
        return 62.92 + 0.32217 * t + 0.005589 * t**2
    if year < 2150:
        return -20 + 32 * ((y - 1820) / 100)**2 - 0.5628 * (2150 - y)

    u = (year - 1820) / 100
    return -20 + 32 * u**2

