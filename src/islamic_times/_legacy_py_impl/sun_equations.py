"""Deprecated pure-Python sun-equation implementations."""

from __future__ import annotations

import math
from dataclasses import replace
from datetime import datetime, timedelta
from typing import List, Tuple

import numpy as np

from islamic_times import calculation_equations as ce
from islamic_times import time_equations as te
from islamic_times.it_dataclasses import Angle, DateTimeInfo, ObserverInfo, RightAscension
from islamic_times import sun_equations as active


def oblique_eq(jde: float) -> Angle:
    u = ((jde - te.J2000) / te.JULIAN_CENTURY) / 100
    eps = 23 + 26 / 60 + (21.448 / 3600)

    powers = np.power(u, np.arange(1, len(active.__OBLIQUITY_TERMS) + 1))
    eps += np.dot(active.__OBLIQUITY_TERMS, powers) / 3600
    return Angle(eps)


def sun_nutation(jde: float) -> Tuple[Angle, Angle]:
    t = (jde - te.J2000) / te.JULIAN_CENTURY
    t2 = t * t
    t3 = t * t2

    ta = np.array([
        297.850363 + 445267.11148 * t - 0.0019142 * t2 + t3 / 189474.0,
        357.52772 + 35999.05034 * t - 0.0001603 * t2 - t3 / 300000.0,
        134.96298 + 477198.867398 * t + 0.0086972 * t2 + t3 / 56250.0,
        93.27191 + 483202.017538 * t - 0.0036825 * t2 + t3 / 327270,
        125.04452 - 1934.136261 * t + 0.0020708 * t2 + t3 / 450000.0,
    ]) % 360
    ta = np.radians(ta)

    sun_args = np.array(active.__SUN_NUTATION_ARGUMENTS).reshape(-1, 5)
    sun_coeff = np.array(active.__SUN_NUTATION_COEFFICIENTS).reshape(-1, 4)

    ang = np.dot(sun_args, ta)
    dp = np.sum((sun_coeff[:, 0] + sun_coeff[:, 1] * t) * np.sin(ang))
    de = np.sum((sun_coeff[:, 2] + sun_coeff[:, 3] * t) * np.cos(ang))

    deltaPsi = dp / (3600.0 * 10000.0)
    deltaEpsilon = de / (3600.0 * 10000.0)
    return (Angle(deltaPsi), Angle(deltaEpsilon))


def sunrise_or_sunset(observer_date: DateTimeInfo, observer: ObserverInfo, rise_or_set: str, angle: Angle = Angle(5 / 6)) -> datetime | float:
    if rise_or_set not in ['rise', 'set', 'sunrise', 'sunset']:
        raise ValueError("Invalid value for rise_or_set. Please use 'rise' or 'set'.")

    ymd = datetime(observer_date.date.year, observer_date.date.month, observer_date.date.day)
    new_jd = te.gregorian_to_jd(observer_date.date) - te.fraction_of_day(observer_date.date)
    new_deltaT = te.delta_t_approx(ymd.year, ymd.month)

    sun_params: List[active.Sun] = []
    for i in range(3):
        ymd_temp = te.jd_to_gregorian(new_jd + i - 1, observer_date.utc_offset)
        delT_temp = te.delta_t_approx(ymd_temp.year, ymd_temp.month)
        sun_params.append(active.sunpos(replace(observer_date, date=ymd_temp, jd=(new_jd + i - 1), deltaT=delT_temp), observer))

    h_zero: Angle = Angle(-angle.decimal)
    cosH_zero: float = (math.sin(h_zero.radians) - math.sin(observer.latitude.radians) * math.sin(sun_params[1].apparent_declination.radians)) / (
        math.cos(observer.latitude.radians) * math.cos(sun_params[1].apparent_declination.radians)
    )

    if abs(cosH_zero) < 1:
        H_zero = Angle(math.degrees(math.acos(cosH_zero)))
    else:
        return math.inf

    sidereal_time: Angle = te.greenwich_mean_sidereal_time(new_jd)
    m0: float = (sun_params[1].apparent_right_ascension.decimal_degrees.decimal - observer.longitude.decimal - sidereal_time.decimal) / 360

    event = rise_or_set.lower()
    if event in ['rise', 'sunrise']:
        m_event: float = m0 - H_zero.decimal / 360
    else:
        m_event = m0 + H_zero.decimal / 360

    for _ in range(3):
        theta_event: Angle = Angle((sidereal_time.decimal + 360.985647 * m_event) % 360)
        n_event: float = m_event + new_deltaT / 86400
        interp_dec_event = Angle(
            ce.interpolation(
                n_event,
                sun_params[0].apparent_declination.decimal,
                sun_params[1].apparent_declination.decimal,
                sun_params[2].apparent_declination.decimal,
            )
        )

        interp_ra_event = RightAscension(
            ce.interpolation(
                n_event,
                sun_params[0].apparent_right_ascension.decimal_degrees.decimal,
                sun_params[1].apparent_right_ascension.decimal_degrees.decimal,
                sun_params[2].apparent_right_ascension.decimal_degrees.decimal,
            )
            / 15
        )

        local_hour_angle_event = Angle((theta_event.decimal - (-observer.longitude.decimal) - interp_ra_event.decimal_degrees.decimal) % 360)
        sun_alt = Angle(
            math.degrees(
                math.asin(
                    math.sin(observer.latitude.radians) * math.sin(interp_dec_event.radians)
                    + math.cos(observer.latitude.radians) * math.cos(interp_dec_event.radians) * math.cos(local_hour_angle_event.radians)
                )
            )
        )

        deltaM = (sun_alt.decimal - h_zero.decimal) / (
            360 * math.cos(interp_dec_event.radians) * math.cos(observer.latitude.radians) * math.sin(local_hour_angle_event.radians)
        )
        m_event += deltaM

    return datetime(ymd.year, ymd.month, ymd.day) + timedelta(days=m_event) - timedelta(hours=observer_date.utc_offset)

