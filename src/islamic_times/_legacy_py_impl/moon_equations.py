"""Deprecated pure-Python moon-equation implementations."""

from __future__ import annotations

import math
from dataclasses import replace
from datetime import datetime, timedelta
from typing import List, Tuple

import numpy as np

from islamic_times import calculation_equations as ce
from islamic_times import moon_equations as active
from islamic_times import sun_equations as se
from islamic_times import time_equations as te
from islamic_times.it_dataclasses import Angle, DateTimeInfo, Moon, ObserverInfo, RightAscension


def moon_nutation(jde: float) -> Tuple[np.ndarray, float, float, float]:
    t = (jde - te.J2000) / te.JULIAN_CENTURY
    t2 = t ** 2
    t3 = t ** 3
    t4 = t ** 4

    fundamental_arguments = np.array(
        [
            297.8501921 + 445267.1114034 * t - 0.0018819 * t2 + t3 / 545868 - t4 / 113065000,
            357.5291092 + 35999.0502909 * t - 0.0001536 * t2 + t3 / 24490000,
            134.9633964 + 477198.8675055 * t + 0.0087414 * t2 + t3 / 69699 - t4 / 14712000,
            93.2720950 + 483202.0175233 * t - 0.0036539 * t2 - t3 / 3526000 + t4 / 863310000,
            218.3164477 + 481267.88123421 * t - 0.0015786 * t2 + t3 / 538841 - t4 / 65194000,
        ]
    ) % 360

    a = np.array([119.75 + 131.849 * t, 53.09 + 479264.290 * t, 313.45 + 481266.484 * t]) % 360
    eccentricity = 1 - 0.002516 * t - 0.0000074 * t2

    args_lr = np.array(active.__MOON_NUTATION_ARGUMENTS_LR).reshape(-1, 4)
    coeff_lr = np.array(active.__MOON_NUTATION_COEFF_LR).reshape(-1, 2)
    args_b = np.array(active.__MOON_NUTATION_ARGUMENTS_B).reshape(-1, 4)
    coeff_b = np.array(active.__MOON_NUTATION_COEFF_B)

    temp_lr = np.dot(args_lr, fundamental_arguments[:4])
    temp_b = np.dot(args_b, fundamental_arguments[:4])

    eccentricity_comp_lr = np.where(np.abs(args_lr[:, 1]) == 0, 1, eccentricity ** np.abs(args_lr[:, 1]))
    eccentricity_comp_b = np.where(np.abs(args_b[:, 1]) == 0, 1, eccentricity ** np.abs(args_b[:, 1]))

    sum_l = np.sum(eccentricity_comp_lr * coeff_lr[:, 0] * ce.sin(temp_lr))
    sum_r = np.sum(eccentricity_comp_lr * coeff_lr[:, 1] * np.cos(np.radians(temp_lr)))
    sum_b = np.sum(eccentricity_comp_b * coeff_b * ce.sin(temp_b))

    sum_l += 3958 * ce.sin(a[0]) + 1962 * ce.sin(fundamental_arguments[4] - fundamental_arguments[3]) + 318 * ce.sin(a[1])
    sum_b += (
        -2235 * ce.sin(fundamental_arguments[4])
        + 382 * ce.sin(a[2])
        + 175 * ce.sin(a[0] - fundamental_arguments[3])
        + 175 * ce.sin(a[0] + fundamental_arguments[3])
        + 127 * ce.sin(fundamental_arguments[4] - fundamental_arguments[2])
        - 115 * ce.sin(fundamental_arguments[4] + fundamental_arguments[2])
    )

    return (fundamental_arguments, sum_l, sum_b, sum_r)


def moonrise_or_moonset(observer_date: DateTimeInfo, observer: ObserverInfo, rise_or_set: str = 'set') -> datetime | float:
    if rise_or_set not in ['rise', 'set', 'moonrise', 'moonset']:
        raise ValueError("Invalid value for rise_or_set. Please use 'rise' or 'set'.")

    ymd = datetime(observer_date.date.year, observer_date.date.month, observer_date.date.day)
    new_jd = te.gregorian_to_jd(observer_date.date) - te.fraction_of_day(observer_date.date)
    new_deltaT = te.delta_t_approx(ymd.year, ymd.month)

    moon_params: List[Moon] = []
    for i in range(3):
        ymd_temp = te.jd_to_gregorian(new_jd + i - 1, observer_date.utc_offset)
        delT_temp = te.delta_t_approx(ymd_temp.year, ymd_temp.month)
        sun_params = se.sunpos(replace(observer_date, date=ymd_temp, jd=(new_jd + i - 1), deltaT=delT_temp), observer)
        delPsi = sun_params.nutation[0]
        moon_params.append(active.moonpos(replace(observer_date, date=ymd_temp, jd=(new_jd + i - 1), deltaT=delT_temp), observer, delPsi, sun_params.true_obliquity))

    h_zero: Angle = Angle(0.7275 * moon_params[1].eh_parallax.decimal - 0.566667)
    cosH_zero: float = (math.sin(h_zero.radians) - math.sin(observer.latitude.radians) * math.sin(moon_params[1].declination.radians)) / (
        math.cos(observer.latitude.radians) * math.cos(moon_params[1].declination.radians)
    )
    if abs(cosH_zero) < 1:
        H_zero = Angle(math.degrees(math.acos(cosH_zero)))
    else:
        return math.inf

    if math.isnan(H_zero.decimal):
        return math.inf

    sidereal_time: Angle = te.greenwich_mean_sidereal_time(new_jd)
    m0: float = (moon_params[1].right_ascension.decimal_degrees.decimal - observer.longitude.decimal - sidereal_time.decimal) / 360

    if rise_or_set in ['rise', 'moonrise']:
        m_event = m0 - H_zero.decimal / 360
    else:
        m_event = m0 + H_zero.decimal / 360

    for _ in range(3):
        theta_event: Angle = Angle((sidereal_time.decimal + 360.985647 * m_event) % 360)
        n_event: float = m_event + new_deltaT / 86400
        interp_dec: Angle = Angle(ce.interpolation(n_event, moon_params[0].declination.decimal, moon_params[1].declination.decimal, moon_params[2].declination.decimal))
        interp_ra = RightAscension(
            ce.interpolation(
                n_event,
                moon_params[0].right_ascension.decimal_degrees.decimal,
                moon_params[1].right_ascension.decimal_degrees.decimal,
                moon_params[2].right_ascension.decimal_degrees.decimal,
            )
            / 15
        )

        local_hour_angle = Angle((theta_event.decimal - (-observer.longitude.decimal) - interp_ra.decimal_degrees.decimal) % 360)
        moon_alt = Angle(
            math.degrees(
                math.asin(
                    math.sin(observer.latitude.radians) * math.sin(interp_dec.radians)
                    + math.cos(observer.latitude.radians) * math.cos(interp_dec.radians) * math.cos(local_hour_angle.radians)
                )
            )
        )

        deltaM = (moon_alt.decimal - h_zero.decimal) / (360 * math.cos(interp_dec.radians) * math.cos(observer.latitude.radians) * math.sin(local_hour_angle.radians))
        m_event += deltaM

    return datetime(ymd.year, ymd.month, ymd.day) + timedelta(days=m_event) - timedelta(hours=observer_date.utc_offset)


def calculate_visibility(sun_az: Angle, sun_alt: Angle, moon_az: Angle, moon_alt: Angle, moon_pi: Angle, criterion: int = 0) -> float:
    daz = Angle(abs(sun_az.decimal - moon_az.decimal))
    arcv = Angle(abs(sun_alt.decimal - moon_alt.decimal))
    arcl = Angle(math.degrees(math.acos(math.cos(daz.radians) * math.cos(arcv.radians))))

    semi_diameter = 0.27245 * (moon_pi.decimal * 60)
    semi_diameter_prime = semi_diameter * (1 + math.sin(moon_alt.radians) * math.sin(moon_pi.radians))
    w_prime = semi_diameter_prime * (1 - math.cos(arcl.radians))

    if criterion == 0:
        return arcv.decimal - (-0.1018 * w_prime ** 3 + 0.7319 * w_prime ** 2 - 6.3226 * w_prime + 7.1651)
    return (arcv.decimal - (11.8371 - 6.3226 * w_prime + 0.7319 * w_prime ** 2 - 0.1018 * w_prime ** 3)) / 10


def calculate_visibility_shaukat(sun_az: Angle, sun_long: Angle, moon_long: Angle, moon_lat: Angle, moon_az: Angle, moon_pi: Angle, moon_illumin: float) -> float:
    semi_diameter: float = 0.27245 * (moon_pi.decimal * 60)
    width: float = 2 * semi_diameter * moon_illumin

    arcv: Angle = Angle(
        math.degrees(
            math.acos((math.cos(moon_long.radians - sun_long.radians) * math.cos(moon_lat.radians)) / math.cos(sun_az.radians - moon_az.radians))
        )
    )
    return (arcv.decimal - (11.8371 - 6.3226 * width + 0.7319 * width ** 2 - 0.1018 * width ** 3)) / 10


def classify_visibility(q: float, criterion: int = 1) -> str:
    if q == -999:
        return "Moonset before the new moon."
    if q == -998:
        return "Moonset before sunset."
    if q == -997:
        return "Moonset & Sunset don't exist."
    if q == -996:
        return "Sunset doesn't exist."
    if q == -995:
        return "Moonset doesn't exist."

    if criterion == 0:
        if q >= 5.65:
            return "A: Crescent is visible by naked eyes."
        if 5.65 > q >= 2:
            return "B: Crescent is visible by optical aid, and it could be seen by naked eyes."
        if 2 > q >= -0.96:
            return "C: Crescent is visible by optical aid only."
        if -0.96 > q:
            return "D: Crescent is not visible even by optical aid."
        raise ValueError("Invalid q value. Must be a float.")

    if q > 0.216:
        return "A: Easily visible."
    if 0.216 >= q > -0.014:
        return "B: Visible under perfect conditions."
    if -0.014 >= q > -0.160:
        return "C: May need optical aid to find crescent."
    if -0.160 >= q > -0.232:
        return "D: Will need optical aid to find crescent."
    if -0.232 >= q > -0.293:
        return "E: Not visible with a [conventional] telescope."
    if -0.293 >= q:
        return "F: Not visible; below the Danjon limit."
    raise ValueError("Invalid q value. Must be a float.")

