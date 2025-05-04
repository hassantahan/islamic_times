'''
Module for Solar astronomical calculations.

This module provides functions to:
  - Calculate the Suns's position (declination, right ascension, altitude, azimuth, etc.).
  - Calculate the Sun's nutation.
  - Calculate the obliquity of the ecliptic.
  - Calculate the equation of time offset.

This module should not be used directly unless the user is familiar with the underlying calculations.

References:
  - Jean Meeus, *Astronomical Algorithms*, 2nd Edition, Willmann-Bell, Inc., 1998.
'''

import math
import numpy as np
from typing import List, Tuple
from datetime import datetime, timedelta
from dataclasses import dataclass, replace
from islamic_times import time_equations as te
from islamic_times import calculation_equations as ce
from islamic_times.it_dataclasses import *
from warnings import warn

__OBLIQUITY_TERMS = [
    -4680.93,
    -1.55,
    1999.25,
    -51.38,
    -249.67,
    -39.05,
    7.12,
    27.87,
    5.79,
    2.45
]

__SUN_NUTATION_ARGUMENTS = [
     0,  0,  0,  0,  1,
    -2,  0,  0,  2,  2,
     0,  0,  0,  2,  2,
     0,  0,  0,  0,  2,
     0,  1,  0,  0,  0,
     0,  0,  1,  0,  0,
    -2,  1,  0,  2,  2,
     0,  0,  0,  2,  1,
     0,  0,  1,  2,  2,
    -2, -1,  0,  2,  2,
    -2,  0,  1,  0,  0,
    -2,  0,  0,  2,  1,
     0,  0, -1,  2,  2,
     2,  0,  0,  0,  0,
     0,  0,  1,  0,  1,
     2,  0, -1,  2,  2,
     0,  0, -1,  0,  1,
     0,  0,  1,  2,  1,
    -2,  0,  2,  0,  0,
     0,  0, -2,  2,  1,
     2,  0,  0,  2,  2,
     0,  0,  2,  2,  2,
     0,  0,  2,  0,  0,
    -2,  0,  1,  2,  2,
     0,  0,  0,  2,  0,
    -2,  0,  0,  2,  0,
     0,  0, -1,  2,  1,
     0,  2,  0,  0,  0,
     2,  0, -1,  0,  1,
    -2,  2,  0,  2,  2,
     0,  1,  0,  0,  1,
    -2,  0,  1,  0,  1,
     0, -1,  0,  0,  1,
     0,  0,  2, -2,  0,
     2,  0, -1,  2,  1,
     2,  0,  1,  2,  2,
     0,  1,  0,  2,  2,
    -2,  1,  1,  0,  0,
     0, -1,  0,  2,  2,
     2,  0,  0,  2,  1,
     2,  0,  1,  0,  0,
    -2,  0,  2,  2,  2,
    -2,  0,  1,  2,  1,
     2,  0, -2,  0,  1,
     2,  0,  0,  0,  1,
     0, -1,  1,  0,  0,
    -2, -1,  0,  2,  1,
    -2,  0,  0,  0,  1,
     0,  0,  2,  2,  1,
    -2,  0,  2,  0,  1,
    -2,  1,  0,  2,  1,
     0,  0,  1, -2,  0,
    -1,  0,  1,  0,  0,
    -2,  1,  0,  0,  0,
     1,  0,  0,  0,  0,
     0,  0,  1,  2,  0,
     0,  0, -2,  2,  2,
    -1, -1,  1,  0,  0,
     0,  1,  1,  0,  0,
     0, -1,  1,  2,  2,
     2, -1, -1,  2,  2,
     0,  0,  3,  2,  2,
     2, -1,  0,  2,  2
]

__SUN_NUTATION_COEFFICIENTS = [
    -171996,  -174.2,   92025,     8.9,          #  0,  0,  0,  0,  1 
     -13187,    -1.6,    5736,    -3.1,          # -2,  0,  0,  2,  2 
      -2274,     -.2,     977,     -.5,          #  0,  0,  0,  2,  2 
       2062,      .2,    -895,      .5,          #  0,  0,  0,  0,  2 
       1426,    -3.4,      54,     -.1,          #  0,  1,  0,  0,  0 
        712,      .1,      -7,       0,          #  0,  0,  1,  0,  0 
       -517,     1.2,     224,     -.6,          # -2,  1,  0,  2,  2 
       -386,     -.4,     200,       0,          #  0,  0,  0,  2,  1 
       -301,      .0,     129,     -.1,          #  0,  0,  1,  2,  2 
        217,     -.5,     -95,      .3,          # -2, -1,  0,  2,  2 
       -158,       0,       0,       0,          # -2,  0,  1,  0,  0 
        129,      .1,     -70,       0,          # -2,  0,  0,  2,  1 
        123,       0,     -53,       0,          #  0,  0, -1,  2,  2 
         63,       0,       0,       0,          #  2,  0,  0,  0,  0 
         63,      .1,     -33,       0,          #  0,  0,  1,  0,  1 
        -59,       0,      26,       0,          #  2,  0, -1,  2,  2 
        -58,     -.1,      32,       0,          #  0,  0, -1,  0,  1 
        -51,       0,      27,       0,          #  0,  0,  1,  2,  1 
         48,       0,       0,       0,          # -2,  0,  2,  0,  0 
         46,       0,     -24,       0,          #  0,  0, -2,  2,  1 
        -38,       0,      16,       0,          #  2,  0,  0,  2,  2 
        -31,       0,      13,       0,          #  0,  0,  2,  2,  2 
         29,       0,       0,       0,          #  0,  0,  2,  0,  0 
         29,       0,     -12,       0,          # -2,  0,  1,  2,  2 
         26,       0,       0,       0,          #  0,  0,  0,  2,  0 
        -22,       0,       0,       0,          # -2,  0,  0,  2,  0 
         21,       0,     -10,       0,          #  0,  0, -1,  2,  1 
         17,     -.1,       0,       0,          #  0,  2,  0,  0,  0 
         16,       0,      -8,       0,          #  2,  0, -1,  0,  1 
        -16,      .1,       7,       0,          # -2,  2,  0,  2,  2 
        -15,       0,       9,       0,          #  0,  1,  0,  0,  1 
        -13,       0,       7,       0,          # -2,  0,  1,  0,  1 
        -12,       0,       6,       0,          #  0, -1,  0,  0,  1 
         11,       0,       0,       0,          #  0,  0,  2, -2,  0 
        -10,       0,       5,       0,          #  2,  0, -1,  2,  1 
         -8,       0,       3,       0,          #  2,  0,  1,  2,  2 
          7,       0,      -3,       0,          #  0,  1,  0,  2,  2 
         -7,       0,       0,       0,          # -2,  1,  1,  0,  0 
         -7,       0,       3,       0,          #  0, -1,  0,  2,  2 
         -7,       0,       3,       0,          #  2,  0,  0,  2,  1 
          6,       0,       0,       0,          #  2,  0,  1,  0,  0 
          6,       0,      -3,       0,          # -2,  0,  2,  2,  2 
          6,       0,      -3,       0,          # -2,  0,  1,  2,  1 
         -6,       0,       3,       0,          #  2,  0, -2,  0,  1 
         -6,       0,       3,       0,          #  2,  0,  0,  0,  1 
          5,       0,       0,       0,          #  0, -1,  1,  0,  0 
         -5,       0,       3,       0,          # -2, -1,  0,  2,  1 
         -5,       0,       3,       0,          # -2,  0,  0,  0,  1 
         -5,       0,       3,       0,          #  0,  0,  2,  2,  1 
          4,       0,       0,       0,          # -2,  0,  2,  0,  1 
          4,       0,       0,       0,          # -2,  1,  0,  2,  1 
          4,       0,       0,       0,          #  0,  0,  1, -2,  0 
         -4,       0,       0,       0,          # -1,  0,  1,  0,  0 
         -4,       0,       0,       0,          # -2,  1,  0,  0,  0 
         -4,       0,       0,       0,          #  1,  0,  0,  0,  0 
          3,       0,       0,       0,          #  0,  0,  1,  2,  0 
         -3,       0,       0,       0,          # -1, -1,  1,  0,  0 
         -3,       0,       0,       0,          #  0,  1,  1,  0,  0 
         -3,       0,       0,       0,          #  0, -1,  1,  2,  2 
         -3,       0,       0,       0,          #  2, -1, -1,  2,  2 
         -3,       0,       0,       0,          #  0,  0, -2,  2,  2 
         -3,       0,       0,       0,          #  0,  0,  3,  2,  2 
         -3,       0,       0,       0           #  2, -1,  0,  2,  2 
]

@dataclass(frozen=True, slots=True)
class Sun:
    '''
    A class to compute the position of the Sky in the sky based on given astronomical parameters.
    '''

    # Orbital elements
    mean_longitude: Angle
    mean_anomaly: Angle
    earth_orbit_eccentricity: float
    sun_centre: Angle
    true_longitude: Angle
    true_anomaly: Angle
    geocentric_distance: Distance

    # Nutation and obliquity
    omega: Angle
    apparent_longitude: Angle
    nutation: Tuple[Angle, Angle]
    delta_obliquity: Angle
    mean_obliquity: Angle
    true_obliquity: Angle

    # Apparent coordinates
    true_right_ascension: RightAscension
    true_declination: Angle
    apparent_right_ascension: RightAscension
    apparent_declination: Angle

    # Hour angles
    greenwich_hour_angle: Angle
    local_hour_angle: Angle

    # Topocentric quantities
    eh_parallax: Angle
    topocentric_ascension: RightAscension
    topocentric_declination: Angle
    topocentric_local_hour_angle: Angle

    # Horizontal coordinates
    true_altitude: Angle
    true_azimuth: Angle
    apparent_altitude: Angle

# Chapter 22
def oblique_eq(jde: float) -> Angle:
    '''
    Calculate the obliquity of the ecliptic for a given Julian Ephemeris Day. See Chapter 22 of *Astronomical Algorthims* for more information.

    Parameters:
        jde (float): The Julian Ephemeris Day.

    Returns:
        Angle: The obliquity of the ecliptic.
    '''
    warn('This particular function will no longer be supported in python. The proper function is in the C extension but cannot be called from python. Use "islamic_times.astro_core.compute_sun()" or "islamic_times.sun_equations.sunpos()" instead.', DeprecationWarning)

    u = ((jde - te.J2000) / te.JULIAN_CENTURY) / 100

    eps = 23 + 26 / 60 + (21.448 / 3600)

    powers = np.power(u, np.arange(1, len(__OBLIQUITY_TERMS) + 1))
    eps += np.dot(__OBLIQUITY_TERMS, powers) / 3600

    return Angle(eps)

# Chapter 22
def sun_nutation(jde: float) -> Tuple[Angle, Angle]:
    '''
    Calculate the sun's nutation for a given Julian Ephemeris Day. See Chapter 22 of *Astronomical Algorthims* for more information.

    Parameters:
        jde (float): The Julian Ephemeris Day.

    Returns:
        tuple (Tuple[Angle, Angle]): The nutation in longitude and obliquity.
    '''
    warn('This particular function will no longer be supported in python. The proper function is in the C extension but cannot be called from python. Use "islamic_times.astro_core.compute_sun()" or "islamic_times.sun_equations.sunpos()" instead.', DeprecationWarning)

    # Precompute time variables
    t = (jde - te.J2000) / te.JULIAN_CENTURY
    t2 = t * t
    t3 = t * t2

    ta = np.array([
        297.850363 + 445267.11148 * t - 0.0019142 * t2 + t3 / 189474.0,
        357.52772 + 35999.05034 * t - 0.0001603 * t2 - t3 / 300000.0,
        134.96298 + 477198.867398 * t + 0.0086972 * t2 + t3 / 56250.0,
        93.27191 + 483202.017538 * t - 0.0036825 * t2 + t3 / 327270,
        125.04452 - 1934.136261 * t + 0.0020708 * t2 + t3 / 450000.0
    ]) % 360  

    ta = np.radians(ta)

    sun_args = np.array(__SUN_NUTATION_ARGUMENTS).reshape(-1, 5)
    sun_coeff = np.array(__SUN_NUTATION_COEFFICIENTS).reshape(-1, 4)

    ang = np.dot(sun_args, ta)

    dp = np.sum((sun_coeff[:, 0] + sun_coeff[:, 1] * t) * np.sin(ang))
    de = np.sum((sun_coeff[:, 2] + sun_coeff[:, 3] * t) * np.cos(ang))

    deltaPsi = dp / (3600.0 * 10000.0)
    deltaEpsilon = de / (3600.0 * 10000.0)

    return (Angle(deltaPsi), Angle(deltaEpsilon))

# Chapter 25
def sunpos(observer_date: DateTimeInfo, observer: ObserverInfo) -> Sun:
    '''
    Calculate the various solar positional parameters for a given Julian Ephemeris Day, ΔT, and observer coordinates. See Chapter 25 of the *Astronomical Algorthims* for more information.

    Parameters:
        observer_date (DateTimeInfo): The date and time of the observer.
        observer (ObserverInfo): The observer's coordinates and elevation.

    Returns:
        Sun (obj): A `Sun` object that contains various attributes that describe its position and parameters. 

    Notes: 
    - The temperature and pressure are used for atmospheric refraction calculations. Currently, this feature is disabled.
    '''
    import islamic_times.astro_core as fast_astro
    the_sun: Sun = fast_astro.compute_sun(observer_date.jde, observer_date.deltaT, observer.latitude.decimal, observer.longitude.decimal, observer.elevation.value, observer.temperature, observer.pressure)

    return the_sun

def equation_of_time(deltaPsi: float, L0: float, epsilon: float, alpha: float) -> float:
    '''
    Calculate the equation of time offset for a given certain solar parameters. See Chapter 28 of *Astronomical Algorthims* for more information.

    Parameters:
        deltaPsi (float): The sun's nutation in longitude.
        L0 (float): The sun's mean longitude.
        epsilon (float): The obliquity of the ecliptic.
        alpha (float): The sun's right ascension.

    Returns:
        float: The equation of time offset in minutes.
    '''

    # Only adjust alpha when it appears to have wrapped
    if L0 > 300 and alpha < 50:
        alpha += 360
    elif L0 < 50 and alpha > 300:
        alpha -= 360

    E = (L0 - 0.0057183 - alpha + deltaPsi * ce.cos(epsilon)) * 4  # convert degrees to minutes

    return E

def find_sun_transit(observer_date: DateTimeInfo, observer: ObserverInfo) -> datetime:
    """
    Calculate the transit of the sun (specifically culmination) for a given date and observer coordinates. See Chapter 15 of *Astronomical Algorithms* for more information.

    Parameters:
        observer_date (DateTimeInfo): The date and time of the observer.
        observer (ObserverInfo): The observer's coordinates and elevation.

    Returns:
        datetime: The time of moonset.
    """
    import islamic_times.astro_core as fast_astro
    sun_transit_dt: datetime = fast_astro.find_sun_transit(observer_date.jd, observer_date.deltaT, 
                                       observer.latitude.decimal, observer.longitude.decimal, 
                                       observer.elevation.in_unit(DistanceUnits.METRE), observer.temperature, observer.pressure, 
                                       observer_date.utc_offset)
                                                                                      
    return sun_transit_dt.replace(tzinfo=observer_date.date.tzinfo)

def sunrise_or_sunset(observer_date: DateTimeInfo, observer: ObserverInfo, rise_or_set: str, angle: Angle = Angle(5 / 6)) -> datetime:
    """
    Calculate either the sunrise or the sunset time for a given date and observer coordinates.
    This function computes only the requested event (sunrise or sunset) without calculating both.
    
    Parameters:
        observer_date (DateTimeInfo): The date and time of the observer.
        observer (ObserverInfo): The observer's coordinates and elevation.
        rise_or_set (str): 'rise' or 'sunrise' to calculate sunrise; 'set' or 'sunset' to calculate sunset.
        angle (Angle): The standard altitude of the sun (default is 5/6°).
        
    Returns:
        datetime: The computed time of the requested event (sunrise or sunset).
        
    Raises:
        ValueError: If the value of rise_or_set is not recognized.
    """
    warn('This particular function will no longer be supported in python. The proper function is in the C extension but cannot be called from python. Use "islamic_times.astro_core.find_proper_suntime()" or "islamic_times.sun_equations.find_proper_suntime()" instead.', DeprecationWarning)

    if rise_or_set not in ['rise', 'set', 'sunrise', 'sunset']:
        raise ValueError("Invalid value for rise_or_set. Please use 'rise' or 'set'.")

    # First find the Year Month Day at UT 0h from JDE
    ymd = datetime(observer_date.date.year, observer_date.date.month, observer_date.date.day)
    new_jd = te.gregorian_to_jd(observer_date.date) - te.fraction_of_day(observer_date.date)
    new_deltaT = te.delta_t_approx(ymd.year, ymd.month)

    # Calculate new sun params with the new_jd
    sun_params: List[Sun] = []
    for i in range(3):
        ymd_temp = te.jd_to_gregorian(new_jd + i - 1, observer_date.utc_offset)
        delT_temp = te.delta_t_approx(ymd_temp.year, ymd_temp.month)
        sun_params.append(
                        sunpos(
                            replace(observer_date, date=ymd_temp, jd=(new_jd + i - 1), deltaT=delT_temp), 
                            observer
                        )
                    )

    # Compute H0: the hour angle corresponding to the desired altitude.
    h_zero: Angle = Angle(-angle)
    cosH_zero: float = (math.sin(h_zero.radians) - math.sin(observer.latitude.radians) * math.sin(sun_params[1].apparent_declination.radians)) / \
                (math.cos(observer.latitude.radians) * math.cos(sun_params[1].apparent_declination.radians))
    
    if abs(cosH_zero) < 1:
        H_zero = Angle(math.degrees(math.acos(cosH_zero)))
    else:
        return math.inf

    # Compute Greenwich Mean Sidereal Time (GMST) at new_jd.
    sidereal_time: Angle = te.greenwich_mean_sidereal_time(new_jd)

    # Compute the transit estimate m0.
    m0: float = (sun_params[1].apparent_right_ascension.decimal_degrees.decimal - observer.longitude.decimal - sidereal_time.decimal) / 360

    # Choose which event to compute.
    event = rise_or_set.lower()
    if event in ['rise', 'sunrise']:
        # Initial estimate for sunrise.
        m_event: float = m0 - H_zero.decimal / 360
    else:
        # Initial estimate for sunset.
        m_event: float = m0 + H_zero.decimal / 360

    # Iteratively refine m0 (transit) and m_event (rise or set).
    # We use three iterations which are typically sufficient.
    for _ in range(3):
        # --- Update the event (sunrise or sunset) estimate ---
        theta_event: Angle = Angle((sidereal_time.decimal + 360.985647 * m_event) % 360)
        n_event: float = m_event + new_deltaT / 86400
        interp_dec_event = Angle(ce.interpolation(n_event,
                                            sun_params[0].apparent_declination.decimal,
                                            sun_params[1].apparent_declination.decimal,
                                            sun_params[2].apparent_declination.decimal))
        
        interp_ra_event = RightAscension(ce.interpolation(n_event,
                                            sun_params[0].apparent_right_ascension.decimal_degrees.decimal,
                                            sun_params[1].apparent_right_ascension.decimal_degrees.decimal,
                                            sun_params[2].apparent_right_ascension.decimal_degrees.decimal) / 15
                                        )
        
        local_hour_angle_event = Angle((theta_event.decimal - (-observer.longitude.decimal) - interp_ra_event.decimal_degrees.decimal) % 360)
        sun_alt = Angle(math.degrees(math.asin(math.sin(observer.latitude.radians) * math.sin(interp_dec_event.radians) +
                                        math.cos(observer.latitude.radians) * math.cos(interp_dec_event.radians) *
                                       math.cos(local_hour_angle_event.radians))))
        
        # Correct m_event using the difference between computed altitude and desired h_zero.
        deltaM = (sun_alt.decimal - h_zero.decimal) / (360 * math.cos(interp_dec_event.radians) * math.cos(observer.latitude.radians) * math.sin(local_hour_angle_event.radians))
        m_event += deltaM

    # Convert the fractional day to a datetime object, adjusting for the UTC offset.
    event_dt = datetime(ymd.year, ymd.month, ymd.day) + timedelta(days=m_event) - timedelta(hours=observer_date.utc_offset)
    return event_dt

def find_proper_suntime(observer_date: DateTimeInfo, observer: ObserverInfo, rise_or_set: str, angle: Angle = Angle(5 / 6)) -> datetime:
    """
    Determines the proper local time for a setting or rising moon. It finds the time that corresponds to the reference date given.

    Parameters:
        observer_date (DateTimeInfo): The date and time of the observer.
        observer (ObserverInfo): The observer's coordinates and elevation.
        rise_or_set (str): Find either the setting or rising option.
        angle (Angle): Find the time when the sun reaches the given angle above or below the observer's horizon (controlled by `rise_or_set`). Default is 50 arcminutes (0.8333°) for visible sunset/sunrise.

    Returns:
        datetime: The date and time of the sun event. If the sun event is not found (does not set or rise), returns `math.inf`.

    Raises:
        ValueError: If `rise_or_set` is not set correctly to either 'rise' or 'set'.
        ArithmeticError: If the sun event does not exist for the given location at the given date & time.
    """
    import islamic_times.astro_core as fast_astro

    if rise_or_set not in ['rise', 'set', 'sunrise', 'sunset']:
        raise ValueError("Invalid value for rise_or_set. Please use 'rise' or 'set'.")

    if rise_or_set in ["set", "sunset"]:
        event = 's'
    else:
        event = 'r'

    try:
        suntime: datetime = fast_astro.find_proper_suntime(observer_date.jd, 
                                       observer.latitude.decimal, observer.longitude.decimal, 
                                       observer.elevation.in_unit(DistanceUnits.METRE), observer.temperature, observer.pressure, 
                                       observer_date.utc_offset, angle.decimal, event)
    except:
        raise ArithmeticError("Sun event does not exist for the given location at the given date & time.")
    
    return suntime.replace(tzinfo=observer_date.date.tzinfo)