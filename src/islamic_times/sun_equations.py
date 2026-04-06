"""
Solar ephemeris and event-time calculations.

This module provides low-level solar calculations used by the higher-level prayer
and calendar APIs, including:

- solar position and topocentric coordinates,
- solar nutation and obliquity terms,
- equation-of-time offsets,
- transit and rise/set event-time solvers.

Most heavy numerical work is delegated to the native C extension. Prefer using
the higher-level public APIs unless direct control is required.
"""

import math
import numpy as np
from typing import List, Tuple
from datetime import datetime, timedelta
from dataclasses import dataclass, replace
from islamic_times import time_equations as te
from islamic_times import calculation_equations as ce
from islamic_times.it_dataclasses import Angle, DateTimeInfo, Distance, DistanceUnits, ObserverInfo, RightAscension
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
    """Container for solar ephemeris values for one observer and instant."""

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
    """Compute mean ecliptic obliquity for a Julian Ephemeris Day.

    Parameters
    ----------
    jde : float
        Julian Ephemeris Day.

    Returns
    -------
    Angle
        Mean obliquity in degrees.

    Notes
    -----
    This Python implementation is deprecated in favor of the native C
    implementation exposed via ``islamic_times.astro_core``.
    """
    warn('This particular function will no longer be supported in python. The proper function is in the C extension but cannot be called from python. Use "islamic_times.astro_core.compute_sun()" or "islamic_times.sun_equations.sunpos()" instead.', DeprecationWarning)
    from islamic_times._legacy_py_impl import sun_equations as legacy_se
    return legacy_se.oblique_eq(jde)

# Chapter 22
def sun_nutation(jde: float) -> Tuple[Angle, Angle]:
    """Compute solar nutation in longitude and obliquity for a JDE value.

    Parameters
    ----------
    jde : float
        Julian Ephemeris Day.

    Returns
    -------
    tuple[Angle, Angle]
        Nutation in longitude and nutation in obliquity, both in degrees.

    Notes
    -----
    This Python implementation is deprecated in favor of the native C
    implementation exposed via ``islamic_times.astro_core``.
    """
    warn('This particular function will no longer be supported in python. The proper function is in the C extension but cannot be called from python. Use "islamic_times.astro_core.compute_sun()" or "islamic_times.sun_equations.sunpos()" instead.', DeprecationWarning)
    from islamic_times._legacy_py_impl import sun_equations as legacy_se
    return legacy_se.sun_nutation(jde)

# Chapter 25
def sunpos(observer_date: DateTimeInfo, observer: ObserverInfo) -> Sun:
    """Compute the solar ephemeris for one observer instant.

    Parameters
    ----------
    observer_date : DateTimeInfo
        Observer date/time metadata, including JD/JDE and Delta T.
    observer : ObserverInfo
        Observer latitude, longitude, elevation, and weather context.

    Returns
    -------
    Sun
        Solar position container with orbital, equatorial, topocentric, and
        horizontal fields.

    Notes
    -----
    The implementation delegates to the native extension
    ``islamic_times.astro_core.compute_sun``.
    Temperature and pressure are currently forwarded for API compatibility;
    the active solar apparent-altitude path does not apply weather-based
    refraction correction.
    """
    import islamic_times.astro_core as fast_astro
    the_sun: Sun = fast_astro.compute_sun(
        observer_date.jde_tt,
        observer_date.jd_ut1,
        observer.latitude.decimal,
        observer.longitude.decimal,
        observer.elevation.value,
        observer.temperature,
        observer.pressure,
    )

    return the_sun

def equation_of_time(deltaPsi: float, L0: float, epsilon: float, alpha: float) -> float:
    """Compute equation-of-time offset from solar terms.

    Parameters
    ----------
    deltaPsi : float
        Nutation in longitude, in degrees.
    L0 : float
        Mean solar longitude, in degrees.
    epsilon : float
        True obliquity, in degrees.
    alpha : float
        Apparent solar right ascension, in degrees.

    Returns
    -------
    float
        Equation-of-time offset in minutes.
    """

    # Only adjust alpha when it appears to have wrapped
    if L0 > 300 and alpha < 50:
        alpha += 360
    elif L0 < 50 and alpha > 300:
        alpha -= 360

    E = (L0 - 0.0057183 - alpha + deltaPsi * ce.cos(epsilon)) * 4  # convert degrees to minutes

    return E

def find_sun_transit(observer_date: DateTimeInfo, observer: ObserverInfo) -> datetime:
    """Compute the Sun's meridian transit (culmination) time.

    Parameters
    ----------
    observer_date : DateTimeInfo
        Observer date/time metadata.
    observer : ObserverInfo
        Observer location and environment values.

    Returns
    -------
    datetime
        Local datetime of solar transit on the reference civil date.

    Notes
    -----
    Weather fields are accepted for cross-body API parity but are not currently
    used by the active solar apparent-altitude correction path.
    """
    import islamic_times.astro_core as fast_astro
    sun_transit_dt: datetime = fast_astro.find_sun_transit(observer_date.jd, observer_date.deltaT, 
                                       observer.latitude.decimal, observer.longitude.decimal, 
                                       observer.elevation.in_unit(DistanceUnits.METRE), observer.temperature, observer.pressure, 
                                       observer_date.utc_offset)
                                                                                      
    return sun_transit_dt.replace(tzinfo=observer_date.date.tzinfo)

def sunrise_or_sunset(observer_date: DateTimeInfo, observer: ObserverInfo, rise_or_set: str, angle: Angle = Angle(5 / 6)) -> datetime | float:
    """Compute either sunrise or sunset for a reference date.

    Parameters
    ----------
    observer_date : DateTimeInfo
        Observer date/time metadata.
    observer : ObserverInfo
        Observer location and environment values.
    rise_or_set : str
        Event selector: ``"rise"``/``"sunrise"`` or ``"set"``/``"sunset"``.
    angle : Angle, default Angle(5 / 6)
        Target altitude magnitude in degrees (typically 0.8333 degrees for the
        visible horizon event).

    Returns
    -------
    datetime | float
        Event datetime if found; ``math.inf`` when the requested event does not
        occur at the location/date.

    Raises
    ------
    ValueError
        If ``rise_or_set`` is not one of the supported values.

    Notes
    -----
    This pure-Python routine is retained for compatibility and is deprecated in
    favor of native extension solvers.
    """
    warn('This particular function will no longer be supported in python. The proper function is in the C extension but cannot be called from python. Use "islamic_times.astro_core.find_proper_suntime()" or "islamic_times.sun_equations.find_proper_suntime()" instead.', DeprecationWarning)
    from islamic_times._legacy_py_impl import sun_equations as legacy_se
    return legacy_se.sunrise_or_sunset(observer_date, observer, rise_or_set, angle)

def find_proper_suntime(observer_date: DateTimeInfo, observer: ObserverInfo, rise_or_set: str, angle: Angle = Angle(5 / 6)) -> datetime:
    """Return the local sunrise/sunset aligned to the reference civil day.

    Parameters
    ----------
    observer_date : DateTimeInfo
        Observer date/time metadata.
    observer : ObserverInfo
        Observer location and environment values.
    rise_or_set : str
        Event selector: ``"rise"``/``"sunrise"`` or ``"set"``/``"sunset"``.
    angle : Angle, default Angle(5 / 6)
        Target altitude magnitude in degrees.

    Returns
    -------
    datetime
        Local datetime of the requested solar event.

    Raises
    ------
    ValueError
        If ``rise_or_set`` is not recognized.
    ArithmeticError
        If the event does not exist for the given location/date.

    Notes
    -----
    Weather fields are currently accepted for API consistency; solar event
    solving does not yet apply weather-based apparent-altitude correction.
    """
    import islamic_times.astro_core as fast_astro

    if rise_or_set not in ['rise', 'set', 'sunrise', 'sunset']:
        raise ValueError("Invalid value for rise_or_set. Please use 'rise' or 'set'.")

    if rise_or_set in ["set", "sunset"]:
        event = 's'
    else:
        event = 'r'

    try:
        suntime: datetime = fast_astro.find_proper_suntime(
            observer_date.jd,
            observer.latitude.decimal,
            observer.longitude.decimal,
            observer.elevation.in_unit(DistanceUnits.METRE),
            observer.temperature,
            observer.pressure,
            observer_date.utc_offset,
            angle.decimal,
            event,
        )
    except ValueError as exc:
        raise ArithmeticError("Sun event does not exist for the given location at the given date & time.") from exc
    
    return suntime.replace(tzinfo=observer_date.date.tzinfo)
