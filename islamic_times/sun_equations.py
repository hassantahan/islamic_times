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
from dataclasses import dataclass
from islamic_times import calculation_equations as ce
from islamic_times import time_equations as te

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

@dataclass
class Sun:
    jde: float
    deltaT: float
    local_latitude: float
    local_longitude: float
    temperature: float = 10
    pressure: float = 101

    def calculate(self):
        t = (self.jde - te.J2000) / te.JULIAN_MILLENNIUM
        t2 = t ** 2
        t3 = t ** 3
        t4 = t ** 4
        t5 = t ** 5

        # Coordinates
        self.mean_longitude: float = (280.4664567 + (360007.6982779 * t) + (0.03032028 * t2) + (t3 / 49931) - (t4 / 15300) - (t5 / 2000000)) % 360
        self.mean_anomaly: float = (357.52911 + 359990.50340 * t - 0.001603 * t2 - t3 / 30000) % 360

        self.earth_oribit_eccentricity: float = 0.016708634 - 0.00042037 * t - 0.000001267 * t2
        
        self.sun_centre: float = (1.914602 - 0.04817 * t - 0.000014 * t2) * ce.sin(self.mean_anomaly) + \
                            (0.019993 - 0.000101 * t) * ce.sin(2 * self.mean_anomaly) + \
                            0.000289 * ce.sin(3 * self.mean_anomaly)

        self.true_longitude: float = self.mean_longitude + self.sun_centre
        self.true_anomaly: float = self.mean_anomaly + self.sun_centre

        self.geocentric_distance: float = (1.000001018 * (1 - self.earth_oribit_eccentricity ** 2)) / (1 + (self.earth_oribit_eccentricity * ce.cos(self.true_anomaly)))
        self.omega: float = 125.04452 - 19341.36261 * t + 0.020708 * t2 + t3 / 45000

        self.apparent_longitude: float = self.true_longitude - 0.00569 - 0.00478 * ce.sin(self.omega)

        # TODO: Chapter 26 for referencing to J2000.0 standard

        self.nutation: Tuple[float, float] = sun_nutation(self.jde)
        self.delta_obliquity: float = self.nutation[1]
        self.mean_obliquity: float = oblique_eq(self.jde)
        self.true_obliquity: float = self.mean_obliquity + self.delta_obliquity

        # True Right ascension (RA) & declination calculations
        self.true_right_ascension: float = (np.rad2deg(math.atan2(ce.cos(self.mean_obliquity) * ce.sin(self.true_longitude), ce.cos(self.true_longitude)))) % 360
        self.true_declination: float = np.rad2deg(math.asin(ce.sin(self.mean_obliquity) * ce.sin(self.true_longitude)))

        # Adjust RA and declination to find their apparents
        self.apparent_right_ascension: float = (np.rad2deg(math.atan2(ce.cos(self.true_obliquity + 0.00256 * ce.cos(self.omega)) * ce.sin(self.apparent_longitude), ce.cos(self.apparent_longitude)))) % 360
        self.apparent_declination: float = np.rad2deg(math.asin(ce.sin(self.true_obliquity + 0.00256 * ce.cos(self.omega)) * ce.sin(self.apparent_longitude)))

        greenwich_hour_angle: float = te.greenwich_mean_sidereal_time(self.jde - self.deltaT / 86400) % 360
        nutation_in_longitude_dms = ce.decimal_to_dms(self.nutation[0])
        
        # Make correction for the apparent sidereal time according to pg. 88
        st_correction: float = (nutation_in_longitude_dms[2] + nutation_in_longitude_dms[1] * 60 + nutation_in_longitude_dms[0] * 3600) * ce.cos(self.true_obliquity) / 15
        greenwich_hour_angle += (st_correction / 240)

        # Local Hour angle is then simply the GHA minus the apparent RA adjusted for the local longitude
        self.local_hour_angle: float = (greenwich_hour_angle + self.local_longitude - self.apparent_right_ascension) % 360

        # Topocentric R.A., Declination, and LHA
        self.eh_parallax = np.rad2deg(np.arcsin(ce.sin(8.794 / 3600) / self.geocentric_distance))
        self.topocentric_ascension, self.topocentric_declination = ce.correct_ra_dec(self.apparent_right_ascension, self.apparent_declination, self.local_hour_angle, self.eh_parallax, self.local_latitude, 76)
        self.topocentric_local_hour_angle = (greenwich_hour_angle + self.local_longitude - self.topocentric_ascension) % 360

        self.altitude = np.rad2deg(math.asin(ce.sin(self.local_latitude) * ce.sin(self.topocentric_declination) + ce.cos(self.local_latitude) * ce.cos(self.topocentric_declination) * ce.cos(self.topocentric_local_hour_angle))) 
        self.azimuth = np.rad2deg(np.arctan2(-1 * ce.cos(self.topocentric_declination) * ce.sin(self.topocentric_local_hour_angle), ce.sin(self.topocentric_declination) * ce.cos(self.local_latitude) \
                                          - ce.cos(self.topocentric_declination) * ce.sin(self.local_latitude) * ce.cos(self.topocentric_local_hour_angle))) % 360

        # # Geocentric Altitude & Azimuth calculations (disabled)
        # self.altitude: float = np.rad2deg(math.asin(ce.sin(self.local_latitude) * ce.sin(self.true_declination) + ce.cos(self.local_latitude) * ce.cos(self.true_declination) * ce.cos(self.local_hour_angle))) 
        # self.azimuth: float = np.rad2deg(np.arctan2(-1 * ce.cos(self.true_declination) * ce.sin(self.local_hour_angle), ce.sin(self.true_declination) * ce.cos(self.local_latitude) - ce.cos(self.true_declination) \
        #                                             * ce.sin(self.local_latitude) * ce.cos(self.local_hour_angle))) % 360
        
        # # Correct for atmospheric refraction (taken from https://en.wikipedia.org/wiki/Atmospheric_refraction)
        # # Currently disabled
        # refraction = 1.02 / ce.tan(altitude + 10.3 / (altitude + 5.11)) * pressure / 101 * 283 / (273 + temperature)
        # self.altitude += refraction / 60

# Chapter 22
def oblique_eq(jde: float) -> float:
    '''
    Calculate the obliquity of the ecliptic for a given Julian Ephemeris Day. See Chapter 22 of *Astronomical Algorthims* for more information.

    Parameters:
        jde (float): The Julian Ephemeris Day.

    Returns:
        float: The obliquity of the ecliptic in degrees.
    '''
    u = ((jde - te.J2000) / te.JULIAN_CENTURY) / 100

    eps = 23 + 26 / 60 + (21.448 / 3600)

    powers = np.power(u, np.arange(1, len(__OBLIQUITY_TERMS) + 1))
    eps += np.dot(__OBLIQUITY_TERMS, powers) / 3600

    return eps

# Chapter 22
def sun_nutation(jde: float) -> Tuple[float, float]:
    '''
    Calculate the sun's nutation for a given Julian Ephemeris Day. See Chapter 22 of *Astronomical Algorthims* for more information.

    Parameters:
        jde (float): The Julian Ephemeris Day.

    Returns:
        tuple (Tuple[float, float]): The nutation in longitude and obliquity.
    '''

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

    return (deltaPsi, deltaEpsilon)

# Chapter 25
def sunpos(jde: float, deltaT: float, local_latitude: float, local_longitude: float, temperature: float = 10, pressure: float = 101) -> Sun:
    '''
    Calculate the various solar positional parameters for a given Julian Ephemeris Day, ΔT, and observer coordinates. See Chapter 25 of the *Astronomical Algorthims* for more information.

    Parameters:
        jde (float): The Julian Ephemeris Day.
        deltaT (float): The difference between Terrestrial Dynamical Time and Universal Time (ΔT = TT - UT1) in seconds.
        local_latitude (float): The observer's latitude in degrees.
        local_longitude (float): The observer's longitude in degrees.
        temperature (float): The observer's temperature in degrees Celsius.
        pressure (float): The observer's pressure in kPa.

    Returns:
        Sun (obj): A `Sun` object that contains various attributes that describe its position. 

    Notes: 
    - The temperature and pressure are used for atmospheric refraction calculations. Currently, this feature is disabled.
    '''
    
    the_sun = Sun(
        jde,
        deltaT,
        local_latitude,
        local_longitude
    )

    the_sun.calculate()
    
    return the_sun

def equation_of_time(deltaPsi: float, L0: float, epsilon: float, alpha: float) -> float:
    '''
    Calculate the equation of time offset for a given Julian Ephemeris Day, ΔT, and observer coordinates. See Chapter 28 of _Astronomical Algorthims_ for more information.

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

def solar_hour_angle(latitude: float, declination: float, angle: float = 0.8333) -> float:
    '''
    Calculate the solar hour angle for a given latitude, solar declination, and angle. The angle is the zenith angle of the sun at the horizon. The default value is 0.8333° representing visible sunset.

    Parameters:
        latitude (float): The observer's latitude in degrees.
        declination (float): The solar declination in degrees.
        angle (float): The zenith angle of the sun at the horizon. Default is 0.8333° representing visible sunset.

    Returns:
        float: The solar hour angle in degrees, but 'np.inf' if the angle is not possible to achieve.
    '''
    dec = np.deg2rad(declination)
    lat = np.deg2rad(latitude)

    num = -1 * np.sin(np.deg2rad(angle)) - np.sin(lat) * np.sin(dec)
    denom = np.cos(lat) * np.cos(dec)

    # At extreme latitudes at different times of the year, some angles are not possible to achieve
    # TODO: Warnings
    ratio = np.clip(num / denom, -1, 1)

    if ratio == 1 or ratio == -1:
        return np.inf
    else:
        solar_angle = np.arccos(ratio)

    return np.rad2deg(solar_angle)

# -1 for Sunrise
# 1 for Sunset
def sunrise_sunset(set_or_rise: int, hour_angle: float) -> float:
    '''
    Calculate the time of sunrise or sunset for a given hour angle. The hour angle is the angle between the observer's meridian and the sun's position. The hour angle is positive for sunset and negative for sunrise.

    Parameters:
        set_or_rise (int): The value -1 for sunset and 1 for sunrise.
        hour_angle (float): The hour angle in degrees.

    Raises:
        ValueError: If the 'set_or_rise' is not -1 or 1.

    Returns:
        float: The time of sunrise or sunset in hours from noon. Returns 'np.inf' if the angle is not possible to achieve
    '''
    
    if set_or_rise not in [1, -1]:
        raise ValueError("'set_or_rise' from sun_equations.sunrise_sunset() accepts only -1 or 1 for sunset or sunrise respectively.")
    
    if hour_angle == np.inf:
        return np.inf
    
    hours_offset_from_noon = hour_angle / 15

    offset = 12 + set_or_rise * hours_offset_from_noon
    return  offset