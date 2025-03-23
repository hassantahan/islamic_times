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
from datetime import datetime, timedelta
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
    '''
    A class to compute the position of the Sky in the sky based on given astronomical parameters.

    Attributes:
        jde (float): The Julian Ephemeris Day.
        deltaT (float): The difference between Terrestrial Time and Universal Time (ΔT).
        local_latitude (float): The observer's latitude in degrees.
        local_longitude (float): The observer's longitude in degrees.
        temperature (float): The observer's ambient temperature in degrees Celcius.
        pressure (float): The observer's ambient pressure in kPa.
        mean_longitude (float): The mean longitude of the Sun in degrees.
        mean_anomaly (float): The mean anomaly of the Sun in degrees.
        earth_orbit_eccentricity (float): The eccentricity of Earth's orbit.
        sun_centre (float): The Sun's equation of center in degrees.
        true_longitude (float): The Sun's true longitude in degrees.
        true_anomaly (float): The Sun's true anomaly in degrees.
        geocentric_distance (float): The distance from the Earth to the Sun in Astronomical Units (AU).
        omega (float): The longitude of the ascending node of the Moon's orbit in degrees.
        apparent_longitude (float): The apparent longitude of the Sun in degrees.
        nutation (Tuple[float, float]): Nutation in longitude and obliquity in degrees.
        delta_obliquity (float): Change in Earth's obliquity due to nutation.
        mean_obliquity (float): The mean obliquity of the ecliptic in degrees.
        true_obliquity (float): The true obliquity of the ecliptic in degrees.
        true_right_ascension (float): The true right ascension of the Sun in degrees.
        true_declination (float): The true declination of the Sun in degrees.
        apparent_right_ascension (float): The apparent right ascension of the Sun in degrees.
        apparent_declination (float): The apparent declination of the Sun in degrees.
        local_hour_angle (float): The local hour angle of the Sun in degrees.
        eh_parallax (float): The equatorial horizontal parallax of the Sun in degrees.
        topocentric_ascension (float): The topocentric right ascension of the Sun in degrees.
        topocentric_declination (float): The topocentric declination of the Sun in degrees.
        topocentric_local_hour_angle (float): The topocentric local hour angle in degrees.
        altitude (float): The altitude of the Sun above the horizon in degrees.
        azimuth (float): The azimuth angle of the Sun in degrees.

    Methods:
        `calculate()`: Computes various positional attributes of the Sun, including its right ascension, declination, altitude, and azimuth.

    Notes:
        - All attributes after `elev` are only computed and available after the `calculate()` method is called.
        - The `calculate()` method is called automatically by the `sunpos()` method.
    '''
        
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
        self.true_right_ascension: float = (np.rad2deg(np.arctan2(ce.cos(self.mean_obliquity) * ce.sin(self.true_longitude), ce.cos(self.true_longitude)))) % 360
        self.true_declination: float = np.rad2deg(np.arcsin(ce.sin(self.mean_obliquity) * ce.sin(self.true_longitude)))

        # Adjust RA and declination to find their apparents
        self.apparent_right_ascension: float = (np.rad2deg(np.arctan2(ce.cos(self.true_obliquity + 0.00256 * ce.cos(self.omega)) * ce.sin(self.apparent_longitude), ce.cos(self.apparent_longitude)))) % 360
        self.apparent_declination: float = np.rad2deg(np.arcsin(ce.sin(self.true_obliquity + 0.00256 * ce.cos(self.omega)) * ce.sin(self.apparent_longitude)))

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

        self.altitude = np.rad2deg(np.arcsin(ce.sin(self.local_latitude) * ce.sin(self.topocentric_declination) + ce.cos(self.local_latitude) * ce.cos(self.topocentric_declination) * ce.cos(self.topocentric_local_hour_angle))) 
        self.azimuth = np.rad2deg(np.arctan2(-1 * ce.cos(self.topocentric_declination) * ce.sin(self.topocentric_local_hour_angle), ce.sin(self.topocentric_declination) * ce.cos(self.local_latitude) \
                                          - ce.cos(self.topocentric_declination) * ce.sin(self.local_latitude) * ce.cos(self.topocentric_local_hour_angle))) % 360

        # # Geocentric Altitude & Azimuth calculations (disabled)
        # self.altitude: float = np.rad2deg(np.arcsin(ce.sin(self.local_latitude) * ce.sin(self.true_declination) + ce.cos(self.local_latitude) * ce.cos(self.true_declination) * ce.cos(self.local_hour_angle))) 
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

def solar_hour_angle(latitude: float, declination: float, angle: float = 5 / 6) -> float:
    '''
    Calculate the solar hour angle for a given latitude, solar declination, and angle. The angle is the zenith angle of the sun at the horizon. The default value is 0.8333° representing visible sunset.

    Parameters:
        latitude (float): The observer's latitude in degrees.
        declination (float): The solar declination in degrees.
        angle (float): The zenith angle of the sun at the horizon. Default is 0.8333° (5/6) representing visible sunset/sunrise.

    Returns:
        float: The solar hour angle in degrees, but 'np.inf' if the angle is not possible to achieve.
    '''
    dec = np.deg2rad(declination)
    lat = np.deg2rad(latitude)
    ang  =np.deg2rad(angle)

    num =  np.sin(-1 * ang) - np.sin(lat) * np.sin(dec)
    denom = np.cos(lat) * np.cos(dec)

    # At extreme latitudes at different times of the year, some angles are not possible to achieve
    # TODO: Warnings
    ratio = np.clip(num / denom, -1, 1)

    if ratio == 1 or ratio == -1:
        return np.inf
    else:
        return np.rad2deg(np.arccos(ratio))

def find_sun_transit(date: datetime, lat: float, long: float, elev: float, utc_offset: float, angle: float = 5 / 6) -> Tuple[datetime, datetime, datetime]:
    """
    Calculate the times of sunrise, sunset, and its transit for a given date and observer coordinates. See Chapter 15 of *Astronomical Algorithms* for more information.

    Parameters:
        date (datetime): The date to calculate the moonset for.
        lat (float): The observer's latitude (°).
        long (float): The observer's longitude (°).
        elev (float): The observer's elevation above sea level (m).
        utc_offset (float): The observer's difference from UTC (hours).

    Returns:
        datetime: The time of moonset.
    """

    # First find the Year Month Day at UT 0h from JDE
    ymd = datetime(date.year, date.month, date.day)
    new_jd = te.gregorian_to_jd(date) - te.fraction_of_day(date)
    new_deltaT = te.delta_t_approx(ymd.year, ymd.month)
    new_jde = new_jd + new_deltaT / 86400

    # Calculate new sun params with the new_jd
    sun_params: List[Sun] = []
    for i in range(3):
        ymd_temp = te.jd_to_gregorian(new_jd + i - 1, utc_offset)
        delT_temp = te.delta_t_approx(ymd_temp.year, ymd_temp.month)
        sun_params.append(sunpos(new_jde + i - 1, delT_temp, lat, long))

    # GMST
    sidereal_time = te.greenwich_mean_sidereal_time(new_jd)

    # Compute m0 and m2 without wrapping
    # Transit
    m0 = (sun_params[1].topocentric_ascension - long - sidereal_time) / 360

    # Minor corrective steps thru iteration
    for _ in range(3):
        little_theta_zero = (sidereal_time + 360.985647 * m0) % 360
        n = m0 + new_deltaT / 86400
        interpolated_sun_ra = ce.interpolation(n, sun_params[0].topocentric_ascension, 
                                                sun_params[1].topocentric_ascension, 
                                                sun_params[2].topocentric_ascension)

        solar_local_hour_angle = (little_theta_zero - (-long) - interpolated_sun_ra) % 360
        m0 -= solar_local_hour_angle / 360

    # Compute final rise, transit, and set time by adding days to base date
    m0 %= 1
    sun_transit_dt = datetime(ymd.year, ymd.month, ymd.day) + timedelta(days=m0) - timedelta(hours=utc_offset)
                                                                                        
    return sun_transit_dt

def sunrise_or_sunset(date: datetime, lat: float, long: float, elev: float, utc_offset: float, rise_or_set: str, angle: float = 5/6) -> datetime:
    """
    Calculate either the sunrise or the sunset time for a given date and observer coordinates.
    This function computes only the requested event (sunrise or sunset) without calculating both.
    
    Parameters:
        date (datetime): The date for which to calculate the event.
        lat (float): Observer's latitude in degrees.
        long (float): Observer's longitude in degrees.
        elev (float): Observer's elevation above sea level in meters.
        utc_offset (float): The observer's time difference from UTC in hours.
        rise_or_set (str): 'rise' or 'sunrise' to calculate sunrise; 'set' or 'sunset' to calculate sunset.
        angle (float): The standard altitude of the sun (default is 5/6°).
        
    Returns:
        datetime: The computed time of the requested event (sunrise or sunset).
        
    Raises:
        ValueError: If the value of rise_or_set is not recognized.
    """

    if rise_or_set not in ['rise', 'set', 'sunrise', 'sunset']:
        raise ValueError("Invalid value for rise_or_set. Please use 'rise' or 'set'.")

    # Base date at UT 0h and associated Julian day calculations.
    ymd = datetime(date.year, date.month, date.day)
    new_jd = te.gregorian_to_jd(date) - te.fraction_of_day(date)
    new_deltaT = te.delta_t_approx(ymd.year, ymd.month)
    new_jde = new_jd + new_deltaT / 86400

    # Calculate sun parameters for three consecutive days.
    sun_params: List[Sun] = []
    for i in range(3):
        ymd_temp = te.jd_to_gregorian(new_jd + i - 1, utc_offset)
        delT_temp = te.delta_t_approx(ymd_temp.year, ymd_temp.month)
        sun_params.append(sunpos(new_jde + i - 1, delT_temp, lat, long))

    # Compute H0: the hour angle corresponding to the desired altitude.
    h_zero = -angle
    cosH_zero = (ce.sin(h_zero) - ce.sin(lat) * ce.sin(sun_params[1].topocentric_declination)) / (
                ce.cos(lat) * ce.cos(sun_params[1].topocentric_declination))
    H_zero = np.rad2deg(np.arccos(cosH_zero))

    # If H_zero is not defined (NaN), the event does not occur (e.g., polar day/night).
    if np.isnan(H_zero):
        return np.inf

    # Compute Greenwich Mean Sidereal Time (GMST) at new_jd.
    sidereal_time = te.greenwich_mean_sidereal_time(new_jd)

    # Compute the transit estimate m0.
    m0 = (sun_params[1].topocentric_ascension - long - sidereal_time) / 360

    # Choose which event to compute.
    event = rise_or_set.lower()
    if event in ['rise', 'sunrise']:
        # Initial estimate for sunrise.
        m_event = m0 - H_zero / 360
    else:
        # Initial estimate for sunset.
        m_event = m0 + H_zero / 360

    # Iteratively refine m0 (transit) and m_event (rise or set).
    # We use three iterations which are typically sufficient.
    for _ in range(3):
        # --- Update the event (sunrise or sunset) estimate ---
        theta_event = (sidereal_time + 360.985647 * m_event) % 360
        n_event = m_event + new_deltaT / 86400
        interp_dec_event = ce.interpolation(n_event,
                                            sun_params[0].topocentric_declination,
                                            sun_params[1].topocentric_declination,
                                            sun_params[2].topocentric_declination)
        interp_ra_event = ce.interpolation(n_event,
                                           sun_params[0].topocentric_ascension,
                                           sun_params[1].topocentric_ascension,
                                           sun_params[2].topocentric_ascension)
        local_hour_angle_event = (theta_event - (-long) - interp_ra_event) % 360
        sun_alt = np.rad2deg(np.arcsin(ce.sin(lat) * ce.sin(interp_dec_event) +
                                       ce.cos(lat) * ce.cos(interp_dec_event) *
                                       ce.cos(local_hour_angle_event)))
        
        # Correct m_event using the difference between computed altitude and desired h_zero.
        deltaM = (sun_alt - h_zero) / (360 * ce.cos(interp_dec_event) * ce.cos(lat) * ce.sin(local_hour_angle_event))
        m_event += deltaM

    # Convert the fractional day to a datetime object, adjusting for the UTC offset.
    event_dt = datetime(ymd.year, ymd.month, ymd.day) + timedelta(days=m_event) - timedelta(hours=utc_offset)
    return event_dt

def find_proper_suntime(true_date: datetime, latitude, longitude, elevation, utc_offset, rise_or_set: str, angle: float = 5 / 6) -> datetime:
        """
        Determines the proper local sunset time.

        Adjusts the calculated sunset time to account for local UTC differences.

        Parameters:
            date (datetime): The reference date.

        Returns:
            datetime: Adjusted sunset time. If sunset is not found, returns `datetime.min`.
        """

        temp_utc_offset = np.floor(longitude / 15) - 1
        temp_suntime = sunrise_or_sunset(true_date, latitude, longitude, elevation, utc_offset, rise_or_set, angle)
        date_doy = true_date.timetuple().tm_yday
        if temp_suntime == np.inf:
            return np.inf

        i = 1
        while(True):
            temp_suntime_doy = (temp_suntime + timedelta(hours=temp_utc_offset, minutes=-20)).timetuple().tm_yday
            if (temp_suntime_doy < date_doy and temp_suntime.year == true_date.year) or ((temp_suntime + timedelta(hours=temp_utc_offset)).year < true_date.year):
                temp_suntime = sunrise_or_sunset(true_date + timedelta(days=i), latitude, longitude, elevation, utc_offset, rise_or_set, angle)
                i += 1
            else: 
                return temp_suntime