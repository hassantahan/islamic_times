import math
import numpy as np
from typing import List
from islamic_times import calculation_equations as ce
from islamic_times import time_equations as te

obliquity_terms = [
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

sun_nutation_arguments = [
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

sun_nutation_coefficients = [
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

# Chapter 22
def oblique_eq(jde: float) -> float:
    u = ((jde - te.J2000) / te.JULIAN_CENTURY) / 100

    eps = 23 + 26 / 60 + (21.448 / 3600)

    for i in range(10):
        eps += (obliquity_terms[i] / 3600) * (u ** (i + 1))

    return eps

# Chapter 22
def sun_nutation(jde: float) -> float:
    t = (jde - te.J2000) / te.JULIAN_CENTURY
    t2 = t ** 2
    t3 = t ** 3

    ta = [np.deg2rad(297.850363 + 445267.11148 * t - 0.0019142 * t2 + t3 / 189474.0),
          np.deg2rad(357.52772 + 35999.05034 * t - 0.0001603 * t2 - t3 / 300000.0),
          np.deg2rad(134.96298 + 477198.867398 * t + 0.0086972 * t2 + t3 / 56250.0),
          np.deg2rad(93.27191 + 483202.017538 * t - 0.0036825 * t2 + t3 / 327270),
          np.deg2rad(125.04452 - 1934.136261 * t + 0.0020708 * t2 + t3 / 450000.0)]

    for i in range(5):
        ta[i] = ce.bound_angle_rad(ta[i])

    dp, de = 0, 0

    for i in range(63):
        ang = 0
        for j in range(5):
            if sun_nutation_arguments[(i * 5) + j] != 0:
                ang += sun_nutation_arguments[(i * 5) + j] * ta[j]
        dp += (sun_nutation_coefficients[(i * 4) + 0] + sun_nutation_coefficients[(i * 4) + 1] * t) * math.sin(ang)
        de += (sun_nutation_coefficients[(i * 4) + 2] + sun_nutation_coefficients[(i * 4) + 3] * t) * math.cos(ang)

    deltaPsi, deltaEpsilon = dp / (3600.0 * 10000.0), de / (3600.0 * 10000.0)

    return deltaPsi, deltaEpsilon

# Chapter 25
def sunpos(jde: float, deltaT: float, local_latitude: float, local_longitude: float, temperature: float = 10, pressure: float = 101) -> List[float]:
    T = (jde - te.J2000) / te.JULIAN_MILLENNIUM
    T2 = T ** 2
    T3 = T ** 3
    T4 = T ** 4
    T5 = T ** 5

    L0 = 280.4664567 + (360007.6982779 * T) + (0.03032028 * T2) + (T3 / 49931) - (T4 / 15300) - (T5 / 2000000)
    L0 = ce.bound_angle_deg(L0)

    M = 357.52911 + 359990.50340 * T - 0.001603 * T2 - T3 / 30000
    M = ce.bound_angle_deg(M)

    e = 0.016708634 - 0.00042037 * T - 0.000001267 * T2

    C = (1.914602 - 0.04817 * T - 0.000014 * T2) * ce.sin(M) + \
        (0.019993 - 0.000101 * T) * ce.sin(2 * M) + \
        0.000289 * ce.sin(3 * M)

    sunLong = L0 + C
    sunAnomaly = M + C

    sunR = (1.000001018 * (1 - e ** 2)) / (1 + (e * ce.cos(sunAnomaly)))

    omega = 125.04452 - 19341.36261 * T + 0.020708 * T2 + T3 / 45000
    Lambda = sunLong - 0.00569 - 0.00478 * ce.sin(omega)

    nut = sun_nutation(jde)
    delta_epsilon = nut[1]
    epsilon0 = oblique_eq(jde)
    epsilon = epsilon0 + delta_epsilon
    

    # Right ascension (RA) & declination calculations
    alpha = ce.bound_angle_deg(np.rad2deg(math.atan2(ce.cos(epsilon0) * ce.sin(sunLong), ce.cos(sunLong))))
    delta = np.rad2deg(math.asin(ce.sin(epsilon0) * ce.sin(sunLong)))

    # Adjust RA and declination to find their apparents
    alphaApp = ce.bound_angle_deg(np.rad2deg(math.atan2(ce.cos(epsilon) * ce.sin(Lambda), ce.cos(Lambda))))
    deltaApp = np.rad2deg(math.asin(ce.sin(epsilon) * ce.sin(Lambda)))

    # Local Hour Angle calculation
    # Start by calculating Mean Greenwich Sidereal Time
    greenwich_hour_angle = ce.bound_angle_deg(te.greenwich_mean_sidereal_time(jde - deltaT / 86400))

    # Attain the sun's nutation in the longitude in DMS
    nut_long_dms = ce.decimal_to_dms(nut[0])
    
    # Make correction for the apparent sidereal time according to pg. 88
    st_correction = (nut_long_dms[2] + nut_long_dms[1] * 60 + nut_long_dms[0] * 3600) * ce.cos(epsilon) / 15
    greenwich_hour_angle += (st_correction / 240)

    # Local Hour angle is then simply the GHA minus the apparent RA adjusted for the local longitude
    local_hour_angle = ce.bound_angle_deg(greenwich_hour_angle + local_longitude - alphaApp)

    # Altitude & Azimuth calculations
    altitude = np.rad2deg(math.asin(ce.sin(local_latitude) * ce.sin(delta) + ce.cos(local_latitude) * ce.cos(delta) * ce.cos(local_hour_angle))) 
    azimuth = np.rad2deg(np.arctan2(-1 * ce.cos(delta) * ce.sin(local_hour_angle), ce.sin(delta) * ce.cos(local_latitude) - ce.cos(delta) * ce.sin(local_latitude) * ce.cos(local_hour_angle))) % 360

    # Correct for atmospheric refraction (taken from https://en.wikipedia.org/wiki/Atmospheric_refraction)
    refraction = 1.02 / ce.tan(altitude + 10.3 / (altitude + 5.11)) * pressure / 101 * 283 / (273 + temperature)
    #altitude += refraction / 60

    # Correct for parallax 
    eh_parallax = np.rad2deg(np.arcsin(te.EARTH_RADIUS_KM / (sunR * te.ASTRONOMICAL_UNIT) ))
    parallax_correction = -np.rad2deg(np.deg2rad(eh_parallax) * ce.cos(altitude))
    altitude += parallax_correction


    return [
        L0,                 # 0; mean longitude
        M,                  # 1; mean anomaly
        e,                  # 2; eccentricity
        C,                  # 3; sun's centre
        sunLong,            # 4; true longitude
        sunAnomaly,         # 5; true anomaly
        sunR,               # 6; distance to sun (in AU)
        Lambda,             # 7; apparent longitude
        alpha,              # 8; right ascension
        delta,              # 9; declination
        alphaApp,           # 10; apparent right ascension
        deltaApp,           # 11; apparent declination
        epsilon0,           # 12; mean obliquity of the ecliptic
        epsilon,            # 13; true obliquity of the ecliptic
        omega,              # 14; Longitude of the ascending node of the Moon’s mean orbit on the ecliptic
        altitude,           # 15; apparent
        azimuth,            # 16; apparent
        local_hour_angle    # 17
    ]

def equation_of_time(jde: float, deltaT: float, local_latitude: float, local_longitude: float) -> float:
    sun_factors = sunpos(jde, deltaT, local_latitude, local_longitude)
    L0 = sun_factors[0]
    nut = sun_nutation(jde)
    epsilon = sun_factors[13]
    deltaPsi = nut[0]
    alpha = sun_factors[10]

    # Only adjust alpha when it appears to have wrapped
    if L0 > 300 and alpha < 50:
        alpha += 360
    elif L0 < 50 and alpha > 300:
        alpha -= 360

    E = (L0 - 0.0057183 - alpha + deltaPsi * ce.cos(epsilon)) * 4  # convert degrees to minutes

    return E

def solar_hour_angle(latitude: float, declination: float, angle: float = 0.8333):
    dec = np.deg2rad(declination)
    lat = np.deg2rad(latitude)

    num = -1 * np.sin(np.deg2rad(angle)) - np.sin(lat) * np.sin(dec)
    denom = np.cos(lat) * np.cos(dec)

    # At extreme latitudes at different times of the year, some angles are not possible to achieve
    # TODO: Warnings
    ratio = np.clip(num / denom, -1, 1)

    solar_angle = np.arccos(ratio)

    return np.rad2deg(solar_angle)

# -1 for Sunrise
# 1 for Sunset
def sunrise_sunset(set_or_rise: int, hour_angle: float) -> float:
    if set_or_rise not in [1, -1]:
        raise ValueError("set_or_rise from sun_equations.sunrise_sunset() accepts only -1 or 1 for sunset or sunrise respectively.")
    
    hours_offset_from_noon = hour_angle / 15

    offset = 12 + set_or_rise * hours_offset_from_noon
    return  offset