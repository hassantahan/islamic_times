import numpy as np
import math
from time_equations import *

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

def oblique_eq(julian_day):
    u = ((julian_day - J2000) / JULIAN_CENTURY) / 100

    eps = 23 + 26 / 60 + (21.448 / 3600)

    for i in range(10):
        eps += (obliquity_terms[i] / 3600) * (u ** (i + 1))

    return eps

def sun_nutation(julian_day):
    t = (julian_day - J2000) / JULIAN_CENTURY
    t2 = t ** 2
    t3 = t ** 3

    ta = [np.deg2rad(297.850363 + 445267.11148 * t - 0.0019142 * t2 + t3 / 189474.0),
          np.deg2rad(357.52772 + 35999.05034 * t - 0.0001603 * t2 - t3 / 300000.0),
          np.deg2rad(134.96298 + 477198.867398 * t + 0.0086972 * t2 + t3 / 56250.0),
          np.deg2rad(93.27191 + 483202.017538 * t - 0.0036825 * t2 + t3 / 327270),
          np.deg2rad(125.04452 - 1934.136261 * t + 0.0020708 * t2 + t3 / 450000.0)]

    for i in range(5):
        ta[i] = bound_angle_rad(ta[i])

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

def sunpos(julian_day, local_latitude, local_longitude):
    T = (julian_day - J2000) / JULIAN_MILLENNIUM
    T2 = T ** 2
    T3 = T ** 3
    T4 = T ** 4
    T5 = T ** 5

    L0 = 280.4664567 + (360007.6982779 * T) + (0.03032028 * T2) + (T3 / 49931) - (T4 / 15300) - (T5 / 2000000)
    L0 = bound_angle_deg(L0)

    M = 357.52911 + 359990.50340 * T - 0.001603 * T2 - T3 / 30000
    M = bound_angle_deg(M)

    e = 0.016708634 - 0.00042037 * T - 0.000001267 * T2

    C = (1.914602 - 0.04817 * T - 0.000014 * T2) * sin(M) + \
        (0.019993 - 0.000101 * T) * sin(2 * M) + \
        0.000289 * sin(3 * M)

    sunLong = L0 + C
    sunAnomaly = M + C

    sunR = (1.000001018 * (1 - e ** 2)) / (1 + (e * cos(sunAnomaly)))

    omega = 125.04452 - 19341.36261 * T + 0.020708 * T2 + T3 / 45000
    Lambda = sunLong - 0.00569 - 0.00478 * sin(omega)

    nut = sun_nutation(julian_day)
    delta_epsilon = nut[1]
    epsilon0 = oblique_eq(julian_day)
    epsilon = epsilon0 + delta_epsilon
    

    alpha = bound_angle_deg(np.rad2deg(math.atan2(cos(epsilon0) * sin(sunLong), cos(sunLong))))
    delta = np.rad2deg(math.asin(sin(epsilon0) * sin(sunLong)))

    alphaApp = bound_angle_deg(np.rad2deg(math.atan2(cos(epsilon) * sin(Lambda), cos(Lambda))))
    deltaApp = np.rad2deg(math.asin(sin(epsilon) * sin(Lambda)))

    greenwich_hour_angle = bound_angle_deg(siderial_time(julian_day))
    st_correction = decimal_to_dms(nut[0])[2] * cos(epsilon) / 15
    greenwich_hour_angle += (st_correction / 240)
    local_hour_angle = bound_angle_deg(greenwich_hour_angle + local_longitude - alpha)

    altitude = np.rad2deg(math.asin(sin(local_latitude) * sin(delta) + cos(local_latitude) * cos(delta) * cos(local_hour_angle))) 

    azimuth = np.rad2deg(np.arccos((sin(delta) * cos(local_latitude) - cos(delta) * sin(local_latitude) * cos(local_hour_angle)) / cos(altitude)))
    azimuth = bound_angle_deg(azimuth)

    if local_hour_angle >= 0:
        azimuth = 360 - azimuth

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
        altitude,           # 15
        azimuth             # 16
    ]

def equation_of_time(julian_day, local_latitude, local_longitude):
    sun_factors = sunpos(julian_day, local_latitude, local_longitude)
    L0 = sun_factors[0]
    M = sun_factors[1]
    eccentricity = sun_factors[2]
 
    nut = sun_nutation(julian_day)
    
    epsilon = oblique_eq(julian_day) + nut[1]
    y = (math.tan(np.deg2rad(epsilon) / 2)) ** 2

    alpha = sun_factors[10]
    deltaPsi = nut[0]
    E = L0 - 0.0057183 - alpha + deltaPsi * cos(epsilon)
    E *= 4

    return E

def solar_hour_angle(latitude, declination, angle = 0.833):
    dec = np.deg2rad(declination)
    lat = np.deg2rad(latitude)

    num = -1 * np.sin(np.deg2rad(angle)) - np.sin(lat) * np.sin(dec)
    denom = np.cos(lat) * np.cos(dec)

    angle = np.arccos(num / denom)
    return np.rad2deg(angle)

def sunrise_sunset(set_or_rise, hour_angle):
    hours_offset_from_noon = hour_angle / 15

    offset = 12 + set_or_rise * hours_offset_from_noon
    return  offset