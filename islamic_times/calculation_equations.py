"""
Module for various calculations used in the Islamic Times package.

This module provides functions for angle normalization, trigonometric operations in degrees,
and conversions between different angle representations.
"""

import math
import numpy as np
from typing import Tuple
from islamic_times.time_equations import EARTH_RADIUS_KM
from islamic_times.it_dataclasses import Angle, RightAscension, Distance, DistanceUnits

def sin(angle: float) -> float:
    '''sin(x) function that auto-converts the argument from degrees to radians before evaluating.
    Parameters:
        angle (float): The angle in degrees.

    Returns:
        float: The tangent of the angle.
    '''
    return np.sin(np.radians(angle))

def cos(angle: float) -> float:
    '''cos(x) function that auto-converts the argument from degrees to radians before evaluating.

    Parameters:
        angle (float): The angle in degrees.

    Returns:
        float: The cosine of the angle.
    '''
    return np.cos(np.radians(angle))

def calculate_angle_diff(azimuth1: float, altitude1: float, azimuth2: float, altitude2: float) -> float:
    '''Difference between two angles in a radial coordinate system using the haversine formula. 

    Parameters:
        azimuth1 (float): The azimuth angle of the first point (°).
        altitude1 (float): The altitude angle of the first point (°).
        azimuth2 (float): The azimuth angle of the second point (°).
        altitude2 (float): The altitude angle of the second point (°).

    Returns:
        float: The angle difference between the two points (°).
    '''
    # Convert degrees to radians
    azimuth1, altitude1, azimuth2, altitude2 = map(np.radians, [azimuth1, altitude1, azimuth2, altitude2])

    # Apply the haversine formula
    delta_azimuth = azimuth2 - azimuth1
    delta_altitude = altitude2 - altitude1

    a = math.sin(delta_altitude / 2)**2 + math.cos(altitude1) * math.cos(altitude2) * math.sin(delta_azimuth / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

    # Convert the angle difference from radians to degrees
    angle_diff = math.degrees(c)
    
    return angle_diff

# Modified calculate_angle_diff for finding course angle
def haversine(lat1: float, lon1: float, lat2: float, lon2: float) -> Tuple[float, float]:
    '''
    Calculate the great-circle distance and initial bearing between two points using the haversine formula.
    
    Parameters:
        lat1 (float): Latitude of point 1 (°).
        lon1 (float): Longitude of point 1 (°).
        lat2 (float): Latitude of point 2 (°).
        lon2 (float): Longitude of point 2 (°).

    Returns:
        tuple (Tuple[float, float]): The great-circle distance (km) and initial bearing (°).
    '''
    c = math.radians(calculate_angle_diff(lon1, lat1, lon2, lat2))
    distance = EARTH_RADIUS_KM * c  # Earth radius in km

    # Calculate course angle
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    y = math.sin(lon2 - lon1) * math.cos(lat2)
    x = math.cos(lat1) * math.sin(lat2) - math.sin(lat1) * math.cos(lat2) * math.cos(lon2 - lon1)
    course_angle = math.degrees(math.atan2(y, x))

    return distance, course_angle

# Used for finding direction to Mecca 
def get_cardinal_direction(degree: float) -> str:
    '''Determine the cardinal direction corresponding to an angle. Input is in degrees and 0 corresponds to North.

    Parameters:
        degree (float): The angle in degrees.

    Returns:
        str: The cardinal direction corresponding to the
    '''
    cardinals = [
        'N', 'NNE', 'NE', 'ENE',
        'E', 'ESE', 'SE', 'SSE',
        'S', 'SSW', 'SW', 'WSW',
        'W', 'WNW', 'NW', 'NNW'
    ]
    idx = int(((degree + 11.25) % 360) / 22.5)
    return cardinals[idx]

# Returns the signed difference (b - a) in the range (-180, +180].
def angle_diff(a: float, b: float) -> float:
    '''Compute the signed difference between two angles, normalized to (-180°, 180°].

    Parameters:
        a (float): The first angle in degrees.
        b (float): The second angle in degrees.
    
    Returns:
        float: The signed difference between the two angles
    '''
    d = (b - a) % 360.0      # now d is in [0, 360)
    if d > 180.0:
        d -= 360.0           # shift to (-180, +180]
    return d

# As found in Chapter 3 of AA
# Only for angles
def interpolation(n: float, y1: float, y2: float, y3: float) -> float:
    '''Interpolate between angles with the method described in Jean Meeus' *Astronomical Algorithms*, Chapter 3.

    Parameters:
        n (float): The fraction of the interpolation.
        y1 (float): The angle at the first point.
        y2 (float): The angle at the second point.
        y3 (float): The angle at the third point.
    
    Returns:
        float: The interpolated angle.

    Notes:
    - The function has been modified so that it works for angles in degrees and proper wrapping occurs.
    '''
    a = angle_diff(y1, y2)
    b = angle_diff(y2, y3)
    c = b - a

    value = y2 + n / 2 * (a + b + n * c)

    return (value + 360) % 360

# Fixing RA and Dec for apparency pg. 279
def correct_ra_dec(ra: RightAscension, dec: Angle, lha: Angle, parallax: Angle, lat: Angle, elev: Distance, dist: Distance = Distance(EARTH_RADIUS_KM, DistanceUnits.KILOMETRE)) -> Tuple[RightAscension, Angle]:
	'''
	Correct the Moon's Right Ascension and Declination for apparent position. See Chapter 40 of *Astronomical Algorthims* for more information.

	Parameters:
		ra (RightAscension): The Moon's Right Ascension.
		dec (Angle): The Moon's Declination.
		lha (Angle): The Local Hour Angle.
		parallax (Angle): The Moon's parallax.
		lat (Angle): The observer's latitude.
		elev (Distance): The observer's elevation above sea level.
		dist (float): The observer's distance from the celestial body, usually the Moon or the Sun.

	Returns:
		tuple (Tuple[RightAscension, Angle]): The topocentric Right Ascension and Declination.
	'''
	a: float = dist.value
	f: float = 1 / 298.257223563
	b: float = a * (1 - f)

	u: Angle = Angle(math.degrees(math.atan2(b / a * math.tan(lat.radians), 1)))
	p_sin_phi_prime: float = b / a * math.sin(u.radians) + elev.in_unit(DistanceUnits.KILOMETRE) / dist.in_unit(DistanceUnits.KILOMETRE) * math.sin(lat.radians)
	p_cos_phi_prime: float = math.cos(u.radians) + elev.in_unit(DistanceUnits.KILOMETRE) / dist.in_unit(DistanceUnits.KILOMETRE) * math.cos(lat.radians)

	temp_num: float = -1 * p_cos_phi_prime * math.sin(parallax.radians) * math.sin(lha.radians)
	temp_denom: float = math.cos(dec.radians) - p_cos_phi_prime * math.sin(parallax.radians) * math.cos(lha.radians)
	deltaA: Angle = Angle(math.degrees(math.atan2(temp_num, temp_denom)))

	temp_num: float = (math.sin(dec.radians) - p_sin_phi_prime * math.sin(parallax.radians)) * math.cos(deltaA.radians)

	ascension_prime: RightAscension = RightAscension((ra.decimal_degrees.decimal + deltaA.decimal) / 15)
	declination_prime: Angle = Angle(math.degrees(math.atan2(temp_num, temp_denom)))

	return (ascension_prime, declination_prime)

def geocentric_horizontal_coordinates(observer_latitude: Angle, body_declination: Angle, body_lha: Angle) -> Tuple[Angle, Angle]:
    """
    Calculate the geocentric horizontal coordinates (altitude and azimuth) of a celestial body.

    Parameters:
        observer_latitude (Angle): The latitude of the observer.
        body_declination (Angle): The declination of the celestial body.
        body_lha (Angle): The local hour angle of the celestial body.

    Returns:
        tuple (Tuple[Angle, Angle]): The geocentric altitude and azimuth of the celestial body.
    """
    f = 1 / 298.257223563
    tan_phi = (1 - f) ** 2 * math.tan(observer_latitude.radians)
    observer_geocentric_latitude = Angle(math.degrees(math.atan2(tan_phi, 1)))

    # Geocentric Alt 
    body_geo_alt: Angle = Angle(math.degrees(math.asin(
                    math.sin(observer_geocentric_latitude.radians) * math.sin(body_declination.radians) + 
                    math.cos(observer_geocentric_latitude.radians) * math.cos(body_declination.radians) * math.cos(body_lha.radians)
                )))
    
    # Geocentric Az
    cos_az = (math.sin(body_declination.radians) - math.sin(body_geo_alt.radians) * math.sin(observer_geocentric_latitude.radians)) / \
                (math.cos(body_geo_alt.radians) * math.cos(observer_geocentric_latitude.radians))
    sin_az = (-1 * math.cos(body_declination.radians) * math.sin(body_lha.radians)) / (math.cos(body_geo_alt.radians))
    body_geo_az: Angle = Angle(math.degrees(math.atan2(sin_az, cos_az)) % 360)

    return (body_geo_alt, body_geo_az)