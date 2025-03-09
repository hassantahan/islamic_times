"""
Module for various calculations used in the Islamic Times package.

This module provides functions for angle normalization, trigonometric operations in degrees,
and conversions between different angle representations.
"""

import math
import numpy as np
from typing import Tuple
from islamic_times.time_equations import EARTH_RADIUS_KM

def sin(angle: float) -> float:
    '''sin(x) function that auto-converts the argument from degrees to radians before evaluating.
    Parameters:
        angle (float): The angle in degrees.

    Returns:
        float: The tangent of the angle.
    '''
    return math.sin(np.deg2rad(angle))

def cos(angle: float) -> float:
    '''cos(x) function that auto-converts the argument from degrees to radians before evaluating.

    Parameters:
        angle (float): The angle in degrees.

    Returns:
        float: The cosine of the angle.
    '''
    return math.cos(np.deg2rad(angle))

def tan(angle: float) -> float:
    '''tan(x) function that auto-converts the argument from degrees to radians before evaluating.
    
    Parameters:
        angle (float): The angle in degrees.

    Returns:
        float: The tangent of the angle.
    '''
    return math.tan(np.deg2rad(angle))

def decimal_to_dms(decimal_deg: float) -> Tuple[int, int, float]:
    '''Convert a degree value from degrees decimal to degrees-minutes-seconds.
    
    Parameters:
        decimal_deg (float): The angle in degrees.

    Returns:
        tuple (Tuple[int, int, float]): The angle in degrees-minutes-seconds.
    '''
    degrees = int(decimal_deg)
    minutes = int((decimal_deg - degrees) * 60)
    seconds = np.round((decimal_deg - degrees - minutes / 60) * 3600, 2)

    return (degrees, minutes, seconds)

def decimal_to_hms(decimal_degrees: float) -> Tuple[int, int, float]:
    '''Convert a degree value from degrees decimal to hours-minutes-seconds.
    
    Parameters:
        decimal_degrees (float): The angle in degrees.

    Returns:
        tuple (Tuple[int, int, float]): The angle in hours-minutes-seconds.
    '''
    hours = int(decimal_degrees / 15)
    minutes = int((decimal_degrees / 15 - hours) * 60)
    seconds = (decimal_degrees / 15 - hours - minutes / 60) * 3600
    return (hours , round(minutes), seconds)

def hms_to_decimal(hour_angle: Tuple[int, int, float]) -> float:
    '''Convert a degree value from hours-minutes-seconds to degrees decimal.

    Parameters:
        hour_angle (Tuple[int, int, float]): The angle in hours-minutes-seconds.

    Returns:
        float: The angle in degrees.
    '''
    degree = hour_angle[0] + hour_angle[1] / 60 + hour_angle[2] / 3600
    degree *= 15
    return degree

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
    azimuth1, altitude1, azimuth2, altitude2 = map(math.radians, [azimuth1, altitude1, azimuth2, altitude2])

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
    c = np.deg2rad(calculate_angle_diff(lon1, lat1, lon2, lat2))
    distance = EARTH_RADIUS_KM * c  # Earth radius in km

    # Calculate course angle
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
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