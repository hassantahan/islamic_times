"""
Module for calculating Islamic prayer times.

This module provides functions to compute Islamic prayer times based on various methodologies.
It accounts for differences in ʿAṣr calculations, Islamic midnight determination, and prayer 
adjustments for extreme latitudes.

### Features:
- Computes **all major prayer times**: Fajr, Sunrise, Solar Noon (Ẓuhr), ʿAṣr, Sunset, Maghrib, ʿIshāʾ, and Midnight.
- Supports different **Asr calculation methods** (standard & Hanafi).
- Handles **Makkah-based** Isha timing in Ramadan.
- Computes **Islamic midnight** based on the next day's Fajr or Solar Midnight.
- Accounts for **extreme latitude conditions** (work in progress).

### References:
- Jean Meeus, *Astronomical Algorithms*, 2nd Edition, Willmann-Bell, Inc., 1998.
- [PrayTimes.org - Prayer Times Calculation Methods](http://praytimes.org/wiki/Calculation_Methods)
"""

from islamic_times import sun_equations as se
from islamic_times import time_equations as te
from islamic_times import calculation_equations as ce
from datetime import datetime, timedelta
import numpy as np
from typing import Tuple, Dict


def calculate_prayer_times(date: datetime, latitude: float, longitude: float, elevation: float, utc_offset: float, sun_declination: float,
                           prayer_times_parameters: Tuple [float, float, float, int, int], 
                           is_ramadan = False
                           ) -> Dict[str, datetime | str]:
    """
    Computes Islamic prayer times for a given date and location.

    This function calculates the standard prayer times based on astronomical data. It accounts for 
    different methodologies for Fajr, Maghrib, Isha, ʿAṣr, and Islamic midnight.

    Parameters:
    	jde (float): Julian Ephemeris Date (JDE) for the calculation.
    	deltaT (float): Difference between Terrestrial Time (TT) and Universal Time (UT).
    	latitude (float): Observer's latitude in decimal degrees.
    	longitude (float): Observer's longitude in decimal degrees.
        elevation (float): Observer's elevation above sea level in metres.
    	utc_offset (float): Difference between local time and UTC (in hours).
    	prayer_times_parameters (Tuple[float, float, float, int, int]):
            - Fajr angle (degrees).
            - Maghrib angle (degrees).
            - Isha angle (degrees, or `inf for Makkah method).
            - Asr method (0 = Standard, 1 = Hanafi).
            - Midnight method (0 = Sunset to Sunrise, 1 = Sunset to Fajr).
    	is_ramadan (bool, optional): If `True`, applies the special **Makkah ʿishāʾ method** for Ramadan. Defaults to `False`.

    Returns:
        dict: A dictionary of computed prayer times, where each value is either:\n- A `datetime` object (valid prayer time).\n- A `str` message (`"Does not exist."` or `"ʿAṣr time cannot be calculated."`).

    Notes:
    - If the sun does not satisfy the required angle for a prayer time, the function returns `"Does not exist."`.
    - The **Makkah method** sets **Isha time** as **one hour after Maghrib (two hours in Ramadan)**.
    - **Islamic midnight** is determined based on either:
        	midnight_type=0`: Midpoint between sunset and sunrise.
        	midnight_type=1`: Midpoint between sunset and the next day''s Fajr.
    """

    # Extract parameters
    fajr_angle = prayer_times_parameters[0]
    maghrib_angle = prayer_times_parameters[1]
    isha_angle = prayer_times_parameters[2]
    midnight_type = prayer_times_parameters[3]
    asr_type = prayer_times_parameters[4]

    # Standard sun times calculations
    standard_fajr = se.find_proper_suntime(date, latitude, longitude, elevation, utc_offset, 'rise', fajr_angle)
    standard_sunrise = se.find_proper_suntime(date, latitude, longitude, elevation, utc_offset, 'rise')
    standard_noon = se.find_sun_transit(date, latitude, longitude, elevation, utc_offset)
    standard_sunset = se.find_proper_suntime(date, latitude, longitude, elevation, utc_offset, 'set')
    
    # Calculate ʿAṣr time
    try:
        standard_asr = asr_time(standard_noon, latitude, sun_declination, asr_type)
    except ArithmeticError as e:
        standard_asr = str(e)
    
    # Calculate Maghrib time
    # Only if maghrib is not at sunset
    if maghrib_angle > 0:
        standard_maghrib = se.sunrise_or_sunset(date, latitude, longitude, elevation, utc_offset, 'set', maghrib_angle)
    else:
        # Otherwise Maghrib is sunset
        standard_maghrib = standard_sunset

    # Calculate ʿishāʾ time
    # If NOT makkah method
    if isha_angle is not np.inf:
        standard_isha = se.sunrise_or_sunset(date, latitude, longitude, elevation, utc_offset, 'set', isha_angle)
    # Makkah method is special
    else:
        # Ramadan has isha as two hours after maghrib
        if is_ramadan:
            standard_isha = standard_maghrib + timedelta(hours=2)
        # Otherwise it is one hour
        else:
            standard_isha = standard_maghrib + timedelta(hours=1)
    
    # Calculate Midnight time
    try:
        if midnight_type:
            standard_midnight = te.time_midpoint(standard_sunset, find_tomorrow_time(date, latitude, longitude, elevation, utc_offset, fajr_angle))
        else:
            standard_midnight = te.time_midpoint(standard_sunset, find_tomorrow_time(date, latitude, longitude, elevation, utc_offset))
    except ValueError as e:
        standard_midnight = np.inf

    prayer_dict = {
        'fajr': standard_fajr,
        'sunrise': standard_sunrise,
        'noon': standard_noon,
        'asr': standard_asr,
        'sunset': standard_sunset,
        'maghrib': standard_maghrib,
        'isha': standard_isha,
        'midnight': standard_midnight
    }

    # Deal with locations at extreme latitudes
    if any(dt == np.inf for dt in prayer_dict.values()):
        prayer_dict = extreme_latitudes(prayer_dict)

    return prayer_dict

# Used to calculate islamic midnight
def find_tomorrow_time(date: datetime, lat: float, long: float, elev: float, utc_offset: float, angle: float = 5 / 6) -> datetime:
    """
    Calculates the solar time for the next day to determine **Islamic midnight**.

    This function is used to compute the **Fajr time for the next day**, which is needed 
    for the **midnight calculation** in the **Sunset-to-Fajr method**.

    Parameters:
    	date (datetime): Date of the observer.
    	lat (float): Latitude of the observer.
    	long (float): Longitude of the observer.
    	elevation (float): Observer's elevation above sea level in metres.
    	utc_change (float): Observer's UTC offset.

    Returns:
    	datetime: The calculated **fajr time** for the next day in **local standard time**.
    """
    tomorrow_date = date + timedelta(days=1)
    tomorrow_standard_time = se.sunrise_or_sunset(tomorrow_date, lat, long, elev, utc_offset, 'rise', angle)

    return tomorrow_standard_time

# ʿAṣr time
# Type 0: Most schools
# Type 1: Ḥanafī definition
def asr_time(noon: datetime, lat: float, dec: float, ts: int = 1) -> datetime:
    """
    Computes the **ʿaṣr prayer time** based on the observer's **shadow ratio**.

    The ʿaṣr prayer time is determined based on the length of an object's shadow 
    relative to its height. Two methodologies exist:
    - **Standard Method (`t=1`)**: Shadow is equal to the object's height.
    - **Ḥanafī Method (`t=2`)**: Shadow is **twice** the object's height.

    Parameters:
        noon (datetime): Time of ẓuhr/solar noon/sun transit/sun culmination.
        lat (float): Observer's latitude in decimal degrees.
        dec (float): Sun's **declination** at the given time.
        ts (int, optional): ʿaṣr calculation method:
            - 1 (Standard) → Shadow ratio of 1:1.
            - 2 (Hanafi) → Shadow ratio of 2:1.  

    Returns:
        float/np.inf: Number of hours after solar noon when ʿaṣr occurs, or ʿaṣr cannot be calculated due to extreme solar geometry.

    Notes:
    - **If the sun never reaches the required shadow ratio**, the function returns **infinity (`np.inf`)**.
    - This condition occurs in **polar regions** during certain seasons.
    """

    temp_num = ce.sin(np.rad2deg(np.arctan2(1, ts + ce.tan(lat - dec)))) - ce.sin(lat) * ce.sin(dec)
    temp_denom = ce.cos(lat) * ce.cos(dec)
    
    # Sometimes, ʿaṣr cannot be calculated because the sun's geometry at the given date and coordintes does not satisfy the shadow ratio
    # In such a scenario, it will just return that message.
    # TODO: Another fix might be needed.
    if temp_num > temp_denom:
        raise ArithmeticError("ʿAṣr time cannot be calculated because the sun's geometry at the given date and coordintes does not satisfy the shadow ratio.")
    else:
        asr_hours = 1 / 15 * np.rad2deg(np.arccos(temp_num / temp_denom))
    
    return noon + timedelta(hours=asr_hours)

def extreme_latitudes(prayer_times: Dict[str, datetime]):
    """
    Placeholder function for handling **prayer time adjustments in extreme latitudes**.

    ### Purpose:
    - **In progress:** This function is intended to implement methods for adjusting prayer times 
      in locations where the sun **never rises or sets for extended periods** (e.g., near the poles).
    - Methods that may be implemented:
        - Nearest latitude method.
        - Middle of the night method.
        - 1/7th night division method.
        - Angle-based approximation.

    ### Returns:
    - (To be implemented)
    """

    for key, value in prayer_times.items():
        if value == np.inf:
            prayer_times[key] = "The observer at a latitude such that the sun does not reach the given angle for the prayer."

    return prayer_times