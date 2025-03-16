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


def calculate_prayer_times(jde: float, deltaT: float, latitude: float, longitude: float, utc_diff: float, solar_parameters: Tuple[float, float], 
                           prayer_times_parameters: Tuple [float, float, float, int, int], is_ramadan = False) -> Dict[str, datetime]:
    """
    Computes Islamic prayer times for a given date and location.

    This function calculates the standard prayer times based on astronomical data. It accounts for 
    different methodologies for Fajr, Maghrib, Isha, ʿAṣr, and Islamic midnight.

    Parameters:
    	jde (float): Julian Ephemeris Date (JDE) for the calculation.
    	deltaT (float): Difference between Terrestrial Time (TT) and Universal Time (UT).
    	latitude (float): Observer's latitude in decimal degrees.
    	longitude (float): Observer's longitude in decimal degrees.
    	utc_diff (float): Difference between local time and UTC (in hours).
    	solar_parameters (Tuple[float, float, float, float, float, float]):
            - Sun declination (degrees).
            - Solar angle for sunrise/sunset.
            - Delta Psi (nutation in longitude).
            - Mean longitude of the Sun.
            - Obliquity of the ecliptic.
            - Right ascension of the Sun.
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

    sun_declination = solar_parameters[0]
    solar_angle = solar_parameters[1]
    deltaPsi = solar_parameters[2]
    L0 = solar_parameters[3]
    epsilon = solar_parameters[4]
    alpha = solar_parameters[5]

    jd = jde - deltaT / 86400

    # Equation of Time
    eq_of_time = se.equation_of_time(deltaPsi, L0, epsilon, alpha)

    # Calculate prayer times in solar time
    solar_fajr = se.sunrise_sunset(-1, se.solar_hour_angle(latitude, sun_declination, fajr_angle))
    solar_sunrise = se.sunrise_sunset(-1, solar_angle)
    solar_sunset = se.sunrise_sunset(1, solar_angle)

    # Convert prayer times from solar to standard time
    standard_fajr = te.solar2standard(jd, solar_fajr, utc_diff, longitude, eq_of_time)
    standard_sunrise = te.solar2standard(jd, solar_sunrise, utc_diff, longitude, eq_of_time)
    standard_noon = te.solar2standard(jd, 12.0, utc_diff, longitude, eq_of_time)
    standard_sunset = te.solar2standard(jd, solar_sunset, utc_diff, longitude, eq_of_time)
    
    # Calculate ʿAṣr time
    asr_hours = asr_time(latitude, sun_declination, ts=asr_type + 1)
    
    if asr_hours != np.inf:
        standard_asr = standard_noon + timedelta(hours=asr_hours)
    else:
        standard_asr = "ʿAṣr time cannot be calculated because the sun's geometry at the given date and coordintes does not satisfy the shadow ratio."
    
    # Calculate Maghrib time
    # Only if maghrib is not at sunset
    if maghrib_angle > 0:
        solar_maghrib = se.sunrise_sunset(1, se.solar_hour_angle(latitude, sun_declination, maghrib_angle))
        standard_maghrib = te.solar2standard(jd, solar_maghrib, utc_diff, longitude, eq_of_time)
    else:
        # Otherwise Maghrib is sunset
        standard_maghrib = standard_sunset

    # Calculate ʿishāʾ time
    # If NOT makkah method
    if isha_angle is not np.inf:
        solar_isha = se.sunrise_sunset(1, se.solar_hour_angle(latitude, sun_declination, isha_angle))
        standard_isha = te.solar2standard(jd, solar_isha, utc_diff, longitude, eq_of_time)
    # Makkah method is special
    else:
        # Ramadan has isha as two hours after maghrib
        if is_ramadan:
            standard_isha = standard_maghrib + timedelta(hours=2)
        # Otherwise it is one hour
        else:
            standard_isha = standard_maghrib + timedelta(hours=2)
    
    # Calculate Midnight time
    if midnight_type:
        standard_midnight = te.time_midpoint(standard_sunset, find_tomorrow_time(jde, deltaT, utc_diff, latitude, longitude, eq_of_time, fajr_angle))
    else:
        standard_midnight = te.time_midpoint(standard_sunset, find_tomorrow_time(jde, deltaT, utc_diff, latitude, longitude, eq_of_time, 0))

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

    for key, value in prayer_dict.items():
        if value == datetime.min:
            prayer_dict[key] = "Does not exist."

    return prayer_dict

# Used to calculate islamic midnight
def find_tomorrow_time(jde: float, deltaT: float, utc_change: float, lat: float, long: float, eq_of_time_minutes: float, angle: float) -> datetime:
    """
    Calculates the solar time for the next day to determine **Islamic midnight**.

    This function is used to compute the **Fajr time for the next day**, which is needed 
    for the **midnight calculation** in the **Sunset-to-Fajr method**.

    Parameters:
    	jde (float): Julian Ephemeris Date (JDE) for the current day.
    	deltaT (float): Difference between Terrestrial Time and Universal Time.
    	utc_change (float): Observer's UTC offset.
    	lat (float): Latitude of the observer.
    	long (float): Longitude of the observer.
    	eq_of_time_minutes (float): Equation of Time (minutes) for correcting solar time.
    	angle (float): Solar hour angle for computing Fajr time.

    Returns:
    	datetime: The calculated **fajr time** for the next day in **local standard time**.
    """
    jde_tomorrow = jde + 1

    sun_params = se.sunpos(jde_tomorrow, deltaT, lat, long)
    tomorrow_sun_declination = sun_params.apparent_declination
    tomorrow_time_solar_angle = se.solar_hour_angle(lat, tomorrow_sun_declination, angle)
    tomorrow_solar_time = se.sunrise_sunset(-1, tomorrow_time_solar_angle)
    tomorrow_standard_time = te.solar2standard(jde_tomorrow - deltaT / 86400, tomorrow_solar_time, utc_change, long, eq_of_time_minutes)

    return tomorrow_standard_time

# ʿAṣr time
# Type 0: Most schools
# Type 1: Ḥanafī definition
def asr_time(lat: float, dec: float, ts: int = 1) -> float:
    """
    Computes the **ʿaṣr prayer time** based on the observer's **shadow ratio**.

    The ʿaṣr prayer time is determined based on the length of an object's shadow 
    relative to its height. Two methodologies exist:
    - **Standard Method (`t=1`)**: Shadow is equal to the object's height.
    - **Ḥanafī Method (`t=2`)**: Shadow is **twice** the object's height.

    Parameters:
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
        return np.inf
    else:
        return 1 / 15 * np.rad2deg(np.arccos(temp_num / temp_denom))
    
def extreme_latitutdes():
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
    ...