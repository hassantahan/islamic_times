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


def calculate_prayer_times(date: datetime, latitude: float, longitude: float, elevation: float, utc_diff: float, sun_declination: float,
                           prayer_times_parameters: Tuple [float, float, float, int, int], 
                           is_ramadan = False
                           ) -> Dict[str, datetime]:
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
    	utc_diff (float): Difference between local time and UTC (in hours).
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

    standard_fajr = find_proper_suntime(date, latitude, longitude, elevation, utc_diff, 'rise', fajr_angle)
    standard_sunrise = find_proper_suntime(date, latitude, longitude, elevation, utc_diff, 'rise')
    _, standard_noon, _ = se.rise_transit_set(date, latitude, longitude, 0, utc_diff)
    standard_sunset = find_proper_suntime(date, latitude, longitude, elevation, utc_diff, 'set')
    
    # Calculate ʿAṣr time
    asr_hours = asr_time(latitude, sun_declination, ts=asr_type + 1)
    
    if asr_hours != np.inf:
        standard_asr = standard_noon + timedelta(hours=asr_hours)
    else:
        standard_asr = "ʿAṣr time cannot be calculated because the sun's geometry at the given date and coordintes does not satisfy the shadow ratio."
    
    # Calculate Maghrib time
    # Only if maghrib is not at sunset
    if maghrib_angle > 0:
        standard_maghrib = se.sunrise_or_sunset(date, latitude, longitude, elevation, utc_diff, 'set', maghrib_angle)
    else:
        # Otherwise Maghrib is sunset
        standard_maghrib = standard_sunset

    # Calculate ʿishāʾ time
    # If NOT makkah method
    if isha_angle is not np.inf:
        standard_isha = se.sunrise_or_sunset(date, latitude, longitude, elevation, utc_diff, 'set', isha_angle)
    # Makkah method is special
    else:
        # Ramadan has isha as two hours after maghrib
        if is_ramadan:
            standard_isha = standard_maghrib + timedelta(hours=2)
        # Otherwise it is one hour
        else:
            standard_isha = standard_maghrib + timedelta(hours=1)
    
    # Calculate Midnight time
    if midnight_type:
        standard_midnight = te.time_midpoint(standard_sunset, find_tomorrow_time(date, latitude, longitude, elevation, utc_diff, fajr_angle))
    else:
        standard_midnight = te.time_midpoint(standard_sunset, find_tomorrow_time(date, latitude, longitude, elevation, utc_diff))

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
def find_tomorrow_time(date: datetime, lat: float, long: float, elev: float, utc_diff: float, angle: float = 5 / 6) -> datetime:
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
    tomorrow_standard_time = se.sunrise_or_sunset(tomorrow_date, lat, long, elev, utc_diff, 'rise', angle)

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
        temp_sunset = se.sunrise_or_sunset(true_date, latitude, longitude, elevation, utc_offset, rise_or_set, angle)
        date_doy = true_date.timetuple().tm_yday
        if temp_sunset == np.inf:
            return datetime.min

        i = 1
        while(True):
            temp_sunset_doy = (temp_sunset + timedelta(hours=temp_utc_offset, minutes=-20)).timetuple().tm_yday
            if (temp_sunset_doy < date_doy and temp_sunset.year == true_date.year) or ((temp_sunset + timedelta(hours=temp_utc_offset)).year < true_date.year):
                temp_sunset = se.sunrise_or_sunset(true_date + timedelta(days=i), latitude, longitude, elevation, utc_offset, rise_or_set, angle)
                i += 1
            else: 
                return temp_sunset

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