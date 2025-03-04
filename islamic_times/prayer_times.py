"""
Special module for specific Islamic prayer times calculations.

This module provides functions to compute specific prayer times such as:
  - ʿAṣr time based on different methodologies.
  - The prayer time for the next day needed to determine Islamic midnight.
"""

from islamic_times import sun_equations as se
from islamic_times import time_equations as te
from islamic_times import calculation_equations as ce
from datetime import datetime, timedelta
import numpy as np
from typing import Tuple, Dict


def calculate_prayer_times(jde: float, deltaT: float, latitude: float, longitude: float, utc_diff: float, solar_parameters: Tuple[float, float], 
                           prayer_times_parameters: Tuple [float, float, float, int, int], is_ramadan = False) -> Dict[str, datetime]:
    '''
    ### Description:
    Calculate the prayer times for a given date and location.

    ### Args:
        - `jde`: Julian Ephemeris Date.
        - `deltaT`: Difference between Terrestrial Time and Universal Time.
        - `latitude`: Latitude of the observer.
        - `longitude`: Longitude of the observer.
        - `utc_diff`: Difference between local time and UTC.
        - `solar_parameters`: Tuple containing the solar declination and the solar angle.
        - `prayer_times_parameters`: Tuple containing the angles for Fajr, Maghrib, and Isha, the ʿAṣr method, and the midnight method
        - `is_ramadan`: Boolean to determine if the calculation is for Ramadan.

    ### Returns:
        - Dictionary containing the prayer times in standard time.
    '''

    # Extract parameters
    fajr_angle = prayer_times_parameters[0]
    maghrib_angle = prayer_times_parameters[1]
    isha_angle = prayer_times_parameters[2]
    midnight_type = prayer_times_parameters[3]
    asr_type = prayer_times_parameters[4]

    sun_declination = solar_parameters[0]
    solar_angle = solar_parameters[1]

    jd = jde - deltaT / 86400

    # Equation of Time
    eq_of_time = se.equation_of_time(jde, deltaT, latitude, longitude)

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
    asr_hours = asr_time(latitude, sun_declination, t=asr_type + 1)
    
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

    # Calculate ʿIshāʾ time
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
    '''
    Calculate the solar time for the next day to determine Islamic midnight.
    '''
    jde_tomorrow = jde + 1

    tomorrow_sun_declination = se.sunpos(jde_tomorrow, deltaT, lat, long)[11]
    tomorrow_time_solar_angle = se.solar_hour_angle(lat, tomorrow_sun_declination, angle)
    tomorrow_solar_time = se.sunrise_sunset(-1, tomorrow_time_solar_angle)
    tomorrow_standard_time = te.solar2standard(jde_tomorrow - deltaT / 86400, tomorrow_solar_time, utc_change, long, eq_of_time_minutes)

    return tomorrow_standard_time

# ʿAṣr time
# Type 0: Most schools
# Type 1: Ḥanafī definition
def asr_time(lat: float, dec: float, t: int = 1) -> float:
    '''
    Calculate the ʿAṣr time based on the shadow ratio.
    '''
    temp_num = ce.sin(np.rad2deg(np.arctan2(1, t + ce.tan(lat - dec)))) - ce.sin(lat) * ce.sin(dec)
    temp_denom = ce.cos(lat) * ce.cos(dec)
    
    # Sometimes, ʿaṣr cannot be calculated because the sun's geometry at the given date and coordintes does not satisfy the shadow ratio
    # In such a scenario, it will just return that message.
    # TODO: Another fix might be needed.
    if temp_num > temp_denom:
        return np.inf
    else:
        return 1 / 15 * np.rad2deg(np.arccos(temp_num / temp_denom))
    
def extreme_latitutdes():
    
    return