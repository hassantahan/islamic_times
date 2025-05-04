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

import warnings
from typing import List
import islamic_times.astro_core as fast_astro
from islamic_times.it_dataclasses import *
from islamic_times import sun_equations as se
from islamic_times import time_equations as te
from datetime import datetime, timedelta
from dataclasses import replace

##### Prayer Methods #####
# Source: http://praytimes.org/wiki/Calculation_Methods
DEFAULT_PRAYER_METHODS: List[PrayerMethod] = [
    PrayerMethod(
        name="Muslim World League (MWL)",
        keys=('MWL', 'MUSLIM WORLD LEAGUE'), 
        fajr_angle=Angle(18), isha_angle=Angle(17)
    ),
    PrayerMethod(
        name="Islamic Society of North America (ISNA)",
        keys=('ISNA', 'ISLAMIC SOCIETY OF NORTH AMERICA'),
        fajr_angle=Angle(15), isha_angle=Angle(15)
    ),
    PrayerMethod(
        name="Egyptian General Authority of Survey (Egypt)",
        keys=('EGYPT', 'EGYPTIAN', 'EGYPTIAN GENERAL AUTHORITY OF SURVEY', 'EGAS'), 
        fajr_angle=Angle(19.5), isha_angle=Angle(17.5)
    ),
    PrayerMethod(
        name="Umm al-Qura University (Makkah)", 
        keys=('MAKKAH', 'MECCA', 'MEKKAH', 'UMM AL-QURA UNIVERSITY', 'UMM AL-QURA', 'UQU'), 
        fajr_angle=Angle(18.5), isha_angle=Angle(999)
    ),
    PrayerMethod(
        name="University of Islamic Sciences, (Karachi)", 
        keys=('KARACHI', 'UNIVERSITY OF ISLAMIC SCIENCES', 'UIS'), 
        fajr_angle=Angle(18), isha_angle=Angle(18)
    ),
    PrayerMethod(
        name="Institute of Geophysics, University of Tehran (Tehran)", 
        keys=('TEHRAN', 'UNIVERSITY OF TEHRAN', 'INSTITUTE OF GEOPHYSICS', 'UOT', 'IOG', 'IOGUOT', 'UOTIOG', 'IGUT', 'UTIG'), 
        fajr_angle=Angle(17.7), isha_angle=Angle(14), maghrib_angle=Angle(4.5), midnight_type=1
    ),
    PrayerMethod(name="Shia Ithna Ashari, Leva Research Institute, Qom (Jafari)",
        keys=('JAFARI', 'JAAFARI', 'SHIA', 'SHIA ITHNA ASHARI', 'LEVA', 'LEVA RESEARCH INSTITUTE', 'QOM', 'QUM', 'LRI', 'SIA', 'SIALRI'), 
        fajr_angle=Angle(16), isha_angle=Angle(14), maghrib_angle=Angle(4), midnight_type=1
    )
]
"""
Reference: http://praytimes.org/wiki/Calculation_Methods
"""

def safe_sun_time(observer_date: DateTimeInfo, observer: ObserverInfo, event: str, angle: float) -> datetime:
    try:
        return se.find_proper_suntime(observer_date, observer, event, angle)
    except ArithmeticError:
        return math.inf
    
# Used to calculate islamic midnight
def find_tomorrow_time(observer_date: DateTimeInfo, observer: ObserverInfo, angle: Angle = Angle(5 / 6)) -> datetime:
    """
    Calculates the solar time for the next day to determine **Islamic midnight**.

    This function is used to compute the **Fajr time for the next day**, which is needed 
    for the **midnight calculation** in the **Sunset-to-Fajr method**.

    Parameters:
    	observer_date (DateTimeInfo): Date information of the observer.
    	observer (ObserverInfo): Position information of the observer.
    	angle (Angle): Fajr angle.

    Returns:
    	datetime: The calculated fajr time for the next day in local time.
    """
    new_date: datetime = observer_date.date + timedelta(days=1) 
    tomorrow_date: DateTimeInfo = replace(observer_date, 
                                            date=new_date,
                                            jd=observer_date.jd + 1,
                                            deltaT=fast_astro.delta_t_approx(new_date.year, new_date.month))
    tomorrow_standard_time = safe_sun_time(tomorrow_date, observer, 'rise', angle)

    return tomorrow_standard_time

# ʿAṣr time
# Type 0: Most schools
# Type 1: Ḥanafī definition
def asr_time(noon: datetime, lat: Angle, dec: Angle, ts: int = 0) -> datetime:
    """
    Computes the ʿaṣr prayer time based on the observer's shadow ratio.

    The ʿaṣr prayer time is determined based on the length of an object's shadow 
    relative to its height. Two methodologies exist:
    - Standard Method (`t=1`): Shadow is equal to the object's height.
    - Ḥanafī Method (`t=2`): Shadow is **twice** the object's height.

    Parameters:
        noon (datetime): Time of ẓuhr/solar noon/sun transit/sun culmination.
        lat (Angle): Observer's latitude.
        dec (Angle): Sun's declination.
        ts (int, optional): ʿaṣr calculation method:
            - 0 (Standard) → Shadow ratio of 1:1.
            - 1 (Ḥanafī) → Shadow ratio of 2:1.  

    Returns:
        float/math.inf: Number of hours after solar noon when ʿaṣr occurs, or ʿaṣr cannot be calculated due to extreme solar geometry.

    Notes:
        - If the sun never reaches the required shadow ratio, the function returns `math.inf`.
    """

    temp_num = math.sin(math.atan2(1, 1 + ts + math.tan(lat.radians - dec.radians))) - math.sin(lat.radians) * math.sin(dec.radians)
    temp_denom = math.cos(lat.radians) * math.cos(dec.radians)
    
    # Sometimes, ʿaṣr cannot be calculated because the sun's geometry at the given date and coordintes does not satisfy the shadow ratio
    # In such a scenario, it will just return that message.
    # TODO: Another fix might be needed.
    if temp_num > temp_denom:
        raise ArithmeticError("ʿAṣr time cannot be calculated because the sun's geometry at the given date and coordintes does not satisfy the shadow ratio.")
    else:
        asr_hours = 1 / 15 * math.degrees(math.acos(temp_num / temp_denom))
    
    return noon + timedelta(hours=asr_hours)

def extreme_latitudes(prayer_list: List[datetime]):
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

    for prayer in prayer_list:
        if prayer == math.inf:
            warnings.warn(f"Extreme latitude warning. Prayer times at this latitude are not well established.")
            prayer = "The observer is at a latitude such that the sun does not reach the given angle for the prayer."

    return prayer_list

def calculate_prayer_times(observer_date: DateTimeInfo, observer: ObserverInfo, sun_info: SunInfo, method: PrayerMethod) -> PrayerTimes:
    """
    Computes Islamic prayer times for a given date and location.

    This function calculates the standard prayer times based on astronomical data. It accounts for 
    different methodologies for Fajr, Maghrib, Isha, ʿAṣr, and Islamic midnight.

    Parameters:
    	observer_date (DateTimeInfo): The date and time of the observer.
        observer (ObserverInfo): The observer's coordinates and elevation.
        sun_info (SunInfo): The sun's position data (sunrise, sunset, solar noon).
        method (PrayerMethod): The prayer calculation method to use.

    Returns:
        PrayerTimes: An object containing the calculated prayer times.

    Notes:
    - If the sun does not satisfy the required angle for a prayer time, the function returns `"Does not exist."`. Later updates will allow for different approaches.
    """

    # Calculate fajr
    fajr_dt = safe_sun_time(observer_date, observer, 'rise', method.fajr_angle)
    
    # Calculate ʿAṣr time
    try:
        asr_dt = asr_time(sun_info.sun_transit, observer.latitude, sun_info.apparent_declination, method.asr_type)
    except ArithmeticError as e:
        asr_dt = str(e)
    
    # Calculate Maghrib time
    # Only if maghrib is not at sunset
    if method.maghrib_angle.decimal > 0:
        maghrib_dt = safe_sun_time(observer_date, observer, 'set', method.maghrib_angle)
    else:
        # Otherwise Maghrib is sunset
        maghrib_dt = sun_info.sunset

    # Calculate ʿishāʾ time
    # If NOT makkah method
    if "Makkah" not in method.name:
        isha_dt = safe_sun_time(observer_date, observer,  'set', method.isha_angle)
    # Makkah method is special and will be the only exception
    else:
        # During Ramadan, ʿishāʾ is set to a flat two hours after maghrib
        if observer_date.hijri.hijri_month == 9:
            isha_dt = maghrib_dt + timedelta(hours=2)
        # Otherwise, it is set to a flat one hour
        else:
            isha_dt = maghrib_dt + timedelta(hours=1)
    
    # Calculate Midnight time
    try:
        if method.midnight_type:
            midnight_dt = te.time_midpoint(sun_info.sunset, find_tomorrow_time(observer_date, observer, method.fajr_angle))
        else:
            midnight_dt = te.time_midpoint(sun_info.sunset, find_tomorrow_time(observer_date, observer))
    except:
        midnight_dt = math.inf
    

    prayer_list = [
        fajr_dt,
        sun_info.sunrise,
        sun_info.sun_transit,
        asr_dt,
        sun_info.sunset,
        maghrib_dt,
        isha_dt,
        midnight_dt
    ]

    # Deal with locations at extreme latitudes
    if any(p == math.inf for p in prayer_list):
        prayer_list = extreme_latitudes(prayer_list)

    return PrayerTimes(
        method=method,
        fajr=Prayer("Fajr", fajr_dt, method),
        sunrise=Prayer("Sunrise", sun_info.sunrise, method),
        zuhr=Prayer("Ẓuhr", sun_info.sun_transit, method),
        asr=Prayer("ʿAṣr", asr_dt, method),
        sunset=Prayer("Sunset", sun_info.sunset, method),
        maghrib=Prayer("Maghrib", maghrib_dt, method),
        isha=Prayer("ʿIshāʾ", isha_dt, method),
        midnight=Prayer("Midnight", midnight_dt, method),
    )