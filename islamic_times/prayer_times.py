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
- [PrayTimes.org - Prayer Times Calculation Methods](https://praytimes.org/docs/calculation)
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
# Source: https://praytimes.org/docs/calculation
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
Reference: https://praytimes.org/docs/calculation
"""

def safe_sun_time(observer_date: DateTimeInfo, observer: ObserverInfo, event: str, angle: Angle) -> datetime:
    try:
        return se.find_proper_suntime(observer_date, observer, event, angle)
    except ArithmeticError:
        raise ArithmeticError
    
# Used to calculate islamic midnight
def find_tomorrow_time(observer_date: DateTimeInfo, observer: ObserverInfo, angle: Angle = Angle(5 / 6), rise_or_set = 'rise') -> datetime:
    """
    Calculates the "rise" time of the sun at a given angle the day after the observer date. This is used to determine Islamic midnight, either from sunset to Fajr, or sunset to sunrise.

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
    tomorrow_standard_time = safe_sun_time(tomorrow_date, observer, rise_or_set, angle)

    return tomorrow_standard_time

# ʿAṣr time
# Type 0: Most schools
# Type 1: Ḥanafī definition
def asr_time(noon: datetime, lat: Angle, dec: Angle, ts: int = 0) -> datetime:
    """
    Computes the ʿaṣr prayer time based on the observer's shadow ratio.

    The ʿaṣr prayer time is determined based on the length of an object's shadow 
    relative to its height. Two methodologies exist:
    - Standard Method (`t=1`): Shadow is equal to the object's height + 
    the height of it's shadow at noon.
    - Ḥanafī Method (`t=2`): Shadow is **twice** the object's height + 
    the height of it's shadow at noon.

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
    # Convert shadow factor from 0/1 to 1/2 for actual calculation
    shadow_factor = ts + 1
    
    # Calculate noon shadow length
    noon_shadow = abs(math.tan(math.radians(lat.decimal - dec.decimal)))
    
    # Total shadow length = object height + (shadow_factor * noon shadow)
    total_shadow = 1 + (shadow_factor * noon_shadow)
    
    # Calculate angle using inverse tangent
    angle_rad = math.atan(1.0 / total_shadow)
    
    # Check if ʿAṣr is possible at this location and date
    temp_num = (math.sin(angle_rad) - math.sin(lat.radians) * math.sin(dec.radians))
    temp_denom = (math.cos(lat.radians) * math.cos(dec.radians))
    
    if abs(temp_num) > abs(temp_denom):
        raise ArithmeticError("ʿAṣr prayer time cannot be calculated - sun never reaches required shadow length")
    
    # Convert to hours from noon
    asr_hours = 1/15.0 * math.degrees(math.acos(temp_num / temp_denom))
    
    return noon + timedelta(hours=asr_hours)

def extreme_latitudes(observer_date: DateTimeInfo, observer: ObserverInfo, prayer_list: List[Prayer]):
    """
    Placeholder function for handling **prayer time adjustments in extreme latitudes**.

    ### Purpose:
    - **In progress:** This function is intended to implement methods for adjusting prayer times 
      in locations where the sun **never rises or sets for extended periods** (e.g., near the poles).
    - Methods that may be implemented:
        - None; just accept the facts.
        - Nearest latitude method.
        - Middle of the night method.
        - 1/7th night division method.
        - Angle-based approximation.

    ### Returns:
    - (To be implemented)
    """

    method: PrayerMethod = prayer_list[0].method
    observer_new_lat: ObserverInfo = observer
    
    # Nearest latitude method
    if method.extreme_lats == 'NearestLat':

        lowest_angle: float = max(method.fajr_angle.decimal, method.isha_angle.decimal, 18)
        eps = se.oblique_eq(observer_date.jde)
        nearest_lat: float = 90 - eps.decimal - lowest_angle - 0.01
        observer_new_lat: ObserverInfo = replace(observer_new_lat, latitude=Angle(nearest_lat))

        for i, prayer in enumerate(prayer_list):
            # Night prayer times
            if prayer.name in ['Maghrib', 'ʿIshāʾ']:
                # Find angle
                ang = Angle(0)
                if prayer.name == 'Maghrib':
                    ang = method.maghrib_angle
                elif prayer.name == 'ʿIshāʾ':
                    ang = method.isha_angle

                prayer_list[i] = replace(prayer, time=safe_sun_time(observer_date, observer_new_lat, 'set', ang))

            # Midnight is exceptional
            elif prayer.name == 'Midnight':
                # Either next fajr or sunrise
                if method.midnight_type:
                    next_prayer = find_tomorrow_time(observer_date, observer_new_lat, method.fajr_angle)
                else:
                    next_prayer = find_tomorrow_time(observer_date, observer_new_lat)

                prayer_list[i] = replace(prayer, time=te.time_midpoint(prayer_list[4].time, next_prayer)) # type: ignore

            # Fajr and Sunset
            elif prayer.name == 'Fajr':
                prayer_list[i] = replace(prayer, time=safe_sun_time(observer_date, observer_new_lat, 'rise', method.fajr_angle))
                
    elif method.extreme_lats in ['MiddleNight', 'OneSeventh', 'AngleBased']:
        EXTREME_LATITUDES_FACTORS = {
            'MiddleNight' : 1 / 2,
            'OneSeventh' : 1 / 7,
            'AngleBased' : 1 / 60
        }

        midnight_float = (timedelta(hours=24) + prayer_list[1].time - prayer_list[4].time).total_seconds() / 3600 # type: ignore
        portion = EXTREME_LATITUDES_FACTORS[method.extreme_lats] * midnight_float

        for i, prayer in enumerate(prayer_list):
            if prayer.name in ['Fajr', 'Maghrib', 'ʿIshāʾ']:
                if prayer.name == 'Fajr':
                    sign = -1
                    base = prayer_list[1] # sunrise
                else:
                    sign = 1
                    base = prayer_list[4] # sunset

                if method.extreme_lats == 'AngleBased':
                    if prayer.name == 'Fajr': ang = method.fajr_angle
                    elif prayer.name == 'Maghrib': ang = method.maghrib_angle
                    else: ang = method.isha_angle
                    portion_mod = portion * ang.decimal
                else:
                    portion_mod = portion

                try:
                    delta: timedelta = prayer.time - base.time # type: ignore
                except:
                    delta = math.inf # type: ignore

                if delta == math.inf or delta.total_seconds() / 3600 > portion_mod:
                    prayer_list[i] = replace(prayer, time=base.time + timedelta(hours=sign*portion_mod)) # type: ignore
        
        # Midnight needs readjustment
        next_sunrise: datetime = find_tomorrow_time(observer_date, observer)
        next_sunset: datetime = find_tomorrow_time(observer_date, observer, rise_or_set='set')
        if method.midnight_type:
            # Find next fajr thru the same method
            try:
                next_fajr = find_tomorrow_time(observer_date, observer, method.fajr_angle)
            except:
                next_fajr = math.inf

            midnight_float = (timedelta(hours=24) + next_sunrise - next_sunset).total_seconds() / 3600
            portion = EXTREME_LATITUDES_FACTORS[method.extreme_lats] * midnight_float
            if method.extreme_lats == 'AngleBased': portion *= method.fajr_angle.decimal

            try:
                delta: timedelta = next_fajr - next_sunrise # type: ignore
            except:
                delta = math.inf # type: ignore

            if delta == math.inf or delta.total_seconds() / 3600 > portion:
                next_fajr = next_sunrise - timedelta(hours=portion)
            
            prayer_list[-1] = replace(prayer, time=te.time_midpoint(prayer_list[4].time, next_fajr)) # type: ignore
        else:
            prayer_list[-1] = replace(prayer, time=te.time_midpoint(prayer_list[4].time, next_sunrise)) # type: ignore

    else:
        # warnings.warn(f"Extreme latitude warning. Prayer times at this latitude are not well established.")
        for i, prayer in enumerate(prayer_list):
            prayer_list[i] = replace(prayer, time="The observer is at a latitude such that the sun does not reach the given angle for the prayer.")

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
    try:
        fajr_dt = safe_sun_time(observer_date, observer, 'rise', method.fajr_angle)
    except ArithmeticError:
        fajr_dt = math.inf
    
    # Calculate ʿAṣr time
    try:
        asr_dt = asr_time(sun_info.sun_transit, observer.latitude, sun_info.apparent_declination, method.asr_type)
    except ArithmeticError as e:
        asr_dt = str(e)
    
    # Calculate Maghrib time
    # Only if maghrib is not at sunset
    if method.maghrib_angle.decimal > 0:
        try:
            maghrib_dt = safe_sun_time(observer_date, observer, 'set', method.maghrib_angle)
        except:
            maghrib_dt = math.inf
    else:
        # Otherwise Maghrib is sunset
        maghrib_dt = sun_info.sunset + timedelta(minutes=1)

    # Calculate ʿishāʾ time
    # If NOT makkah method
    if "Makkah" not in method.name:
        try:
            isha_dt = safe_sun_time(observer_date, observer,  'set', method.isha_angle)
            if isha_dt.time() < sun_info.sun_transit.time():
                isha_dt = find_tomorrow_time(observer_date, observer, method.isha_angle, 'set')
        except ArithmeticError:
            isha_dt = math.inf
    # Makkah method is special and will be the only exception
    else:
        # During Ramadan, ʿishāʾ is set to a flat two hours after maghrib
        if observer_date.hijri is not None and observer_date.hijri.hijri_month == 9:
            isha_dt = maghrib_dt + timedelta(hours=2) # type: ignore

        # Otherwise, it is set to a flat 1.5 hours
        else:
            isha_dt = maghrib_dt + timedelta(hours=1.5) # type: ignore
    
    # Calculate Midnight time
    try:
        if method.midnight_type:
            midnight_dt = te.time_midpoint(sun_info.sunset, find_tomorrow_time(observer_date, observer, method.fajr_angle))
        else:
            midnight_dt = te.time_midpoint(sun_info.sunset, find_tomorrow_time(observer_date, observer))
    except ArithmeticError:
        midnight_dt = math.inf
    

    prayer_list: List[Prayer] = [
        Prayer("Fajr", fajr_dt, method),
        Prayer("Sunrise", sun_info.sunrise, method),
        Prayer("Ẓuhr", sun_info.sun_transit, method),
        Prayer("ʿAṣr", asr_dt, method),
        Prayer("Sunset", sun_info.sunset, method),
        Prayer("Maghrib", maghrib_dt, method),
        Prayer("ʿIshāʾ", isha_dt, method),
        Prayer("Midnight", midnight_dt, method)
    ]

    # Deal with locations at extreme latitudes
    if any(p.time == math.inf for p in prayer_list):
        warnings.warn(f"Extreme latitude warning. Prayer times at this latitude are not well established.")
        prayer_list = extreme_latitudes(observer_date, observer, prayer_list)

    return PrayerTimes(
        method=method,
        fajr=prayer_list[0],
        sunrise=prayer_list[1],
        zuhr=prayer_list[2],
        asr=prayer_list[3],
        sunset=prayer_list[4],
        maghrib=prayer_list[5],
        isha=prayer_list[6],
        midnight=prayer_list[7],
    )