"""Prayer-time calculation algorithms and extreme-latitude adjustments.

This module computes the canonical prayer-time set used by :class:`ITLocation`:
Fajr, Sunrise, Ẓuhr, ʿAṣr, Sunset, Maghrib, ʿIshāʾ, and Islamic midnight.
It includes method profiles from PrayTimes, handling for Makkah-specific ʿIshāʾ
rules, and fallback strategies for high-latitude locations where standard
sun-angle events may be unavailable.
"""

import warnings
import math
from typing import List
import islamic_times.astro_core as fast_astro
from islamic_times.it_dataclasses import Angle, DateTimeInfo, ObserverInfo, Prayer, PrayerMethod, PrayerTimes, SunInfo
from islamic_times import sun_equations as se
from islamic_times import time_equations as te
from datetime import datetime, timedelta
from dataclasses import replace

# Prayer method profiles based on:
# https://praytimes.org/docs/calculation
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
def safe_sun_time(observer_date: DateTimeInfo, observer: ObserverInfo, event: str, angle: Angle) -> datetime:
    """Return a solar event datetime or re-raise an astronomical failure.

    Parameters
    ----------
    observer_date : DateTimeInfo
        Date and time context for the observer.
    observer : ObserverInfo
        Observer geodetic and atmospheric configuration.
    event : str
        Solar event selector accepted by ``sun_equations.find_proper_suntime``.
    angle : Angle
        Solar altitude angle used for event solving.

    Returns
    -------
    datetime
        The resolved event time in observer-local timezone.

    Raises
    ------
    ArithmeticError
        Raised when the requested event does not exist for the given geometry.
    """
    try:
        return se.find_proper_suntime(observer_date, observer, event, angle)
    except ArithmeticError:
        raise ArithmeticError
    
# Used by midnight calculations to sample the next/previous day.
def find_tomorrow_time(observer_date: DateTimeInfo, observer: ObserverInfo, angle: Angle = Angle(5 / 6), rise_or_set = 'rise', num_days: int = 1) -> datetime:
    """Resolve a sun event for a date offset relative to ``observer_date``.

    Parameters
    ----------
    observer_date : DateTimeInfo
        Date/time state for the current observer context.
    observer : ObserverInfo
        Observer geodetic and atmospheric configuration.
    angle : Angle, default=Angle(5 / 6)
        Target solar altitude angle.
    rise_or_set : str, default='rise'
        Event kind forwarded to ``safe_sun_time``.
    num_days : int, default=1
        Day offset from ``observer_date.date``. Negative values are supported.

    Returns
    -------
    datetime
        Event datetime in observer-local timezone.
    """
    new_date: datetime = observer_date.date + timedelta(days=num_days) 
    tomorrow_date: DateTimeInfo = replace(observer_date, 
                                            date=new_date,
                                            jd=observer_date.jd + num_days,
                                            deltaT=fast_astro.delta_t_approx(new_date.year, new_date.month))
    tomorrow_standard_time = safe_sun_time(tomorrow_date, observer, rise_or_set, angle)

    return tomorrow_standard_time

# ʿAṣr time
# Type 0: majority schools
# Type 1: Ḥanafī definition
def asr_time(noon: datetime, lat: Angle, dec: Angle, ts: int = 0) -> datetime:
    """Compute ʿAṣr time from solar-noon geometry and shadow-ratio rule.

    Parameters
    ----------
    noon : datetime
        Local solar-noon datetime.
    lat : Angle
        Observer latitude.
    dec : Angle
        Apparent solar declination.
    ts : int, default=0
        ʿAṣr method selector: ``0`` for standard (shadow factor 1),
        ``1`` for Ḥanafī (shadow factor 2).

    Returns
    -------
    datetime
        Computed ʿAṣr datetime.

    Raises
    ------
    ArithmeticError
        Raised if the required shadow condition is not reachable.
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

def extreme_latitudes(observer_date: DateTimeInfo, observer: ObserverInfo, prayer_list: List[Prayer], sun_dec: Angle) -> List[Prayer]:
    """Adjust prayer outputs for locations with extreme-latitude constraints.

    Parameters
    ----------
    observer_date : DateTimeInfo
        Date/time context used for event solving.
    observer : ObserverInfo
        Observer geodetic and atmospheric configuration.
    prayer_list : list[Prayer]
        Baseline prayer outputs before extreme-latitude corrections.
    sun_dec : Angle
        Apparent solar declination for the working date.

    Returns
    -------
    list[Prayer]
        Adjusted prayer sequence preserving original ordering.
    """

    if not prayer_list:
        return prayer_list

    # Local constants and small helpers kept inside per request
    IDX_SUNRISE = 1
    IDX_SUNSET = 4
    IDX_MAGHRIB = 5
    IDX_MIDNIGHT = -1

    STANDARD_SUN_ANGLE = Angle(5/6)

    FACTOR = {
        'MIDDLENIGHT': 1 / 2,
        'ONESEVENTH': 1 / 7,
        'ANGLEBASED': 1 / 60,
    }

    def _angle_for(name: str, method: PrayerMethod) -> Angle:
        if name == 'Fajr':
            return method.fajr_angle
        if name == 'Maghrib':
            return method.maghrib_angle
        return method.isha_angle  # for ʿIshāʾ
    
    method: PrayerMethod = prayer_list[0].method
    out = list(prayer_list)
    lat_sign: int = 1 if observer.latitude.decimal > 0 else -1

    lowest_angle = method.fajr_angle.decimal # in practice, always the lowest
    eps = se.oblique_eq(observer_date.jde)  # axial tilt at JDE
    nearest_lat = 90.0 - eps.decimal - lowest_angle - 0.01
    adjusted_obs = replace(observer, latitude=Angle(lat_sign*nearest_lat))

    # If only sunrise/sunset are undefined, repair those first and preserve
    # already valid non-solar prayers.
    polar_lat_ex = 90.0 - eps.decimal - 0.01

    try: 
        only_sunrise_sunset_undef = (
            all((prayer_list[i].time == math.inf or 'does not exist' in prayer_list[i].time) for i in [IDX_SUNRISE, IDX_SUNSET, IDX_MAGHRIB, IDX_MIDNIGHT]) # type: ignore[arg-type]
            and all((prayer_list[j].time != math.inf or not ('does not exist' in prayer_list[j].time)) for j in range(3) if j not in [IDX_SUNRISE, IDX_SUNSET, IDX_MAGHRIB, IDX_MIDNIGHT]) # type: ignore[arg-type]
        )

        if only_sunrise_sunset_undef:
            polar_obs_ex = replace(observer, latitude=Angle(lat_sign*abs(polar_lat_ex)))
            out[IDX_SUNRISE] = replace(out[IDX_SUNRISE], time=safe_sun_time(observer_date, polar_obs_ex, 'rise', STANDARD_SUN_ANGLE))
            out[IDX_SUNSET] = replace(out[IDX_SUNSET], time=safe_sun_time(observer_date, polar_obs_ex, 'set', STANDARD_SUN_ANGLE))
            out[IDX_MAGHRIB] = replace(out[IDX_MAGHRIB], time=safe_sun_time(observer_date, polar_obs_ex, 'set', method.maghrib_angle))
            out[IDX_MAGHRIB+1] = replace(out[IDX_MAGHRIB+1], time=safe_sun_time(observer_date, polar_obs_ex, 'set', method.isha_angle))

            if method.midnight_type:
                out[IDX_MIDNIGHT] = replace(out[IDX_MIDNIGHT], time=te.time_midpoint(out[IDX_SUNSET].time, find_tomorrow_time(observer_date, polar_obs_ex, method.fajr_angle))) # type: ignore[arg-type]
            else:
                out[IDX_MIDNIGHT] = replace(out[IDX_MIDNIGHT], time=te.time_midpoint(out[IDX_SUNSET].time, find_tomorrow_time(observer_date, polar_obs_ex))) # type: ignore[arg-type]

            if isinstance(out[IDX_SUNRISE+2].time, str):
                out[IDX_SUNRISE + 2] = replace(out[IDX_SUNRISE + 2], time=asr_time(out[IDX_SUNRISE + 1].time, Angle(lat_sign*abs(polar_lat_ex)), sun_dec, method.asr_type)) # type: ignore[arg-type]
            
            return out
    except (TypeError, AttributeError, ArithmeticError, ValueError):
        pass

    # Calculate latitude where minimum night/day length
    H_zero: Angle = Angle(15 / 2 * (24 - 2)) # 2 hour minimum time between sunset and sunrise.

    A = math.sin(sun_dec.radians)
    B = math.cos(sun_dec.radians) * math.cos(H_zero.radians)
    C = math.sin(Angle(-5/6).radians)
    R = math.hypot(A, B)

    if R == 0:
        raise ValueError("Degenerate case: δ=0 and H0=90°.")
    if abs(C) > R:
        raise ValueError("Target unattainable for these inputs.")

    u = max(-1.0, min(1.0, C / R))

    alpha = math.atan2(A, B)
    phi_plus = (alpha + math.acos(u)) * 180 / math.pi
    phi_minus = (alpha - math.acos(u)) * 180 / math.pi

    polar_lat = -999
    for phi in (phi_plus, phi_minus):
        if -90.0 <= phi <= 90.0:
            polar_lat = phi
    if polar_lat == -999: raise RuntimeError("No physical latitude in range.")
    
    polar_obs = replace(observer, latitude=Angle(lat_sign*abs(polar_lat)))

    # Recalculate sunset and sunrise times if latitude is above that limit
    if abs(observer.latitude.decimal) > abs(polar_lat):
        out[IDX_SUNRISE] = replace(out[IDX_SUNRISE], time=safe_sun_time(observer_date, polar_obs, 'rise', STANDARD_SUN_ANGLE))
        out[IDX_SUNSET] = replace(out[IDX_SUNSET], time=safe_sun_time(observer_date, polar_obs, 'set', STANDARD_SUN_ANGLE))

        out[IDX_SUNRISE + 2] = replace(out[IDX_SUNRISE + 2], time=asr_time(out[IDX_SUNRISE + 1].time, Angle(lat_sign*abs(polar_lat)), sun_dec, method.asr_type)) # type: ignore[arg-type]
        
    # Repair day-boundary ordering when polar adjustments invert sunrise/sunset.
    if out[IDX_SUNSET].time < out[IDX_SUNRISE].time: # type: ignore[arg-type]
        if observer_date.date.day == out[IDX_SUNRISE].time.day: # type: ignore[arg-type]
            if abs(observer.latitude.decimal) > abs(polar_lat):
                out[IDX_SUNSET] = replace(out[IDX_SUNSET], time=find_tomorrow_time(observer_date, polar_obs, STANDARD_SUN_ANGLE, 'set'))
            elif method.extreme_lats == 'NEARESTLAT':
                out[IDX_SUNSET] = replace(out[IDX_SUNSET], time=find_tomorrow_time(observer_date, adjusted_obs, STANDARD_SUN_ANGLE, 'set'))
        elif observer_date.date.day == out[IDX_SUNSET].time.day: # type: ignore[arg-type]
            if abs(observer.latitude.decimal) > abs(polar_lat):
                out[IDX_SUNRISE] = replace(out[IDX_SUNRISE], time=find_tomorrow_time(observer_date, polar_obs, STANDARD_SUN_ANGLE, 'rise', -1))
            elif method.extreme_lats == 'NEARESTLAT':
                out[IDX_SUNRISE] = replace(out[IDX_SUNRISE], time=find_tomorrow_time(observer_date, adjusted_obs, STANDARD_SUN_ANGLE, 'rise', -1))

    # Nearest latitude method
    if method.extreme_lats == 'NEARESTLAT':
        for i, prayer in enumerate(out):
            if prayer.name in ['Maghrib', 'ʿIshāʾ']:
                ang = _angle_for(prayer.name, method)
                out[i] = replace(prayer, time=safe_sun_time(observer_date, adjusted_obs, 'set', ang)) # type: ignore[arg-type]

                # Change sunset and sunrise times so they make more sense.
                if prayer.name == 'Maghrib':
                    if out[IDX_SUNRISE].time < out[0].time: # type: ignore[arg-type]
                        out[IDX_SUNRISE] = replace(out[IDX_SUNSET], time=safe_sun_time(observer_date, adjusted_obs, 'rise', STANDARD_SUN_ANGLE))
                    
                    if out[IDX_SUNSET].time > out[IDX_MAGHRIB].time: # type: ignore[arg-type]
                        out[IDX_SUNSET] = replace(out[IDX_SUNSET], time=safe_sun_time(observer_date, adjusted_obs, 'set', STANDARD_SUN_ANGLE))

            elif prayer.name == 'Fajr':
                out[i] = replace(prayer, time=safe_sun_time(observer_date, adjusted_obs, 'rise', method.fajr_angle)) # type: ignore[arg-type]
                if out[i].time.day != out[IDX_SUNRISE].time.day: # type: ignore[arg-type]
                    out[i] = replace(prayer, time=find_tomorrow_time(observer_date, adjusted_obs, method.fajr_angle, num_days=-1))
            
            elif prayer.name == 'Midnight':
                if method.midnight_type:
                    next_target = find_tomorrow_time(observer_date, adjusted_obs, method.fajr_angle) # type: ignore[arg-type]
                else:
                    next_target = find_tomorrow_time(observer_date, adjusted_obs) # type: ignore[arg-type]
                
                diff = -999
                if observer_date.date.day + 1 != next_target.day:
                    diff = observer_date.date.day + 2 - next_target.day
                    if method.midnight_type:
                        next_target = find_tomorrow_time(observer_date, adjusted_obs, method.fajr_angle, num_days=diff) # type: ignore[arg-type]
                    else:
                        next_target = find_tomorrow_time(observer_date, adjusted_obs, num_days=diff) # type: ignore[arg-type]

                out[i] = replace(prayer, time=te.time_midpoint(out[IDX_SUNSET].time, next_target))  # type: ignore[arg-type]

                if out[i].time < out[6].time: # type: ignore[arg-type]
                    next_target = find_tomorrow_time(observer_date, adjusted_obs)
                    out[i] = replace(prayer, time=te.time_midpoint(out[IDX_SUNSET].time, next_target)) # type: ignore[arg-type]

        return out

    # Portion-of-night clamping methods
    if method.extreme_lats in ['MIDDLENIGHT', 'ONESEVENTH', 'ANGLEBASED']:
        # Use sunrise and sunset already present in prayer_list like the original
        midnight_float = (timedelta(hours=24) + out[IDX_SUNRISE].time - out[IDX_SUNSET].time).total_seconds() / 3600.0 # type: ignore[arg-type]
        base_portion = FACTOR[method.extreme_lats] * midnight_float

        for i, prayer in enumerate(out):
            if prayer.name not in ['Fajr', 'Maghrib', 'ʿIshāʾ']:
                continue

            if prayer.name == 'Fajr':
                sign = -1
                base = out[IDX_SUNRISE]  # sunrise
            else:
                sign = 1
                base = out[IDX_SUNSET]   # sunset

            if method.extreme_lats == 'ANGLEBASED':
                ang = _angle_for(prayer.name, method)
                portion_mod = base_portion * ang.decimal
            else:
                portion_mod = base_portion

            # Match original: compare raw delta, not absolute value
            try:
                delta_h = (prayer.time - base.time).total_seconds() / 3600.0  # type: ignore[attr-defined]
            except Exception:
                delta_h = math.inf

            if delta_h == math.inf or delta_h > portion_mod or abs(observer.latitude.decimal) > abs(polar_lat):
                out[i] = replace(prayer, time=base.time + timedelta(hours=sign * portion_mod))  # type: ignore[attr-defined]

        # Midnight recompute exactly like original
        if abs(observer.latitude.decimal) > abs(polar_lat):
            next_sunrise: datetime = find_tomorrow_time(observer_date, polar_obs)
            next_sunset: datetime = find_tomorrow_time(observer_date, polar_obs, rise_or_set='set')
        else:
            next_sunrise: datetime = find_tomorrow_time(observer_date, adjusted_obs)
            next_sunset: datetime = find_tomorrow_time(observer_date, adjusted_obs, rise_or_set='set')
        
        if next_sunset < next_sunrise:
            if observer_date.date.day + 1 == next_sunrise.day:
                if abs(observer.latitude.decimal) > abs(polar_lat):
                    next_sunset = find_tomorrow_time(observer_date, polar_obs, rise_or_set='set', num_days=2)
                else:
                    next_sunset = find_tomorrow_time(observer_date, adjusted_obs, rise_or_set='set', num_days=2)
            if observer_date.date.day + 1 == next_sunset.day:
                if abs(observer.latitude.decimal) > abs(polar_lat):
                    next_sunrise = safe_sun_time(observer_date, polar_obs, 'rise', STANDARD_SUN_ANGLE)
                else:
                    next_sunrise = safe_sun_time(observer_date, polar_obs, 'rise', STANDARD_SUN_ANGLE)

        if getattr(method, 'midnight_type', False):
            try:
                next_fajr = find_tomorrow_time(observer_date, observer, method.fajr_angle)
            except Exception:
                next_fajr = math.inf

            midnight_float_next = (timedelta(hours=24) + next_sunrise - next_sunset).total_seconds() / 3600.0

            portion_next = FACTOR[method.extreme_lats] * midnight_float_next
            if method.extreme_lats == 'ANGLEBASED':
                portion_next *= method.fajr_angle.decimal

            try:
                delta_h_next = (next_fajr - next_sunrise).total_seconds() / 3600.0  # type: ignore[operator]
            except Exception:
                delta_h_next = math.inf

            if delta_h_next == math.inf or delta_h_next > portion_next:
                next_fajr = next_sunrise - timedelta(hours=portion_next)

            out[IDX_MIDNIGHT] = replace(out[IDX_MIDNIGHT], time=te.time_midpoint(out[IDX_SUNSET].time, next_fajr)) # type: ignore[arg-type]
        else:
            out[IDX_MIDNIGHT] = replace(out[IDX_MIDNIGHT], time=te.time_midpoint(out[IDX_SUNSET].time, next_sunrise)) # type: ignore[arg-type]

        return out

    # Fallback message
    msg = "The observer is at a latitude such that the sun does not reach the given angle for the prayer."
    return [replace(p, time=msg) for p in out]

def calculate_prayer_times(observer_date: DateTimeInfo, observer: ObserverInfo, sun_info: SunInfo, method: PrayerMethod) -> PrayerTimes:
    """Compute daily prayer times for a configured observer and method profile.

    Parameters
    ----------
    observer_date : DateTimeInfo
        Date/time context for astronomical calculations.
    observer : ObserverInfo
        Observer geodetic and atmospheric configuration.
    sun_info : SunInfo
        Precomputed solar events and coordinates for the working date.
    method : PrayerMethod
        Prayer-angle and rule configuration.

    Returns
    -------
    PrayerTimes
        Structured prayer output dataclass.

    Notes
    -----
    When a required solar event cannot be solved directly, this function delegates
    to ``extreme_latitudes`` for fallback behavior.
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
        except (ArithmeticError, TypeError, ValueError):
            maghrib_dt = math.inf
    else:
        # Otherwise Maghrib is sunset
        try:
            maghrib_dt = sun_info.sunset + timedelta(minutes=1)
        except (TypeError, ValueError):
            maghrib_dt = math.inf

    # Calculate ʿishāʾ time
    # If NOT makkah method
    if "Makkah" not in method.name:
        try:
            isha_dt = safe_sun_time(observer_date, observer,  'set', method.isha_angle)
            if isha_dt.time() < sun_info.sun_transit.time() and isha_dt.day == sun_info.sun_transit.day:
                isha_dt = find_tomorrow_time(observer_date, observer, method.isha_angle, 'set')
        except ArithmeticError:
            isha_dt = math.inf
    # Makkah method is special and will be the only exception
    else:
        # During Ramadan, ʿishāʾ is set to a flat two hours after maghrib
        if observer_date.hijri is not None and observer_date.hijri.hijri_month == 9:
            isha_dt = maghrib_dt + timedelta(hours=2) # type: ignore[arg-type]

        # Otherwise, it is set to a flat 1.5 hours
        else:
            isha_dt = maghrib_dt + timedelta(hours=1.5) # type: ignore[arg-type]
    
    # Calculate Midnight time
    try:
        if method.midnight_type:
            midnight_dt = te.time_midpoint(sun_info.sunset, find_tomorrow_time(observer_date, observer, method.fajr_angle))
        else:
            midnight_dt = te.time_midpoint(sun_info.sunset, find_tomorrow_time(observer_date, observer))
    except Exception:
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
    if abs(observer.latitude.decimal) > 46.5: # 90 - obliquity (23.5) - max fajr angle (20)
        if any(p.time == math.inf for p in prayer_list):
            warnings.warn("Extreme latitude warning. Prayer times at this latitude are not well established.")
            prayer_list = extreme_latitudes(observer_date, observer, prayer_list, sun_info.apparent_declination)

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
