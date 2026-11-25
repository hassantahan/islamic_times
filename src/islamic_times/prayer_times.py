"""Prayer-time calculation algorithms and extreme-latitude adjustments.

This module computes the canonical prayer-time set used by :class:`ITLocation`:
Fajr, Sunrise, Ẓuhr, ʿAṣr, Sunset, Maghrib, ʿIshāʾ, and Islamic midnight.
It includes method profiles from PrayTimes, handling for Makkah-specific ʿIshāʾ
rules, and fallback strategies for high-latitude locations where standard
sun-angle events may be unavailable.
"""

import math
import warnings
from dataclasses import dataclass, replace
from datetime import datetime, time, timedelta
from typing import List

from islamic_times import astro_core
from islamic_times import sun_equations as se
from islamic_times import time_equations as te
from islamic_times.it_dataclasses import (
    Angle, AsrConvention, DateTimeInfo, ObserverInfo, Prayer, PrayerMethod, PrayerTimes, SunInfo,
)
from islamic_times.prayer_fadilah import calculate_sistani_daylight_fadilah

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
        fajr_angle=Angle(18.5), isha_angle=Angle(999),
        isha_offset_minutes=90, ramadan_isha_offset_minutes=120,
    ),
    PrayerMethod(
        name="University of Islamic Sciences, (Karachi)", 
        keys=('KARACHI', 'UNIVERSITY OF ISLAMIC SCIENCES', 'UIS'), 
        fajr_angle=Angle(18), isha_angle=Angle(18)
    ),
    PrayerMethod(
        name="Institute of Geophysics, University of Tehran (Tehran)", 
        keys=('TEHRAN', 'UNIVERSITY OF TEHRAN', 'INSTITUTE OF GEOPHYSICS', 'UOT', 'IOG', 'IOGUOT', 'UOTIOG', 'IGUT', 'UTIG'), 
        fajr_angle=Angle(17.7), isha_angle=Angle(14), maghrib_angle=Angle(4.5), midnight_type=1,
        asr_type=AsrConvention.SISTANI_PREFERRED,
    ),
    PrayerMethod(name="Shia Ithna Ashari, Leva Research Institute, Qom (Jafari)",
        keys=('JAFARI', 'JAAFARI', 'SHIA', 'SHIA ITHNA ASHARI', 'LEVA', 'LEVA RESEARCH INSTITUTE', 'QOM', 'QUM', 'LRI', 'SIA', 'SIALRI'), 
        fajr_angle=Angle(16), isha_angle=Angle(14), maghrib_angle=Angle(4), midnight_type=1,
        asr_type=AsrConvention.SISTANI_PREFERRED,
    ),
    PrayerMethod(
        name="Muslims of France", 
        keys=('FRANCE',), 
        fajr_angle=Angle(12), isha_angle=Angle(12)
    ),
    PrayerMethod(
        name="Spiritual Administration of Muslims of Russia", 
        keys=('RUSSIA',), 
        fajr_angle=Angle(16), isha_angle=Angle(15)
    ),
    PrayerMethod(
        name="Islamic Religious Council of Singapore", 
        keys=('SINGAPORE',), 
        fajr_angle=Angle(20), isha_angle=Angle(18)
    ),
]

HIGH_LATITUDE_RULES = ("NONE", "NEARESTLAT", "MIDDLENIGHT", "ONESEVENTH", "ANGLEBASED")


def prayer_method_catalog() -> tuple[PrayerMethod, ...]:
    """Return the supported presets in display order, including their accepted aliases."""
    return tuple(DEFAULT_PRAYER_METHODS)


def _build_prayer_method_alias_map() -> dict[str, PrayerMethod]:
    """Index the preset catalog by case-insensitive aliases."""
    aliases: dict[str, PrayerMethod] = {}
    for method in prayer_method_catalog():
        for key in method.keys or ():
            aliases.setdefault(key.upper(), method)
    return aliases


_PRAYER_METHOD_ALIASES = _build_prayer_method_alias_map()


def get_prayer_method(method_key: str) -> PrayerMethod:
    """Return the immutable preset matching a primary key or alias."""
    if not isinstance(method_key, str):
        raise TypeError("'method_key' must be a string.")
    key = method_key.strip().upper()
    method = _PRAYER_METHOD_ALIASES.get(key)
    if method is None:
        valid_options = ", ".join(method.keys[0] for method in prayer_method_catalog() if method.keys)
        raise ValueError(f"Invalid prayer method '{key}'. Valid options are: {valid_options}")
    return method


@dataclass(frozen=True, slots=True)
class ExtremeLatitudeContext:
    """Shared context passed to extreme-latitude strategy functions."""

    observer_date: DateTimeInfo
    observer: ObserverInfo
    method: PrayerMethod
    standard_sun_angle: Angle
    idx_sunrise: int
    idx_sunset: int
    idx_maghrib: int
    idx_midnight: int
    polar_lat: float
    adjusted_obs: ObserverInfo
    polar_obs: ObserverInfo


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
    return se.find_proper_suntime(observer_date, observer, event, angle)
    
# Used by midnight calculations to sample the next/previous day.
def find_tomorrow_time(
    observer_date: DateTimeInfo,
    observer: ObserverInfo,
    angle: Angle = Angle(5 / 6),
    rise_or_set: str = 'rise',
    num_days: int = 1,
) -> datetime:
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
    day = observer_date.date.date() + timedelta(days=num_days)
    noon = datetime.combine(day, time(12))
    local_noon = te.localize_or_convert_datetime(noon, observer_date.date.tzinfo)
    return safe_sun_time(te.observer_datetime_info(local_noon), observer, rise_or_set, angle)

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
    if ts not in (AsrConvention.STANDARD, AsrConvention.HANAFI):
        raise ValueError("The shadow-factor calculation requires the standard or Ḥanafī convention.")
    shadow_factor = ts + 1
    
    # Calculate noon shadow length
    noon_shadow = abs(math.tan(math.radians(lat.decimal - dec.decimal)))
    
    # Shadow lengths are expressed in units of the object's height.
    total_shadow = shadow_factor + noon_shadow
    
    # Calculate angle using inverse tangent
    angle_rad = math.atan(1.0 / total_shadow)
    
    # Check if ʿAṣr is possible at this location and date
    temp_num = (math.sin(angle_rad) - math.sin(lat.radians) * math.sin(dec.radians))
    temp_denom = (math.cos(lat.radians) * math.cos(dec.radians))
    
    if abs(temp_num) > abs(temp_denom):
        raise ArithmeticError("ʿAṣr prayer time cannot be calculated - sun never reaches required shadow length")
    
    # Convert to hours from noon
    asr_hours = 1/15.0 * math.degrees(math.acos(temp_num / temp_denom))
    
    return te.add_elapsed_time(noon, timedelta(hours=asr_hours))


def _angle_for(name: str, method: PrayerMethod) -> Angle:
    """Return the configured twilight depression for a prayer."""
    if name == 'Fajr':
        return method.fajr_angle
    if name == 'Maghrib':
        return method.maghrib_angle
    return method.isha_angle


def _isha_after_maghrib(maghrib: datetime | str | float, minutes: int) -> datetime | float:
    """Calculate a fixed ʿIshāʾ interval only when Maghrib has a known time."""
    if isinstance(maghrib, datetime):
        return te.add_elapsed_time(maghrib, timedelta(minutes=minutes))
    return math.inf


def _require_time(prayer: Prayer) -> datetime:
    """Return a fallback anchor or report why the rule cannot proceed."""
    when = prayer.when
    if when is None:
        raise ArithmeticError(f"{prayer.name} is unavailable as a reference event.")
    return when


def _strategy_nearest_lat(ctx: ExtremeLatitudeContext, out: List[Prayer]) -> List[Prayer]:
    """Solve twilight events using the reference latitude."""
    date = ctx.observer_date
    reference = ctx.adjusted_obs
    method = ctx.method
    sunrise = _require_time(out[ctx.idx_sunrise])
    sunset = _require_time(out[ctx.idx_sunset])

    fajr = safe_sun_time(date, reference, 'rise', method.fajr_angle)
    if fajr.day != sunrise.day:
        fajr = find_tomorrow_time(date, reference, method.fajr_angle, num_days=-1)
    out[0] = replace(out[0], time=fajr)

    maghrib = safe_sun_time(date, reference, 'set', method.maghrib_angle)
    out[ctx.idx_maghrib] = replace(out[ctx.idx_maghrib], time=maghrib)
    if sunrise < fajr:
        sunrise = safe_sun_time(date, reference, 'rise', ctx.standard_sun_angle)
        out[ctx.idx_sunrise] = replace(out[ctx.idx_sunrise], time=sunrise)
    if sunset > maghrib:
        sunset = safe_sun_time(date, reference, 'set', ctx.standard_sun_angle)
        out[ctx.idx_sunset] = replace(out[ctx.idx_sunset], time=sunset)

    isha_index = ctx.idx_maghrib + 1
    if method.isha_offset_minutes is None:
        isha = safe_sun_time(date, reference, 'set', method.isha_angle)
        out[isha_index] = replace(out[isha_index], time=isha)

    midnight_angle = method.fajr_angle if method.midnight_type else ctx.standard_sun_angle
    next_target = find_tomorrow_time(date, reference, midnight_angle)
    midnight = te.time_midpoint(sunset, next_target)
    isha = out[isha_index].when
    if isha is not None and midnight < isha:
        midnight = te.time_midpoint(sunset, find_tomorrow_time(date, reference))
    out[ctx.idx_midnight] = replace(out[ctx.idx_midnight], time=midnight)
    return out


def _strategy_portion_of_night(ctx: ExtremeLatitudeContext, out: List[Prayer]) -> List[Prayer]:
    """Estimate twilight events from a fraction of the reference night."""
    date = ctx.observer_date
    method = ctx.method
    fraction = {'MIDDLENIGHT': 1 / 2, 'ONESEVENTH': 1 / 7, 'ANGLEBASED': 1 / 60}[method.extreme_lats]
    above_polar_lat = abs(ctx.observer.latitude.decimal) > abs(ctx.polar_lat)
    sunrise = _require_time(out[ctx.idx_sunrise])
    sunset = _require_time(out[ctx.idx_sunset])
    # The rule estimates this night's duration from rise and set on one civil day.
    day = timedelta(days=1)
    night_hours = (day + sunrise - sunset).total_seconds() / 3600

    for i, prayer in enumerate(out):
        if prayer.name not in ('Fajr', 'Maghrib', 'ʿIshāʾ'):
            continue
        if prayer.name == 'ʿIshāʾ' and method.isha_offset_minutes is not None:
            continue

        base, direction = (sunrise, -1) if prayer.name == 'Fajr' else (sunset, 1)
        portion = fraction * night_hours
        if method.extreme_lats == 'ANGLEBASED':
            portion *= _angle_for(prayer.name, method).decimal
        when = prayer.when
        if when is None or above_polar_lat or (when - base).total_seconds() / 3600 > portion:
            out[i] = replace(prayer, time=base + timedelta(hours=direction * portion))

    reference = ctx.polar_obs if above_polar_lat else ctx.adjusted_obs
    next_sunrise = find_tomorrow_time(date, reference)
    next_sunset = find_tomorrow_time(date, reference, rise_or_set='set')
    if next_sunset < next_sunrise:
        # Handle reference events that cross the observer's day boundary.
        if date.date.day + 1 == next_sunrise.day:
            next_sunset = find_tomorrow_time(date, reference, rise_or_set='set', num_days=2)
        if date.date.day + 1 == next_sunset.day:
            next_sunrise = safe_sun_time(date, ctx.polar_obs, 'rise', ctx.standard_sun_angle)

    next_target = next_sunrise
    if method.midnight_type:
        try:
            next_fajr = find_tomorrow_time(date, ctx.observer, method.fajr_angle)
        except ArithmeticError:
            next_fajr = None

        next_night_hours = (day + next_sunrise - next_sunset).total_seconds() / 3600
        portion = fraction * next_night_hours
        if method.extreme_lats == 'ANGLEBASED':
            portion *= method.fajr_angle.decimal
        if next_fajr is None or (next_fajr - next_sunrise).total_seconds() / 3600 > portion:
            next_fajr = next_sunrise - timedelta(hours=portion)
        next_target = next_fajr

    out[ctx.idx_midnight] = replace(out[ctx.idx_midnight], time=te.time_midpoint(sunset, next_target))
    return out

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

    IDX_ZUHR = 2
    IDX_ASR = 3
    IDX_SUNRISE = 1
    IDX_SUNSET = 4
    IDX_MAGHRIB = 5
    IDX_MIDNIGHT = -1

    STANDARD_SUN_ANGLE = Angle(5/6)

    method: PrayerMethod = prayer_list[0].method
    if method.extreme_lats == 'NONE':
        return prayer_list
    out = list(prayer_list)
    lat_sign: int = 1 if observer.latitude.decimal > 0 else -1

    # The nearest-latitude rule uses the selected Fajr depression as its reference.
    fajr_angle = method.fajr_angle.decimal
    obliquity = astro_core.mean_obliquity(observer_date.jde)
    nearest_lat = 90.0 - obliquity - fajr_angle - 0.01
    adjusted_obs = replace(observer, latitude=Angle(lat_sign*nearest_lat))

    # If only sunrise/sunset are undefined, repair those first and preserve
    # already valid non-solar prayers.
    polar_lat_ex = 90.0 - obliquity - 0.01

    try: 
        only_sunrise_sunset_undef = (
            all(prayer_list[i].when is None for i in (IDX_SUNRISE, IDX_SUNSET, IDX_MAGHRIB, IDX_MIDNIGHT))
            and all(prayer_list[i].when is not None for i in (0, IDX_ZUHR))
        )

        if only_sunrise_sunset_undef:
            polar_obs_ex = replace(observer, latitude=Angle(lat_sign*abs(polar_lat_ex)))
            out[IDX_SUNRISE] = replace(out[IDX_SUNRISE], time=safe_sun_time(observer_date, polar_obs_ex, 'rise', STANDARD_SUN_ANGLE))
            sunset = safe_sun_time(observer_date, polar_obs_ex, 'set', STANDARD_SUN_ANGLE)
            out[IDX_SUNSET] = replace(out[IDX_SUNSET], time=sunset)
            out[IDX_MAGHRIB] = replace(out[IDX_MAGHRIB], time=safe_sun_time(observer_date, polar_obs_ex, 'set', method.maghrib_angle))
            if method.isha_offset_minutes is None:
                out[IDX_MAGHRIB+1] = replace(out[IDX_MAGHRIB+1], time=safe_sun_time(observer_date, polar_obs_ex, 'set', method.isha_angle))

            midnight_angle = method.fajr_angle if method.midnight_type else STANDARD_SUN_ANGLE
            next_target = find_tomorrow_time(observer_date, polar_obs_ex, midnight_angle)
            out[IDX_MIDNIGHT] = replace(out[IDX_MIDNIGHT], time=te.time_midpoint(sunset, next_target))

            if method.asr_type != AsrConvention.SISTANI_PREFERRED and out[IDX_ASR].when is None:
                asr = asr_time(_require_time(out[IDX_ZUHR]), polar_obs_ex.latitude, sun_dec, method.asr_type)
                out[IDX_ASR] = replace(out[IDX_ASR], time=asr)
            
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

        if method.asr_type != AsrConvention.SISTANI_PREFERRED:
            asr = asr_time(_require_time(out[IDX_ZUHR]), polar_obs.latitude, sun_dec, method.asr_type)
            out[IDX_ASR] = replace(out[IDX_ASR], time=asr)
        
    # Repair day-boundary ordering when polar adjustments invert sunrise/sunset.
    sunset = _require_time(out[IDX_SUNSET])
    sunrise = _require_time(out[IDX_SUNRISE])
    if sunset < sunrise:
        if observer_date.date.day == sunrise.day:
            if abs(observer.latitude.decimal) > abs(polar_lat):
                out[IDX_SUNSET] = replace(out[IDX_SUNSET], time=find_tomorrow_time(observer_date, polar_obs, STANDARD_SUN_ANGLE, 'set'))
            elif method.extreme_lats == 'NEARESTLAT':
                out[IDX_SUNSET] = replace(out[IDX_SUNSET], time=find_tomorrow_time(observer_date, adjusted_obs, STANDARD_SUN_ANGLE, 'set'))
        elif observer_date.date.day == sunset.day:
            if abs(observer.latitude.decimal) > abs(polar_lat):
                out[IDX_SUNRISE] = replace(out[IDX_SUNRISE], time=find_tomorrow_time(observer_date, polar_obs, STANDARD_SUN_ANGLE, 'rise', -1))
            elif method.extreme_lats == 'NEARESTLAT':
                out[IDX_SUNRISE] = replace(out[IDX_SUNRISE], time=find_tomorrow_time(observer_date, adjusted_obs, STANDARD_SUN_ANGLE, 'rise', -1))

    ctx = ExtremeLatitudeContext(
        observer_date=observer_date,
        observer=observer,
        method=method,
        standard_sun_angle=STANDARD_SUN_ANGLE,
        idx_sunrise=IDX_SUNRISE,
        idx_sunset=IDX_SUNSET,
        idx_maghrib=IDX_MAGHRIB,
        idx_midnight=IDX_MIDNIGHT,
        polar_lat=polar_lat,
        adjusted_obs=adjusted_obs,
        polar_obs=polar_obs,
    )

    if method.extreme_lats == 'NEARESTLAT':
        return _strategy_nearest_lat(ctx, out)

    if method.extreme_lats in ['MIDDLENIGHT', 'ONESEVENTH', 'ANGLEBASED']:
        return _strategy_portion_of_night(ctx, out)

    return out

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
        Structured prayer output dataclass. When a high-latitude fallback path
        is used, ``extreme_latitude_applied``/``extreme_latitude_rule``/
        ``extreme_latitude_reason`` are populated.

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
    
    asr_fadilah = None
    uses_fadilah = method.asr_type == AsrConvention.SISTANI_PREFERRED
    if uses_fadilah:
        asr_fadilah = calculate_sistani_daylight_fadilah(observer, sun_info).get("asr")
        asr_dt = asr_fadilah.start if asr_fadilah else "Preferred ʿAṣr window unavailable."
    else:
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

    # Fixed-interval methods derive ʿIshāʾ from Maghrib after latitude adjustment.
    isha_minutes = method.isha_offset_minutes
    if observer_date.hijri is not None and observer_date.hijri.hijri_month == 9:
        if method.ramadan_isha_offset_minutes is not None:
            isha_minutes = method.ramadan_isha_offset_minutes
    if isha_minutes is None:
        try:
            isha_dt = safe_sun_time(observer_date, observer,  'set', method.isha_angle)
            if isha_dt.time() < sun_info.sun_transit.time() and isha_dt.day == sun_info.sun_transit.day:
                isha_dt = find_tomorrow_time(observer_date, observer, method.isha_angle, 'set')
        except ArithmeticError:
            isha_dt = math.inf
    else:
        isha_dt = _isha_after_maghrib(maghrib_dt, isha_minutes)
    
    # Calculate Midnight time
    try:
        if method.midnight_type:
            midnight_dt = te.time_midpoint(sun_info.sunset, find_tomorrow_time(observer_date, observer, method.fajr_angle))
        else:
            midnight_dt = te.time_midpoint(sun_info.sunset, find_tomorrow_time(observer_date, observer))
    except (ArithmeticError, TypeError, ValueError):
        midnight_dt = math.inf
    

    prayer_list = [
        Prayer("Fajr", fajr_dt, method),
        Prayer("Sunrise", sun_info.sunrise, method),
        Prayer("Ẓuhr", sun_info.sun_transit, method),
        Prayer("ʿAṣr", asr_dt, method, fadilah=asr_fadilah),
        Prayer("Sunset", sun_info.sunset, method),
        Prayer("Maghrib", maghrib_dt, method),
        Prayer("ʿIshāʾ", isha_dt, method),
        Prayer("Midnight", midnight_dt, method),
    ]
    extreme_latitude_rule: str | None = None
    extreme_latitude_reason: str | None = None

    # Deal with locations at extreme latitudes
    fallback_prayers = prayer_list[:3] + prayer_list[4:] if uses_fadilah else prayer_list
    if abs(observer.latitude.decimal) > 46.5 and any(prayer.when is None for prayer in fallback_prayers):
        warnings.warn("Extreme latitude warning. Prayer times at this latitude are not well established.", stacklevel=2)
        try:
            adjusted = extreme_latitudes(observer_date, observer, prayer_list, sun_info.apparent_declination)
        except ArithmeticError:
            # A requested reference solar event can also be unavailable.
            extreme_latitude_reason = "The selected high-latitude rule could not resolve a reference solar event."
        else:
            prayer_list = [
                replace(prayer, fallback_rule=method.extreme_lats)
                if prayer is not before and prayer.when is not None else prayer
                for prayer, before in zip(adjusted, prayer_list)
            ]

    if isha_minutes is not None:
        sunset = prayer_list[4]
        maghrib = prayer_list[5]
        if method.maghrib_angle.decimal == 0 and sunset.when is not None:
            maghrib = replace(
                maghrib, time=te.add_elapsed_time(sunset.when, timedelta(minutes=1)), fallback_rule=sunset.fallback_rule
            )
            prayer_list[5] = maghrib
        prayer_list[6] = replace(
            prayer_list[6],
            time=_isha_after_maghrib(maghrib.time, isha_minutes),
            fallback_rule=maghrib.fallback_rule,
        )

    extreme_latitude_applied = any(prayer.fallback_rule is not None for prayer in prayer_list)
    if extreme_latitude_applied:
        extreme_latitude_rule = method.extreme_lats
        extreme_latitude_reason = (
            "One or more required solar events were unavailable for standard-angle solving at this latitude/date."
        )

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
        extreme_latitude_applied=extreme_latitude_applied,
        extreme_latitude_rule=extreme_latitude_rule,
        extreme_latitude_reason=extreme_latitude_reason,
    )
