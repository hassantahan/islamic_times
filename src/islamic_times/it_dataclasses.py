"""
Contents:
    - Angle: Represents angles in decimal degrees with DMS and radians conversion.
    - RightAscension: Handles right ascension values and their conversions.
    - DistanceUnit & DistanceUnits: Define distance units and conversion utilities.
    - Distance: Represents physical distances with conversion support.
    - ObserverInfo: Stores observer's geographical and environmental parameters.
    - IslamicDateInfo: Stores and formats Islamic (Hijri) calendar information.
    - DateTimeInfo: Combines Gregorian, Islamic, Julian date, and delta-T data.
    - PrayerMethod: Encapsulates parameters for prayer time calculation methods.
    - Prayer: Represents a single prayer time and its calculation context.
    - PrayerTimes: Aggregates all prayer times for a given day.
    - MeccaInfo: Provides Qibla direction and distance.
    - SunInfo: Stores information about solar positions and events.
    - MoonInfo: Stores information about lunar positions, events, and illumination.
    - Visibilities: Stores information about new moon visibility calculations.
"""


import math
from typing import Tuple, ClassVar, Any
from datetime import datetime
from dataclasses import dataclass


PUBLIC_SCHEMA_VERSION: str = "1.0"


def _serialize_time_value(value: datetime | str | float) -> dict[str, Any]:
    """Serialize prayer/sun/moon time values with explicit kind metadata."""
    if isinstance(value, datetime):
        return {
            "kind": "datetime",
            "iso": value.isoformat(),
            "display": value.strftime("%X %d-%m-%Y"),
        }
    if isinstance(value, str):
        return {
            "kind": "message",
            "message": value,
        }
    if value == math.inf:
        return {
            "kind": "missing",
            "message": "Does not exist.",
        }
    return {
        "kind": "numeric",
        "value": float(value),
    }

# === Helper Classes ===

@dataclass(frozen=True, slots=True)
class Angle:
    """
    Represents an angle in decimal degrees.

    Provides utilities for conversion to degrees-minutes-seconds (DMS) format and radians.

    Attributes:
        decimal (float): The angle in decimal degrees.
        dms (Tuple[int, int, float]): The angle in degrees, minutes, and seconds.
        radians (float): The angle in radians.
        dms_str (str): The angle formatted as a DMS string.
    """

    decimal: float

    @property
    def dms(self) -> Tuple[int, int, float]:
        """
        Returns:
            Tuple[int, int, float]: The angle in (degrees, minutes, seconds).
        """
        degrees = int(self.decimal)
        abs_decimal = abs(self.decimal)
        minutes_full = (abs_decimal - abs(degrees)) * 60
        minutes = int(minutes_full)
        seconds = (minutes_full - minutes) * 60
        return (degrees, minutes, seconds)

    @property
    def radians(self) -> float:
        """
        Returns:
            float: The angle in radians.
        """
        return math.radians(self.decimal)

    def dms_str(self) -> str:
        """
        Returns:
            str: The angle formatted as a DMS string.
        """
        d, m, s = self.dms
        return f"{d:+04}\u00b0 {m:02}\u2032 {s:05.2f}\u2033"

    def to_dict(self) -> dict[str, Any]:
        """Serialize angle with stable numeric and display fields."""
        d, m, s = self.dms
        return {
            "decimal_degrees": self.decimal,
            "radians": self.radians,
            "dms": {
                "degrees": d,
                "minutes": m,
                "seconds": s,
                "display": self.dms_str(),
            },
        }

    def __str__(self):
        return f"{self.decimal:+08.3f}\u00b0\t\t({self.dms_str()})"

@dataclass(frozen=True, slots=True)
class RightAscension:
    """
    Represents Right Ascension (RA) in decimal hours.

    Provides utilities to convert RA to HMS, DMS, degrees, and radians.

    Attributes:
        decimal_hours (float): Right ascension in decimal hours.
        hms (Tuple[int, int, float]): Right ascension in hours, minutes, and seconds.
        decimal_degrees (Angle): Right ascension in degrees.
        radians (float): Right ascension in radians.
        dms (Tuple[int, int, float]): Right ascension in degrees, minutes, and seconds.
        dms_str (str): Right ascension formatted as a DMS string.
    """
    decimal_hours: float

    @property
    def hms(self) -> Tuple[int, int, float]:
        """
        Returns:
            Tuple[int, int, float]: RA in (hours, minutes, seconds).
        """
        hours = int(self.decimal_hours)
        minutes_full = (self.decimal_hours - hours) * 60
        minutes = int(minutes_full)
        seconds = (minutes_full - minutes) * 60
        return (hours, minutes, seconds)

    def hms_str(self) -> str:
        """
        Returns
            str: Pretty formatted RA in HMS format
        """
        h, m, s = self.hms
        return f"{h:02}h {m:02}m {s:05.2f}s"

    @property
    def decimal_degrees(self) -> Angle:
        """
        Returns:
            Angle: RA converted to degrees.
        """
        return Angle(self.decimal_hours * 15)
    
    @property
    def radians(self) -> float:
        """
        Returns:
            float: RA converted to radians.
        """
        return self.decimal_degrees.radians

    @property
    def dms(self) -> Tuple[int, int, float]:
        deg = self.decimal_degrees.decimal
        degrees = int(deg)
        decimal_minutes = (abs(deg) - degrees) * 60
        minutes = int(decimal_minutes)
        seconds = (decimal_minutes - minutes) * 60
        return (degrees, minutes, seconds)

    def dms_str(self) -> str:
        d, m, s = self.dms
        return f"{d:+04}\u00b0 {m:02}\u2032 {s:05.2f}\u2033"

    def to_dict(self) -> dict[str, Any]:
        """Serialize right ascension in both hour and degree forms."""
        h, m, s = self.hms
        return {
            "decimal_hours": self.decimal_hours,
            "decimal_degrees": self.decimal_degrees.decimal,
            "radians": self.radians,
            "hms": {
                "hours": h,
                "minutes": m,
                "seconds": s,
                "display": self.hms_str(),
            },
            "dms": {
                "display": self.dms_str(),
            },
        }

    def __str__(self):
        return f"{self.hms_str()}\t\t({self.decimal_degrees})"

@dataclass(frozen=True, slots=True)
class DistanceUnit:
    """
    Represents a unit of distance and provides conversion to other distance units.

    Attributes:
        name (str): The unit's full name.
        symbol (str): The unit's symbol.
        to_m_factor (float): The conversion factor to meters.
    """

    name: str
    symbol: str
    to_m_factor: float

    def convert_to(self, value: float, target_unit: 'DistanceUnit') -> float:
        """
        Converts a value from this unit to another.

        Args:
            value (float): Value in this unit.
            target_unit (DistanceUnit): Target unit.

        Returns:
            float: Converted value in target unit.
        """

        value_in_m = value * self.to_m_factor
        return value_in_m / target_unit.to_m_factor

    def __str__(self):
        return self.symbol

class DistanceUnits:
    """Standard distance units used throughout the application.
    
    Attributes:
        METRE (DistanceUnit): Meter unit.
        KILOMETRE (DistanceUnit): Kilometer unit.
        AU (DistanceUnit): Astronomical unit.
        NAUTICAL_MILE (DistanceUnit): Nautical mile unit.
        MILE (DistanceUnit): Mile unit.
        FOOT (DistanceUnit): Foot unit.
    """
    METRE = DistanceUnit("meter", "m", 1)
    KILOMETRE = DistanceUnit("kilometer", "km", 1_000)
    AU = DistanceUnit("astronomical unit", "AU", 149_597_870_691)
    NAUTICAL_MILE = DistanceUnit("nautical mile", "nmi", 1_852)
    MILE = DistanceUnit("mile", "mi", 1_609.3445)
    FOOT = DistanceUnit("foot", "ft", 0.3048)

@dataclass(frozen=True, slots=True)
class Distance:
    """
    Represents a physical distance with unit conversion utilities.

    Attributes:
        value (float): The numerical value.
        unit (DistanceUnit): The distance unit.
        in_unit (float): Converts the distance to the target unit.
        to (Distance): Returns a new Distance object in the target unit.
    """
    value: float
    unit: DistanceUnit = DistanceUnits.METRE

    def in_unit(self, target_unit: DistanceUnit) -> float:
        """
        Converts the distance into another unit.

        Args:
            target_unit (DistanceUnit): Target unit.

        Returns:
            float: Converted distance value.
        """
        return self.unit.convert_to(self.value, target_unit)

    def to(self, target_unit: DistanceUnit) -> 'Distance':
        """
        Converts and returns a new Distance object with the specified unit.

        Args:
            target_unit (DistanceUnit): Target unit.

        Returns:
            Distance: New Distance object.
        """
        return Distance(self.in_unit(target_unit), target_unit)

    def to_dict(self) -> dict[str, Any]:
        """Serialize distance preserving explicit unit metadata."""
        return {
            "value": self.value,
            "unit": {
                "name": self.unit.name,
                "symbol": self.unit.symbol,
                "to_m_factor": self.unit.to_m_factor,
            },
            "meters": self.in_unit(DistanceUnits.METRE),
        }

    def __str__(self):
        if self.unit == DistanceUnits.METRE:
            return f"{self.value:.2f} {self.unit}"
        else:
            if self.unit == DistanceUnits.AU:
                km_equiv = self.in_unit(DistanceUnits.KILOMETRE)
                return f"{self.value:.6f} {self.unit} ({km_equiv:.2} km)"
            else:
                m_equiv = self.in_unit(DistanceUnits.METRE)
                return f"{self.value:.2f} {self.unit} ({m_equiv:.2f} m)"


# === Data Classes ===

@dataclass(frozen=True, slots=True)
class ObserverInfo:
    """
    Contains geographical and environmental information of the observer.

    Attributes:
        latitude (Angle): Observer's latitude.
        longitude (Angle): Observer's longitude.
        elevation (Distance): Observer's elevation.
        pressure (float): Atmospheric pressure in kPa.
        temperature (float): Temperature in Celsius.
    """

    latitude: Angle
    longitude: Angle
    elevation: Distance
    pressure: float = 101.325
    temperature: float = 10

    def to_dict(self) -> dict[str, Any]:
        """Serialize observer configuration with stable field names."""
        return {
            "schema_version": PUBLIC_SCHEMA_VERSION,
            "type": "ObserverInfo",
            "latitude": self.latitude.to_dict(),
            "longitude": self.longitude.to_dict(),
            "elevation": self.elevation.to_dict(),
            "pressure_kpa": self.pressure,
            "temperature_c": self.temperature,
        }

    def __str__(self):
        return ("Observer Parameters\n"
                f"\tLatitude:\t\t{self.latitude.decimal:09.6f}\u00b0\n"
                f"\tLongitude:\t\t{self.longitude.decimal:09.6f}\u00b0\n"
                f"\tElevation:\t\t{self.elevation}\n"
                f"\tPressure:\t\t{self.pressure} kPa\n"
                f"\tTemperature:\t\t{self.temperature:+06.2f}°C")

@dataclass(frozen=True, slots=True)
class IslamicDateInfo:
    """
    Stores information about the Islamic (Hijri) date.

    Attributes:
        hijri_year (int): Hijri year.
        hijri_month (int): Hijri month.
        hijri_day (int): Hijri day.
    """

    ISLAMIC_MONTHS: ClassVar[dict[int, str]] = {
        1: 'Muḥarram',
        2: 'Ṣaffar',
        3: 'Rabīʿ al-Awwal',
        4: 'Rabīʿ al-Thānī',
        5: 'Jumādā al-Ūlā',
        6: 'Jumādā al-Thāniyah',
        7: 'Rajab',
        8: 'Shaʿbān',
        9: 'Ramaḍān',
        10: 'Shawwāl',
        11: 'Dhū al-Qaʿdah',
        12: 'Dhū al-Ḥijjah',
    }
    ISLAMIC_DAYS: ClassVar[dict[str, str]] = {
            'Sunday': 'al-Aḥad',
            'Monday': 'al-Ithnayn',
            'Tuesday': 'al-Thulāthāʾ',
            'Wednesday': 'al-Arbiʿāʾ',
            'Thursday': 'al-Khamīs',
            'Friday': 'al-Jumuʿah',
            'Saturday': 'al-Sabt',
        }

    hijri_year: int
    hijri_month: int
    hijri_day: int

    @property
    def hijri_month_name(self) -> str:
        return self.ISLAMIC_MONTHS[self.hijri_month]
    
    def hijri_day_of_week_name(self, day_of_week: str) -> str:
        return self.ISLAMIC_DAYS[day_of_week]
    
    def full_date(self, day_of_week: str) -> str:
        return f"{self.hijri_day_of_week_name(day_of_week)}, {self.hijri_day} {self.hijri_month_name}, {self.hijri_year}"

    def to_dict(self) -> dict[str, Any]:
        """Serialize Hijri date with canonical numeric/date-name fields."""
        return {
            "schema_version": PUBLIC_SCHEMA_VERSION,
            "type": "IslamicDateInfo",
            "hijri_year": self.hijri_year,
            "hijri_month": self.hijri_month,
            "hijri_day": self.hijri_day,
            "hijri_month_name": self.hijri_month_name,
        }

@dataclass(frozen=True, slots=True)
class DateTimeInfo:
    """
    Represents Gregorian and Islamic date-time information along with JD and ΔT.

    Attributes:
        date (datetime): Gregorian date and time.
        jd (float): Julian Date.
        deltaT (float): Difference TT - UT in seconds.
        hijri (IslamicDateInfo | None): Islamic date information.
    """

    date: datetime
    jd: float
    deltaT: float
    hijri: IslamicDateInfo | None = None

    @property
    def gregorian_date(self) -> str:
        return self.date.strftime("%A, %d %B, %Y")
    
    @property
    def clock(self) -> str:
        return self.date.strftime("%X")
    
    @property
    def timezone(self) -> str | None:
        return self.date.timetz().tzname()
    
    @property
    def utc_offset(self) -> float:
        tz_delta = self.date.utcoffset()
        return 0.0 if tz_delta is None else tz_delta.total_seconds() / 3600
    
    @property
    def jde(self) -> float:
        return self.jd + self.deltaT / 86400
    
    def format_utc_offset(self) -> str:
        '''String formatting for UTC Offsets.'''

        hours = int(self.utc_offset)
        minutes = abs(int((self.utc_offset - hours) * 60))
        return f"UTC{hours:+03d}:{minutes:02d}"

    def __str__(self):
        hijri_date = self.hijri.full_date(self.date.strftime("%A")) if self.hijri is not None else "N/A"
        return ("Time & Date\n"
                f"\tGregorian Date:\t\t{self.gregorian_date}\n"
                f"\tIslamic Date:\t\t{hijri_date}\n"
                f"\t24h-Time:\t\t{self.clock}\n"
                f"\tTime Zone:\t\t{self.timezone} {self.format_utc_offset()}\n"
                f"\tJulian Day:\t\t{self.jd}\n"
                f"\tEstimated ΔT:\t\t{self.deltaT:.2f} s")

    def to_dict(self) -> dict[str, Any]:
        """Serialize date-time info with stable machine-readable schema."""
        return {
            "schema_version": PUBLIC_SCHEMA_VERSION,
            "type": "DateTimeInfo",
            "date_iso": self.date.isoformat(),
            "gregorian_date": self.gregorian_date,
            "clock_24h": self.clock,
            "timezone": self.timezone,
            "utc_offset_hours": self.utc_offset,
            "utc_offset_display": self.format_utc_offset(),
            "julian_day": self.jd,
            "delta_t_seconds": self.deltaT,
            "julian_ephemeris_day": self.jde,
            "hijri": None if self.hijri is None else self.hijri.to_dict(),
        }

@dataclass(frozen=True, slots=True)
class PrayerMethod:
    """
    Describes a prayer time calculation method.

    Attributes:
        name (str): Method name.
        fajr_angle (Angle): Fajr twilight angle.
        isha_angle (Angle): Isha twilight angle.
        keys (Tuple [str, ...] | None): Optional method identifiers.
        maghrib_angle (Angle): Maghrib angle (if applicable).
        asr_type (int): Asr shadow factor type.
        midnight_type (int): Midnight definition type.
        extreme_lats (str): The method on how to deal with prayer times at extreme latitudes. Default is Nearest Latitude method.
    """
    name: str
    fajr_angle: Angle
    isha_angle: Angle
    keys: Tuple[str, ...] | None = None
    maghrib_angle: Angle = Angle(0)
    asr_type: int = 0
    midnight_type: int = 0
    extreme_lats: str = 'ANGLEBASED'

    def to_dict(self) -> dict[str, Any]:
        """Serialize prayer-method configuration with explicit units."""
        return {
            "schema_version": PUBLIC_SCHEMA_VERSION,
            "type": "PrayerMethod",
            "name": self.name,
            "keys": list(self.keys) if self.keys is not None else None,
            "fajr_angle_degrees": self.fajr_angle.decimal,
            "maghrib_angle_degrees": self.maghrib_angle.decimal,
            "isha_angle_degrees": self.isha_angle.decimal,
            "asr_type": self.asr_type,
            "midnight_type": self.midnight_type,
            "extreme_lats": self.extreme_lats,
        }

@dataclass(frozen=True, slots=True)
class Prayer:
    """
    Represents a single prayer time.

    Attributes:
        name (str): Name of the prayer.
        time (datetime | str | float): Time of the prayer or status.
        method (PrayerMethod): The method used for calculation.
    """
    name: str
    time: datetime | str | float
    method: PrayerMethod

    @property
    def time_str(self) -> str:
        if self.time == math.inf:
            return "Does not exist."
        elif isinstance(self.time, datetime):
            return self.time.strftime("%X %d-%m-%Y")
        elif isinstance(self.time, str):
            return self.time
        else:
            return str(self.time)

    def to_dict(self) -> dict[str, Any]:
        """Serialize a prayer event and its underlying method context."""
        return {
            "schema_version": PUBLIC_SCHEMA_VERSION,
            "type": "Prayer",
            "name": self.name,
            "time": _serialize_time_value(self.time),
            "time_display": self.time_str,
            "method_name": self.method.name,
        }

@dataclass(frozen=True, slots=True)
class PrayerTimes:
    """
    Holds all calculated prayer times for a day.

    Attributes:
        method (PrayerMethod): Prayer calculation method.
        fajr (Prayer): Fajr prayer.
        sunrise (Prayer): Sunrise.
        zuhr (Prayer): Ẓuhr prayer.
        asr (Prayer): ʿAṣr prayer.
        sunset (Prayer): Sunset.
        maghrib (Prayer): Maghrib prayer.
        isha (Prayer): ʿishāʾ prayer.
        midnight (Prayer): Midnight marker.
        extreme_latitude_applied (bool): Whether an extreme-latitude fallback was applied.
        extreme_latitude_rule (str | None): Rule used when fallback was applied.
        extreme_latitude_reason (str | None): Human-readable reason for fallback activation.
    """

    method: PrayerMethod
    fajr: Prayer
    sunrise: Prayer
    zuhr: Prayer
    asr: Prayer
    sunset: Prayer
    maghrib: Prayer
    isha: Prayer
    midnight: Prayer
    extreme_latitude_applied: bool = False
    extreme_latitude_rule: str | None = None
    extreme_latitude_reason: str | None = None

    def __str__(self):
        extreme_line = ""
        if self.extreme_latitude_applied:
            reason = self.extreme_latitude_reason or "Extreme latitude fallback applied."
            extreme_line = f"\tExtreme Latitude:\t{self.extreme_latitude_rule} ({reason})\n"
        return ("Prayer Times at Observer Timezone\n"
                f"\tMethod:\t\t\t{self.method.name}\n"
                f"{extreme_line}"
                f"\t{self.fajr.name}:\t\t\t{self.fajr.time_str}\n"
                f"\t{self.sunrise.name}:\t\t{self.sunrise.time_str}\n"
                f"\t{self.zuhr.name}:\t\t\t{self.zuhr.time_str}\n"
                f"\t{self.asr.name}:\t\t\t{self.asr.time_str}\n"
                f"\t{self.sunset.name}:\t\t\t{self.sunset.time_str}\n"
                f"\t{self.maghrib.name}:\t\t{self.maghrib.time_str}\n"
                f"\t{self.isha.name}:\t\t\t{self.isha.time_str}\n"
                f"\t{self.midnight.name}:\t\t{self.midnight.time_str}")

    def to_dict(self) -> dict[str, Any]:
        """Serialize full prayer-time output as a stable API payload."""
        return {
            "schema_version": PUBLIC_SCHEMA_VERSION,
            "type": "PrayerTimes",
            "method": self.method.to_dict(),
            "prayers": {
                "fajr": self.fajr.to_dict(),
                "sunrise": self.sunrise.to_dict(),
                "zuhr": self.zuhr.to_dict(),
                "asr": self.asr.to_dict(),
                "sunset": self.sunset.to_dict(),
                "maghrib": self.maghrib.to_dict(),
                "isha": self.isha.to_dict(),
                "midnight": self.midnight.to_dict(),
            },
            "extreme_latitude": {
                "applied": self.extreme_latitude_applied,
                "rule": self.extreme_latitude_rule,
                "reason": self.extreme_latitude_reason,
            },
        }

@dataclass(frozen=True, slots=True)
class MeccaInfo:
    """
    Represents the Qibla direction and distance towards Mecca.

    Attributes:
        distance (Distance): Distance to Mecca.
        angle (Angle): Azimuth angle towards Mecca.
        cardinal (str): Cardinal direction.
    """

    distance: Distance
    angle: Angle
    cardinal: str

    def __str__(self):
        return ("Mecca\n"
                f"\tDistance:\t\t{self.distance.value:,.0f} {self.distance.unit}\t\t({self.distance.in_unit(DistanceUnits.MILE):,.0f} mi)\n"
                f"\tDirection:\t\t{self.cardinal}\t\t\t({self.angle})")

    def to_dict(self) -> dict[str, Any]:
        """Serialize Qibla distance/direction output."""
        return {
            "schema_version": PUBLIC_SCHEMA_VERSION,
            "type": "MeccaInfo",
            "distance": self.distance.to_dict(),
            "angle": self.angle.to_dict(),
            "cardinal": self.cardinal,
        }

@dataclass(frozen=True, slots=True)
class SunInfo:
    """
    Contains Sun position and event information.

    Attributes:
        sunrise (datetime): Sunrise time.
        sun_transit (datetime): Solar transit (culmination) time.
        sunset (datetime): Sunset time.
        apparent_altitude (Angle): Apparent altitude.
        true_azimuth (Angle): True azimuth.
        geocentric_distance (Distance): Geocentric distance from the Sun to the centre of the Earth.
        apparent_declination (Angle): Apparent declination.
        apparent_right_ascension (RightAscension): Apparent right ascension.
        greenwich_hour_angle (Angle): Greenwich hour angle.
        local_hour_angle (Angle): Local hour angle.
    """

    sunrise: datetime | str
    sun_transit: datetime
    sunset: datetime | str
    apparent_altitude: Angle
    true_azimuth: Angle
    geocentric_distance: Distance
    apparent_declination: Angle
    apparent_right_ascension: RightAscension
    greenwich_hour_angle: Angle
    local_hour_angle: Angle

    def sun_time_str(self, time: datetime | str) -> str:
        if isinstance(time, datetime):
            return time.strftime("%X %d-%m-%Y")
        elif isinstance(time, str):
            return time
        else:
            return str(time)

    def __str__(self):
        return ("The Sun\n"
                f"\tSunrise:\t\t{self.sun_time_str(self.sunrise)}\n"
                f"\tSun Transit:\t\t{self.sun_transit.strftime('%X %d-%m-%Y')}\n"
                f"\tSunset:\t\t\t{self.sun_time_str(self.sunset)}\n"
                f"\tApp. Altitude:\t\t{self.apparent_altitude}\n"
                f"\tApp. Azimuth:\t\t{self.true_azimuth}\n"
                f"\tDistance:\t\t{self.geocentric_distance}\n"
                f"\tApp. Declination:\t{self.apparent_declination}\n"
                f"\tApp. Right Ascension:\t{self.apparent_right_ascension}\n"
                f"\tGreenwich Hour Angle:\t{self.greenwich_hour_angle}\n"
                f"\tLocal Hour Angle:\t{self.local_hour_angle}")

    def to_dict(self) -> dict[str, Any]:
        """Serialize solar event and coordinate outputs."""
        return {
            "schema_version": PUBLIC_SCHEMA_VERSION,
            "type": "SunInfo",
            "sunrise": _serialize_time_value(self.sunrise),
            "sun_transit": _serialize_time_value(self.sun_transit),
            "sunset": _serialize_time_value(self.sunset),
            "apparent_altitude": self.apparent_altitude.to_dict(),
            "true_azimuth": self.true_azimuth.to_dict(),
            "geocentric_distance": self.geocentric_distance.to_dict(),
            "apparent_declination": self.apparent_declination.to_dict(),
            "apparent_right_ascension": self.apparent_right_ascension.to_dict(),
            "greenwich_hour_angle": self.greenwich_hour_angle.to_dict(),
            "local_hour_angle": self.local_hour_angle.to_dict(),
        }

@dataclass(frozen=True, slots=True)
class MoonInfo:
    """
    Contains Moon position, transit, and illumination information.

    Attributes:
        moonrise (datetime): Moonrise time.
        moon_transit (datetime): Moon's upper transit time.
        moonset (datetime): Moonset time.
        illumination (float): Illuminated fraction (0 to 1).
        apparent_altitude (Angle): Apparent altitude.
        true_azimuth (Angle): True azimuth.
        geocentric_distance (Distance): Geocentric distance from the Moon to the centre of the Earth.
        parallax (Angle): Lunar parallax.
        topocentric_declination (Angle): Topocentric declination.
        topocentric_right_ascension (RightAscension): Topocentric right ascension.
        greenwich_hour_angle (Angle): Greenwich hour angle.
        local_hour_angle (Angle): Local hour angle.
    """

    moonrise: datetime | str
    moon_transit: datetime
    moonset: datetime | str
    illumination: float
    apparent_altitude: Angle
    true_azimuth: Angle
    geocentric_distance: Distance
    parallax: Angle
    topocentric_declination: Angle
    topocentric_right_ascension: RightAscension
    greenwich_hour_angle: Angle
    local_hour_angle: Angle

    def moon_time_str(self, time: datetime | str) -> str:
        if isinstance(time, datetime):
            return time.strftime("%X %d-%m-%Y")
        elif isinstance(time, str):
            return time
        else:
            return str(time)

    def __str__(self):
        return ("The Moon\n"
                f"\tMoonrise:\t\t{self.moon_time_str(self.moonrise)}\n"
                f"\tMoon Transit:\t\t{self.moon_transit.strftime('%X %d-%m-%Y')}\n"
                f"\tMoonset:\t\t{self.moon_time_str(self.moonset)}\n"
                f"\tIllumination:\t\t{self.illumination * 100:.2f}%\n"
                f"\tApp. Altitude:\t\t{self.apparent_altitude}\n"
                f"\tAzimuth:\t\t{self.true_azimuth}\n"
                f"\tDistance:\t\t{self.geocentric_distance.value:,.0f} {self.geocentric_distance.unit}\t\t({self.geocentric_distance.in_unit(DistanceUnits.MILE):,.0f} mi)\n"
                f"\tParallax:\t\t{self.parallax}\n"
                f"\tTop. Declination:\t{self.topocentric_declination}\n"
                f"\tTop. Right Ascension:\t{self.topocentric_right_ascension}\n"
                f"\tGreenwich Hour Angle:\t{self.greenwich_hour_angle}\n"
                f"\tLocal Hour Angle:\t{self.local_hour_angle}")

    def to_dict(self) -> dict[str, Any]:
        """Serialize lunar event and coordinate outputs."""
        return {
            "schema_version": PUBLIC_SCHEMA_VERSION,
            "type": "MoonInfo",
            "moonrise": _serialize_time_value(self.moonrise),
            "moon_transit": _serialize_time_value(self.moon_transit),
            "moonset": _serialize_time_value(self.moonset),
            "illumination_fraction": self.illumination,
            "apparent_altitude": self.apparent_altitude.to_dict(),
            "true_azimuth": self.true_azimuth.to_dict(),
            "geocentric_distance": self.geocentric_distance.to_dict(),
            "parallax": self.parallax.to_dict(),
            "topocentric_declination": self.topocentric_declination.to_dict(),
            "topocentric_right_ascension": self.topocentric_right_ascension.to_dict(),
            "greenwich_hour_angle": self.greenwich_hour_angle.to_dict(),
            "local_hour_angle": self.local_hour_angle.to_dict(),
        }

@dataclass(frozen=True, slots=True)
class Visibilities:
    """
    Contains visibility information for the new moon crescent.

    Attributes:
        criterion (str): Visibility criterion used.
        dates (Tuple[datetime]): Dates of visibility calculations.
        q_values (Tuple[float]): Q values for each date.
        classifications (Tuple[str]): Classifications for each date.
    """
    criterion: str
    dates: Tuple[datetime, ...]
    q_values: Tuple[float, ...]
    classifications: Tuple[str, ...]

    def __str__(self):
        base: str = f"Visibility of New Moon Crescent:\n\tCriterion:\t\t{self.criterion}\n"
        for i, q in enumerate(self.q_values):
            if q == 0:
                formatted_q = "+0.000"
            else:
                int_digits = int(math.log10(abs(q))) + 1
                decimal_digits = max(4 - int_digits, 0)
                formatted_q = f"{q:+.{decimal_digits}f}"

            base += f"\t{self.dates[i].strftime('%X %d-%m-%Y')}:\t{formatted_q}\t{self.classifications[i]}\n"
        return base

    def to_dict(self) -> dict[str, Any]:
        """Serialize visibility sequence with stable ordered entries."""
        entries = [
            {
                "date_iso": date.isoformat(),
                "q_value": q,
                "classification": classification,
            }
            for date, q, classification in zip(self.dates, self.q_values, self.classifications)
        ]
        return {
            "schema_version": PUBLIC_SCHEMA_VERSION,
            "type": "Visibilities",
            "criterion": self.criterion,
            "entries": entries,
        }


@dataclass(frozen=True, slots=True)
class BatchVisibilities:
    """Batch visibility output for multiple coordinates and days.

    Attributes:
        criterion (str): Criterion name used for the run.
        output (str): Output mode: ``raw``, ``classification``, or ``code``.
        date (datetime): Input reference datetime.
        days (int): Number of days generated per location.
        latitudes (Tuple[float, ...]): Input latitude sequence in order.
        longitudes (Tuple[float, ...]): Input longitude sequence in order.
        values (Tuple[Tuple[float | str | int, ...], ...]): Matrix indexed as
            ``[location_index][day_index]``.
    """

    criterion: str
    output: str
    date: datetime
    days: int
    latitudes: Tuple[float, ...]
    longitudes: Tuple[float, ...]
    values: Tuple[Tuple[float | str | int, ...], ...]

    @property
    def shape(self) -> tuple[int, int]:
        """Return matrix shape as (locations, days)."""
        location_count = len(self.values)
        if location_count == 0:
            return (0, 0)
        return (location_count, len(self.values[0]))

    def to_dict(self) -> dict[str, Any]:
        """Serialize batch output with stable schema metadata."""
        return {
            "schema_version": PUBLIC_SCHEMA_VERSION,
            "type": "BatchVisibilities",
            "criterion": self.criterion,
            "output": self.output,
            "date_iso": self.date.isoformat(),
            "days": self.days,
            "shape": list(self.shape),
            "latitudes": list(self.latitudes),
            "longitudes": list(self.longitudes),
            "values": [list(row) for row in self.values],
        }
