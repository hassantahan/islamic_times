import math
from typing import Tuple, ClassVar, List
from datetime import datetime
from dataclasses import dataclass

# === Helper Classes ===

@dataclass(frozen=True, slots=True)
class Angle:
    decimal: float

    @property
    def dms(self) -> Tuple[int, int, float]:
        degrees = int(self.decimal)
        abs_decimal = abs(self.decimal)
        minutes_full = (abs_decimal - abs(degrees)) * 60
        minutes = int(minutes_full)
        seconds = (minutes_full - minutes) * 60
        return (degrees, minutes, seconds)

    @property
    def radians(self) -> float:
        return math.radians(self.decimal)

    def dms_str(self) -> str:
        d, m, s = self.dms
        return f"{d:+04}\u00b0 {m:02}\u2032 {s:05.2f}\u2033"

    def __str__(self):
        return f"{self.decimal:+08.3f}\u00b0\t\t({self.dms_str()})"

@dataclass(frozen=True, slots=True)
class RightAscension:
    decimal_hours: float

    @property
    def hms(self) -> Tuple[int, int, float]:
        hours = int(self.decimal_hours)
        minutes_full = (self.decimal_hours - hours) * 60
        minutes = int(minutes_full)
        seconds = (minutes_full - minutes) * 60
        return (hours, minutes, seconds)

    def hms_str(self) -> str:
        h, m, s = self.hms
        return f"{h:02}h {m:02}m {s:05.2f}s"

    @property
    def decimal_degrees(self) -> Angle:
        return Angle(self.decimal_hours * 15)
    
    @property
    def radians(self) -> float:
        return self.decimal_degrees.radians

    @property
    def dms(self) -> Tuple[int, int, float]:
        deg = self.decimal_degrees
        degrees = int(deg)
        decimal_minutes = (abs(deg) - degrees) * 60
        minutes = int(decimal_minutes)
        seconds = (decimal_minutes - minutes) * 60
        return (degrees, minutes, seconds)

    def dms_str(self) -> str:
        d, m, s = self.dms
        return f"{d:+04}\u00b0 {m:02}\u2032 {s:05.2f}\u2033"

    def __str__(self):
        return f"{self.hms_str()}\t\t({self.decimal_degrees})"

@dataclass(frozen=True, slots=True)
class DistanceUnit:
    name: str
    symbol: str
    to_m_factor: float

    def convert_to(self, value: float, target_unit: 'DistanceUnit') -> float:
        value_in_m = value * self.to_m_factor
        return value_in_m / target_unit.to_m_factor

    def __str__(self):
        return self.symbol

class DistanceUnits:
    METRE = DistanceUnit("meter", "m", 1)
    KILOMETRE = DistanceUnit("kilometer", "km", 1_000)
    AU = DistanceUnit("astronomical unit", "AU", 149_597_870_691)
    NAUTICAL_MILE = DistanceUnit("nautical mile", "nmi", 1_852)
    MILE = DistanceUnit("mile", "mi", 1_609.3445)
    FOOT = DistanceUnit("foot", "ft", 0.3048)

@dataclass(frozen=True, slots=True)
class Distance:
    value: float
    unit: DistanceUnit = DistanceUnits.METRE

    def in_unit(self, target_unit: DistanceUnit) -> float:
        return self.unit.convert_to(self.value, target_unit)

    def to(self, target_unit: DistanceUnit) -> 'Distance':
        return Distance(self.in_unit(target_unit), target_unit)

    def __str__(self):
        if self.unit == DistanceUnits.METRE:
            return f"{self.value:.2f} {self.unit}"
        else:
            m_equiv = self.in_unit(DistanceUnits.METRE)
            return f"{self.value:.2f} {self.unit} ({m_equiv:.2f} m)"


# === Data Classes ===

@dataclass(frozen=True, slots=True)
class ObserverInfo:
    latitude: Angle
    longitude: Angle
    elevation: Distance
    pressure: float = 10
    temperature: float = 101.325

    def __str__(self):
        return ("Observer Parameters\n"
                f"\tLatitude:\t\t{self.latitude.decimal:09.6f}\u00b0\n"
                f"\tLongitude:\t\t{self.longitude.decimal:09.6f}\u00b0\n"
                f"\tElevation:\t\t{self.elevation}\n"
                f"\tPressure:\t\t{self.pressure} kPa\n"
                f"\tTemperature:\t\t{self.temperature:+06.2f}°C")

@dataclass(frozen=True, slots=True)
class IslamicDateInfo:
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

@dataclass(frozen=True, slots=True)
class DateTimeInfo:
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
    def timezone(self) -> str:
        return self.date.timetz().tzname()
    
    @property
    def utc_offset(self) -> float:
        return self.date.utcoffset().total_seconds() / 3600 * -1
    
    @property
    def jde(self) -> float:
        return self.jd + self.deltaT / 86400
    
    def format_utc_offset(self) -> str:
        '''String formatting for UTC Offsets.'''

        hours = int(self.utc_offset)
        minutes = abs(int((self.utc_offset - hours) * 60))
        return f"UTC{hours:+03d}:{minutes:02d}"

    def __str__(self):
        return ("Time & Date\n"
                f"\tGregorian Date:\t\t{self.gregorian_date}\n"
                f"\tIslamic Date:\t\t{self.hijri.full_date(self.date.strftime("%A"))}\n"
                f"\t24h-Time:\t\t{self.clock}\n"
                f"\tTime Zone:\t\t{self.timezone} {self.format_utc_offset()}\n"
                f"\tLocal JD:\t\t{self.jd}\n"
                f"\tEstimated ΔT:\t\t{self.deltaT:.2f} s")

@dataclass(frozen=True, slots=True)
class PrayerMethod:
    name: str
    fajr_angle: Angle
    isha_angle: Angle
    keys: Tuple[str, ...] | None = None
    maghrib_angle: Angle = Angle(0)
    asr_type: int = 0
    midnight_type: int = 0

@dataclass(frozen=True, slots=True)
class Prayer:
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

@dataclass(frozen=True, slots=True)
class PrayerTimes:
    method: PrayerMethod
    fajr: Prayer
    sunrise: Prayer
    zuhr: Prayer
    asr: Prayer
    sunset: Prayer
    maghrib: Prayer
    isha: Prayer
    midnight: Prayer

    def __str__(self):
        return ("Prayer Times at Observer Timezone\n"
                f"\tMethod:\t\t\t{self.method.name}\n"
                f"\t{self.fajr.name}:\t\t\t{self.fajr.time_str}\n"
                f"\t{self.sunrise.name}:\t\t{self.sunrise.time_str}\n"
                f"\t{self.zuhr.name}:\t\t\t{self.zuhr.time_str}\n"
                f"\t{self.asr.name}:\t\t\t{self.asr.time_str}\n"
                f"\t{self.sunset.name}:\t\t\t{self.sunset.time_str}\n"
                f"\t{self.maghrib.name}:\t\t{self.maghrib.time_str}\n"
                f"\t{self.isha.name}:\t\t\t{self.isha.time_str}\n"
                f"\t{self.midnight.name}:\t\t{self.midnight.time_str}")

@dataclass(frozen=True, slots=True)
class MeccaInfo:
    distance: Distance
    angle: Angle
    cardinal: str

    def __str__(self):
        return ("Mecca\n"
                f"\tDistance:\t\t{self.distance.value:.0f} {self.distance.unit}\n"
                f"\tDirection:\t\t{self.cardinal}\t\t\t({self.angle})")

@dataclass(frozen=True, slots=True)
class SunInfo:
    sunrise: datetime
    sun_transit: datetime
    sunset: datetime
    apparent_declination: Angle
    apparent_right_ascension: RightAscension
    apparent_altitude: Angle
    true_azimuth: Angle

    def __str__(self):
        return ("The Sun\n"
                f"\tSunrise:\t\t{self.sunrise.strftime("%X %d-%m-%Y")}\n"
                f"\tSun Transit:\t\t{self.sun_transit.strftime("%X %d-%m-%Y")}\n"
                f"\tSunset:\t\t\t{self.sunset.strftime("%X %d-%m-%Y")}\n"
                f"\tApp. Declination:\t{self.apparent_declination}\n"
                f"\tApp. Right Ascension:\t{self.apparent_right_ascension}\n"
                f"\tApp. Altitude:\t\t{self.apparent_altitude}\n"
                f"\tApp. Azimuth:\t\t{self.true_azimuth}")

@dataclass(frozen=True, slots=True)
class MoonInfo:
    moonrise: datetime
    moon_transit: datetime
    moonset: datetime
    topocentric_declination: Angle
    topocentric_right_ascension: RightAscension
    apparent_altitude: Angle
    true_azimuth: Angle
    parallax: Angle
    illumination: float

    def __str__(self):
        return ("The Moon\n"
                f"\tMoonrise:\t\t{self.moonrise.strftime("%X %d-%m-%Y")}\n"
                f"\tMoon Transit:\t\t{self.moon_transit.strftime("%X %d-%m-%Y")}\n"
                f"\tMoonset:\t\t{self.moonset.strftime("%X %d-%m-%Y")}\n"
                f"\tTop. Declination:\t{self.topocentric_declination}\n"
                f"\tTop. Right Ascension:\t{self.topocentric_right_ascension}\n"
                f"\tApp. Altitude:\t\t{self.apparent_altitude}\n"
                f"\tAzimuth:\t\t{self.true_azimuth}\n"
                f"\tParallax:\t\t{self.parallax}\n"
                f"\tIllumination:\t\t{self.illumination * 100:.2f}%")