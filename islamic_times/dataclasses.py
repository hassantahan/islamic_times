from typing import Tuple
from datetime import datetime
from dataclasses import dataclass

# === Helper Classes ===

@dataclass(frozen=True)
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

    def dms_str(self) -> str:
        d, m, s = self.dms
        return f"{d:+04}\u00b0 {m:02}\u2032 {s:05.2f}\u2033"

    def __str__(self):
        return f"{self.decimal:+07.3f}\u00b0 ({self.dms_str()})"


@dataclass(frozen=True)
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
    def decimal_degrees(self) -> float:
        return self.decimal_hours * 15

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
        return f"{self.hms_str()} ({self.dms_str()})"


@dataclass(frozen=True)
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


@dataclass(frozen=True)
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

@dataclass(frozen=True)
class ObserverInfo:
    latitude: float
    longitude: float
    elevation: float
    pressure: float
    temperature: float

    def __str__(self):
        return ("Observer Parameters\n"
                f"\tLatitude:\t\t\t{self.latitude}\n"
                f"\tLongitude:\t\t\t{self.longitude}\n"
                f"\tElevation:\t\t\t{self.elevation}\n"
                f"\tPressure:\t\t\t{self.pressure}\n"
                f"\tTemperature:\t\t{self.temperature}")


@dataclass(frozen=True)
class DateTimeInfo:
    date: datetime
    hijri: str
    jd: float
    eq_of_time: float
    deltaT: float

    @property
    def gregorian_date(self) -> str:
        return self.date.strftime("%A, %d %B, %Y")
    
    @property
    def clock(self) -> str:
        return self.date.strftime("%X")
    
    @property
    def timezone(self) -> str:
        return self.date.timetz()
    
    @property
    def utc_offset(self) -> float:
        return self.date.utcoffset().total_seconds() / 3600
    
    def format_utc_offset(self) -> str:
        '''String formatting for UTC Offsets.'''

        hours = int(self.utc_offset)
        minutes = abs(int((self.utc_offset - hours) * 60))
        return f"UTC{hours:+03d}:{minutes:02d}"

    def __str__(self):
        return ("Time & Date\n"
                f"\tGregorian Date:\t\t{self.gregorian_date}\n"
                f"\tIslamic Date:\t\t{self.hijri}\n"
                f"\t24h-Time:\t\t\t{self.clock}\n"
                f"\tTime Zone:\t\t\t{self.timezone} {self.format_utc_offset()}\n"
                f"\tLocal JD:\t\t\t{self.jd}\n"
                f"\tEquation of Time:\t{self.eq_of_time} minutes\n"
                f"\tEstimated ΔT:\t\t{self.deltaT}s")


@dataclass(frozen=True)
class PrayerTimes:
    method: str
    fajr: str
    sunrise: str
    zuhr: str
    asr: str
    sunset: str
    maghrib: str
    isha: str
    midnight: str

    def __str__(self):
        return ("Prayer Times at Observer Timezone\n"
                f"\tMethod:\t\t\t{self.method}\n"
                f"\tFajr:\t\t\t{self.fajr}\n"
                f"\tSunrise:\t\t\t{self.sunrise}\n"
                f"\tZ\u0323uhr:\t\t\t {self.zuhr}\n"
                f"\tʿAṣr:\t\t\t{self.asr}\n"
                f"\tSunset:\t\t\t{self.sunset}\n"
                f"\tMaghrib:\t\t\t{self.maghrib}\n"
                f"\tʿIsh\u0101ʾ:\t\t\t{self.isha}\n"
                f"\tMidnight:\t\t\t{self.midnight}")


@dataclass(frozen=True)
class MeccaInfo:
    distance: Distance
    angle: Angle
    cardinal: str

    def __str__(self):
        return ("Mecca\n"
                f"\tDistance:\t\t\t{self.distance}\n"
                f"\tDirection:\t\t\t{self.cardinal} ({self.angle.decimal:.2f}\u00b0)")


@dataclass(frozen=True)
class SunInfo:
    sunrise: str
    sun_transit: str
    sunset: str
    declination: Angle
    right_ascension: RightAscension
    altitude: Angle
    azimuth: Angle

    def __str__(self):
        return ("The Sun\n"
                f"\tSunrise:\t\t\t{self.sunrise}\n"
                f"\tSun Transit:\t\t{self.sun_transit}\n"
                f"\tSunset:\t\t\t{self.sunset}\n"
                f"\tApp. Declination:\t{self.declination}\n"
                f"\tApp. Right Ascension:\t{self.right_ascension}\n"
                f"\tAltitude:\t\t\t{self.altitude}\n"
                f"\tAzimuth:\t\t\t{self.azimuth}")


@dataclass(frozen=True)
class MoonInfo:
    moonrise: datetime
    moon_transit: datetime
    moonset: datetime
    declination: Angle
    right_ascension: RightAscension
    altitude: Angle
    azimuth: Angle
    parallax: Angle
    illumination: float

    def __str__(self):
        return ("The Moon\n"
                f"\tMoonrise:\t\t\t{self.moonrise.strftime("%X %d-%m-%Y")}\n"
                f"\tMoon Transit:\t\t{self.moon_transit.strftime("%X %d-%m-%Y")}\n"
                f"\tMoonset:\t\t\t{self.moonset.strftime("%X %d-%m-%Y")}\n"
                f"\tDeclination:\t\t{self.declination}\n"
                f"\tRight Ascension:\t{self.right_ascension}\n"
                f"\tAltitude:\t\t\t{self.altitude}\n"
                f"\tAzimuth:\t\t\t{self.azimuth}\n"
                f"\tParallax:\t\t\t{self.parallax}\n"
                f"\tIllumination:\t\t{self.illumination}%")