import math
import datetime
import numpy as np
import pytz
from timezonefinder import TimezoneFinder

J2000              	= 2451545.0
JULIAN_CENTURY    	= 36525.0
JULIAN_MILLENNIUM 	= JULIAN_CENTURY * 10
ASTRONOMICAL_UNIT 	= 149597870.7
TROPICAL_YEAR     	= 365.24219878
TO_LAT 				= 43.74533
TO_LONG             = -79.30945

def bound_angle_deg(a):
    return a - 360.0 * (math.floor(a / 360.0))

def bound_angle_rad(a):
    return a - (2 * math.pi) * (math.floor(a / (2 * math.pi)))

def sin(a):
    return math.sin(np.deg2rad(a))

def cos(a):
    return math.cos(np.deg2rad(a))

def tan(a):
    return math.tan(np.deg2rad(a))

def decimal_to_dms(decimal_deg):
    degrees = int(decimal_deg)
    minutes = int((decimal_deg - degrees) * 60)
    seconds = (decimal_deg - degrees - minutes / 60) * 3600

    return [degrees, minutes, seconds]

def decimal_to_hms(decimal_degrees):
    hours = int(decimal_degrees / 15)
    minutes = int((decimal_degrees / 15 - hours) * 60)
    seconds = (decimal_degrees / 15 - hours - minutes / 60) * 3600
    return [hours , round(minutes), seconds]


def fraction_of_day(date):
    return (date - datetime.datetime.combine(date.date(), datetime.datetime.time(date))).total_seconds() / (3600 * 24)

def leap_gregorian(year):
    if year % 4 == 0 and (year % 100 != 0 or year % 400 == 0):
        return True
    else:
        return False

def siderial_time(julian_day):
    t = (julian_day - J2000) / JULIAN_CENTURY
    t2 = t ** 2
    t3 = t ** 3

    theta_zero = 280.46061837 + 360.98564736629 * (julian_day - J2000) + \
                 0.000387933 * t2 - t3 / 38710000

    return theta_zero

# Look to Jean Meeus' "Astronomical Algorithms"
def gregorian_to_jd(year, month, day, zone = 0):
    y = year
    m = month

    if month <= 2:
        y -= 1
        m += 12
    
    a = int(y / 100)
    b = 2 - a + int(a / 4)

    jd = int(365.25 * (year + 4716)) + int(30.6001 * (m + 1)) + day + b - 1524.5 - zone / 24
    return jd

# Look to Jean Meeus' "Astronomical Algorithms"
def jd_to_gregorian(jd):
    jd = jd + 0.5
    z = int(jd)
    f = jd - z

    if z < 2299161:
        a = z
    else:
        alpha = int((z - 1867216.25) / 36524.25)
        a = z + 1 + alpha - int(alpha / 4)

    b = a + 1524
    c = int((b - 122.1) / 365.25)
    d = int(365.25 * c)
    e = int((b - d) / 30.6001)

    day = b - d - int(30.6001 * e) + f
    if e < 14:
        month = e - 1
    else:
        month = e - 13
    if month > 2:
        year = c - 4716
    else:
        year = c - 4715

    # calculate hours, minutes, and seconds
    f_day = day - int(day)
    hour = f_day * 24
    minute = (hour - int(hour)) * 60
    second = (minute - int(minute)) * 60
    microsec = (second - int(second)) * 1000

    return datetime.datetime(int(year), int(month), int(day), int(hour), int(minute), int(second), int(microsec))

# Based off of https://eclipse.gsfc.nasa.gov/LEcat5/deltatpoly.html
def delta_t_approx(year):
    t = (year - 2000) / 100

    if year < -500:
        u = (year - 1820) / 100
        return -20 + 32 * u**2
    elif year < 500:
        return 10583.6 - 1014.41 * t + 33.78311 * t**2 - 5.952053 * t**3 - 0.1798452 * t**4 + 0.022174192 * t**5 + 0.0090316521 * t**6
    elif year < 1600:
        u = (year - 1000) / 100
        return 1574.2 - 556.01 * u + 71.23472 * u**2 + 0.319781 * u**3 - 0.8503463 * u**4 - 0.005050998 * u**5 + 0.0083572073 * u**6
    elif year < 1700:
        t = year - 1600
        return 120 - 0.9808 * t - 0.01532 * t**2 + t**3 / 7129
    elif year < 1800:
        t = year - 1700
        return 8.83 + 0.1603 * t - 0.0059285 * t**2 + 0.00013336 * t**3 - t**4 / 1174000
    elif year < 1860:
        t = year - 1800
        return 13.72 - 0.332447 * t + 0.0068612 * t**2 + 0.0041116 * t**3 - 0.00037436 * t**4 + 0.0000121272 * t**5 - 0.0000001699 * t**6 + 0.000000000875 * t**7
    elif year < 1900:
        t = year - 1860
        return 7.62 + 0.5737 * t - 0.251754 * t**2 + 0.01680668 * t**3 -0.0004473624 * t**4 + t**5 / 233174
    elif year < 1920:
        t = year - 1900
        return -2.79 + 1.494119 * t - 0.0598939 * t**2 + 0.0061966 * t**3 - 0.000197 * t**4
    elif year < 1941:
        t = year - 1920
        return 21.20 + 0.84493*t - 0.076100 * t**2 + 0.0020936 * t**3
    elif year < 1961:
        t = year - 1950
        return 29.07 + 0.407*t - t**2/233 + t**3 / 2547
    elif year < 1986:
        t = year - 1975
        return 45.45 + 1.067*t - t**2/260 - t**3/718
    elif year < 2005:
        t = year - 2000
        return 63.86 + 0.3345 * t - 0.060374 * t**2 + 0.0017275 * t**3 + 0.000651814 * t**4 + 0.00002373599 * t**5
    elif year < 2050:
        return 62.92 + 0.32217 * t + 0.005589 * t**2
    elif year < 2150:
        return -20 + 32 * (year-1820)/100**2
    else:
        u = (year - 1820) / 100
        return -20 + 32 * u**2 - 0.5628 * (2150 - year)

def float_to_24time(_hour_):
    # Convert the float to integer hours and minutes
    _inthour_ = int(_hour_)
    _minutes_ = int((_hour_ - _inthour_) * 60)

    time_24hr = "{:02d}:{:02d}".format(_inthour_, _minutes_)
    return time_24hr

def solar2standard(solar_time, utc_diff, longitude, eq_of_time):
    local_standard_meridian = utc_diff * 15
    error_in_minutes = 4 * (local_standard_meridian + longitude) + eq_of_time
    standard_time = solar_time - error_in_minutes / 60
    return standard_time

def find_utc_offset(lat, long, day):
    # Create a TimezoneFinder object
    tf = TimezoneFinder()

    # Get the timezone name for the coordinates
    timezone_str = tf.certain_timezone_at(lat=lat, lng=long)

    # Get the timezone object for the given name
    timezone = pytz.timezone(timezone_str)

    # Localize the given date with the timezone
    localized_datetime = timezone.localize(datetime.datetime.combine(day, datetime.datetime.min.time()))

    # Get the UTC offset for the given timezone and date
    utc_offset = localized_datetime.utcoffset()

    # Calculate the UTC offset in hours
    return utc_offset.total_seconds() / 3600

def time_midpoint(time1, time2):
    # Convert hours to datetime objects
    time1 = datetime.datetime.strptime(f'{int(time1):02}:{int((time1*60) % 60):02}', "%H:%M")
    time2 = datetime.datetime.strptime(f'{int(time2):02}:{int((time2*60) % 60):02}', "%H:%M")

    # If time2 is earlier than time1, assume it's from the next day
    if time2 < time1:
        time2 += datetime.timedelta(days=1)

    # Calculate middle time
    middle_time = time1 + (time2 - time1) / 2

    # Convert the middle time back to hours as a float
    middle_time_in_hours = middle_time.hour + middle_time.minute / 60.0

    # If the middle time is on the next day, subtract 24 to keep the output in the range [0, 24)
    if middle_time_in_hours >= 24:
        middle_time_in_hours -= 24

    return middle_time_in_hours