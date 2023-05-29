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

def time_midpoint(t1, t2):
    # Calculate the absolute difference between the times
    dt = abs(t2 - t1)
    
    # If the difference is greater than 12 hours, we need to adjust the midpoint
    if dt > 12:
        if t2 > t1:
            t1 += 24
        else:
            t2 += 24
    
    # Calculate the average of the two times
    t_mid = (t1 + t2) / 2
    
    # If the midpoint is greater than 24 hours, subtract 24 to get the correct time
    if t_mid >= 24:
        t_mid -= 24
    
    return t_mid