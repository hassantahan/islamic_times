import datetime
import pytz
from timezonefinder import TimezoneFinder
from islamic_times import calculation_equations as ce
import numpy as np

J2000              	= 2451545.0
HIJRI_JD            = 1948439.389675
HIJRI_STANDARD_LNG  = 39.71658
JULIAN_CENTURY    	= 36525.0
JULIAN_MILLENNIUM 	= JULIAN_CENTURY * 10
ASTRONOMICAL_UNIT 	= 149597870.7
TROPICAL_YEAR     	= 365.24219878
EARTH_RADIUS_KM     = 6378.14

def fraction_of_day(date):
    return (date - datetime.datetime.combine(date.date(), datetime.time(0))).total_seconds() / (3600 * 24)

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

    ce.bound_angle_deg(theta_zero)

    return theta_zero

# Look to Jean Meeus' "Astronomical Algorithms"
def gregorian_to_jd(year, month, day, zone = 0):
    y = year
    m = month

    if month <= 2:
        y -= 1
        m += 12
    
    a = int(np.floor(y / 100))
    b = 2 - a + int(np.floor(a / 4))

    jd = int(np.floor(365.25 * (y + 4716))) + int(np.floor(30.6001 * (m + 1))) + day + b - 1524.5 - zone / 24
    return jd

# Look to Jean Meeus' "Astronomical Algorithms"
def jd_to_gregorian(jd):
    jd = jd + 0.5
    z = int(np.floor(jd))
    f = jd - z

    if z < 2299161:
        a = z
    else:
        alpha = int(np.floor((z - 1867216.25) / 36524.25))
        a = z + 1 + alpha - int(np.floor(alpha / 4))

    b = a + 1524
    c = int(np.floor((b - 122.1) / 365.25))
    d = int(np.floor(365.25 * c))
    e = int(np.floor((b - d) / 30.6001))

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
    f_day = day - int(np.floor(day))
    hour = f_day * 24
    minute = (hour - int(np.floor(hour))) * 60
    second = (minute - int(np.floor(minute))) * 60
    microsec = (second - int(np.floor(second))) * 1000

    return datetime.datetime(int(year), int(month), int(day), int(hour), int(minute), int(second), int(microsec))

# Look to Jean Meeus' "Astronomical Algorithms"
def gregorian_to_hijri(year, month, day, zone = 0):
    x = year
    m = month

    if m < 3:
        x -= 1
        m += 12
    
    alpha = int(x / 100)
    beta = 2 - alpha + int(alpha / 4)
 
    day += beta
    if day == 0:
        m -= 1
        if m in [1, 3, 5, 7, 8, 10, 12]:
            day = 31
        elif m == 2:
            day = 28
        else:
            day = 30

    b = int(365.25 * x) + int(30.6001 * (m + 1)) + day + 1722519 + beta
    c = int((b - 122.1) / 365.25)
    d = int(365.25 * c)
    e = int((b - d) / 30.6001)

    if x % 4 == 0:
        w = 1
    else:
        w = 2
    
    n = int(275 * m / 9) - w * int((m + 9) / 12) + day - 30
    A = x - 623
    B = int(A / 4)
    C = A % 4
    C_1 = 365.2501 * C
    C_2 = int(C_1)

    if C_1 - C_2 > 0.5:
        C_2 += 1

    D_prime = 1461 * B + 170 + C_2
    q = int(D_prime / 10631)
    r = D_prime % 10631
    j = int(r / 354)
    k = r % 354
    o = int((11 * j + 14) / 30)
    h = 30 * q + j + 1
    jj = k - o + n - 1

    if jj > 354:
        cl = h % 30
        dl = (11 * cl + 3) % 30
        if dl < 19:
            jj -= 354
            h += 1
        elif dl > 18:
            jj -= 355
            h += 1
    elif jj == 0:
        jj = 355
        h -= 1
    
    s = int((jj - 1) / 29.5)
    if jj != 355:
        i_month = 1 + s
        i_day = int(jj - 29.5 * s)
    else:
        i_month = 12
        i_day = 30
    i_year = h

    return [i_year, i_month, i_day]

def get_islamic_month(month):
    islamic_months = {
        1: 'Muḥarram',
        2: 'Ṣaffar',
        3: 'Rabīʿ al-’Awwal',
        4: 'Rabīʿ al-Thānī',
        5: 'Jumādā al-’Ūlā',
        6: 'Jumādā al-Thāniyah',
        7: 'Rajab',
        8: 'Sha‘bān',
        9: 'Ramaḍān',
        10: 'Shawwāl',
        11: 'Ḏū al-Qa‘dah',
        12: 'Ḏū al-Ḥijjah',
    }
    return islamic_months.get(month, "Invalid month number")

def get_islamic_day(day):
    islamic_days = {
        'Sunday': 'al-’Aḥad',
        'Monday': 'al-Ithnayn',
        'Tuesday': 'al-Thulāthā’',
        'Wednesday': 'al-’Arbi‘ā’',
        'Thursday': 'al-Khamīs',
        'Friday': 'al-Jum‘ah',
        'Saturday': 'al-Sabt',
    }
    return islamic_days.get(day, "Invalid day")

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
    return timezone_str, utc_offset.total_seconds() / 3600

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

def format_utc_offset(utc_offset):
    hours = int(utc_offset)
    minutes = int((utc_offset - hours) * 60)

    # Since minutes are absolute, take the absolute value
    minutes = abs(minutes)

    # Use the built-in `format` function to convert to the required format
    offset_str = "UTC{:+03d}:{:02d}".format(hours, minutes)
    return offset_str