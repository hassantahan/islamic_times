import numpy as np
import datetime
from datetime import timedelta
from sun_equations import *
from moon_equations import *
from time_equations import *

##### Definitions #####
FAJR_ANGLE = 16
MAGHRIB_ANGLE = 4
#ISHA_ANGLE = 
MECCA_LAT = 21.420164986 
MECCA_LONG = 39.822330044


##### Inputs #####
today = datetime.datetime.now() #datetime.datetime(2023, 3, 17, 12, 0 ,0)
latitude = 43.74533  #45.508888 #39.0 
longitude = -79.30945  #-73.561668 #-76.8 

##### Functions #####
def find_tomorrow_fajr(utc_change, long, eq_of_time_minutes):
    tomorrow = today + timedelta(days=1)
    jd_tomorrow = gregorian_to_jd(tomorrow.year, tomorrow.month, tomorrow.day)

    tomorrow_sun_declination = sunpos(jd_tomorrow, latitude, longitude)[11]
    tomorrow_fajr_solar_angle = solar_hour_angle(tomorrow_sun_declination, FAJR_ANGLE)
    tomorrow_solar_fajr = sunrise_sunset(-1, tomorrow_fajr_solar_angle)
    tomorrow_standard_fajr = solar2standard(tomorrow_solar_fajr, utc_change, long, eq_of_time_minutes)

    return tomorrow_standard_fajr

def haversine(lat1, lon1, lat2 = MECCA_LAT, lon2 = MECCA_LONG):
    # Convert latitude and longitude to radians
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])

    # Haversine formula
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = math.sin(dlat / 2) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2) ** 2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    distance = 6371 * c  # Earth radius in km

    # Calculate course angle
    y = math.sin(dlon) * math.cos(lat2)
    x = math.cos(lat1) * math.sin(lat2) - math.sin(lat1) * math.cos(lat2) * math.cos(dlon)
    course_angle = math.degrees(math.atan2(y, x))

    return distance, course_angle

def get_cardinal_direction(degree):
    cardinals = [
        'N', 'NNE', 'NE', 'ENE',
        'E', 'ESE', 'SE', 'SSE',
        'S', 'SSW', 'SW', 'WSW',
        'W', 'WNW', 'NW', 'NNW'
    ]
    idx = int(((degree + 11.25) % 360) / 22.5)
    return cardinals[idx]


##### Calculations #####
### Find UTC Offset According to Lat/Long
utc_diff = -1 * find_utc_offset(latitude, longitude, today)


### Calculate Julian Date
jd = gregorian_to_jd(today.year, today.month, today.day + fraction_of_day(today), -1 * utc_diff)


### Sun & Moon Properties Calculations
# Get factors arrays
sun_factors = sunpos(jd, latitude, longitude)
moon_factors = moonpos(jd, latitude, longitude)

# Important Sun Factors placed into local variables
sun_declination = sun_factors[11]
alpha_hour = sun_factors[10] // 15
alpha_minute = int((sun_factors[10] / 15 - alpha_hour) * 60)
alpha_second = (sun_factors[10] / 15 * 60 - alpha_hour * 60 - alpha_minute) * 60
solar_angle = solar_hour_angle(latitude, sun_declination)

# Important Moon Factors placed into local variables



### Equation of Time
eq_of_time = equation_of_time(jd, latitude, longitude)


### Prayer Time Calculations
# Calculate prayer times in solar time
solar_fajr = sunrise_sunset(-1, solar_hour_angle(latitude, sun_declination, FAJR_ANGLE))
solar_sunrise = sunrise_sunset(-1, solar_angle)
solar_sunset = sunrise_sunset(1, solar_angle)
solar_maghrib = sunrise_sunset(1, solar_hour_angle(latitude, sun_declination, MAGHRIB_ANGLE))

# Convert prayer times from solar to standard time
standard_fajr = solar2standard(solar_fajr, utc_diff, longitude, eq_of_time)
standard_sunset = solar2standard(solar_sunset, utc_diff, longitude, eq_of_time)
standard_noon = solar2standard(12.0, utc_diff, longitude, eq_of_time)
#standard_asr = solar2standard()
standard_sunrise = solar2standard(solar_sunrise, utc_diff, longitude, eq_of_time)
standard_maghrib = solar2standard(solar_maghrib, utc_diff, longitude, eq_of_time)
#standard_isha = solar2standard(solar_isha, utc_diff, longitude, eq_of_time)
standard_midnight = time_midpoint(standard_maghrib, find_tomorrow_fajr(utc_diff, longitude, eq_of_time))

# TODO: Find Next New Moon (and the rest of the phases)

# TODO: Find Day of Clear Visibility of Next New Moon

# TODO: Calculate Moon Illumination

# TODO: Calculate Current Islamic Date


### Misc
# Distance from Given Coordinates to Mecca + direction
mecca_distance, mecca_direction = haversine(latitude, longitude)


##### Outputs #####
# Date & Time
#TODO: Islamic Calendar
print("Time & Date\n\tGregorian Date:\t\t{}".format(today.strftime("%A, %d %B, %Y")))
print("\t24h-Time:\t\t{}\n\tTime Zone:\t\tUTC{}".format(today.strftime("%X"),float_to_24time(utc_diff * -1)))

# Prayer Times
print("Prayer Times\n\tFajr:\t\t\t{}".format(float_to_24time(standard_fajr)))
print("\tSunrise:\t\t{}".format(float_to_24time(standard_sunrise)))
print("\tẒuhr:\t\t\t {}".format(float_to_24time(standard_noon)))
#print("\t`Asr:\t\t\t {}".format(float_to_24time(standard_asr)))
print("\tSunset: \t\t{}".format(float_to_24time(standard_sunset)))
print("\tMaghrib: \t\t{}".format(float_to_24time(standard_maghrib)))
#print("\t`Isha: \t\t{}".format(float_to_24time(standard_isha)))
print("\tMidnight: \t\t{}".format(float_to_24time(standard_midnight)))

# Mecca
print("Mecca\n\tDistance: \t\t{:.2f} km".format(mecca_distance))
print("\tDirection: \t\t{} ({:.2f}°)".format(get_cardinal_direction(np.round(mecca_direction)), mecca_direction))

# The Sun
print("The Sun\n\tApparent Declination:\t{:.3f}°".format(sun_declination))
print("\tRight Ascenscion\t{}h {}m {:.2f}s".format(alpha_hour, alpha_minute, alpha_second))
print("\tSolar Hour Angle:\t{:.2f}°".format(solar_angle))
print("\tEquation of time:\t{:.2f} mins".format(equation_of_time(jd, latitude, longitude)))

# The Moon pg 347- 349
print("The Moon\n\tDeclination:\t{:.3f}°".format(sun_declination))