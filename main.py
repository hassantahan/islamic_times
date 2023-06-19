from moon_equations import *

##### Definitions #####
FAJR_ANGLE = 16
MAGHRIB_ANGLE = 4
ISHA_ANGLE = 14
MECCA_LAT = 21.420164986 
MECCA_LONG = 39.822330044

##### Inputs #####
today = datetime.datetime.now() #datetime.datetime(2023, 3, 17, 12, 0 ,0)
latitude = 43.74533 #27.7172 #  #45.508888 #39.0 
longitude = -79.30945 #85.3240 #  #-73.561668 #-76.8 

##### Functions #####
def find_tomorrow_fajr(jd, utc_change, long, eq_of_time_minutes):
    jd_tomorrow = jd + 1

    tomorrow_sun_declination = sunpos(jd_tomorrow, latitude, longitude)[11]
    tomorrow_fajr_solar_angle = solar_hour_angle(latitude, tomorrow_sun_declination, FAJR_ANGLE)
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

def asr_time(lat, dec, t = 1):
    temp_num = sin(np.rad2deg(np.arctan2(1, t + tan(lat - dec)))) - sin(lat) * sin(dec)
    temp_denom = cos(lat) * cos(dec)
    return 1 / 15 * np.rad2deg(np.arccos(temp_num / temp_denom))

##### Calculations #####
### Find UTC Offset According to Lat/Long
tz_name, utc_diff = find_utc_offset(latitude, longitude, today)
utc_diff *= -1

# Fix time for timezone
today = today + datetime.timedelta(hours=-utc_diff + 4)


### Calculate Julian Date
jd = gregorian_to_jd(today.year, today.month, today.day + fraction_of_day(today), -1 * utc_diff)


### Sun & Moon Properties Calculations
# Get factors arrays
sun_factors = sunpos(jd, latitude, longitude)
delPsi, delEps = sun_nutation(jd)
moon_factors = moonpos(jd, latitude, longitude, delPsi, sun_factors[13], TO_ELEV)

# Important Sun Factors placed into local variables
sun_declination = sun_factors[11]
alpha_hour = sun_factors[10] // 15
alpha_minute = int((sun_factors[10] / 15 - alpha_hour) * 60)
alpha_second = (sun_factors[10] / 15 * 60 - alpha_hour * 60 - alpha_minute) * 60
solar_angle = solar_hour_angle(latitude, sun_declination)
sun_alt = sun_factors[15]
sun_az = sun_factors[16]

# Important Moon Factors placed into local variables
moon_declination = moon_factors[7]
moon_alpha = decimal_to_hms(moon_factors[6])
moon_alt = moon_factors[9]
moon_az = moon_factors[10]


### Equation of Time
eq_of_time = equation_of_time(jd, latitude, longitude)


### Prayer Time Calculations
# Calculate prayer times in solar time
solar_fajr = sunrise_sunset(-1, solar_hour_angle(latitude, sun_declination, FAJR_ANGLE))
solar_sunrise = sunrise_sunset(-1, solar_angle)
solar_sunset = sunrise_sunset(1, solar_angle)
solar_maghrib = sunrise_sunset(1, solar_hour_angle(latitude, sun_declination, MAGHRIB_ANGLE))
solar_isha = sunrise_sunset(1, solar_hour_angle(latitude, sun_declination, ISHA_ANGLE))

# Convert prayer times from solar to standard time
standard_fajr = solar2standard(solar_fajr, utc_diff, longitude, eq_of_time)
standard_sunrise = solar2standard(solar_sunrise, utc_diff, longitude, eq_of_time)
standard_noon = solar2standard(12.0, utc_diff, longitude, eq_of_time)
standard_asr = standard_noon + asr_time(latitude, sun_declination)
standard_sunset = solar2standard(solar_sunset, utc_diff, longitude, eq_of_time)
standard_maghrib = solar2standard(solar_maghrib, utc_diff, longitude, eq_of_time)
standard_isha = solar2standard(solar_isha, utc_diff, longitude, eq_of_time)
standard_midnight = time_midpoint(standard_sunset, find_tomorrow_fajr(jd, utc_diff, longitude, eq_of_time))


### Find Next New Moon (and the rest of the phases)
moon_phases = next_phases_of_moon_utc(today)
for i, phase in enumerate(moon_phases):
    phase_str = ""
    if i == 0:
        phase_str = "New Moon"
    elif i == 1:
        phase_str = "First Quarter"
    elif i == 2:
        phase_str = "Full Moon"
    else:
        phase_str = "Last Quarter"
    
    moon_phases[i] = {"phase": phase_str, "datetime": phase - datetime.timedelta(hours=utc_diff)}

moon_phases = sorted(moon_phases, key = lambda item: item["datetime"])

# TODO: Calculate Moon Illumination
moon_illumin = moon_illumination(sun_declination, sun_factors[10], sun_factors[4], 
                                 moon_declination, moon_factors[6], moon_factors[1], 
                                 moon_factors[0], sun_factors[6], moon_factors[2] / ASTRONOMICAL_UNIT)

# TODO: Find Day of Clear Visibility of Next New Moon
q_values = [
    [1.0, "Easily Visible"],
    [2.0, "Easily Visible"],
    [3.0, "Easily Visible"],
]


### TODO: Calculate Current Islamic Date


### Misc
# Distance from Given Coordinates to Mecca + direction
mecca_distance, mecca_direction = haversine(latitude, longitude)
mecca_direction = bound_angle_deg(mecca_direction)


##### Outputs #####
# Date & Time
#TODO: Islamic Calendar
print("Time & Date\n\tGregorian Date:\t\t{}".format(today.strftime("%A, %d %B, %Y")))
print("\t24h-Time:\t\t{}\n\tTime Zone:\t\t{} {}".format(today.strftime("%X"), tz_name, format_utc_offset(utc_diff * -1)))
print("\tEquation of time:\t{:.2f} mins".format(equation_of_time(jd, latitude, longitude)))

# Prayer Times
print("Prayer Times\n\tFajr:\t\t\t{}".format(float_to_24time(standard_fajr)))
print("\tSunrise:\t\t{}".format(float_to_24time(standard_sunrise)))
print("\tẒuhr:\t\t\t {}".format(float_to_24time(standard_noon)))
print("\t`Asr:\t\t\t{}".format(float_to_24time(standard_asr)))
print("\tSunset: \t\t{}".format(float_to_24time(standard_sunset)))
print("\tMaghrib: \t\t{}".format(float_to_24time(standard_maghrib)))
print("\t`Isha: \t\t\t{}".format(float_to_24time(standard_isha)))
print("\tMidnight: \t\t{}".format(float_to_24time(standard_midnight)))

# Mecca
print("Mecca\n\tDistance: \t\t{:.2f} km".format(mecca_distance))
print("\tDirection: \t\t{} ({:.2f}°)".format(get_cardinal_direction(np.round(mecca_direction)), mecca_direction))

# The Sun
print("The Sun\n\tApp. Declination:\t{:.3f}°".format(sun_declination))
print("\tApp. Right Ascenscion:\t{}h {}m {:.2f}s".format(alpha_hour, alpha_minute, alpha_second))
print("\tAltitude:\t\t{:.2f}°".format(sun_alt))
print("\tAzimuth:\t\t{:.2f}°".format(sun_az))
# print("\tSolar Hour Angle:\t{:.2f}°".format(solar_angle))

# The Moon
print("The Moon\n\tApp. Declination:\t{:.3f}°".format(moon_declination))
print("\tApp. Right Ascenscion:\t{}h {}m {:.2f}s".format(moon_alpha[0], moon_alpha[1], moon_alpha[2]))
print("\tAltitude:\t\t{:.2f}°".format(moon_alt))
print("\tAzimuth:\t\t{:.2f}°".format(moon_az))
print("\tIllumination:\t\t{:.2f}%".format(moon_illumin * 100))

# Moon Phases
print("Moon Phases\n\t{}:\t\t{}".format(moon_phases[0]["phase"], moon_phases[0]["datetime"].strftime("%H:%M:%S %A, %d %B, %Y")))
print("\t{}:\t\t{}".format(moon_phases[1]["phase"], moon_phases[1]["datetime"].strftime("%H:%M:%S %A, %d %B, %Y")))
print("\t{}:\t\t{}".format(moon_phases[2]["phase"], moon_phases[2]["datetime"].strftime("%H:%M:%S %A, %d %B, %Y")))
print("\t{}:\t\t{}".format(moon_phases[3]["phase"], moon_phases[3]["datetime"].strftime("%H:%M:%S %A, %d %B, %Y")))

# TODO: New Moon First Visibility
print("Visibility (Q) of New Moon\n\t0 days after:\t\t{:.3f} ({})".format(q_values[0][0], q_values[0][1]))
print("\t1 days after:\t\t{:.3f} ({})".format(q_values[1][0], q_values[1][1]))
print("\t2 days after:\t\t{:.3f} ({})".format(q_values[2][0], q_values[2][1]))