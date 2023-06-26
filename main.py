from prayer_times import *
from hijri_converter import convert

##### Definitions #####
FAJR_ANGLE          = 16
MAGHRIB_ANGLE       = 4
ISHA_ANGLE          = 14
MECCA_LAT           = 21.420164986 
MECCA_LONG          = 39.822330044
TO_LAT 				= 43.74506
TO_LONG             = -79.30947
TO_ELEV             = 170.5

##### Inputs #####
latitude = TO_LAT #37.336111
longitude = TO_LONG #-121.890556
elev = TO_ELEV #25

##### Calculations #####
### Find UTC Offset According to Lat/Long & adjust datetime
today = datetime.datetime.utcnow()
tz_name, utc_diff = find_utc_offset(latitude, longitude, today)
today += datetime.timedelta(hours=utc_diff)
utc_diff *= -1


### Calculate Julian Date
jd = gregorian_to_jd(today.year, today.month, today.day + fraction_of_day(today), -1 * utc_diff)


### Sun & Moon Properties Calculations
# Get factors arrays
sun_factors = sunpos(jd, latitude, longitude)
delPsi, delEps = sun_nutation(jd)
moon_factors = moonpos(jd, latitude, longitude, delPsi, sun_factors[13], elev)

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
standard_midnight = time_midpoint(standard_sunset, find_tomorrow_fajr(jd, utc_diff, latitude, longitude, eq_of_time, FAJR_ANGLE))


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

### Calculate Moon Illumination
moon_illumin = moon_illumination(sun_declination, sun_factors[10], sun_factors[4], 
                                 moon_declination, moon_factors[6], moon_factors[1], 
                                 moon_factors[0], sun_factors[6], moon_factors[2] / ASTRONOMICAL_UNIT)

### Calculate New Moon Visibilities

# Get New Moon Date from moon_phases list
for item in moon_phases:
        if item['phase'] == "New Moon":
            new_moon = item['datetime']

# Find JD for the given date; adjust day for difference in UTC and local timezone
jd_new_moon = gregorian_to_jd(new_moon.year, new_moon.month, new_moon.day + fraction_of_day(new_moon), -1 * utc_diff)
if new_moon.day != jd_to_gregorian(jd_new_moon).day:
    if new_moon.day < jd_to_gregorian(jd_new_moon).day:
        new_moon += datetime.timedelta(days=1)
    else:
        new_moon -= datetime.timedelta(days=1)
    jd_new_moon = gregorian_to_jd(new_moon.year, new_moon.month, new_moon.day, -1 * utc_diff)

# Find local sunset as visibilities are calculated from then
nm_sun_factors = sunpos(jd_new_moon, latitude, longitude)
nm_sunset = solar2standard(sunrise_sunset(1, solar_hour_angle(latitude, nm_sun_factors[11])), utc_diff, longitude, equation_of_time(jd_new_moon, latitude, longitude))
jd_new_moon = gregorian_to_jd(new_moon.year, new_moon.month, new_moon.day + nm_sunset / 24, -1 * utc_diff)

# Find visibilities for the three days
visibilities = []
i = 0
while (i < 3):
    nm_sun_factors = sunpos(jd_new_moon + i, latitude, longitude)
    delPsi, delEps = sun_nutation(jd_new_moon + i)
    nm_moon_factors = moonpos(jd_new_moon + i, latitude, longitude, delPsi, nm_sun_factors[13], elev)

    #print(jd_to_gregorian(jd_new_moon + i))
    visibilities.append(calculate_visibility(nm_sun_factors[16], nm_sun_factors[15], nm_moon_factors[10], nm_moon_factors[9], np.deg2rad(nm_moon_factors[8])))
    i += 1

# Arrange and classify visibilties
q_values = [
    [visibilities[0], classify_visibility(visibilities[0])],
    [visibilities[1], classify_visibility(visibilities[1])],
    [visibilities[2], classify_visibility(visibilities[2])],
]

### Calculate Current Islamic Date (estimate)
# TODO: Look into newer versions of this, see if it can be corrected.
islamic_date = gregorian_to_hijri(today.year, today.month, today.day)
hijri_date = convert.Gregorian(today.year, today.month, today.day).to_hijri()
#print(hijri_date)

### Misc
# Distance from Given Coordinates to Mecca + direction
mecca_distance, mecca_direction = haversine(latitude, longitude, MECCA_LAT, MECCA_LONG)
mecca_direction = bound_angle_deg(mecca_direction)


##### Outputs #####
# Date & Time
print("Time & Date\n\tGregorian Date:\t\t{}".format(today.strftime("%A, %d %B, %Y")))
print(f"\tIslamic Date:\t\t{get_islamic_day(today.strftime('%A'))}, {islamic_date[2]} {get_islamic_month(islamic_date[1])}, {islamic_date[0]}")
print("\t24h-Time:\t\t{}\n\tTime Zone:\t\t{} {}".format(today.strftime("%X"), tz_name, format_utc_offset(utc_diff * -1)))
print("\tEquation of time:\t{:.2f} minutes".format(equation_of_time(jd, latitude, longitude)))

# Prayer Times
print("Prayer Times\n\tFajr:\t\t\t{}".format(float_to_24time(standard_fajr)))
print("\tSunrise:\t\t{}".format(float_to_24time(standard_sunrise)))
print("\tẒuhr:\t\t\t {}".format(float_to_24time(standard_noon)))
print("\t‘Asr:\t\t\t{}".format(float_to_24time(standard_asr)))
print("\tSunset: \t\t{}".format(float_to_24time(standard_sunset)))
print("\tMaghrib: \t\t{}".format(float_to_24time(standard_maghrib)))
print("\t‘Isha: \t\t\t{}".format(float_to_24time(standard_isha)))
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
print("\tIllumination:\t\t{:.1f}%".format(moon_illumin * 100))

# Moon Phases
print("Moon Phases\n\t{}:\t\t{}".format(moon_phases[0]["phase"], moon_phases[0]["datetime"].strftime("%H:%M:%S %A, %d %B, %Y")))
print("\t{}:\t\t{}".format(moon_phases[1]["phase"], moon_phases[1]["datetime"].strftime("%H:%M:%S %A, %d %B, %Y")))
print("\t{}:\t\t{}".format(moon_phases[2]["phase"], moon_phases[2]["datetime"].strftime("%H:%M:%S %A, %d %B, %Y")))
print("\t{}:\t\t{}".format(moon_phases[3]["phase"], moon_phases[3]["datetime"].strftime("%H:%M:%S %A, %d %B, %Y")))

# TODO: New Moon First Visibility
print("Visibility Values of New Moon after\n\t0 days:\t\t\t{:.3f} ({})".format(q_values[0][0], q_values[0][1]))
print("\t1 day:\t\t\t{:.3f} ({})".format(q_values[1][0], q_values[1][1]))
print("\t2 days:\t\t\t{:.3f} ({})".format(q_values[2][0], q_values[2][1]))