# Used to test ideas
import math
import ephem
from main import gregorian_to_jd, datetime, fraction_of_day, moonpos, sunpos, TO_LAT, TO_LONG
from main_2 import dms_to_decimal
# from astropy.time import Time
# from astropy.coordinates import get_body, EarthLocation, AltAz
# from astropy import units as u

def calculate_angle_diff(azimuth1, altitude1, azimuth2, altitude2):
    # Convert degrees to radians
    azimuth1, altitude1, azimuth2, altitude2 = map(math.radians, [azimuth1, altitude1, azimuth2, altitude2])

    # Apply the haversine formula
    delta_azimuth = azimuth2 - azimuth1
    delta_altitude = altitude2 - altitude1

    a = math.sin(delta_altitude / 2)**2 + math.cos(altitude1) * math.cos(altitude2) * math.sin(delta_azimuth / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

    # Convert the angle difference from radians to degrees
    angle_diff = math.degrees(c)
    
    return angle_diff

today = datetime.datetime.now()
today_utc = today + datetime.timedelta(hours=4)
jd = gregorian_to_jd(today.year, today.month, today.day + fraction_of_day(today), -4)
moon_factors = moonpos(jd, TO_LAT, TO_LONG)
sun_factors = sunpos(jd, TO_LAT, TO_LONG)

moon_az = moon_factors[8]
moon_alt = moon_factors[7]

sun_az = sun_factors[16]
sun_alt = sun_factors[15]

observer = ephem.Observer()
observer.date = today_utc.strftime('%Y/%m/%d %H:%M:%S')
observer.lat = str(TO_LAT)
observer.long = str(TO_LONG)

# Get the Sun and Moon for the given date and location
sun = ephem.Sun(observer)
moon = ephem.Moon(observer)

# Calculate ARCL
arcl_rad = ephem.separation((sun.az, sun.alt), (moon.az, moon.alt))
arcl = math.degrees(arcl_rad)

# # Choose your location
# location = EarthLocation(lon=f'{str(TO_LONG)}d', lat=f'{str(TO_LAT)}d', height=0*u.m)

# # Choose the time for which you want to calculate the declination
# time = Time(today)  # you can also specify a particular time

# # Get the position of the moon at this time
# moon_2 = get_body("moon", time, location)
# sun_2 = get_body("sun", time, location)

# # Convert to AltAz coordinates
# moon_altaz = moon_2.transform_to(AltAz(obstime=time, location=location))
# sun_altaz = sun_2.transform_to(AltAz(obstime=time, location=location))


#print(today, observer.date, time)
print(sun_factors[11], dms_to_decimal(str(sun.dec)))#, sun_2.dec.degree)
print(moon_factors[5], dms_to_decimal(str(moon.dec)))#, moon_2.dec.degree)
print(sun_az, dms_to_decimal(str(sun.az)))#, sun_altaz.az.degree)
print(sun_alt, dms_to_decimal(str(sun.alt)))#, sun_altaz.alt.degree)
print(moon_az, dms_to_decimal(str(moon.az)))#, moon_altaz.az.degree)
print(moon_alt, dms_to_decimal(str(moon.alt)))#, moon_altaz.alt.degree)
print(calculate_angle_diff(sun_az, sun_alt, moon_az, moon_alt), arcl)