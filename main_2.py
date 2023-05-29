import ephem
import math
import pytz
from datetime import datetime, timedelta, timezone
from timezonefinder import TimezoneFinder

def dms_to_decimal(dms):
    sign = -1 if dms.startswith('-') else 1
    dms = dms.lstrip('-')
    degrees, minutes, seconds = map(float, dms.split(':'))
    decimal_degrees = sign * (abs(degrees) + (minutes / 60) + (seconds / 3600))
    return decimal_degrees

def calculate_values(date, v_time, lat, long):
    # Define the observer
    observer = ephem.Observer()
    observer.date = f'{date} {v_time}'
    observer.lat = str(lat)
    observer.long = str(long)

    # Get the Sun and Moon for the given date and location
    sun = ephem.Sun(observer)
    moon = ephem.Moon(observer)

    # Calculate ARCL
    arcl_rad = ephem.separation((sun.az, sun.alt), (moon.az, moon.alt))
    arcl = math.degrees(arcl_rad)

    # Calculate ARCV
    arcv = abs(dms_to_decimal(str(sun.alt)) - dms_to_decimal(str(moon.alt)))

    # Calculate DAZ
    daz = abs(dms_to_decimal(str(sun.az)) - dms_to_decimal(str(moon.az)))
    if daz > 180:
        daz -= 360
    elif daz < -180:
        daz += 360

    # Calculate Parallax
    # We convert it to kilometers for calculation
    moon_distance_km = moon.earth_distance * 149597870.7 # 1 AU = 149597870.7 KM

    # Earth's radius in kilometers
    earth_radius_km = 6378.1

    # Parallax in radians
    moon_parallax_rad = earth_radius_km / moon_distance_km

    # Convert it to arcminutes for easier interpretation
    moon_pi = moon_parallax_rad * 206265 / 60

    moon_alt = math.degrees(moon.alt)

    return arcl, arcv, daz, moon_pi, moon_alt

def next_new_moon(date):
    next_new_moon_date = ephem.next_new_moon(date)
    return next_new_moon_date

def find_sunset_time(date, lat, long):
    observer = ephem.Observer()
    observer.date = date
    observer.lat = lat
    observer.long = long
    return observer.next_setting(ephem.Sun())

def convert_utc_to_offset(utc_time, offset_hours):
    # Create a UTC timezone object
    utc_timezone = pytz.timezone('UTC')

    # Convert the input UTC time to a datetime object
    utc_datetime = datetime.strptime(utc_time, '%Y/%m/%d %H:%M:%S')

    # Set the UTC timezone to the datetime object
    utc_datetime = utc_timezone.localize(utc_datetime)

    # Create the target timezone object based on the specified offset
    target_timezone = pytz.FixedOffset(offset_hours * 60)

    # Convert the UTC time to the target timezone
    target_datetime = utc_datetime.astimezone(target_timezone)

    # Format the target datetime as a string
    target_time = target_datetime.strftime('%Y/%m/%d %H:%M:%S')

    return target_time

def get_utc_offset(lat, long, date):
    # start_time = time.time()
    # # Create a TimezoneFinder object
    # tf = TimezoneFinder()

    # # Get the timezone name for the coordinates
    # timezone_str = tf.certain_timezone_at(lat=lat, lng=long)

    # end_time = time.time()

    # execution_time = end_time - start_time
    # print(f"Execution Time: {execution_time:.2f} seconds")

    # # Get the timezone object for the given name
    # timezone = pytz.timezone(timezone_str)

    # # Localize the given date with the timezone
    # localized_datetime = timezone.localize(datetime.combine(date, datetime.min.time()))

    # # Get the UTC offset for the given timezone and date
    # utc_offset = localized_datetime.utcoffset()

    # # Calculate the UTC offset in hours
    # offset_hours = int(utc_offset.total_seconds() / 3600)

    offset_hours = long / 15

    return offset_hours

def find_moon_visibility(lat, long, date, time):
    # Convert string input to datetime for manipulation
    current_date = datetime.strptime(f'{date} {time}', '%Y/%m/%d %H:%M:%S')

    # Step 2: Call the next_new_moon function passing the local date as the date for the function
    new_moon_date = next_new_moon(current_date)

    # Define the observer
    observer = ephem.Observer()
    observer.lat = lat
    observer.long = long

    while True:
        # Convert new_moon_date to string format for ephem library
        new_moon_date_str = new_moon_date.strftime('%Y/%m/%d')

        # Step 4: Find the exact time of sunset for the new_moon_date
        observer.date = new_moon_date_str
        sunset_time = ephem.localtime(observer.next_setting(ephem.Sun())).time()

        # Step 5: Call the calculate_values function passing the geographic coordinates, 
        # the new_moon_date, and the time of sunset on that new_moon_date
        arcv, arcl, daz = calculate_values(new_moon_date_str, sunset_time.strftime('%H:%M:%S'), lat, long)

        # Step 6: Check the visibility condition
        q_value = arcv - (10.3743 - 0.0137 * abs(daz) - 0.0097 * daz**2)
        if q_value >= 0:
            # Step 7: If the check is true, terminate and print the message
            print(f'The new moon will be easily visible at your location on {new_moon_date_str} at sunset.')
            print(f'q-value is {q_value:.2f}.')
            break
        else:
            # If the check is false, repeat steps 4 through 6 but use one day after.
            new_moon_date += timedelta(days=1)

# date = datetime(1984, 10, 26)
# observer = ephem.Observer()
# observer.lat = "15.6"
# observer.long = "35.6"
# lat = float(str(observer.lat).split(':', 1)[0]) + float(str(observer.lat).split(':', 2)[1])/60
# long = float(str(observer.long).split(':', 1)[0]) + float(str(observer.long).split(':', 2)[1])/60

# utc_offset = get_utc_offset(lat, long, date)
# print(utc_offset)

# k = date - timedelta(hours=utc_offset)
# observer.date = f'{k:%Y/%m/%d %H:%M:%S}'
# print(observer.date)

# sunset_utc = observer.next_setting(ephem.Sun())
# print(sunset_utc)
# sunset_local = convert_utc_to_offset(str(sunset_utc), utc_offset)
# print(sunset_local)

# new_moon_datetime_utc = str(next_new_moon(date))
# print(new_moon_datetime_utc)
# new_moon_datetime_local = convert_utc_to_offset(new_moon_datetime_utc, utc_offset)
# print(new_moon_datetime_local)

# d1 = datetime.strptime(sunset_local, '%Y/%m/%d %H:%M:%S')
# d2 = datetime.strptime(new_moon_datetime_local, '%Y/%m/%d %H:%M:%S')

# if d1 > d2:
#     print("Sunset happens AFTER the new moon.")
# else:
#     print("Sunset happens BEFORE the new moon.")

# d3 = datetime.strptime(str(sunset_utc), '%Y/%m/%d %H:%M:%S')

# arcl, arcv, daz, pi, alt = calculate_values(f'{d3:%Y/%m/%d}', f'{d3:%H:%M:%S}', lat, long)

# semi_diameter = 0.27245 * pi
# semi_diameter_prime = semi_diameter * (1 + math.sin(math.radians(alt)) * math.sin(math.radians(pi / 60)))

# w_prime = semi_diameter_prime * (1 - math.cos(math.radians(arcl)))

# q_value = (arcv - (11.8371 - 6.3226 * w_prime + 0.7319 * w_prime ** 2 - 0.1018 * w_prime ** 3)) / 10

# print(f'Lat: {lat:.1f}°, Long: {long:.1f}°, ARCL: {arcl:.1f}°, ARCV: {arcv:.1f}°, DAZ: {daz:.1f}°, Moon Pi: {pi:.1f}’, Moon Alt: {alt:.2f}°, W\': {w_prime:.1f}’, \
# q: {q_value:.3f}')