import ephem

def moon_parallax(observer_longitude, observer_latitude, date_time):
    observer = ephem.Observer()
    observer.lon = observer_longitude
    observer.lat = observer_latitude
    observer.date = date_time

    moon = ephem.Moon()
    moon.compute(observer)

    # The moon's distance is stored in the moon.distance attribute and is in Astronomical Units
    # We convert it to kilometers for calculation
    moon_distance_km = moon.earth_distance * 149597870.7 # 1 AU = 149597870.7 KM

    # Earth's radius in kilometers
    earth_radius_km = 6378.1

    # Parallax in radians
    moon_parallax_rad = earth_radius_km / moon_distance_km

    # Convert it to arcseconds for easier interpretation
    moon_parallax_arcsec = moon_parallax_rad * 206265

    return moon_parallax_arcsec

# Test the function
longitude = '15.6'   # Greenwich meridian
latitude = '35.6' # London latitude
date_time = '1984/10/26 15:11:00' # YYYY/MM/DD HH:MM:SS

parallax = moon_parallax(longitude, latitude, date_time) / 60
print(f"The Moon's parallax on {date_time} at coordinates ({longitude}, {latitude}) is {parallax} arcminutes.")