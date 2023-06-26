from prayer_times import *

class ITLocation:
    ##### Definitions #####
    FAJR_ANGLE = 16
    MAGHRIB_ANGLE = 4
    ISHA_ANGLE = 14
    MECCA_LAT = 21.420164986
    MECCA_LONG = 39.822330044

    def __init__(self, latitude, longitude, elev):
        self.latitude = latitude
        self.longitude = longitude
        self.elev = elev
    
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


    def prayertimes(self):
        ### Prayer Time Calculations
        # Calculate prayer times in solar time
        solar_fajr = sunrise_sunset(-1, solar_hour_angle(self.latitude, self.sun_declination, self.FAJR_ANGLE))
        solar_sunrise = sunrise_sunset(-1, self.solar_angle)
        solar_sunset = sunrise_sunset(1, self.solar_angle)
        solar_maghrib = sunrise_sunset(1, solar_hour_angle(self.latitude, self.sun_declination, self.MAGHRIB_ANGLE))
        solar_isha = sunrise_sunset(1, solar_hour_angle(self.latitude, self.sun_declination, self.ISHA_ANGLE))

        # Convert prayer times from solar to standard time
        standard_fajr = solar2standard(solar_fajr, self.utc_diff, self.longitude, self.eq_of_time)
        standard_sunrise = solar2standard(solar_sunrise, self.utc_diff, self.longitude, self.eq_of_time)
        standard_noon = solar2standard(12.0, self.utc_diff, self.longitude, self.eq_of_time)
        standard_asr = standard_noon + asr_time(self.latitude, self.sun_declination)
        standard_sunset = solar2standard(solar_sunset, self.utc_diff, self.longitude, self.eq_of_time)
        standard_maghrib = solar2standard(solar_maghrib, self.utc_diff, self.longitude, self.eq_of_time)
        standard_isha = solar2standard(solar_isha, self.utc_diff, self.longitude, self.eq_of_time)
        standard_midnight = time_midpoint(standard_sunset, find_tomorrow_fajr(self.jd, self.utc_diff, self.latitude, self.longitude, self.eq_of_time, self.FAJR_ANGLE))

        return [standard_fajr, standard_sunrise, standard_noon, standard_asr, standard_sunset, standard_maghrib, standard_isha, standard_midnight]
            

    