from moon_equations import *

##### Functions #####
def find_tomorrow_fajr(jd, utc_change, lat, long, eq_of_time_minutes, angle):
    jd_tomorrow = jd + 1

    tomorrow_sun_declination = sunpos(jd_tomorrow, lat, long)[11]
    tomorrow_fajr_solar_angle = solar_hour_angle(lat, tomorrow_sun_declination, angle)
    tomorrow_solar_fajr = sunrise_sunset(-1, tomorrow_fajr_solar_angle)
    tomorrow_standard_fajr = solar2standard(tomorrow_solar_fajr, utc_change, long, eq_of_time_minutes)

    return tomorrow_standard_fajr

def asr_time(lat, dec, t = 1):
    temp_num = sin(np.rad2deg(np.arctan2(1, t + tan(lat - dec)))) - sin(lat) * sin(dec)
    temp_denom = cos(lat) * cos(dec)
    return 1 / 15 * np.rad2deg(np.arccos(temp_num / temp_denom))