from islamic_times import sun_equations as se
from islamic_times import time_equations as te
from islamic_times import calculation_equations as ce
import numpy as np

# Used to calculate islamic midnight
def find_tomorrow_time(jd, deltaT, utc_change, lat, long, eq_of_time_minutes, angle):
    jd_tomorrow = jd + 1

    tomorrow_sun_declination = se.sunpos(jd_tomorrow, deltaT, lat, long)[11]
    tomorrow_time_solar_angle = se.solar_hour_angle(lat, tomorrow_sun_declination, angle)
    tomorrow_solar_time = se.sunrise_sunset(-1, tomorrow_time_solar_angle)
    tomorrow_standard_time = te.solar2standard(tomorrow_solar_time, utc_change, long, eq_of_time_minutes)

    return tomorrow_standard_time

# ʿAṣr time
# Type 0: Most schools
# Type 1: Ḥanafī definition
def asr_time(lat, dec, t = 1):
    temp_num = ce.sin(np.rad2deg(np.arctan2(1, t + ce.tan(lat - dec)))) - ce.sin(lat) * ce.sin(dec)
    temp_denom = ce.cos(lat) * ce.cos(dec)
    return 1 / 15 * np.rad2deg(np.arccos(temp_num / temp_denom))