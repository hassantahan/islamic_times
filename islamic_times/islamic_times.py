import datetime
import numpy as np
from islamic_times import prayer_times as pt
from islamic_times import sun_equations as se
from islamic_times import moon_equations as me
from islamic_times import time_equations as te
from islamic_times import calculation_equations as ce

class ITLocation:
    ##### Mecca Constants #####
    MECCA_LAT = 21.420164986
    MECCA_LONG = 39.822330044

    ##### Prayer Method Names #####
    mwl = ['mwl', 'muslim world league']
    isna = ['isna', 'islamic society of north america']
    egypt = ['egypt', 'egyptian', 'egyptian general authority of survey', 'egas']
    makkah = ['makkah', 'mecca', 'mekkah', 'umm al-qura university', 'umm al-qura', 'uqu']
    karachi = ['karachi', 'university of islamic sciences', 'uis']
    tehran = ['tehran', 'university of tehran', 'institute of geophysics', 'uot', 'iog', 'ioguot', 'uotiog', 'igut', 'utig']
    jafari = ['jafari', 'jaafari', 'shia', 'shia ithna ashari', 'leva', 'leva research institute', 'qom', 'qum', 'lri', 'sia', 'sialri']

    # fajr, isha, maghrib, midnight (0 --> sunset to sunrise, 1 --> sunset to fajr)
    mwl_vals = (18, 17, 0, 0)
    isna_vals = (15, 15, 0, 0)
    egypt_vals = (19.5, 17.5, 0, 0)
    makkah_vals = (18.5, np.inf, 0, 0)
    karachi_vals = (18, 18, 0, 0)
    tehran_vals = (17.7, 14, 4.5, 1)
    jafari_vals = (16, 14, 4, 1)

    def __init__(self, latitude = 51.477928, longitude = -0.001545, elevation = 76, temperature = 10, pressure = 101, today = datetime.datetime.now(datetime.timezone.utc), method = 'jafari', asr_type = 0):
        self.latitude = latitude
        self.longitude = longitude
        self.elevation = elevation
        self.temperature = temperature
        self.pressure = pressure
        self.today = today

        ### Find UTC Offset According to Lat/Long & adjust datetime
        # This is very computationally expensive
        if self.today.tzinfo == None: 
            self.tz_name, self.utc_diff = te.find_utc_offset(self.latitude, self.longitude, self.today)
        else:  
            self.tz_name, self.utc_diff = self.today.tzinfo, self.today.utcoffset().total_seconds() / 3600

        self.today += datetime.timedelta(hours=self.utc_diff)
        self.utc_diff *= -1

        # Calculate the astronomical parameters
        self.calculate()
        
        # Set the default prayer definitions
        self.asr_type = asr_type
        self.method = method
        self.set_prayer_method(method) # Prayer times are calculated in this method

    def update_time(self, date_time = None):
        # Time is only updates to the latest if no argument is passed
        if date_time is None:
            self.today = datetime.datetime.now(self.today.tzinfo)
        else:
            self.today = date_time

    def calculate(self):
        ### Calculate Julian Date
        self.jd = te.gregorian_to_jd(self.today.year, self.today.month, self.today.day + te.fraction_of_day(self.today), -1 * self.utc_diff)
        self.deltaT = te.delta_t_approx(self.today.year, self.today.month)
        self.jde = self.jd + self.deltaT / 86400

        ### Calculate Current Islamic Date (estimate)
        # TODO: Look into newer versions of this, see if it can be corrected. See function for more info.
        self.islamic_date = te.gregorian_to_hijri(self.today.year, self.today.month, self.today.day)

        ### Sun & Moon Properties Calculations
        # Get factors arrays
        self.sun_factors = se.sunpos(self.jde, self.deltaT, self.latitude, self.longitude, self.temperature, self.pressure)
        self.delPsi, self.delEps = se.sun_nutation(self.jde)
        self.moon_factors = me.moonpos(self.jde, self.deltaT ,self.latitude, self.longitude, self.delPsi, self.sun_factors[13], self.elevation)

        # Important Sun Factors placed into local variables
        self.sun_declination = self.sun_factors[11]
        self.sun_alpha = ce.decimal_to_hms(self.sun_factors[10])
        self.solar_angle = se.solar_hour_angle(self.latitude, self.sun_declination)
        self.sun_alt = self.sun_factors[15]
        self.sun_az = self.sun_factors[16]

        # Important Moon Factors placed into local variables
        self.moon_declination = self.moon_factors[5]
        self.moon_alpha = ce.decimal_to_hms(self.moon_factors[4])
        self.moon_pi = self.moon_factors[8]
        self.moon_alt = self.moon_factors[9]
        self.moon_az = self.moon_factors[10]

        # Moon calculations
        self.moon_illumin = me.moon_illumination(self.sun_declination, self.sun_factors[10], self.moon_declination, 
                                            self.moon_factors[4], self.sun_factors[6], self.moon_factors[2] / 149597870.7)
        self.moonset = te.float_to_24time(me.calculate_moonset(self.jd, self.deltaT, self.latitude, self.longitude, self.elevation, self.utc_diff))

    # Prayer Time Calculations
    def calculate_prayer_times(self):
        # Equation of Time
        eq_of_time = se.equation_of_time(self.jde, self.deltaT, self.latitude, self.longitude)

        # Calculate prayer times in solar time
        solar_fajr = se.sunrise_sunset(-1, se.solar_hour_angle(self.latitude, self.sun_declination, self.fajr_angle))
        solar_sunrise = se.sunrise_sunset(-1, self.solar_angle)
        solar_sunset = se.sunrise_sunset(1, self.solar_angle)

        # Only if maghrib is not at sunset
        if self.maghrib_angle > 0:
            solar_maghrib = se.sunrise_sunset(1, se.solar_hour_angle(self.latitude, self.sun_declination, self.maghrib_angle))

        # Convert prayer times from solar to standard time
        self.standard_fajr = te.solar2standard(solar_fajr, self.utc_diff, self.longitude, eq_of_time)
        self.standard_sunrise = te.solar2standard(solar_sunrise, self.utc_diff, self.longitude, eq_of_time)
        self.standard_noon = te.solar2standard(12.0, self.utc_diff, self.longitude, eq_of_time)
        
        # Type 1 is Ḥanafī
        if self.asr_type:
            self.standard_asr = self.standard_noon + pt.asr_time(self.latitude, self.sun_declination, t = 2)
        # Type 0 is for the rest
        else:
            self.standard_asr = self.standard_noon + pt.asr_time(self.latitude, self.sun_declination)

        self.standard_sunset = te.solar2standard(solar_sunset, self.utc_diff, self.longitude, eq_of_time)
        
        # Only if maghrib is not at sunset
        if self.maghrib_angle > 0:
            self.standard_maghrib = te.solar2standard(solar_maghrib, self.utc_diff, self.longitude, eq_of_time)
        else:
            # Otherwise Maghrib is sunset
            self.standard_maghrib = self.standard_sunset

        # If NOT makkah method
        if self.isha_angle is not np.inf:
            solar_isha = se.sunrise_sunset(1, se.solar_hour_angle(self.latitude, self.sun_declination, self.isha_angle))
            self.standard_isha = te.solar2standard(solar_isha, self.utc_diff, self.longitude, eq_of_time)
        # Makkah method is special
        else:
            # Ramadan has isha as two hours after maghrib
            if self.islamic_date[2] == 9:
                self.standard_isha = self.standard_maghrib + 2
            # Otherwise it is one hour
            else:
                self.standard_isha = self.standard_maghrib + 1
        
        if self.midnight_type:
            self.standard_midnight = te.time_midpoint(self.standard_sunset, pt.find_tomorrow_time(self.jde, self.deltaT, self.utc_diff, self.latitude, self.longitude, eq_of_time, self.fajr_angle))
        else:
            self.standard_midnight = te.time_midpoint(self.standard_sunset, pt.find_tomorrow_time(self.jde, self.deltaT, self.utc_diff, self.latitude, self.longitude, eq_of_time, 0))

    def set_prayer_method(self, method = 'jafari'):
        # Mapping methods to their respective values
        method_groups = {
            "Muslim World League (MWL)": (self.mwl, self.mwl_vals),
            "Islamic Society of North America (ISNA)": (self.isna, self.isna_vals),
            "Egyptian General Authority of Survey (Egypt)": (self.egypt, self.egypt_vals),
            "Umm al-Qura University (Makkah)": (self.makkah, self.makkah_vals),
            "University of Islamic Sciences, (Karachi)": (self.karachi, self.karachi_vals),
            "Institute of Geophysics, University of Tehran (Tehran)": (self.tehran, self.tehran_vals),
            "Shia Ithna Ashari, Leva Research Institute, Qom (Jafari)": (self.jafari, self.jafari_vals),
        }

        # Find the matching group
        for group_name, (group, values) in method_groups.items():
            if method.lower() in group:
                self.fajr_angle, self.isha_angle, self.maghrib_angle, self.midnight_type = values
                self.method = group_name
                break
        else:
            # Raise an error if no match is found
            raise ValueError("Not a valid method. See documentation for a list of acceptable method names.")
        
        # Update prayer times
        self.calculate_prayer_times()
    
    def set_custom_prayer_angles(self, fajr_angle=None, maghrib_angle=None, isha_angle=None):
        # Helper function to validate and set angle
        def validate_and_set(attribute_name, value):
            if value is not None:
                if isinstance(value, (int, float)):  # Check if it's a number
                    if value > 0:
                        setattr(self, attribute_name, value)
                    else:
                        ValueError(f"{attribute_name} must be greater than 0. Invalid value: {value}")
                else:
                    raise ValueError(f"{attribute_name} must be a number. Invalid value: {value}")

        # Validate and set each angle
        validate_and_set("fajr_angle", fajr_angle)
        validate_and_set("maghrib_angle", maghrib_angle)
        validate_and_set("isha_angle", isha_angle)

        # Method is now custom
        self.method = 'Custom'

        # Update prayer times
        self.calculate_prayer_times()

    
    def set_asr_type(self, asr_type = 0):
        # Helper function to validate and set type
        def validate_and_set(attribute_name, value):
            if value is not None:
                if value in [0, 1]:  # Check if it's 0 or 1
                    setattr(self, attribute_name, value)
                else:
                    raise ValueError(f"{attribute_name} must be either 0 or 1. Check documentation to understand each type. Invalid value: {value}")

        # Validate and set asr type
        validate_and_set("asr_type", asr_type)

        # Update prayer times
        self.calculate_prayer_times()

    def set_midnight_type(self, midnight_type = 0):
        # Helper function to validate and set type
        def validate_and_set(attribute_name, value):
            if value is not None:
                if value in [0, 1]:  # Check if it's 0 or 1
                    setattr(self, attribute_name, value)
                else:
                    raise ValueError(f"{attribute_name} must be either 0 or 1. Check documentation to understand each type. Invalid value: {value}")

        # Validate and set asr type
        validate_and_set("midnight_type", midnight_type)

        # Method is now custom
        self.method = 'Custom'

        # Update prayer times
        self.calculate_prayer_times()

    # Return date and time information
    def datetime(self):
        return {
                "gregorian" : self.today.strftime("%A, %d %B, %Y"), 
                "hijri" : f"{te.get_islamic_day(self.today.strftime('%A'))}, {self.islamic_date[2]} {te.get_islamic_month(self.islamic_date[1])}, {self.islamic_date[0]}",
                "time" : self.today.strftime("%X"), "timezone" : self.tz_name, "utc_offset" : te.format_utc_offset(self.utc_diff * -1),
                "jd" : self.jd,
                "eq_of_time" : np.round(se.equation_of_time(self.jde, self.deltaT, self.latitude, self.longitude), 2),
                "deltaT" : np.round(te.delta_t_approx(self.today.year, self.today.month), 2)
            }

    # Return prayer times
    def prayertimes(self):
        return {
                "method": self.method,
                "fajr" : te.float_to_24time(self.standard_fajr),
                "sunrise" : te.float_to_24time(self.standard_sunrise),
                "noon" : te.float_to_24time(self.standard_noon),
                "asr" : te.float_to_24time(self.standard_asr),
                "sunset" : te.float_to_24time(self.standard_sunset),
                "maghrib" : te.float_to_24time(self.standard_maghrib),
                "isha" : te.float_to_24time(self.standard_isha),
                "midnight" : te.float_to_24time(self.standard_midnight)
        }
    
    # Return Mecca information
    def mecca(self):
        mecca_distance, mecca_direction = ce.haversine(self.latitude, self.longitude, self.MECCA_LAT, self.MECCA_LONG)
        mecca_direction = ce.bound_angle_deg(mecca_direction)
        
        return {
                "distance" : np.round(mecca_distance, 2),
                "angle" : np.round(mecca_direction, 2),
                "cardinal" : ce.get_cardinal_direction(np.round(mecca_direction))
            }
    
    # Return sun properties and position values
    def sun(self):
        return {
                "declination" : f"{self.sun_declination:.3f}°",
                "right_ascension" : f"{self.sun_alpha[0]}h {self.sun_alpha[1]}m {self.sun_alpha[2]:.2f}s",
                "altitude" : f"{self.sun_alt:.2f}°",
                "azimuth" : f"{self.sun_az:.2f}°"
            }
    
    # Return moon properties and position values
    def moon(self):
        return {
                "moonset" : self.moonset,
                "declination" : f"{self.moon_declination:.3f}°",
                "right_ascension" : f"{self.moon_alpha[0]}h {self.moon_alpha[1]}m {self.moon_alpha[2]:.2f}s",
                "altitude" : f"{self.moon_alt:.2f}°",
                "azimuth" : f"{self.moon_az:.2f}°",
                "illumination" : f"{self.moon_illumin * 100:.1f}%"
            }
    
    def moonphases(self):
        ### Find Next New Moon (and the rest of the phases)
        moon_phases = me.next_phases_of_moon_utc(self.today)
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
            
            moon_phases[i] = {"phase": phase_str, "datetime": phase - datetime.timedelta(hours=self.utc_diff)}

        moon_phases = sorted(moon_phases, key = lambda item: item["datetime"])

        return moon_phases

    # Calculate Next New Moon Visibilities
    def visibilities(self, days = 3):
        # Get New Moon Date from moon_phases list
        moon_phases = self.moonphases()
        for item in moon_phases:
                if item['phase'] == "New Moon":
                    new_moon = item['datetime']

        # Find JD for the given date; adjust day for difference in UTC and local timezone
        jd_new_moon = te.gregorian_to_jd(new_moon.year, new_moon.month, new_moon.day + te.fraction_of_day(new_moon), -1 * self.utc_diff)
        ymd_new_moon = te.jd_to_gregorian(jd_new_moon)
        deltaT_new_moon = te.delta_t_approx(ymd_new_moon.year, ymd_new_moon.month)
        jde_new_moon = jd_new_moon + deltaT_new_moon / 86400

        # Forgot what this does lol. Likely something to do with timezone differences
        if new_moon.day != ymd_new_moon.day:
            if new_moon.day < ymd_new_moon.day:
                new_moon += datetime.timedelta(days=1)
            else:
                new_moon -= datetime.timedelta(days=1)
            jd_new_moon = te.gregorian_to_jd(new_moon.year, new_moon.month, new_moon.day, -1 * self.utc_diff)

        # Find local sunset as visibilities are calculated from then
        nm_sun_factors = se.sunpos(jde_new_moon, deltaT_new_moon, self.latitude, self.longitude)     

        # Find visibilities for the three days
        visibilities = []
        for i in range(days):
            # First, check if the moonset is before the new moon for the first day
            nm_moonset = me.calculate_moonset(jd_new_moon, deltaT_new_moon, self.latitude, self.longitude, self.elevation, self.utc_diff)
            from math import modf
            minute, hour = modf(nm_moonset)
            if i != 0 or datetime.time(int(hour), 60 * int(minute)) > ymd_new_moon.time():
                # Set the day parameters
                test_jd_new_moon = jd_new_moon + i
                test_ymd_new_moon = te.jd_to_gregorian(test_jd_new_moon)
                test_deltaT_new_moon = te.delta_t_approx(test_ymd_new_moon.year, test_ymd_new_moon.month)
                test_jde_new_moon = test_jd_new_moon + test_deltaT_new_moon / 86400

                # Set sun parameters
                nm_sun_factors = se.sunpos(test_jde_new_moon, test_deltaT_new_moon, self.latitude, self.longitude)

                # Sunset & moonset calculations
                test_nm_moonset = me.calculate_moonset(test_jde_new_moon, test_deltaT_new_moon, self.latitude, self.longitude, self.elevation, self.utc_diff)
                test_nm_sunset = te.solar2standard(se.sunrise_sunset(1, se.solar_hour_angle(self.latitude, nm_sun_factors[11])), self.utc_diff, self.longitude, 
                                      se.equation_of_time(test_jde_new_moon, test_deltaT_new_moon, self.latitude, self.longitude))
                
                # Find the best time which is four ninths the moonset-sunset lag after sunset 
                lag = test_nm_moonset - test_nm_sunset
                best_time = test_nm_sunset + 4 / 9 * lag
                best_time_jd = te.gregorian_to_jd(test_ymd_new_moon.year, test_ymd_new_moon.month, test_ymd_new_moon.day + best_time / 24, -1 * self.utc_diff)
                best_time_jde = best_time_jd + test_deltaT_new_moon / 86400

                # Recalculate sun & calculate moon parameters
                nm_sun_factors = se.sunpos(best_time_jde, test_deltaT_new_moon, self.latitude, self.longitude)
                delPsi, delEps = se.sun_nutation(best_time_jde)
                nm_moon_factors = me.moonpos(best_time_jde, test_deltaT_new_moon, self.latitude, self.longitude, delPsi, nm_sun_factors[13], self.elevation)

                # Visibility is now calculated
                v = me.calculate_visibility(nm_sun_factors[16], nm_sun_factors[15], nm_moon_factors[10], nm_moon_factors[9], np.deg2rad(nm_moon_factors[8]))
            else:
                # Moon is not visibile before the new moon
                v = -np.inf

            visibilities.append(v)

        # Arrange and classify visibilties
        q_values = [
            [visibility, me.classify_visibility(visibility)] 
            for visibility in visibilities
        ]

        return {"0" : q_values[0], "1" : q_values[1], "2" : q_values[2]}