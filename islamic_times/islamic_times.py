import datetime
from islamic_times import prayer_times as pt
from islamic_times import moon_equations as me
from islamic_times import sun_equations as se
from islamic_times import time_equations as te
from islamic_times import calculation_equations as ce
import numpy as np

class ITLocation:
    ##### Definitions #####
    FAJR_ANGLE = 16
    MAGHRIB_ANGLE = 4
    ISHA_ANGLE = 14
    MECCA_LAT = 21.420164986
    MECCA_LONG = 39.822330044

    def __init__(self, latitude, longitude, elev, today = None):
        self.latitude = latitude
        self.longitude = longitude
        self.elev = elev
        self.today = today
        self.update()

    def update(self):
        ### Find UTC Offset According to Lat/Long & adjust datetime
        local = True
        if self.today == None:
            self.today = datetime.datetime.utcnow()
            local = False
        self.tz_name, self.utc_diff = te.find_utc_offset(self.latitude, self.longitude, self.today)
        if not(local): self.today += datetime.timedelta(hours=self.utc_diff)
        self.utc_diff *= -1

        ### Calculate Julian Date
        self.jd = te.gregorian_to_jd(self.today.year, self.today.month, self.today.day + te.fraction_of_day(self.today), -1 * self.utc_diff)

        ### Sun & Moon Properties Calculations
        # Get factors arrays
        self.sun_factors = se.sunpos(self.jd, self.latitude, self.longitude)
        self.delPsi, self.delEps = se.sun_nutation(self.jd)
        self.moon_factors = me.moonpos(self.jd, self.latitude, self.longitude, self.delPsi, self.sun_factors[13], self.elev)

        # Important Sun Factors placed into local variables
        self.sun_declination = self.sun_factors[11]
        self.sun_alpha = ce.decimal_to_hms(self.sun_factors[10])
        self.solar_angle = se.solar_hour_angle(self.latitude, self.sun_declination)
        self.sun_alt = self.sun_factors[15]
        self.sun_az = self.sun_factors[16]

        # Important Moon Factors placed into local variables
        self.moon_declination = self.moon_factors[7]
        self.moon_alpha = ce.decimal_to_hms(self.moon_factors[6])
        self.moon_alt = self.moon_factors[9]
        self.moon_az = self.moon_factors[10]   

        return

    def datetime(self):
        ### Calculate Current Islamic Date (estimate)
        # TODO: Look into newer versions of this, see if it can be corrected.
        islamic_date = te.gregorian_to_hijri(self.today.year, self.today.month, self.today.day)

        return {"gregorian" : self.today.strftime("%A, %d %B, %Y"), 
                "hijri" : f"{te.get_islamic_day(self.today.strftime('%A'))}, {islamic_date[2]} {te.get_islamic_month(islamic_date[1])}, {islamic_date[0]}",
                "time" : self.today.strftime("%X"), "timezone" : self.tz_name, "utc_offset" : te.format_utc_offset(self.utc_diff * -1),
                "eq_of_time" : np.round(se.equation_of_time(self.jd, self.latitude, self.longitude), 2)}

    def prayertimes(self):
        ### Prayer Time Calculations
        # Equation of Time
        eq_of_time = se.equation_of_time(self.jd, self.latitude, self.longitude)

        # Calculate prayer times in solar time
        solar_fajr = se.sunrise_sunset(-1, se.solar_hour_angle(self.latitude, self.sun_declination, self.FAJR_ANGLE))
        solar_sunrise = se.sunrise_sunset(-1, self.solar_angle)
        solar_sunset = se.sunrise_sunset(1, self.solar_angle)
        solar_maghrib = se.sunrise_sunset(1, se.solar_hour_angle(self.latitude, self.sun_declination, self.MAGHRIB_ANGLE))
        solar_isha = se.sunrise_sunset(1, se.solar_hour_angle(self.latitude, self.sun_declination, self.ISHA_ANGLE))

        # Convert prayer times from solar to standard time
        standard_fajr = te.solar2standard(solar_fajr, self.utc_diff, self.longitude, eq_of_time)
        standard_sunrise = te.solar2standard(solar_sunrise, self.utc_diff, self.longitude, eq_of_time)
        standard_noon = te.solar2standard(12.0, self.utc_diff, self.longitude, eq_of_time)
        standard_asr = standard_noon + pt.asr_time(self.latitude, self.sun_declination)
        standard_sunset = te.solar2standard(solar_sunset, self.utc_diff, self.longitude, eq_of_time)
        standard_maghrib = te.solar2standard(solar_maghrib, self.utc_diff, self.longitude, eq_of_time)
        standard_isha = te.solar2standard(solar_isha, self.utc_diff, self.longitude, eq_of_time)
        standard_midnight = te.time_midpoint(standard_sunset, pt.find_tomorrow_fajr(self.jd, self.utc_diff, self.latitude, self.longitude, eq_of_time, self.FAJR_ANGLE))

        return {"fajr" : te.float_to_24time(standard_fajr), "sunrise" : te.float_to_24time(standard_sunrise), "noon" : te.float_to_24time(standard_noon), 
                "asr" : te.float_to_24time(standard_asr), "sunset" : te.float_to_24time(standard_sunset), "maghrib" : te.float_to_24time(standard_maghrib), 
                "isha" : te.float_to_24time(standard_isha), "midnight" : te.float_to_24time(standard_midnight)}
    
    def mecca(self):
        mecca_distance, mecca_direction = ce.haversine(self.latitude, self.longitude, self.MECCA_LAT, self.MECCA_LONG)
        mecca_direction = ce.bound_angle_deg(mecca_direction)
        
        return {"distance" : np.round(mecca_distance, 2), "angle" : np.round(mecca_direction, 2), "cardinal" : ce.get_cardinal_direction(np.round(mecca_direction))}
    
    def sun(self):
        return {"declination" : f"{self.sun_declination:.3f}°", "right_ascension" : f"{self.sun_alpha[0]}h {self.sun_alpha[0]}m {self.sun_alpha[0]:.2f}s",
                "altitude" : f"{self.sun_alt:.2f}°", "azimuth" : f"{self.sun_az:.2f}°"}
    
    def moon(self):
        ### Calculate Moon Illumination
        moon_illumin = me.moon_illumination(self.sun_declination, self.sun_factors[10], self.sun_factors[4], 
                                            self.moon_declination, self.moon_factors[6], self.moon_factors[1], 
                                            self.moon_factors[0], self.sun_factors[6], self.moon_factors[2] / 149597870.7)

        return {"declination" : f"{self.moon_declination:.3f}°", "right_ascension" : f"{self.moon_alpha[0]}h {self.moon_alpha[1]}m {self.moon_alpha[2]:.2f}s",
                "altitude" : f"{self.moon_alt:.2f}°", "azimuth" : f"{self.moon_az:.2f}°", "illumination" : f"{moon_illumin * 100:.1f}%"}
    
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

    def visibilities(self):
        ### Calculate New Moon Visibilities
        
        # Get New Moon Date from moon_phases list
        moon_phases = self.moonphases()
        for item in moon_phases:
                if item['phase'] == "New Moon":
                    new_moon = item['datetime']

        # Find JD for the given date; adjust day for difference in UTC and local timezone
        jd_new_moon = te.gregorian_to_jd(new_moon.year, new_moon.month, new_moon.day + te.fraction_of_day(new_moon), -1 * self.utc_diff)
        if new_moon.day != te.jd_to_gregorian(jd_new_moon).day:
            if new_moon.day < te.jd_to_gregorian(jd_new_moon).day:
                new_moon += datetime.timedelta(days=1)
            else:
                new_moon -= datetime.timedelta(days=1)
            jd_new_moon = te.gregorian_to_jd(new_moon.year, new_moon.month, new_moon.day, -1 * self.utc_diff)

        # Find local sunset as visibilities are calculated from then
        nm_sun_factors = se.sunpos(jd_new_moon, self.latitude, self.longitude)
        nm_sunset = te.solar2standard(se.sunrise_sunset(1, se.solar_hour_angle(self.latitude, nm_sun_factors[11])), self.utc_diff, self.longitude, se.equation_of_time(jd_new_moon, self.latitude, self.longitude))
        jd_new_moon = te.gregorian_to_jd(new_moon.year, new_moon.month, new_moon.day + nm_sunset / 24, -1 * self.utc_diff)

        # Find visibilities for the three days
        visibilities = []
        i = 0
        while (i < 3):
            nm_sun_factors = se.sunpos(jd_new_moon + i, self.latitude, self.longitude)
            delPsi, delEps = se.sun_nutation(jd_new_moon + i)
            nm_moon_factors = me.moonpos(jd_new_moon + i, self.latitude, self.longitude, delPsi, nm_sun_factors[13], self.elev)

            #print(jd_to_gregorian(jd_new_moon + i))
            visibilities.append(me.calculate_visibility(nm_sun_factors[16], nm_sun_factors[15], nm_moon_factors[10], nm_moon_factors[9], np.deg2rad(nm_moon_factors[8])))
            i += 1

        # Arrange and classify visibilties
        q_values = [
            [visibilities[0], me.classify_visibility(visibilities[0])],
            [visibilities[1], me.classify_visibility(visibilities[1])],
            [visibilities[2], me.classify_visibility(visibilities[2])],
        ]

        return {"0" : q_values[0], "1" : q_values[1], "2" : q_values[2]}