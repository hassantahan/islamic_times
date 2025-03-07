"""
Module for calculating Islamic times and astronomical data.

This module provides the ITLocation class that encapsulates calculations for:
  - Astronomical parameters (sun and moon positions, Julian Date, etc.)
  - Prayer times based on various calculation methods.
  - Islamic calendar (Hijri) conversion.
  - New moon visibilities

References:
  - Jean Meeus, "Astronomical Algorithms"
  - Prayer times calculation methods (http://praytimes.org/wiki/Calculation_Methods)
"""

import numpy as np
from numbers import Number
from typing import Dict, List
from datetime import datetime, timedelta, timezone
from islamic_times import prayer_times as pt
from islamic_times import sun_equations as se
from islamic_times import moon_equations as me
from islamic_times import time_equations as te
from islamic_times import calculation_equations as ce

class ITLocation:
    """Represents an observer's location and time for Islamic astronomical calculations.
    """

    ##### Mecca Constants #####
    __MECCA_LAT = 21.420164986
    __MECCA_LONG = 39.822330044

    ##### Prayer Methods #####
    # http://praytimes.org/wiki/Calculation_Methods
    # fajr, isha, maghrib, midnight (0 --> sunset to sunrise, 1 --> sunset to fajr)
    __PRAYER_METHODS = {
        "Muslim World League (MWL)": (('MWL', 'MUSLIM WORLD LEAGUE'), (18, 17, 0, 0)),
        "Islamic Society of North America (ISNA)": (('ISNA', 'ISLAMIC SOCIETY OF NORTH AMERICA'),(15, 15, 0, 0)),
        "Egyptian General Authority of Survey (Egypt)": (('EGYPT', 'EGYPTIAN', 'EGYPTIAN GENERAL AUTHORITY OF SURVEY', 'EGAS'), (19.5, 17.5, 0, 0)),
        "Umm al-Qura University (Makkah)": (('MAKKAH', 'MECCA', 'MEKKAH', 'UMM AL-QURA UNIVERSITY', 'UMM AL-QURA', 'UQU'), (18.5, np.inf, 0, 0)),
        "University of Islamic Sciences, (Karachi)": (('KARACHI', 'UNIVERSITY OF ISLAMIC SCIENCES', 'UIS'), (18, 18, 0, 0)),
        "Institute of Geophysics, University of Tehran (Tehran)": (('TEHRAN', 'UNIVERSITY OF TEHRAN', 'INSTITUTE OF GEOPHYSICS', 'UOT', 'IOG', 'IOGUOT', 'UOTIOG', 'IGUT', 'UTIG'), (17.7, 14, 4.5, 1)),
        "Shia Ithna Ashari, Leva Research Institute, Qom (Jafari)": (('JAFARI', 'JAAFARI', 'SHIA', 'SHIA ITHNA ASHARI', 'LEVA', 'LEVA RESEARCH INSTITUTE', 'QOM', 'QUM', 'LRI', 'SIA', 'SIALRI'), (16, 14, 4, 1))
    }

    def __init__(self, latitude: float = 51.477928,
                 longitude: float = -0.001545,
                 elevation: float = 76,
                 temperature: float = 10,
                 pressure: float = 101.325,
                 today: datetime = datetime.now(timezone.utc),
                 method: str = 'JAFARI',
                 asr_type: int = 0,
                 find_local_tz: bool = True):
        '''Instantiates and initalizes the `ITLocation` object. The default location values are for the Royal Greenwich Observatory (https://en.wikipedia.org/wiki/Royal_Observatory,_Greenwich).
        '''
        
        # Check the numerical inputs
        float_inputs = {
            "latitude": latitude,
            "longitude": longitude,
            "elevation": elevation,
            "temperature": temperature,
            "pressure": pressure
        }

        for name, val in float_inputs.items():
            if not isinstance(val, Number):
                raise TypeError(f"'{name}' must be of type `float`, but got `{type(val).__name__}`.")
            if name == "latitude":
                if val > 90 or val < -90:
                    raise ValueError(f"Input for '{name}' is out of range. Latitudes must be between -90° and 90°.")
            elif name == "longitude":
                if val > 180 or val < -180:
                    raise ValueError(f"Input for '{name}' is out of range. Longitudes must be between -180° and 180°.")
        self.latitude = latitude
        self.longitude = longitude
        self.elevation = elevation
        self.temperature = temperature
        self.pressure = pressure

        # Check the date if it is valid
        if not isinstance(today, datetime):
            raise TypeError(f"'{today}' must be of type `datetime`, but got `{type(today).__name__}`.")
        self.today = today

        # Check if find_local_tz is either 0, 1, or a bool value
        if not isinstance(find_local_tz, bool):
            if not isinstance(find_local_tz, Number):
                raise ValueError(f"'find_local_tz' is out of range; it must be of type `bool` or either the numerical value of 0 or 1.")
            else:
                raise TypeError(f"'find_local_tz' must be of type `bool` or either the numerical value of 0 or 1, but got `{type(find_local_tz).__name__}`.")  
        self.find_local_tz = find_local_tz
        
        # Determine UTC Offset
        self.__calculate_utc_offset(find_local_tz)
        self.utc_diff *= -1

        # Calculate the astronomical parameters
        self.calculate_astro()
        
        # Check both asr and prayer method input values
        if asr_type in (0, 1):
            setattr(self, "asr_type", asr_type)
        else:
            raise ValueError(f"'asr_type' must be either 0 or 1. Check documentation to understand each type. Invalid value: {asr_type}.")
        self.asr_type = asr_type

        for group_name, (group, values) in self.__PRAYER_METHODS.items():
            if method.upper() in group:
                self.fajr_angle, self.isha_angle, self.maghrib_angle, self.midnight_type = values
                self.method = group_name.upper()
                break
        else:
            # Raise an error if no match is found
            supported_methods = [value[0][0] for value in self.__PRAYER_METHODS.values()]
            raise ValueError(f"Invalid prayer calculation method: {method}. Supported methods: {list(supported_methods)}")
        
        # Set the default prayer definitions
        self.set_prayer_method(method) # Prayer times are calculated in this method

    def __calculate_utc_offset(self, find_local_tz: bool) -> float:
        """ Determine UTC offset in hours based on location if needed.
        """
        if self.today.tzinfo == None:
            ### Find UTC Offset According to Lat/Long & adjust datetime
            # This is very computationally expensive
            if find_local_tz: 
                self.tz_name, self.utc_diff = te.find_utc_offset(self.latitude, self.longitude, self.today)
            else:
                self.tz_name, self.utc_diff = timezone.utc, 0
        else:  
            self.tz_name, self.utc_diff = self.today.tzinfo, self.today.utcoffset().total_seconds() / 3600

    # Used to change observe date & time
    # By default, updates to datetime.now() if argument is not specified
    def update_time(self, date_time: datetime = None):
        '''This method updates the observer time to either `datetime.now()` or to the specified datetime
        '''

        if not isinstance(date_time, datetime):
            raise TypeError(f"'date_time' must be of type `datetime`, but got `{type(date_time).__name__}`.")

        # Time is only updates to the latest if no argument is passed
        if date_time is None:
            self.today = datetime.now(self.today.tzinfo)
        else:
            self.today = date_time

    # Calculates the astronomical variables for the moon and sun
    def calculate_astro(self):
        '''This method updates the parameters used in astronomical calculations such as the local Julian Date, Delta T, Sun & Moon position (including alt/az), and more.
        '''

        ### Calculate Julian Date
        self.jd = te.gregorian_to_jd(self.today, -1 * self.utc_diff)
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
        self.moonset = self.__find_proper_moonset(self.today)
        self.moon_illumin = me.moon_illumination(
                                self.sun_declination, 
                                self.sun_factors[10], 
                                self.moon_declination, 
                                self.moon_factors[4], 
                                self.sun_factors[6], 
                                self.moon_factors[2] / te.ASTRONOMICAL_UNIT
                            )
        
    # This is necessary because UTC offsets for coords not near UTC, but also not using local TZ.
    def __find_proper_moonset(self, date: datetime) -> datetime:
        '''Determine the proper local moonset time for the observer.

        Adjusts the calculated moonset time to account for local UTC differences.

        Args:
            `date (datetime)`: The reference date.

        Returns:
            `datetime`: Adjusted moonset time.
        '''

        date_doy = date.timetuple().tm_yday
        temp_moonset = me.calculate_moonset(date, self.latitude, self.longitude, self.elevation, self.utc_diff)
        if temp_moonset == np.inf:
            return datetime.min

        temp_utc_diff = np.floor(self.longitude / 15)
        i = 1
        while(True):
            temp_moonset_doy = (temp_moonset + timedelta(hours=temp_utc_diff)).timetuple().tm_yday
            if (temp_moonset_doy < date_doy and temp_moonset.year == date.year) or ((temp_moonset + timedelta(hours=temp_utc_diff)).year < date.year):
                temp_moonset = me.calculate_moonset(date + timedelta(days=i), self.latitude, self.longitude, self.elevation, self.utc_diff)
                i += 1
            else: 
                return temp_moonset

    # Prayer Time Calculations
    def calculate_prayer_times(self):
        '''This method calculates the prayer times.
        '''
        self.times_of_prayer = pt.calculate_prayer_times(self.jde, self.deltaT, self.latitude, self.longitude, self.utc_diff, (self.sun_declination, self.solar_angle),
                                                        (self.fajr_angle, self.maghrib_angle, self.isha_angle, self.midnight_type, self.asr_type), 
                                                        self.islamic_date[2] == 9)
        
        def check_datetime(val: datetime | str) -> datetime | str:
            if type(val) is datetime:
                return val.strftime('%X %d-%m-%Y')
            else:
                return val
            
        self.times_of_prayer = {key: check_datetime(value) for key, value in self.times_of_prayer.items()}

    # Set the method of calculating prayer times among the available default options
    # The default option (from creation) is the Jaʿfarī method.
    def set_prayer_method(self, method: str = 'jafari'):
        '''This method sets the method of calculation for the local prayer times. By default, the method is the "Jaʿfarī" method, a.k.a. Shia Ithna Ashari, Leva Research Institute, Qom. 
        The methods are taken from http://praytimes.org/wiki/Calculation_Methods (archive: https://archive.ph/v98Pv).

        Args:
            `method (str)`: string specifying the method. See: http://praytimes.org/wiki/Calculation_Methods.
        '''

        # Find the matching group
        for group_name, (group, values) in self.__PRAYER_METHODS.items():
            if method.upper() in group:
                self.fajr_angle, self.isha_angle, self.maghrib_angle, self.midnight_type = values
                self.method = group_name
                break
        else:
            # Raise an error if no match is found
            supported_methods = [value[0][0] for value in self.__PRAYER_METHODS.values()]
            raise ValueError(f"Invalid prayer calculation method: {method}. Supported methods: {list(supported_methods)}")
        
        # Update prayer times
        self.calculate_prayer_times()
    
    # Alows user to set their own solar hour angles for prayer time calculations.
    def set_custom_prayer_angles(self, fajr_angle: float = None, maghrib_angle: float = None, isha_angle: float = None):
        '''
        ### Description
            This method allows for complete customization of the solar hour angles which are used to calculate the prayer times for fajr, maghrib, and ʿishāʾ prayers.

            If no argument is given for a parameter, the value stays the same.

            This method automatically calls `calculate_prayer_times()`.

            `self.method` is set to 'custom' after this function is called. 

            `set_prayer_method()` must be called with one of the default methods passed as an argument to reset `self.method` to one of the defaults.

        ### Args:
            `fajr_angle (float)`: solar hour angle for the fajr prayer
            `maghrib_angle (float)`: solar hour angle for the maghrib prayer
            `isha_angle (float)`: solar hour angle for the ʿishāʾ prayer

        ### Returns:
            No return.
        '''
        
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
    
    # Separated from angles since it is defined by shadow ratio
    def set_asr_type(self, asr_type: int = 0):
        '''
        ### Description:
            This method allows for customization of the type of ʿaṣr prayer time calculation to use.

            Type 0: Shadow ratio is 1:1; this is the majority method.

            Type 1: Shadow ratio is 2:1; this is the method used by the Ḥanafī school.

            `self.method` is set to 'custom' after this function is called. 

            `set_prayer_method()` must be called with one of the default methods passed as an argument to reset `self.method` to one of the defaults.

        ### Args:
            `asr_type (int)`: type of ʿaṣr prayer time calculation. 

        ### Returns:
            No return.
        '''

        if asr_type in (0, 1):
            setattr(self, "asr_type", asr_type)
        else:
            raise ValueError(f"'asr_type' must be either 0 or 1. Check documentation to understand each type. Invalid value: {asr_type}")

        # Method is now custom
        self.method = 'Custom'

        # Update prayer times
        self.calculate_prayer_times()

    # Set to either 0 (sunset to sunrise; the majority method) or 1 (sunset to fajr, the 'Jaʿfarī' method)
    def set_midnight_type(self, midnight_type: int = 0):
        '''
        ### Description:
            This method allows for customization of the Islamic midnight calculation.

            `midnight_type = 0`: midnight as the midpoint between the observer's sunset at its current date, and sunrise of the next day; the majority method.

            `midnight_type = 1`: midnight as the midpoint between the observer's sunset at its current date, and fajr of the next day; the Jaʿfarī method.

            If no argument is given for a parameter, the value stays the same.

            This method automatically calls `calculate_prayer_times()`.

            `self.method` is set to 'Custom' after this function is called. 

            `set_prayer_method()` must be called with one of the default methods passed as an argument to reset `self.method` to one of the defaults.

        ### Args:
            `midnight_type (int)`: type of Islamic midnight calculation

        ### Returns:
            No return.
        '''

        if midnight_type in (0, 1):
            setattr(self, "midnight_type", midnight_type)
        else:
            raise ValueError(f"'midnight_type' must be either 0 or 1. Check documentation to understand each type. Invalid value: {midnight_type}")

        # Method is now custom
        self.method = 'Custom'

        # Update prayer times
        self.calculate_prayer_times()

    # Return Observer Parameters
    def observer(self) -> Dict[str, float]:
        '''
        ### Description:
            Returns a dictionary of the observer parameters. 

            All values in the dictionary are the `float` type. 

        ### Keys:
            - 'latitude' (°)
            - 'longitude' (°)
            - 'elevation' (m)
            - 'pressure' (kPa)
            - 'temperature' (°C)

        ### Returns:
            `Dict[str, float]`
        '''
        return {
                "latitude" : np.round(self.latitude, 5),
                "longitude" : np.round(self.longitude, 5),
                "elevation" : np.round(self.elevation, 2),
                "pressure" : np.round(self.pressure, 3),
                "temperature" : np.round(self.temperature, 2)
            }

    # Return date and time information
    def dates_times(self) -> Dict[str, str | float]:
        '''
        ### Description:
            Returns a dictionary of the observer dates and times. 

            All values in the dictionary are the `str` or `float` type.

        ### Keys:
            - 'gregorian' (Gregorian date; %A, %d %B, %Y)
            - 'hijri' (Islamic (Hijrī) date; %A, %d %B, %Y)
            - 'time' (24-Hour time at the observer's timezone; %X)
            - 'timezone' (Observer timezone; `str`)
            - 'utc_offset' (Hours offset from UTC; str %H:%M)
            - 'jd' (Local Julian Day; `float`)
            - 'eq_of_time' (Equation of time in minutes; `float`)
            - 'deltaT' (Delta T (TT-UT) in seconds; `float`)

        ### Returns:
            `Dict[str, float]`
        '''

        return {
                "gregorian" : self.today.strftime("%A, %d %B, %Y"), 
                "hijri" : f"{te.get_islamic_day(self.today.strftime('%A'))}, {self.islamic_date[2]} {te.get_islamic_month(self.islamic_date[1])}, {self.islamic_date[0]}",
                "time" : self.today.strftime("%X"), "timezone" : self.tz_name, "utc_offset" : te.format_utc_offset(self.utc_diff * -1),
                "jd" : np.round(self.jd, 7),
                "eq_of_time" : np.round(se.equation_of_time(self.jde, self.deltaT, self.latitude, self.longitude), 2),
                "deltaT" : np.round(self.deltaT, 2)
            }

    # Return prayer times
    def prayer_times(self) -> Dict[str, str]:
        '''
        ### Description:
            Returns a dictionary of the prayer times at the observer timezone. 

            All values in the dictionary are the `str` type. 

        ### Keys:
            - 'method' (The method of calculating prayers)
            - 'fajr' (Time of fajr`)
            - 'sunrise' (Time of sunrise)
            - 'noon' (Time of ẓuhr or "solar noon" or "culmination of the sun")
            - 'asr' (Time of ʿaṣr)
            - 'sunset' (Time of sunset)
            - 'maghrib' (Time of maghrib)
            - 'isha' (Time of ʿishāʾ)
            - 'midnight' (Time of midnight)

        ### Returns:
            `Dict[str, str]`
        '''

        def check_datetime(val: datetime | str) -> datetime | str:
            if type(val) is datetime:
                return val.strftime('%X %d-%m-%Y')
            else:
                return val
        
        return {
                "method": self.method,
                "fajr" : self.times_of_prayer["fajr"],
                "sunrise" : self.times_of_prayer["sunrise"],
                "noon" : self.times_of_prayer["noon"],
                "asr" : self.times_of_prayer["asr"],
                "sunset" : self.times_of_prayer["sunset"],
                "maghrib" : self.times_of_prayer["maghrib"],
                "isha" : self.times_of_prayer["isha"],
                "midnight" : self.times_of_prayer["midnight"]
        }
    
    # Return Mecca information
    def mecca(self) -> Dict[str, float | str]:
        '''
        ### Description:
            Returns a dictionary of the variables related to Mecca. 

            All values in the dictionary are the `float` or `str` type. 

        ### Keys:
            - 'distance' (Distance from observer to the Kaʿbah in km; `float`)
            - 'angle' (Angle corresponding to shortest path to the Kaʿbah in °; `float`)
            - 'cardinal' (Cardinal Direction corresponding to the angle; `str`)

        ### Returns:
            `Dict[str, float | str]`
        '''

        mecca_distance, mecca_direction = ce.haversine(self.latitude, self.longitude, self.__MECCA_LAT, self.__MECCA_LONG)
        mecca_direction %= 360
        
        return {
                "distance" : np.round(mecca_distance, 2),
                "angle" : np.round(mecca_direction, 2),
                "cardinal" : ce.get_cardinal_direction(np.round(mecca_direction))
            }
    
    # Return sun properties and position values
    def sun(self) -> Dict[str, float | str]:
        '''
        ### Description:
            Returns a dictionary of the variables related to the position of the Sun. 

            All values in the dictionary are the `float` or `str` type. 

        ### Keys:
            - 'declination' (Declination of the sun in °; `float`)
            - 'right_ascension' (Right ascension of the sun in HMS; `str`)
            - 'altitude' (Altitude of the sun in °; `float`)
            - 'azimuth' (Azimuth of the sun in °; `float`)

        ### Args:
            No arguments.

        ### Returns:
            `Dict[str, float | str]`
        '''

        return {
                "declination" : np.round(self.sun_declination, 3),
                "right_ascension" : f"{self.sun_alpha[0]}h {self.sun_alpha[1]}m {self.sun_alpha[2]:.2f}s",
                "altitude" : np.round(self.sun_alt, 2),
                "azimuth" : np.round(self.sun_az, 2)
            }
    
    # Return moon properties and position values
    def moon(self) -> Dict[str, str | float]:
        '''
        ### Description:
            Returns a dictionary of the variables related to the position of the Moon. 

            All values in the dictionary are the `float` or `str` type. 

        ### Keys:
            - 'moonset' (Time of moonset; `str`)
            - 'declination' (Declination of the sun in °; `float`)
            - 'right_ascension' (Right ascension of the sun in HMS; `str`)
            - 'altitude' (Altitude of the sun in °; `float`)
            - 'azimuth' (Azimuth of the sun in °; `float`)
            - 'illumination' (Percentage of lunar illumination; `float` %)

        ### Returns:
            `Dict[str, str | float]`
        '''

        return {
                "moonset" : self.moonset.strftime("%X %d-%m-%Y"),
                "declination" : np.round(self.moon_declination, 3),
                "right_ascension" : f"{self.moon_alpha[0]}h {self.moon_alpha[1]}m {self.moon_alpha[2]:.2f}s",
                "altitude" : np.round(self.moon_alt, 2),
                "azimuth" : np.round(self.moon_az, 2),
                "illumination" : np.round(self.moon_illumin * 100, 2)
            }
    
    def moonphases(self) -> List[Dict[str, datetime]]:
        '''
        # Description:
            Returns a list of dictionaries of the "nearest" phases of the moon (not always the next phases necessarily. 

            The list is ordered in chronological order.

            All values in each of the dictionaries are of the `datetime` type. 

        # Keys:
            - 'phase' (Moon Phase; `str`)
            - 'datetime' (Date & Time of the given moon phase; `datetime`)

        # Returns:
            `List[Dict[str, str | float]]`
        '''

        # Find Next New Moon (and the rest of the phases)
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
            
            moon_phases[i] = {"phase": phase_str, "datetime": phase - timedelta(hours=self.utc_diff)}

        moon_phases = sorted(moon_phases, key = lambda item: item["datetime"])

        return moon_phases

    # Calculate Next New Moon Visibilities
    def visibilities(self, days: int = 3, type: int = 0) -> Dict [datetime, List[str | float]]:
        '''
        ### Description:
            Returns a dictionary describing the visibility of the [nearest in time] new moon for the observer.
        
            The size of the dictionary is controlled by `days` which specifies how many days from the new moon to look at visibilities.

            The key of each item in the dictionary corresponds to the "Best Time" `datetime` at which to look for the new moon crescent.

            The value of each item in the dictionary is a list in which the first element is the raw number output of the visibility. The second element is the classification of the first element.

            The `type` argument specifies which new moon visibility classification method to use:
            - Type 0: Odeh, 2006
            - Type 1: Yallop, 1997; a.k.a. HMNAO TN No. 69

        ### Args:
            `days (int)`: How many days from the new moon to look at visibilities.
            `type (int)`: Which method to classify visibilities.

        ### Returns:
            `List[Dict[str, str | float]]`
        '''

        if not isinstance(days, int):
            raise TypeError(f"'days' must be of type `int`, but got `{type(days).__name__}`.")
        
        if days < 1:
            raise ValueError(f"'days' must be greater than 0. Invalid value: {days}.")
        
        if not isinstance(type, int):
            raise TypeError(f"'type' must be of type `int`, but got `{type(type).__name__}`.")
        
        if type not in (0, 1):
            raise ValueError(f"'type' must be either 0 or 1. Invalid value: {type}.")
        

        # Get New Moon Date from moon_phases list
        moon_phases = self.moonphases()
        for item in moon_phases:
                if item['phase'] == "New Moon":
                    new_moon = item['datetime']

        # Find JD for the given date; adjust day for difference in UTC and local timezone
        jd_new_moon = te.gregorian_to_jd(new_moon, -1 * self.utc_diff)
        ymd_new_moon = te.jd_to_gregorian(jd_new_moon)
        deltaT_new_moon = te.delta_t_approx(ymd_new_moon.year, ymd_new_moon.month)
        jde_new_moon = jd_new_moon + deltaT_new_moon / 86400

        # Forgot what this does lol. Likely something to do with timezone differences
        if new_moon.day != ymd_new_moon.day:
            if new_moon.day < ymd_new_moon.day:
                new_moon += timedelta(days=1)
            else:
                new_moon -= timedelta(days=1)
            jd_new_moon = te.gregorian_to_jd(new_moon, -1 * self.utc_diff) - te.fraction_of_day(new_moon)

        # Find local sunset as visibilities are calculated from then
        nm_sun_factors = se.sunpos(jde_new_moon, deltaT_new_moon, self.latitude, self.longitude)     

        # Find visibilities for the three days
        visibilities = []
        best_jds = []
        for i in range(days):
            # First, check if the moonset is before the new moon for the first day
            if i == 0:
                nm_moonset = self.__find_proper_moonset(ymd_new_moon)

                if nm_moonset == datetime.min:
                    # Moonset doesn't exist
                    v = -997
                    visibilities.append(v)
                    best_jds.append(jd_new_moon)
                    continue
                elif nm_moonset < ymd_new_moon:
                    # Moon is not visibile before the new moon
                    v = -999
                    visibilities.append(v)
                    best_jds.append(te.gregorian_to_jd(nm_moonset))
                    continue

            # Set the day parameters
            test_jd_new_moon = jd_new_moon + i
            test_ymd_new_moon = te.jd_to_gregorian(test_jd_new_moon)
            test_deltaT_new_moon = te.delta_t_approx(test_ymd_new_moon.year, test_ymd_new_moon.month)
            test_jde_new_moon = test_jd_new_moon + test_deltaT_new_moon / 86400

            # Set sun parameters
            nm_sun_factors = se.sunpos(test_jde_new_moon, test_deltaT_new_moon, self.latitude, self.longitude)

            # Sunset & moonset calculations
            test_nm_sunset = te.solar2standard(
                                        test_jd_new_moon, 
                                        se.sunrise_sunset(1, se.solar_hour_angle(self.latitude, nm_sun_factors[11])), 
                                        self.utc_diff, 
                                        self.longitude,
                                        se.equation_of_time(test_jde_new_moon, test_deltaT_new_moon, self.latitude, self.longitude)
                                    )

            if i == 0:
                test_nm_moonset = nm_moonset
            else:
                test_nm_moonset = self.__find_proper_moonset(test_ymd_new_moon)

            # For extreme latitudes where the moonset or sunset don't exist:
            if np.abs(self.latitude) > 62:
                if test_nm_sunset == datetime.min and test_nm_moonset == datetime.min:
                    # Moonset and sunset don't exist
                    v = -997
                    visibilities.append(v)
                    best_jds.append(test_jd_new_moon)
                    continue
                elif test_nm_sunset == datetime.min:
                    # Only sunset doesn't exist
                    v = -996
                    visibilities.append(v)
                    best_jds.append(te.gregorian_to_jd(test_nm_moonset))
                    continue
                elif test_nm_moonset == datetime.min:
                    # Only moonset doesn't exist
                    v = -995
                    visibilities.append(v)
                    best_jds.append(te.gregorian_to_jd(test_nm_sunset))
                    continue
            
            # If moonset is before sunset, continue
            if test_nm_moonset < test_nm_sunset:
                v = -998
                visibilities.append(v)
                best_jds.append(te.gregorian_to_jd(test_nm_moonset))
                continue

            # Find the best time which is four ninths the moonset-sunset lag after sunset 
            lag = (test_nm_moonset - test_nm_sunset).total_seconds() / 3600
            best_time = test_nm_sunset + timedelta(hours=4 / 9 * lag)
            best_time_jd = te.gregorian_to_jd(best_time, -1 * self.utc_diff)
            best_time_jde = best_time_jd + test_deltaT_new_moon / 86400
            best_jds.append(best_time_jd)

            # Recalculate sun & calculate moon parameters
            nm_sun_factors = se.sunpos(best_time_jde, test_deltaT_new_moon, self.latitude, self.longitude)
            delPsi, delEps = se.sun_nutation(best_time_jde)
            nm_moon_factors = me.moonpos(best_time_jde, test_deltaT_new_moon, self.latitude, self.longitude, delPsi, nm_sun_factors[13], self.elevation)

            # Visibility is now calculated
            v = me.calculate_visibility(nm_sun_factors[16], nm_sun_factors[15], nm_moon_factors[10], nm_moon_factors[9], np.deg2rad(nm_moon_factors[8]), type)

            visibilities.append(v)

        # Arrange and classify visibilties
        q_values = [
            [visibility, me.classify_visibility(visibility, type)] 
            for visibility in visibilities
        ]

        # Convert best times from JD to datetime
        best_dates = [
            te.jd_to_gregorian(jd, self.utc_diff)
            for jd in best_jds
        ]

        # Label each q_value to its associated date 
        visibility_dictionary = {
            dt.strftime('%Y-%m-%d %X'): q_values[i]
            for i, dt in enumerate(best_dates)
        }

        return visibility_dictionary