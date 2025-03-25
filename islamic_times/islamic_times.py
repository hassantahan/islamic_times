"""
Module for calculating Islamic times and astronomical data.

This module provides the ITLocation class that encapsulates calculations for:
  - Astronomical parameters (sun and moon positions, Julian Date, etc.)
  - Prayer times based on various calculation methods.
  - Islamic calendar (Hijri) conversion.
  - New moon visibilities

References:
  - Jean Meeus, *Astronomical Algorithms*, 2nd Edition, Willmann-Bell, Inc., 1998.
  - Prayer times calculation methods (http://praytimes.org/wiki/Calculation_Methods)
"""

import numpy as np
from numbers import Number
from typing import Dict, List, Tuple
from datetime import datetime, timedelta, timezone
from islamic_times import prayer_times as pt
from islamic_times import sun_equations as se
from islamic_times import moon_equations as me
from islamic_times import time_equations as te
from islamic_times import calculation_equations as ce

class ITLocation:
    """
    Represents an observer's location and time for Islamic astronomical calculations.

    This class provides methods for computing astronomical parameters and Islamic prayer times 
    based on various calculation methods. It supports Hijri date conversion, new moon visibility 
    calculations, and Qibla direction determination.

    ### Features:
    - Calculates prayer times using multiple predefined methods.
    - Computes astronomical parameters (Sun/Moon position, Julian Date, Delta T, etc.).
    - Supports custom prayer time calculations by modifying solar hour angles.
    - Determines Qibla direction and distance to Mecca.
    - Converts Gregorian dates to Hijri (Islamic calendar).
    - Estimates new moon visibility and moon phases.
    - Supports both automatic and manual calculation modes.

    ### Attributes:
    - `latitude` (float): Observer's latitude (-90 to 90 degrees).
    - `longitude` (float): Observer's longitude (-180 to 180 degrees).
    - `elevation` (float): Elevation above sea level (meters).
    - `temperature` (float): Local temperature (°C).
    - `pressure` (float): Atmospheric pressure (kPa).
    - `today` (datetime): The observer's date and time.
    - `find_local_tz` (bool): Whether to automatically determine the time zone.
    - `utc_offset` (float): UTC offset in hours.
    - `auto_calculate` (bool): Whether astronomical calculations are performed automatically.
    - `asr_type` (int): Type of Asr calculation method (0 = standard, 1 = Hanafi).
    - `method` (str): The selected prayer calculation method.
    - `fajr_angle` (float): Solar angle for Fajr prayer.
    - `isha_angle` (float): Solar angle for Isha prayer.
    - `maghrib_angle` (float): Solar angle for Maghrib prayer.
    - `midnight_type` (int): Midnight calculation method (0 = sunset to sunrise, 1 = sunset to Fajr).
    - `datetime_modified` (bool): Indicates if the date/time has been manually modified.
    - `prayers_modified` (bool): Indicates if prayer times have been manually modified.

    ### Methods:
    - `update_time(date_time)`: Updates the observer's time.
    - `calculate_astro()`: Computes astronomical parameters.
    - `calculate_prayer_times()`: Computes prayer times based on the selected method.
    - `set_prayer_method(method)`: Sets the prayer time calculation method.
    - `set_custom_prayer_angles(fajr_angle, maghrib_angle, isha_angle)`: Customizes solar angles for prayer time calculations.
    - `set_asr_type(asr_type)`: Sets the Asr prayer calculation method.
    - `set_midnight_type(midnight_type)`: Sets the Islamic midnight calculation method.
    - `observer()`: Returns observer location parameters.
    - `dates_times()`: Returns observer’s date and time details.
    - `prayer_times()`: Returns calculated prayer times.
    - `mecca()`: Returns observer's distance and direction to Mecca.
    - `sun()`: Returns properties and position of the Sun.
    - `moon()`: Returns properties and position of the Moon.
    - `moonphases()`: Returns the nearest moon phases.

    ### References:
    - Jean Meeus, "Astronomical Algorithms"
    - Prayer times calculation methods: [PrayTimes Wiki](http://praytimes.org/wiki/Calculation_Methods)
    """

    # TODO
    __slots__ = (
        "observer_latitude", "observer_longitude", "observer_elevation", "temperature", "pressure", "observer_date", "find_local_tz", "auto_calculate", "datetime_modified", "prayers_modified", "tz_name", "utc_offset", 
        "asr_type", "fajr_angle", "isha_angle", "maghrib_angle", "midnight_type", "method", "times_of_prayer", "islamic_date", 
        "jd", "deltaT", "jde", 
        "sun_params", "delPsi", "sun_declination", "sun_right_ascension", "sunrise", "sun_transit", "sunset", "solar_angle", "sun_alt", "sun_az", 
        "moon_params", "moon_declination", "moon_right_ascension", "moon_pi", "moon_alt", "moon_az", "moonrise", "moon_transit", "moonset", "moon_illumin"
    )

    ##### Prayer Methods #####
    # http://praytimes.org/wiki/Calculation_Methods
    # fajr, isha, maghrib, midnight (0 --> sunset to sunrise, 1 --> sunset to fajr)
    PrayerNames = Tuple[str, ...]  # A tuple of one or more prayer names
    PrayerAngles = Tuple[int, int, int, int]  # A fixed tuple of four integers
    __PRAYER_METHODS: Dict[str , Tuple[PrayerNames, PrayerAngles]] = {
        "Muslim World League (MWL)": (('MWL', 'MUSLIM WORLD LEAGUE'), (18, 17, 0, 0)),
        "Islamic Society of North America (ISNA)": (('ISNA', 'ISLAMIC SOCIETY OF NORTH AMERICA'),(15, 15, 0, 0)),
        "Egyptian General Authority of Survey (Egypt)": (('EGYPT', 'EGYPTIAN', 'EGYPTIAN GENERAL AUTHORITY OF SURVEY', 'EGAS'), (19.5, 17.5, 0, 0)),
        "Umm al-Qura University (Makkah)": (('MAKKAH', 'MECCA', 'MEKKAH', 'UMM AL-QURA UNIVERSITY', 'UMM AL-QURA', 'UQU'), (18.5, np.inf, 0, 0)),
        "University of Islamic Sciences, (Karachi)": (('KARACHI', 'UNIVERSITY OF ISLAMIC SCIENCES', 'UIS'), (18, 18, 0, 0)),
        "Institute of Geophysics, University of Tehran (Tehran)": (('TEHRAN', 'UNIVERSITY OF TEHRAN', 'INSTITUTE OF GEOPHYSICS', 'UOT', 'IOG', 'IOGUOT', 'UOTIOG', 'IGUT', 'UTIG'), (17.7, 14, 4.5, 1)),
        "Shia Ithna Ashari, Leva Research Institute, Qom (Jafari)": (('JAFARI', 'JAAFARI', 'SHIA', 'SHIA ITHNA ASHARI', 'LEVA', 'LEVA RESEARCH INSTITUTE', 'QOM', 'QUM', 'LRI', 'SIA', 'SIALRI'), (16, 14, 4, 1))
    }

    """
    ### `__PRAYER_METHODS` Dictionary

    This dictionary maps prayer time calculation methods to their respective identifiers
    and angle values.

    ### Structure:
    `__PRAYER_METHODS` is a dictionary where:
    
    - Key (`str`): The name of the prayer calculation method.
    - Value (`Tuple[PrayerNames, PrayerAngles]`):
        • PrayerNames (`Tuple[str, ...]`): A tuple of one or more method identifiers.
        • PrayerAngles (`Tuple[int, int, int, int]`): A fixed tuple of four numerical values:
            \t1. Fajr angle (degrees)
            \t2. Isha angle (degrees)
            \t3. Maghrib angle (degrees)
            \t4. Midnight type (0 → sunset to sunrise, 1 → sunset to Fajr)

    ### Example:
    ```python
    {
        "Muslim World League (MWL)": (("MWL", "MUSLIM WORLD LEAGUE"), (18, 17, 0, 0)),
        "Islamic Society of North America (ISNA)": (("ISNA", "ISLAMIC SOCIETY OF NORTH AMERICA"), (15, 15, 0, 0)),
    }
    ```

    ### Reference:
        http://praytimes.org/wiki/Calculation_Methods
    """

    def __setattr__(self, key, value):
        if key in self.__slots__:
            raise AttributeError(f"'{key}' is immutable. Modify attributes using class methods.")
        else:
            raise AttributeError(f"Cannot add new attribute '{key}' to immutable class.")

    def __init__(self, latitude: float = 51.477928,
                 longitude: float = -0.001545,
                 elevation: float = 76,
                 temperature: float = 10,
                 pressure: float = 101.325,
                 date: datetime = datetime.now(timezone.utc),
                 method: str = 'JAFARI',
                 asr_type: int = 0,
                 find_local_tz: bool = False,
                 auto_calculate: bool = True
                 ) -> None:
        """
        Initializes an ITLocation instance with geolocation and astronomical parameters.

        This constructor sets up the ITLocation object, which is used to calculate astronomical 
        parameters and Islamic prayer times. The default location is the Royal Greenwich Observatory.

        Parameters:
            latitude (float, optional): Geographical latitude in decimal degrees (-90 to 90). Defaults to 51.477928.
            longitude (float, optional): Geographical longitude in decimal degrees (-180 to 180). Defaults to -0.001545.
            elevation (float, optional): Elevation above sea level in meters. Defaults to 76.
            temperature (float, optional): Temperature in degrees Celsius. Defaults to 10.
            pressure (float, optional): Atmospheric pressure in kPa. Defaults to 101.325.
            today (datetime, optional): Current date and time (UTC). Defaults to `datetime.now(timezone.utc)`.
            method (str, optional): Prayer calculation method (e.g., 'JAFARI', 'ISNA', etc.). Defaults to 'JAFARI'.
            asr_type (int, optional): Asr calculation type (0 for standard, 1 for Hanafi). Defaults to 0.
            find_local_tz (bool, optional): Whether to determine the local time zone automatically. Defaults to True.
            auto_calculate (bool, optional): Whether to compute astronomical parameters upon initialization. Defaults to True.

        Raises:
            TypeError: If latitude, longitude, elevation, temperature, or pressure are not floats.
            ValueError: If latitude or longitude values are out of range.
            TypeError: If `today` is not a `datetime` object.
            TypeError/ValueError: If `find_local_tz` is not a boolean or a valid numerical representation.
            ValueError: If `asr_type` is not 0 or 1.
            ValueError: If the provided prayer calculation method is invalid.
            
        Notes:
            - **DO NOT UNDER ANY CIRCUMSTANCES MANUALLY CHANGE THE `auto_calculate` PARAMETER**
            - If `auto_calculate` is True, `calculate_astro()` is called to compute astronomical parameters. Otherwise, it must be called manually.
            - The selected prayer calculation method must be among the supported ones, otherwise, an error is raised.
        """
        
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
        super().__setattr__('observer_latitude', latitude)
        super().__setattr__('observer_longitude', longitude)
        super().__setattr__('observer_elevation', elevation)
        super().__setattr__('temperature', temperature)
        super().__setattr__('pressure', pressure)

        # Check the date if it is valid
        if not isinstance(date, datetime):
            raise TypeError(f"'{date}' must be of type `datetime`, but got `{type(date).__name__}`.")
        super().__setattr__('observer_date', date)

        # Check if find_local_tz is either 0, 1, or a bool value
        if not isinstance(find_local_tz, bool):
            if not isinstance(find_local_tz, Number):
                raise ValueError(f"'find_local_tz' is out of range; it must be of type `bool` or either the numerical value of 0 or 1.")
            else:
                raise TypeError(f"'find_local_tz' must be of type `bool` or either the numerical value of 0 or 1, but got `{type(find_local_tz).__name__}`.")  
        super().__setattr__('find_local_tz', find_local_tz)
        
        # Determine UTC Offset
        self.__calculate_utc_offset(find_local_tz)
        super().__setattr__('utc_offset',  self.utc_offset* -1)

        # Autocalculation for astronomical parameters
        super().__setattr__('auto_calculate', auto_calculate)
        if self.auto_calculate:
            self.calculate_astro()
            super().__setattr__('datetime_modified', False)
            super().__setattr__('prayers_modified', False)
        else:
            super().__setattr__('datetime_modified', True)
            super().__setattr__('prayers_modified', True)
        
        # Check both ʿasr and prayer method input values
        if asr_type in (0, 1):
            super().__setattr__('asr_type', asr_type)
        else:
            raise ValueError(f"'asr_type' must be either 0 or 1. Check documentation to understand each type. Invalid value: {asr_type}.")
        for group_name, (group, values) in self.__PRAYER_METHODS.items():
            if method.upper() in group:
                fajr_angle, isha_angle, maghrib_angle, midnight_type = map(lambda x: x[0](x[1]), zip([float, float, float, int], values))
                super().__setattr__('fajr_angle', fajr_angle)
                super().__setattr__('isha_angle', isha_angle)
                super().__setattr__('maghrib_angle', maghrib_angle)
                super().__setattr__('midnight_type', midnight_type)
                super().__setattr__('method', group_name.upper())
                break
        else:
            # Raise an error if no match is found
            supported_methods: str = [value[0][0] for value in self.__PRAYER_METHODS.values()]
            raise ValueError(f"Invalid prayer calculation method: {method}. Supported methods: {list(supported_methods)}")
        
        # Set the default prayer definitions
        self.set_prayer_method(method) # This also calculates the prayer times if `auto_calculate` is True

    def __calculate_utc_offset(self, find_local_tz: bool):
        """ Determine UTC offset in hours based on location if needed.

        Parameters:
            find_local_tz (bool): Controls whether or not to use the timezonefinder library to fine the timezone of the observer
        """
        # Find UTC Offset According to Lat/Long datetime
        # This is very computationally expensive
        if find_local_tz: 
            tz_name, utc_offset = te.find_utc_offset(self.observer_latitude, self.observer_longitude, self.observer_date)
            super().__setattr__('tz_name', tz_name)
            super().__setattr__('utc_offset', utc_offset)
        elif self.observer_date.tzinfo == None or self.observer_date.tzinfo == timezone.utc:
            super().__setattr__('tz_name', timezone.utc)
            super().__setattr__('utc_offset', 0)
        else:  
            super().__setattr__('tz_name', self.observer_date.tzinfo)
            super().__setattr__('utc_offset', self.observer_date.utcoffset().total_seconds() / 3600)

    # Used to change observe date & time
    # By default, updates to datetime.now() if argument is not specified
    def update_time(self, date_time: datetime = None) -> None:
        """
        Updates the observer's time.

        This method updates the observer's time to either `datetime.now()` (if no argument is provided) 
        or to a specified `datetime` object. After updating the time, `update_astro()` must be called 
        to recalculate astronomical parameters.

        Parameters:
            date_time (datetime, optional): The new observer time. Defaults to the current time.

        Raises:
            TypeError: If `date_time` is not a `datetime` object.
        """

        if not isinstance(date_time, datetime):
            raise TypeError(f"'date_time' must be of type `datetime`, but got `{type(date_time).__name__}`.")

        # Time is only updates to the latest if no argument is passed
        if date_time is None:
            super().__setattr__('observer_date', datetime.now(self.observer_date.tzinfo))
        else:
            super().__setattr__('observer_date', date_time)
        
        # Set bools for astro_calculation stuff
        if not self.auto_calculate:
            super().__setattr__('datetime_modified', True)
            super().__setattr__('prayers_modified', True)
        else:
            self.calculate_astro()

    # Calculates the astronomical variables for the moon and sun
    def calculate_astro(self) -> None:
        """
        Calculates astronomical parameters.

        This method computes the parameters used in astronomical calculations, such as:
          - Local Julian Date (JD)
          - Delta T (ΔT = TT - UT)
          - Sun and Moon positions (including altitude, azimuth, declination, etc.)
          - Estimated Islamic (Hijri) date
          - Moon illumination percentage

        The computed values are stored in instance attributes for later use in prayer time calculations.

        Notes:
            - If `auto_calculate` is enabled, this method is called automatically when needed. Otherwise, it must be called manually after changing the observer's time or when initializing the object.
            - `datetime_modified` is set to `False` after execution.
        """

        ### Calculate Julian Date
        super().__setattr__('jd', te.gregorian_to_jd(self.observer_date, -1 * self.utc_offset))
        super().__setattr__('deltaT', te.delta_t_approx(self.observer_date.year, self.observer_date.month))
        super().__setattr__('jde', self.jd + self.deltaT / 86400)

        ### Calculate Current Islamic Date (estimate)
        # TODO: Look into newer versions of this, see if it can be corrected. See function for more info.
        super().__setattr__('islamic_date', te.gregorian_to_hijri(self.observer_date.year, self.observer_date.month, self.observer_date.day))

        ### Sun & Moon Properties Calculations
        # Get Sun and Moon Objects with their parameters
        temp_sun: se.Sun = se.sunpos(self.jde, self.deltaT, self.observer_latitude, self.observer_longitude, self.temperature, self.pressure)
        super().__setattr__('sun_params', temp_sun)

        temp_moon: me.Moon = me.moonpos(self.jde, self.deltaT ,self.observer_latitude, self.observer_longitude, temp_sun.nutation[0], temp_sun.true_obliquity, self.observer_elevation)
        super().__setattr__('moon_params', temp_moon)

        # Important Sun Factors placed into local variables
        super().__setattr__('sun_declination', temp_sun.apparent_declination)
        super().__setattr__('sun_right_ascension', ce.decimal_to_hms(temp_sun.apparent_right_ascension))
        super().__setattr__('solar_angle', se.solar_hour_angle(self.observer_latitude, self.sun_declination))
        super().__setattr__('sun_alt', temp_sun.altitude)
        super().__setattr__('sun_az', temp_sun.azimuth)

        # Important sun times
        super().__setattr__('sunrise', se.find_proper_suntime(self.observer_date, self.observer_latitude, self.observer_longitude, self.observer_elevation, self.utc_offset, 'rise'))
        super().__setattr__('sun_transit', se.find_sun_transit(self.observer_date, self.observer_latitude, self.observer_longitude, self.observer_elevation, self.utc_offset))
        super().__setattr__('sunset', se.find_proper_suntime(self.observer_date, self.observer_latitude, self.observer_longitude, self.observer_elevation, self.utc_offset, 'set'))

        # Important Moon Factors placed into local variables
        super().__setattr__('moon_declination', temp_moon.top_declination)
        super().__setattr__('moon_right_ascension', ce.decimal_to_hms(temp_moon.right_ascension))
        super().__setattr__('moon_pi', temp_moon.eh_parallax)
        super().__setattr__('moon_alt', temp_moon.altitude)
        super().__setattr__('moon_az', temp_moon.azimuth)

        # Moon calculations
        super().__setattr__('moonrise', me.find_proper_moontime(self.observer_date, self.observer_latitude, self.observer_longitude, self.observer_elevation, self.utc_offset, 'rise'))
        super().__setattr__('moon_transit', me.find_moon_transit(self.observer_date, self.observer_latitude, self.observer_longitude, self.observer_elevation, self.utc_offset))
        super().__setattr__('moonset', me.find_proper_moontime(self.observer_date, self.observer_latitude, self.observer_longitude, self.observer_elevation, self.utc_offset, 'set'))
        super().__setattr__('moon_illumin', me.moon_illumination(
                                self.sun_declination, 
                                ce.hms_to_decimal(self.sun_right_ascension), 
                                self.moon_declination, 
                                temp_moon.right_ascension, 
                                temp_sun.geocentric_distance, 
                                temp_moon.geocentric_distance / te.ASTRONOMICAL_UNIT
                            ))
        
        # Astronomical parameters have been calculated so the flag is set to False
        super().__setattr__('datetime_modified', False)

    # Prayer Time Calculations
    def calculate_prayer_times(self) -> None:
        """
        Computes prayer times for the observer.

        This method calculates the times for:
          - Fajr
          - Sunrise
          - Ẓuhr (solar noon)
          - ʿAṣr
          - Maghrib
          - ʿIshāʾ
          - Islamic Midnight

        Uses the selected calculation method and astronomical parameters.

        Raises:
            ValueError: If `auto_calculate` is disabled and astronomical parameters have not been computed.
        """
        can_calculate = self.auto_calculate or not self.datetime_modified
        if not can_calculate:
            raise ValueError("Since auto_calculate has been set to false, prayer times cannot be calculated since astronomical parameters have not been calculated. Call 'cFalculate_astro()' first.")

        super().__setattr__('times_of_prayer', pt.calculate_prayer_times(self.observer_date, self.observer_latitude, self.observer_longitude, self.observer_elevation, self.utc_offset, self.sun_declination,
                                                                        prayer_times_parameters=(self.fajr_angle, self.maghrib_angle, self.isha_angle, self.midnight_type, self.asr_type), 
                                                                        sunset=self.sunset, sun_transit=self.sun_transit, sunrise=self.sunrise,
                                                                        is_ramadan=self.islamic_date[2] == 9
                                                                        ))
        
        def check_datetime(val: datetime | str) -> datetime | str:
            if type(val) is datetime:
                return val.strftime('%X %d-%m-%Y')
            else:
                return val
            
        super().__setattr__('times_of_prayer', {key: check_datetime(value) for key, value in self.times_of_prayer.items()})

        # Prayers have been calculated so the flag is set to False
        super().__setattr__('prayers_modified', False)

    # Set the method of calculating prayer times among the available default options
    # The default option (from creation) is the Jaʿfarī method.
    def set_prayer_method(self, method: str = 'jafari') -> None:
        """
        Sets the prayer time calculation method.

        The available methods follow those documented at:
        http://praytimes.org/wiki/Calculation_Methods.

        Parameters:
            method (str): The name of the prayer calculation method (e.g., 'JAFARI', 'ISNA', etc.). Defaults to 'JAFARI'.

        Raises:
            ValueError: If the method is not among the supported options.

        Notes:
            - If `auto_calculate` is enabled, prayer times are recalculated automatically. Otherwise, `calculate_prayer_times()` must be called.
            - To revert to a default method after using `set_custom_prayer_angles()`, this method must be called again.
        """

        # Find the matching group
        for group_name, (group, values) in self.__PRAYER_METHODS.items():
            if method.upper() in group:
                f, i, m, d = values
                super().__setattr__('fajr_angle', f)
                super().__setattr__('isha_angle', i)
                super().__setattr__('maghrib_angle', m) 
                super().__setattr__('midnight_type', d)
                super().__setattr__('method', group_name)
                break
        else:
            # Raise an error if no match is found
            supported_methods = [value[0][0] for value in self.__PRAYER_METHODS.values()]
            raise ValueError(f"Invalid prayer calculation method: {method}. Supported methods: {list(supported_methods)}")
        
        # Update prayer times if auto_calculate enabled
        if self.auto_calculate:
            self.calculate_prayer_times()
        else:
            super().__setattr__('prayers_modified', True)
    
    # Helper function to validate and set angle
    def __validate_and_set(self, attribute_name: str, value: float | int) -> None:
        if value is not None:
            if isinstance(value, (int, float)):  # Check if it's a number
                if value > 0:
                    super().__setattr__(f"{attribute_name}", value)
                else:
                    ValueError(f"{attribute_name} must be greater than 0. Invalid value: {value}")
            else:
                raise ValueError(f"{attribute_name} must be a number. Invalid value: {value}")

    # Alows user to set their own solar hour angles for prayer time calculations.
    def set_custom_prayer_angles(self, fajr_angle: float = None, maghrib_angle: float = None, isha_angle: float = None) -> None:
        """
        Customizes solar hour angles for prayer time calculations.

        This allows users to manually set the angles used to calculate Fajr, Maghrib, and Isha prayers.

        Parameters:
            fajr_angle (float, optional): Solar hour angle for Fajr.
            maghrib_angle (float, optional): Solar hour angle for Maghrib.
            isha_angle (float, optional): Solar hour angle for Isha.

        Raises:
            ValueError: If any provided angle is not a positive number.

        Notes:
            - If `auto_calculate` is enabled, prayer times are recalculated automatically. Otherwise, `calculate_prayer_times()` must be called.
            - `self.method` is set to 'Custom' after calling this method.
            - Call `set_prayer_method()` to reset to a predefined method.
        """

        # Validate and set each angle
        self.__validate_and_set('fajr_angle', fajr_angle)
        self.__validate_and_set('maghrib_angle', maghrib_angle)
        self.__validate_and_set('isha_angle', isha_angle)

        # Method is now custom
        super().__setattr__('method', 'Custom')

        # Update prayer times
        if self.auto_calculate:
            self.calculate_prayer_times()
        else:
            super().__setattr__('prayers_modified', True)
    
    # Separated from angles since it is defined by shadow ratio
    def set_asr_type(self, asr_type: int = 0) -> None:
        """
        Sets the calculation method for Asr prayer.

        Options:
            - `0`: Shadow ratio of 1:1 (majority method).
            - `1`: Shadow ratio of 2:1 (Hanafi method).

        Parameters:
            asr_type (int): The Asr calculation type.

        Raises:
            ValueError: If `asr_type` is not 0 or 1.

        Notes:
            - If `auto_calculate` is enabled, prayer times are recalculated automatically. Otherwise, `calculate_prayer_times()` must be called.
            - `self.method` is set to 'Custom' after calling this method.
            - Call `set_prayer_method()` to reset to a predefined method.
        """

        if asr_type in (0, 1):
            super().__setattr__("asr_type", asr_type)
        else:
            raise ValueError(f"'asr_type' must be either 0 or 1. Check documentation to understand each type. Invalid value: {asr_type}")

        # Method is now custom
        super().__setattr__('method', 'Custom')

        # Update prayer times
        if self.auto_calculate:
            self.calculate_prayer_times()
        else:
            super().__setattr__('prayers_modified', True)

    # Set to either 0 (sunset to sunrise; the majority method) or 1 (sunset to fajr, the 'Jaʿfarī' method)
    def set_midnight_type(self, midnight_type: int = 0) -> None:
        """
        Sets the calculation method for Islamic midnight.

        Options:
            - `0`: Midpoint between sunset and sunrise (majority method).
            - `1`: Midpoint between sunset and Fajr (Jaʿfarī method).

        Parameters:
            midnight_type (int): The midnight calculation type.

        Raises:
            ValueError: If `midnight_type` is not 0 or 1.

        Notes:
            - If `auto_calculate` is enabled, prayer times are recalculated automatically. Otherwise, `calculate_prayer_times()` must be called.
            - `self.method` is set to 'Custom' after calling this method.
            - Call `set_prayer_method()` to reset to a predefined method.
        """

        if midnight_type in (0, 1):
           super().__setattr__("midnight_type", midnight_type)
        else:
            raise ValueError(f"'midnight_type' must be either 0 or 1. Check documentation to understand each type. Invalid value: {midnight_type}")

        # Method is now custom
        super().__setattr__('method', 'Custom')

        # Update prayer times
        if self.auto_calculate:
            self.calculate_prayer_times()
        else:
            super().__setattr__('prayers_modified', True)

    # Return Observer Parameters
    def observer(self) -> Dict[str, float]:
        """
        Returns observer's date and time information.

        Returns:
            dict: 
            Dictionary containing:
                - 'latitude': Observer latitude (°).
                - 'longitude': Observer longitude (°).
                - 'elevation': Observer elevation above sea level (m).
                - 'pressure': Observer ambient pressure (kPa).
                - 'temperature': Observer ambient temperature (℃).
        """
        return {
                "latitude" : round(self.observer_latitude, 5),
                "longitude" : round(self.observer_longitude, 5),
                "elevation" : round(self.observer_elevation, 2),
                "pressure" : round(self.pressure, 3),
                "temperature" : round(self.temperature, 2)
            }

    # Return date and time information
    def dates_times(self) -> Dict[str, str | float]:
        """
        Returns observer's date and time information.

        Returns:
            dict: Dictionary containing:
                - 'gregorian': Gregorian date (str)
                - 'hijri': Islamic (Hijri) date (str)
                - 'time': Local time (str)
                - 'timezone': Time zone (str)
                - 'utc_offset': UTC offset (str)
                - 'jd': Julian Date (float)
                - 'eq_of_time': Equation of time (float)
                - 'deltaT': Delta T in seconds (float)
        """

        can_print = self.auto_calculate or not self.datetime_modified

        if not can_print:
            raise ValueError("Cannot print dates and times without calculating the astronomical parameters. Set `auto_calculate` to `True` or call `calculate_astro()`.")

        return {
                "gregorian" : self.observer_date.strftime("%A, %d %B, %Y"), 
                "hijri" : f"{te.get_islamic_day(self.observer_date.strftime('%A'))}, {self.islamic_date[2]} {te.get_islamic_month(self.islamic_date[1])}, {self.islamic_date[0]}",
                "time" : self.observer_date.strftime("%X"), "timezone" : self.tz_name, "utc_offset" : te.format_utc_offset(self.utc_offset * -1),
                "jd" : round(self.jd, 7),
                "eq_of_time" : round(se.equation_of_time(self.sun_params.delta_obliquity, self.sun_params.mean_longitude, self.sun_params.true_obliquity, self.sun_params.apparent_right_ascension), 2),
                "deltaT" : round(self.deltaT, 2)
            }

    # Return prayer times
    def prayer_times(self) -> Dict[str, str]:
        """
        Returns calculated prayer times.

        Returns:
            dict: Dictionary containing:
                - 'method'
                - 'fajr'
                - 'sunrise'
                - 'zuhr'
                - 'asr'
                - 'sunset'
                - 'maghrib'
                - 'isha'
                - 'midnight'
        """

        can_print: bool = self.auto_calculate or not self.prayers_modified

        if not can_print:
            raise ValueError("`auto_calculate` has been set to False and a change to the prayer methods and/or the date & time may have been made. `calculate_prayer_times()` must be called first.")
        
        return {
                "method": self.method,
                "fajr" : self.times_of_prayer["fajr"],
                "sunrise" : self.times_of_prayer["sunrise"],
                "zuhr" : self.times_of_prayer["zuhr"],
                "asr" : self.times_of_prayer["asr"],
                "sunset" : self.times_of_prayer["sunset"],
                "maghrib" : self.times_of_prayer["maghrib"],
                "isha" : self.times_of_prayer["isha"],
                "midnight" : self.times_of_prayer["midnight"]
        }
    
    # Return Mecca information
    def mecca(self) -> Dict[str, float | str]:
        """
        Returns observer's distance and direction to Mecca.

        Returns:
            dict: Dictionary containing:
                - 'distance': Distance to Mecca (km)
                - 'angle': Qibla direction (°)
                - 'cardinal': Cardinal direction (str)
        """

        mecca_distance, mecca_direction = ce.haversine(self.observer_latitude, self.observer_longitude, te.MECCA_LAT, te.MECCA_LONG)
        mecca_direction %= 360
        
        return {
                "distance" : round(mecca_distance, 2),
                "angle" : round(mecca_direction, 2),
                "cardinal" : ce.get_cardinal_direction(round(mecca_direction))
            }
    
    # Return sun properties and position values
    def sun(self) -> Dict[str, float | str]:
        """
        Returns properties and position of the Sun.

        Returns:
            dict (Dict[str, float | str]): Dictionary containing:
                - 'sunset': Sunrise time (str)
                - 'sun_transit': Sun transit (culmination) time (str)
                - 'sunset': Sunset time (str)
                - 'declination': The angular distance of the Sun perpendicular to the celestial equator using the true equinox of date (°)
                - 'right_ascension': The angular distance of the Sun eastward along the celestial equator from the vernal equinox to the hour circle passing through the Sun using the true equinox of date (HMS format)
                - 'altitude': The angle between the Sun and the observer's local horizon accounting for atmospheric refraction (°)
                - 'azimuth': The angle of the Sun around the horizon measured from true north (°)
        """

        can_print = self.auto_calculate or not self.datetime_modified

        if not can_print:
            raise ValueError("Cannot print dates and times without calculating the astronomical parameters. First call `calculate_astro()`.")

        return {
                "sunrise" : self.sunrise.strftime("%X %d-%m-%Y"),
                "sun_transit" :self.sun_transit.strftime("%X %d-%m-%Y"),
                "sunset" : self.sunset.strftime("%X %d-%m-%Y"),
                "declination" : round(self.sun_declination, 3),
                "right_ascension" : f"{self.sun_right_ascension[0]}h {self.sun_right_ascension[1]}m {self.sun_right_ascension[2]:.2f}s",
                "altitude" : round(self.sun_alt, 3),
                "azimuth" : round(self.sun_az, 3)
            }
    
    # Return moon properties and position values
    def moon(self) -> Dict[str, str | float]:
        """
        Returns properties and position of the Moon.

        Returns:
            dict: Dictionary containing:
                - 'moonset': Moonrise time (str)
                - 'moon_transit': Moon transit (culmination) time (str)
                - 'moonset': Moonset time (str)
                - 'declination': The angular distance of the Moon perpendicular to the celestial equator using the true equinox of date (°)
                - 'right_ascension': The angular distance of the Moon eastward along the celestial equator from the vernal equinox to the hour circle passing through the Moon using the true equinox of date (HMS format)
                - 'altitude': The angle between the Moon and the observer's local horizon accounting for atmospheric refraction (°)
                - 'azimuth': The angle of the Moon around the horizon measured from true north (°)
                - 'parallax': Equitorial horizontal parallax (′)
                - 'illumination': Percentage of illumination (%)
        """

        can_print = self.auto_calculate or not self.datetime_modified

        if not can_print:
            raise ValueError("Cannot print dates and times without calculating the astronomical parameters. First call `calculate_astro()`.")

        return {
                "moonrise" : self.moonrise.strftime("%X %d-%m-%Y"),
                "moon_transit" : self.moon_transit.strftime("%X %d-%m-%Y"),
                "moonset" : self.moonset.strftime("%X %d-%m-%Y"),
                "declination" : round(self.moon_declination, 3),
                "right_ascension" : f"{self.moon_right_ascension[0]}h {self.moon_right_ascension[1]}m {self.moon_right_ascension[2]:.2f}s",
                "altitude" : round(self.moon_alt, 3),
                "azimuth" : round(self.moon_az, 3),
                "parallax" : round(self.moon_pi * 60, 3),
                "illumination" : round(self.moon_illumin * 100, 2)
            }
    
    def moonphases(self) -> List[Dict[str, datetime]]:
        """
        Returns the nearest moon phases.

        Returns:
            list: List of dictionaries containing:
                - 'phase': Moon phase name
                - 'datetime': Date & time of the phase
        """

        # Find Next New Moon (and the rest of the phases)
        moon_phases = me.next_phases_of_moon_utc(self.observer_date)
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
            
            moon_phases[i] = {"phase": phase_str, "datetime": phase - timedelta(hours=self.utc_offset)}

        moon_phases = sorted(moon_phases, key = lambda item: item["datetime"])

        return moon_phases

    # Calculate Next New Moon Visibilities
    def visibilities(self, days: int = 3, criterion: int = 0) -> Dict [datetime, List[str | float]]:
        """
        Returns a dictionary describing the visibility of the nearest new moon for the observer.

        The dictionary contains `datetime` keys and lists of `str` or `float` values.
    
        The size of the dictionary is controlled by `days` which specifies how many days from the new moon to look at visibilities.

        The key of each item in the dictionary corresponds to the "Best Time" `datetime` at which to look for the new moon crescent.

        The value of each item in the dictionary is a list in which the first element is the raw number output of the visibility. The second element is the classification of the first element.

        The `type` argument specifies which new moon visibility classification method to use:
        - Criterion 0: Odeh, 2006
        - Criterion 1: Yallop, 1997; a.k.a. HMNAO TN No. 69

        Parameters:
            days (int, optional): How many days from the new moon to look at visibilities. Defaults to 3 days.
            criterion (int, optional): Which method to classify visibilities. Either 0 (referring to Odeh, 2006) or 1 (referring to Yallop, 1997). Defaults 0.

        Returns:
            list (List[Dict[str, str | float]]): A list of dictionaries containing:
                - 'datetime': The key of each dictionary which represents the best time to look for the new moon crescent.
                - 'visibility': Raw number output of the visibility.
        """

        if not isinstance(days, int):
            raise TypeError(f"'days' must be of type `int`, but got `{criterion(days).__name__}`.")
        
        if days < 1:
            raise ValueError(f"'days' must be greater than 0. Invalid value: {days}.")
        
        if criterion not in (0, 1):
            raise ValueError(f"'type' must be either 0 or 1. Invalid value: {criterion}.")
        
        
        if not isinstance(criterion, int):
            raise TypeError(f"'type' must be of type `int`, but got `{criterion(criterion).__name__}`.")
        
        if criterion not in (0, 1):
            raise ValueError(f"'type' must be either 0 or 1. Invalid value: {criterion}.")
        

        # Get New Moon Date from moon_phases list
        moon_phases = self.moonphases()
        for item in moon_phases:
                if item['phase'] == "New Moon":
                    new_moon = item['datetime']
                    break

        # Find JD for the given date; adjust day for difference in UTC and local timezone
        jd_new_moon = te.gregorian_to_jd(new_moon, -1 * self.utc_offset)
        ymd_new_moon = te.jd_to_gregorian(jd_new_moon)
        deltaT_new_moon = te.delta_t_approx(ymd_new_moon.year, ymd_new_moon.month)
        jde_new_moon = jd_new_moon + deltaT_new_moon / 86400

        # Forgot what this does. Likely something to do with timezone differences.
        if new_moon.day != ymd_new_moon.day:
            if new_moon.day < ymd_new_moon.day:
                new_moon += timedelta(days=1)
            else:
                new_moon -= timedelta(days=1)
            jd_new_moon = te.gregorian_to_jd(new_moon, -1 * self.utc_offset) - te.fraction_of_day(new_moon)

        # Find local sunset as visibilities are calculated from then
        # nm_sun_factors = se.sunpos(jde_new_moon, deltaT_new_moon, self.latitude, self.longitude)     

        # Find visibilities for the three days
        visibilities = []
        best_jds = []
        for day in range(days):
            # First, check if the moonset is before the new moon for the first day
            if day == 0:
                nm_moonset = me.find_proper_moontime(ymd_new_moon, self.observer_latitude, self.observer_longitude, self.observer_elevation, self.utc_offset, 'set')
                if nm_moonset == datetime.min:
                    # Moonset doesn't exist, for extreme latitudes
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
            test_jd_new_moon = jd_new_moon + day
            test_ymd_new_moon = te.jd_to_gregorian(test_jd_new_moon)
            test_deltaT_new_moon = te.delta_t_approx(test_ymd_new_moon.year, test_ymd_new_moon.month)
            test_jde_new_moon = test_jd_new_moon + test_deltaT_new_moon / 86400

            # Set sun parameters
            nm_sun_params = se.sunpos(test_jde_new_moon, test_deltaT_new_moon, self.observer_latitude, self.observer_longitude)

            # Sunset & moonset calculations
            test_nm_sunset = se.find_proper_suntime(test_ymd_new_moon, 
                                                  self.observer_latitude, 
                                                  self.observer_longitude,
                                                  self.observer_elevation,
                                                  self.utc_offset,
                                                  'set')

            if day == 0:
                test_nm_moonset = nm_moonset
            else:
                test_nm_moonset = me.find_proper_moontime(test_ymd_new_moon, self.observer_latitude, self.observer_longitude, self.observer_elevation, self.utc_offset, 'set')

            # For extreme latitudes where the moonset or sunset don't exist:
            if np.abs(self.observer_latitude) > 62:
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
            best_time_jd = te.gregorian_to_jd(best_time, -1 * self.utc_offset)
            best_time_jde = best_time_jd + test_deltaT_new_moon / 86400
            best_jds.append(best_time_jd)

            # Recalculate sun & calculate moon parameters
            nm_sun_params = se.sunpos(best_time_jde, test_deltaT_new_moon, self.observer_latitude, self.observer_longitude)
            nm_moon_params = me.moonpos(best_time_jde, test_deltaT_new_moon, self.observer_latitude, self.observer_longitude, nm_sun_params.delta_obliquity, nm_sun_params.true_obliquity, self.observer_elevation)

            # Visibility is now calculated
            v = me.calculate_visibility(nm_sun_params.azimuth, nm_sun_params.altitude, nm_moon_params.azimuth, nm_moon_params.altitude, nm_moon_params.eh_parallax, criterion)

            visibilities.append(v)

        # Arrange and classify visibilties
        q_values = [[v, me.classify_visibility(v, criterion)] for v in visibilities]

        # Convert best times from JD to datetime
        best_dates = [
            te.jd_to_gregorian(jd, self.utc_offset)
            for jd in best_jds
        ]

        # Label each q_value to its associated date
        visibility_dictionary = {
            dt.strftime('%Y-%m-%d %X'): q_value
            for dt, q_value in zip(best_dates, q_values)
        }

        return visibility_dictionary