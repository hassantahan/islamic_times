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

from numbers import Number
from typing import List
from dataclasses import replace
from datetime import datetime, timedelta, timezone
from islamic_times.it_dataclasses import *
from islamic_times import prayer_times as pt
from islamic_times import sun_equations as se
from islamic_times import moon_equations as me
from islamic_times import time_equations as te
from islamic_times import calculation_equations as ce
import islamic_times.astro_core as fast_astro

class ITLocation:
    """
    The `ITLocation` class is the primary way to interact with the `islamic_times` library. 
    It is designed to contain all the information and methods necessary to determine all Islamic times and astronomical parameters.

    This class provides methods for computing astronomical parameters and Islamic prayer times 
    based on various calculation methods. It also supports Hijri date conversion, new moon visibility 
    calculations, and Qibla direction determination.

    ### Methods:
    - `update_time(date_time)`: Updates the observer's time.
    - `calculate_astro()`: Computes astronomical parameters.
    - `calculate_prayer_times()`: Computes prayer times based on the selected method.
    - `set_prayer_method(method)`: Sets the prayer time calculation method.
    - `set_custom_prayer_angles(fajr_angle, maghrib_angle, isha_angle)`: Customizes solar angles for prayer time calculations.
    - `set_asr_type(asr_type)`: Sets the Asr prayer calculation method.
    - `set_midnight_type(midnight_type)`: Sets the Islamic midnight calculation method.
    - `observer()`: Returns observer location parameters.
    - `dates_times()`: Returns observer's date and time details.
    - `prayer_times()`: Returns calculated prayer times.
    - `mecca()`: Returns observer's distance and direction to Mecca.
    - `sun()`: Returns properties and position of the Sun.
    - `moon()`: Returns properties and position of the Moon.
    - `moonphases()`: Returns the nearest moon phases.
    """

    # TODO: I might want to add slots to save memory and freeze the class
    # __slots__ = ()

    def __init__(self, latitude: float = 51.477928,
                 longitude: float = -0.001545,
                 elevation: float = 76,
                 temperature: float = 10,
                 pressure: float = 101.325,
                 date: datetime = datetime.now().replace(tzinfo=timezone.utc),
                 method: str = 'JAFARI',
                 asr_type: int = 0,
                 find_local_tz: bool = False,
                 auto_calculate: bool = True) -> None:
        """
        `ITLocation` is initialized with the observer's geographical location, date and time, and other parameters. The default location is the Royal Greenwich Observatory.

        Parameters:
            latitude (float, optional): Geographical latitude in decimal degrees (-90 to 90). Defaults to 51.477928°.
            longitude (float, optional): Geographical longitude in decimal degrees (-180 to 180). Defaults to -0.001545°.
            elevation (float, optional): Elevation above sea level in meters. Defaults to 76 m.
            temperature (float, optional): Temperature in degrees Celsius. Defaults to 10°C.
            pressure (float, optional): Atmospheric pressure in kPa. Defaults to 101.325 kPa.
            date (datetime, optional): Current date and time. Defaults to `datetime.now(timezone.utc)`.
            method (str, optional): Prayer calculation method (e.g., 'JAFARI', 'ISNA', etc.). Defaults to 'JAFARI'.
            asr_type (int, optional): ʿAṣr calculation type (0 for standard, 1 for Ḥanafī). Defaults to 0.
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

        #  Check the numerical inputs
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

        self.observer_info = ObserverInfo(
            latitude=Angle(latitude),
            longitude=Angle(longitude),
            elevation=Distance(elevation),
            pressure=pressure,
            temperature=temperature
        )       

        #  Check the date if it is valid
        if not isinstance(date, datetime):
            raise TypeError(f"'{date}' must be of type `datetime`, but got `{type(date).__name__}`.")
        
        # Check if find_local_tz is either 0, 1, or a bool value
        if not isinstance(find_local_tz, bool):
            if not isinstance(find_local_tz, Number):
                raise ValueError(f"'find_local_tz' is out of range; it must be of type `bool` or either the numerical value of 0 or 1.")
            else:
                raise TypeError(f"'find_local_tz' must be of type `bool` or either the numerical value of 0 or 1, but got `{type(find_local_tz).__name__}`.")
        
        # Determine UTC Offset
        tz = self.get_timezone(find_local_tz, date)
        self.utc_offset = tz.utcoffset(date).total_seconds() / 3600
        date = date.replace(tzinfo=tz)

        jd = fast_astro.gregorian_to_jd(date, date.utcoffset().total_seconds() / 3600)
        deltaT = fast_astro.delta_t_approx(date.year, date.month)
        islamic_dates = te.gregorian_to_hijri(date.year, date.month, date.day)
        self.observer_dateinfo: DateTimeInfo = DateTimeInfo(
            date=date,
            hijri=IslamicDateInfo(*islamic_dates),
            jd=jd,
            deltaT=deltaT
        )
        test = fast_astro.jd_to_gregorian(jd, self.observer_dateinfo.utc_offset)

        # Autocalculation for astronomical parameters
        self.auto_calculate = auto_calculate
        if self.auto_calculate:
            self.calculate_astro()
            self.datetime_modified = False
            self.prayers_modified = False
        else:
            self.datetime_modified = True
            self.prayers_modified = True
        
        # Prayer setting
        self.set_prayer_method(method, asr_type) # This also calculates the prayer times if `auto_calculate` is True

    def get_timezone(self, find_local_tz: bool, date: datetime) -> timezone:
        """ Determine UTC offset in hours based on location if needed.

        Parameters:
            find_local_tz (bool): Controls whether or not to use the timezonefinder library to fine the timezone of the observer.
            date (datetime): The date and time to be used for the calculation.

        Returns:
            timezone: The timezone object representing the UTC offset.

        Notes:  
        - `find_local_tz` is set to `False` by default because the timezonefinder library is computationally expensive.
        """
        # Find UTC Offset According to Lat/Long datetime
        # This is very computationally expensive
        if find_local_tz: 
            tz_name, utc_offset = te.find_utc_offset(self.observer_info.latitude.decimal, self.observer_info.longitude.decimal, date)
            return timezone(offset=timedelta(hours=utc_offset), name=tz_name)
        elif date.tzinfo == None or date.tzinfo == timezone.utc:
            return timezone.utc
        else:  
            return date.tzinfo

    # Used to change observe date & time
    # By default, updates to datetime.now() if argument is not specified
    def update_time(self, new_date: datetime = None) -> None:
        """
        Updates the observer's time.

        This method updates the observer's time to either `datetime.now()` (if no argument is provided) 
        or to a specified `datetime` object. After updating the time, `update_astro()` must be called 
        to recalculate astronomical parameters.

        Parameters:
            new_date (datetime, optional): The new observer time. Defaults to the current time.

        Raises:
            TypeError: If `new_date` is not a `datetime` object.
        """

        if not isinstance(new_date, datetime):
            raise TypeError(f"'date_time' must be of type `datetime`, but got `{type(new_date).__name__}`.")


        # If TZ is NOT specified in new datetime, but a TZ *is* specified in old datetime, add the TZ into the new datetime
        if new_date.tzinfo is None and self.observer_dateinfo.date.tzinfo is not None:
            new_date = datetime.replace(new_date, tzinfo=self.observer_dateinfo.date.tzinfo)

        # Calculate DateInfo params
        jd = fast_astro.gregorian_to_jd(new_date, new_date.utcoffset().total_seconds() / 3600)
        deltaT = fast_astro.delta_t_approx(new_date.year, new_date.month)
        islamic_dates = te.gregorian_to_hijri(new_date.year, new_date.month, new_date.day)

        self.observer_dateinfo: DateTimeInfo = DateTimeInfo(
            date=new_date.replace(tzinfo=self.observer_dateinfo.date.tzinfo),
            hijri=IslamicDateInfo(*islamic_dates),
            jd=jd,
            deltaT=deltaT
        )
        
        # Set bools for astro_calculation stuff
        if not self.auto_calculate:
            self.datetime_modified = True
            self.prayers_modified = True
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

        def safe_sun_time(observer_dateinfo, observer_info, event: str) -> datetime | str:
            try:
                return se.find_proper_suntime(observer_dateinfo, observer_info, event)
            except ArithmeticError:
                return f"Sun{event} does not exist."
            
        def safe_moon_time(observer_dateinfo, observer_info, event: str) -> datetime | str:
            try:
                return me.find_proper_moontime(observer_dateinfo, observer_info, event)
            except ArithmeticError:
                return f"Moon{event} does not exist."

        ### Sun & Moon Properties Calculations
        # Get Sun and Moon Objects with their parameters
        self.sun_params: se.Sun = se.sunpos(self.observer_dateinfo, self.observer_info)
        self.moon_params: me.Moon = me.moonpos(self.observer_dateinfo, self.observer_info, self.sun_params.nutation[0], self.sun_params.true_obliquity)

        # Important Sun Factors placed SunInfo
        self.sun_info = SunInfo(
            sunrise=safe_sun_time(self.observer_dateinfo, self.observer_info, 'rise'),
            sun_transit=se.find_sun_transit(self.observer_dateinfo, self.observer_info),
            sunset=safe_sun_time(self.observer_dateinfo, self.observer_info, 'set'),
            apparent_altitude=self.sun_params.true_altitude,
            true_azimuth=self.sun_params.true_azimuth,
            geocentric_distance=self.sun_params.geocentric_distance,
            apparent_declination=self.sun_params.apparent_declination,
            apparent_right_ascension=self.sun_params.apparent_right_ascension,
            greenwich_hour_angle=self.sun_params.greenwich_hour_angle,
            local_hour_angle=self.sun_params.local_hour_angle
        )

        illumination: float = me.moon_illumination(self.sun_params.apparent_declination, 
                                self.sun_params.apparent_right_ascension, 
                                self.moon_params.declination, 
                                self.moon_params.right_ascension, 
                                self.sun_params.geocentric_distance, 
                                self.moon_params.geocentric_distance.to(DistanceUnits.AU)
                            )

        # Important Moon Factors placed into MoonInfo
        self.moon_info = MoonInfo(
            moonrise=safe_moon_time(self.observer_dateinfo, self.observer_info, 'rise'),
            moon_transit=me.find_moon_transit(self.observer_dateinfo, self.observer_info),
            moonset=safe_moon_time(self.observer_dateinfo, self.observer_info, 'set'),
            illumination=illumination,
            apparent_altitude=self.moon_params.apparent_altitude,
            true_azimuth=self.moon_params.true_azimuth,
            geocentric_distance=self.moon_params.geocentric_distance,
            parallax=self.moon_params.eh_parallax,
            topocentric_declination=self.moon_params.top_declination,
            topocentric_right_ascension=self.moon_params.topocentric_ascension,
            greenwich_hour_angle=self.moon_params.greenwich_hour_angle,
            local_hour_angle=self.moon_params.local_hour_angle
        )

        # Astronomical parameters have been calculated so the flag is set to False
        self.datetime_modified = False

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

        self.times_of_prayer: PrayerTimes = pt.calculate_prayer_times(self.observer_dateinfo, self.observer_info, self.sun_info, self.method)

        # Prayers have been calculated so the flag is set to False
        self.prayers_modified = False

    # Set the method of calculating prayer times among the available default options
    # The default option (from creation) is the Jaʿfarī method.
    def set_prayer_method(self, method_key: str = 'JAFARI', asr_type: int = 0) -> None:
        """
        Sets the prayer time calculation method.

        The available methods follow those documented at:
        http://praytimes.org/wiki/Calculation_Methods.

        Parameters:
            method_key (str): The name of the prayer calculation method (e.g., 'JAFARI', 'ISNA', etc.). Defaults to 'JAFARI'.
            asr_type (int): The ʿaṣr calculation type 0 for standard, 1 for Ḥanafī. Default is set to 0.

        Raises:
            ValueError: If the method is not among the supported options.

        Notes:
            - If `auto_calculate` is enabled, prayer times are recalculated automatically. Otherwise, `calculate_prayer_times()` must be called.
            - To revert to a default method after using `set_custom_prayer_angles()`, this method must be called again.
        """

        method_key = method_key.strip().upper()
    
        for method in pt.DEFAULT_PRAYER_METHODS:
            if method_key in (key.upper() for key in method.keys):
                if asr_type not in (0, 1):
                    raise ValueError(f"'asr_type' must be either 0 or 1. Invalid value: {asr_type}")

                if asr_type == 0:
                    self.method = method
                else:
                    self.method = replace(method, asr_type=asr_type)

                if not self.prayers_modified:
                    self.calculate_prayer_times()
                
                return

        valid_options = [method.keys[0] for method in pt.DEFAULT_PRAYER_METHODS]
        raise ValueError(
            f"Invalid prayer method '{method_key}'. Valid options are: {', '.join(valid_options)}")
    
    # Helper function to validate and set angle
    def __validate_and_set(self, attribute_name: str, value: float | int) -> None:
        if value is not None:
            if isinstance(value, (int, float)):  # Check if it's a number
                if value > 0:
                    self.method = replace(self.method, **{attribute_name: Angle(value)})
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
            TypeError: If any provided angle is not a number.

        Notes:
            - If `auto_calculate` is enabled, prayer times are recalculated automatically. Otherwise, `calculate_prayer_times()` must be called.
            - `self.method` is set to 'Custom' after calling this method.
            - Call `set_prayer_method()` to reset to a predefined method.
        """

        # Validate and set each angle
        self.__validate_and_set('fajr_angle', fajr_angle)
        self.__validate_and_set('maghrib_angle', maghrib_angle)
        self.__validate_and_set('isha_angle', isha_angle)

        self.method = replace(self.method, name="Custom")

        # Update prayer times
        if self.auto_calculate:
            self.calculate_prayer_times()
        else:
            self.prayers_modified = True
    
    # Separated from angles since it is defined by shadow ratio
    def set_asr_type(self, asr_type: int = 0) -> None:
        """
        Sets the calculation method for Asr prayer.

        Options:
        - `0`: Shadow ratio of 1:1 (majority method).
        - `1`: Shadow ratio of 2:1 (Hanafi method).

        Parameters:
            asr_type (int): The ʿaṣr calculation type. 0 for standard, 1 for Ḥanafī. Default is set to 0.

        Raises:
            ValueError: If `asr_type` is not 0 or 1.

        Notes:
        - If `auto_calculate` is enabled, prayer times are recalculated automatically. Otherwise, `calculate_prayer_times()` must be called.
        - `self.method` is set to 'Custom' after calling this method.
        - Call `set_prayer_method()` to reset to a predefined method.
        """

        if asr_type in (0, 1):
            self.method = replace(self.method, asr_type=asr_type)
        else:
            raise ValueError(f"'asr_type' must be either 0 or 1. Check documentation to understand each type. Invalid value: {asr_type}")

        # Method is now custom
        self.method = replace(self.method, name="Custom")

        # Update prayer times
        if self.auto_calculate:
            self.calculate_prayer_times()
        else:
            self.prayers_modified = True

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
           self.method = replace(self.method, midnight_type=midnight_type)
        else:
            raise ValueError(f"'midnight_type' must be either 0 or 1. Check documentation to understand each type. Invalid value: {midnight_type}")

        # Method is now custom
        self.method = replace(self.method, name="Custom")

        # Update prayer times
        if self.auto_calculate:
            self.calculate_prayer_times()
        else:
            self.prayers_modified = True

    # Return Observer Parameters
    def observer(self) -> ObserverInfo:
        """
        Returns observer's date and time information.

        Returns:
            ObserverInfo: Information about the observer's location and conditions.

        Notes:
        - The returned object contains:
            - **latitude** (*Angle*): Latitude of the observer
            - **longitude** (*Angle*): Longitude of the observer
            - **elevation** (*Distance*): Elevation of the observer
            - **temperature** (*float*): Temperature of the observer
            - **pressure** (*float*): Pressure of the observer
        """

        return self.observer_info

    # Return date and time information
    def dates_times(self) -> DateTimeInfo:
        """
        Returns observer's date and time information.

        Returns:
            DateTimeInfo: Information about the observer's date and time.

        Notes:
        - The returned object contains:
            - **date** (*datetime*): The observer's date and time
            - **hijri** (*IslamicDateInfo*): Islamic date information
                - **hijri_year** (*int*): Hijri year
                - **hijri_month** (*int*): Hijri month
                - **hijri_day** (*int*): Hijri day
                - **hijri_month_name** (*str*): Hijri month name
                - **hijri_day_name** (*str*): Hijri day name
            - **jd** (*float*): Julian Date
            - **deltaT** (*float*): Delta T value
        """

        can_print = self.auto_calculate or not self.datetime_modified

        if not can_print:
            raise ValueError("Cannot print dates and times without calculating the astronomical parameters. Set `auto_calculate` to `True` or call `calculate_astro()`.")

        return self.observer_dateinfo

    # Return prayer times
    def prayer_times(self) -> PrayerTimes:
        """
        Returns calculated prayer times.

        Returns:
            PrayerTimes: Dictionary containing prayer times.

        Notes:
        - The returned object contains:
            - **fajr** (*datetime*): Fajr prayer time
            - **sunrise** (*datetime*): Sunrise time
            - **dhuhr** (*datetime*): Dhuhr prayer time
            - **asr** (*datetime*): Asr prayer time
            - **maghrib** (*datetime*): Maghrib prayer time
            - **isha** (*datetime*): Isha prayer time
            - **midnight** (*datetime*): Islamic midnight time
        """

        can_print: bool = self.auto_calculate or not self.prayers_modified

        if not can_print:
            raise ValueError("`auto_calculate` has been set to False and a change to the prayer methods and/or the date & time may have been made. `calculate_prayer_times()` must be called first.")
        
        return self.times_of_prayer
    
    # Return Mecca information
    def mecca(self) -> MeccaInfo:
        """
        Returns observer's distance and direction to Mecca.

        Returns:
            MeccaInfo: Information about the distance and direction to Mecca.
        
        Notes:
        - The returned object contains:
            - **distance** (*Distance*): Distance to Mecca
            - **angle** (*Angle*): Angle to Mecca
            - **cardinal** (*str*): Cardinal direction to Mecca
        """

        mecca_distance, mecca_direction = ce.haversine(self.observer_info.latitude.decimal, self.observer_info.longitude.decimal, te.MECCA_LAT, te.MECCA_LONG)
        mecca_direction %= 360
        
        return MeccaInfo(
            distance=Distance(mecca_distance, DistanceUnits.KILOMETRE),
            angle=Angle(mecca_direction),
            cardinal=ce.get_cardinal_direction(mecca_direction)
        )
    
    # Return sun properties and position values
    def sun(self) -> SunInfo:
        """
        Returns properties and position of the Sun.

        Returns:
            SunInfo: Information about the Sun's position and properties.

        Notes:
        - The returned object contains:
            - **sunrise** (*datetime*): Sunrise time.
            - **sun_transit** (*datetime*): Solar noon time.
            - **sunset** (*datetime*): Sunset time.
            - **apparent_declination** (*Angle*): Sun's declination angle.
            - **apparent_right_ascension** (*RigthAscension*): Sun's right ascension.
            - **apparent_altitude** (*Angle*): Sun's altitude angle.
            - **true_azimuth** (*Angle*): Sun's azimuth angle.
        """

        can_print = self.auto_calculate or not self.datetime_modified

        if not can_print:
            raise ValueError("Cannot print dates and times without calculating the astronomical parameters. First call `calculate_astro()`.")
        
        return self.sun_info
    
    # Return moon properties and position values
    def moon(self) -> MoonInfo:
        """
        Returns properties and position of the Moon.

        Returns:
            MoonInfo: Information about the Moon's position and properties.

        Notes:
        - The returned object contains:
            - **moonrise** (*datetime*): Moonrise time.
            - **moon_transit** (*datetime*): Moon transit time.
            - **moonset** (*datetime*): Moonset time.
            - **topocentric_declination** (*Angle*): Moon's declination angle.
            - **topocentric_right_ascension** (*RightAscension*): Moon's right ascension.
            - **apparent_altitude** (*Angle*): Moon's altitude angle.
            - **true_azimuth** (*Angle*): Moon's azimuth angle.
            - **parallax** (*Angle*): Moon's parallax angle.
            - **illumination** (*float*): Moon's illumination percentage.
        """

        can_print = self.auto_calculate or not self.datetime_modified

        if not can_print:
            raise ValueError("Cannot print dates and times without calculating the astronomical parameters. First call `calculate_astro()`.")

        return self.moon_info
    
    def moonphases(self) -> List[Tuple[str, datetime]]:
        """
        Returns the nearest moon phases.

        Returns:
            list: List of dictionaries containing:
                - 'phase': Moon phase name
                - 'datetime': Date & time of the phase
        """

        # Find Next New Moon (and the rest of the phases)
        phases: Tuple[datetime, datetime, datetime, datetime] = me.next_phases_of_moon_utc(self.observer_dateinfo.date)
        phase_names = ["New Moon", "First Quarter", "Full Moon", "Last Quarter"]

        return [(phase_names[i], phase + timedelta(hours=self.utc_offset)) for i, phase in enumerate(phases)]

    # Calculate Next New Moon Visibilities
    def visibilities(self, days: int = 3, criterion: int = 1) -> Visibilities:
        """
        Returns a `Visibilities` dataclass ojbect describing the visibility of the nearest new moon 
        for the observer for the given amount of days according to a selected criterion.

        The `criterion` argument specifies which new moon visibility classification method to use:
        - Criterion 0: Odeh, 2006
        - Criterion 1: Yallop, 1997; a.k.a. HMNAO TN No. 69

        Parameters:
            days (int, optional): How many days from the new moon to look at visibilities. Defaults to 3 days.
            criterion (int, optional): Which method to classify visibilities. Either 0 (referring to Odeh, 2006) or 1 (referring to Yallop, 1997). Defaults 1.

        Returns:
            Visibilities: Dataclass containing visibility information.
        """

        if not isinstance(days, int):
            raise TypeError(f"'days' must be of type `int`, but got `{criterion(days).__name__}`.")
        
        if days < 1:
            raise ValueError(f"'days' must be greater than 0. Invalid value: {days}.")
        
        if criterion not in (0, 1):
            raise ValueError(f"'criterion' must be either 0 or 1. Invalid value: {criterion}.")
        
        
        if not isinstance(criterion, int):
            raise TypeError(f"'criterion' must be of type `int`, but got `{criterion(criterion).__name__}`.")
        
        if criterion not in (0, 1):
            raise ValueError(f"'criterion' must be either 0 or 1. Invalid value: {criterion}.")

        visibilities: Visibilities = fast_astro.compute_visibilities(self.observer_dateinfo.date, self.observer_dateinfo.utc_offset, 
                                              self.observer_info.latitude.decimal, self.observer_info.longitude.decimal, 
                                              self.observer_info.elevation.in_unit(DistanceUnits.METRE), 
                                              self.observer_info.temperature, self.observer_info.pressure, 
                                              days, criterion)

        return visibilities
