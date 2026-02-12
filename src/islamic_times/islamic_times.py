"""
Module for calculating Islamic times and astronomical data.

This module provides the ITLocation class that encapsulates calculations for:
  - Astronomical parameters (sun and moon positions, Julian Date, etc.)
  - Prayer times based on various calculation methods.
  - Islamic calendar (Hijri) conversion.
  - New moon visibilities

References:
  - Jean Meeus, *Astronomical Algorithms*, 2nd Edition, Willmann-Bell, Inc., 1998.
  - Prayer times calculation methods (https://praytimes.org/docs/calculation)
"""

import islamic_times.astro_core as fast_astro
from numbers import Number
from typing import List, Optional, Tuple
from dataclasses import replace
from datetime import datetime, timedelta, timezone, tzinfo
from islamic_times.it_dataclasses import (
    Angle,
    DateTimeInfo,
    Distance,
    DistanceUnits,
    IslamicDateInfo,
    MeccaInfo,
    MoonInfo,
    ObserverInfo,
    PrayerTimes,
    SunInfo,
    Visibilities,
)
from islamic_times import calculation_equations as ce
from islamic_times import moon_equations as me
from islamic_times import prayer_times as pt
from islamic_times import sun_equations as se
from islamic_times import time_equations as te

class ITLocation:
    """Primary public facade for astronomical and prayer-time calculations.

    ``ITLocation`` encapsulates observer configuration (location, environment,
    method settings, date/time context) and exposes a high-level API for:

    - sun and moon coordinates/events,
    - prayer times,
    - Hijri date conversion context,
    - Qibla direction and distance,
    - nearest new-moon visibility outputs.
    """

    # Potential future optimization: introduce __slots__ after API stability is finalized.
    # __slots__ = ()

    def __init__(
        self,
        latitude: float = 51.477928,
        longitude: float = -0.001545,
        elevation: float = 76,
        temperature: float = 10,
        pressure: float = 101.325,
        date: Optional[datetime] = None,
        method: str = "JAFARI",
        asr_type: int = 0,
        find_local_tz: bool = False,
        auto_calculate: bool = True,
    ) -> None:
        """Initialize observer state and optional eager calculations.

        Parameters
        ----------
        latitude : float, default=51.477928
            Geodetic latitude in decimal degrees in ``[-90, 90]``.
        longitude : float, default=-0.001545
            Geodetic longitude in decimal degrees in ``[-180, 180]``.
        elevation : float, default=76
            Elevation above mean sea level in meters.
        temperature : float, default=10
            Ambient temperature in degrees Celsius.
        pressure : float, default=101.325
            Atmospheric pressure in kPa.
        date : datetime, optional
            Observer-local datetime. If ``None``, current UTC time is used.
        method : str, default='JAFARI'
            Prayer method key (for example ``'JAFARI'``, ``'ISNA'``, ``'MWL'``).
        asr_type : int, default=0
            ʿAṣr rule selector (``0`` standard, ``1`` Ḥanafī).
        find_local_tz : bool | int | float, default=False
            If truthy (or numeric 1), infer timezone from coordinates and date.
        auto_calculate : bool, default=True
            If ``True``, compute astronomy/prayer values during initialization.

        Raises
        ------
        TypeError
            Raised when numeric inputs or ``date`` have invalid types.
        ValueError
            Raised for out-of-range latitude/longitude, unsupported method keys,
            or invalid method selectors.
        """

        if date is None:
            date = datetime.now(timezone.utc)

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
        find_local_tz_flag = self._resolve_find_local_tz(find_local_tz)

        # Determine UTC Offset
        tz = self.get_timezone(find_local_tz_flag, date)
        tz_offset = tz.utcoffset(date)
        if tz_offset is None:
            raise ValueError("Could not determine UTC offset for the provided date and timezone.")
        self.utc_offset = tz_offset.total_seconds() / 3600
        date = date.replace(tzinfo=tz)

        self.observer_dateinfo = self._build_observer_dateinfo(date)

        # Autocalculation for astronomical parameters
        self.auto_calculate = auto_calculate
        if self.auto_calculate:
            self.calculate_astro()
            self._mark_astro_dirty(False)
            self._mark_prayers_dirty(False)
        else:
            self._mark_astro_dirty(True)
            self._mark_prayers_dirty(True)

        # Prayer setting
        self.set_prayer_method(method, asr_type) # This also calculates the prayer times if `auto_calculate` is True

    @staticmethod
    def _resolve_find_local_tz(find_local_tz: bool | int | float) -> bool:
        """Normalize timezone-resolution flag from bool/0/1 inputs."""
        if isinstance(find_local_tz, bool):
            return find_local_tz

        if not isinstance(find_local_tz, Number):
            raise ValueError(
                "'find_local_tz' must be of type `bool` or either the numerical value of 0 or 1."
            )

        if find_local_tz not in (0, 1):
            raise TypeError(
                f"'find_local_tz' must be of type `bool` or either the numerical value of 0 or 1, "
                f"but got `{type(find_local_tz).__name__}`."
            )

        return bool(find_local_tz)

    @staticmethod
    def _safe_sun_event(observer_dateinfo: DateTimeInfo, observer_info: ObserverInfo, event: str) -> datetime | str:
        """Return a sun event datetime or a human-readable unavailable marker."""
        try:
            return se.find_proper_suntime(observer_dateinfo, observer_info, event)
        except ArithmeticError:
            return f"Sun{event} does not exist."

    @staticmethod
    def _safe_moon_event(observer_dateinfo: DateTimeInfo, observer_info: ObserverInfo, event: str) -> datetime | str:
        """Return a moon event datetime or a human-readable unavailable marker."""
        try:
            return me.find_proper_moontime(observer_dateinfo, observer_info, event)
        except ArithmeticError:
            return f"Moon{event} does not exist."

    def _mark_astro_dirty(self, is_dirty: bool) -> None:
        """Set astronomy cache dirty flag."""
        self.datetime_modified = is_dirty

    def _mark_prayers_dirty(self, is_dirty: bool) -> None:
        """Set prayer-times cache dirty flag."""
        self.prayers_modified = is_dirty

    def _refresh_prayers(self) -> None:
        """Recompute or mark prayer outputs dirty based on auto-calc mode."""
        if self.auto_calculate:
            self.calculate_prayer_times()
        else:
            self._mark_prayers_dirty(True)

    def _build_observer_dateinfo(self, date: datetime) -> DateTimeInfo:
        """Build DateTimeInfo from a timezone-aware datetime."""
        utc_offset = date.utcoffset()
        if utc_offset is None:
            raise ValueError("Input datetime must include timezone information.")

        jd = fast_astro.gregorian_to_jd(date, utc_offset.total_seconds() / 3600)
        delta_t = fast_astro.delta_t_approx(date.year, date.month)
        islamic_dates = te.gregorian_to_hijri(date.year, date.month, date.day)
        return DateTimeInfo(
            date=date,
            hijri=IslamicDateInfo(*islamic_dates),
            jd=jd,
            deltaT=delta_t,
        )

    def get_timezone(self, find_local_tz: bool, date: datetime) -> timezone | tzinfo:
        """Resolve timezone from either coordinates or supplied datetime context.

        Parameters
        ----------
        find_local_tz : bool
            If ``True``, infer timezone from observer coordinates and date.
        date : datetime
            Datetime used for timezone interpretation.

        Returns
        -------
        timezone | tzinfo
            Timezone object used by the instance.
        """
        # Find UTC Offset According to Lat/Long datetime
        # This is very computationally expensive
        if find_local_tz:
            tz_name, utc_offset = te.find_utc_offset(self.observer_info.latitude.decimal, self.observer_info.longitude.decimal, date)
            return timezone(offset=timedelta(hours=utc_offset), name=tz_name)
        if date.tzinfo is None or date.tzinfo == timezone.utc:
            return timezone.utc
        return date.tzinfo

    # Used to change observe date & time
    # By default, updates to current datetime in the observer's timezone when no argument is specified.
    def update_time(self, new_date: Optional[datetime] = None) -> None:
        """Update the observer datetime and mark/recompute dependent state.

        Parameters
        ----------
        new_date : datetime, optional
            Replacement datetime. If omitted, current time in the instance
            timezone is used.

        Raises
        ------
        TypeError
            Raised when ``new_date`` is provided with a non-datetime type.
        """

        if new_date is None:
            current_tz = self.observer_dateinfo.date.tzinfo or timezone.utc
            new_date = datetime.now(current_tz)
        elif not isinstance(new_date, datetime):
            raise TypeError(f"'new_date' must be of type `datetime`, but got `{type(new_date).__name__}`.")

        # If TZ is not specified in new datetime, but a TZ is specified in old datetime, preserve the previous TZ.
        if new_date.tzinfo is None and self.observer_dateinfo.date.tzinfo is not None:
            new_date = new_date.replace(tzinfo=self.observer_dateinfo.date.tzinfo)

        self.observer_dateinfo = self._build_observer_dateinfo(new_date)

        if self.auto_calculate:
            self.calculate_astro()
        else:
            self._mark_astro_dirty(True)
            self._mark_prayers_dirty(True)

    # Calculates the astronomical variables for the moon and sun
    def calculate_astro(self) -> None:
        """Compute sun/moon astronomy and refresh cached dataclass outputs.

        Notes
        -----
        This updates internal ``sun_params``/``moon_params`` and public
        ``sun_info``/``moon_info`` fields. It also clears the astronomy-dirty
        flag once complete.
        """

        # Sun & Moon properties calculations
        # Get Sun and Moon Objects with their parameters
        self.sun_params: se.Sun = se.sunpos(self.observer_dateinfo, self.observer_info)
        self.moon_params: me.Moon = me.moonpos(self.observer_dateinfo, self.observer_info, self.sun_params.nutation[0], self.sun_params.true_obliquity)

        # Important Sun Factors placed SunInfo
        self.sun_info = SunInfo(
            sunrise=self._safe_sun_event(self.observer_dateinfo, self.observer_info, "rise"),
            sun_transit=se.find_sun_transit(self.observer_dateinfo, self.observer_info),
            sunset=self._safe_sun_event(self.observer_dateinfo, self.observer_info, "set"),
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
            moonrise=self._safe_moon_event(self.observer_dateinfo, self.observer_info, "rise"),
            moon_transit=me.find_moon_transit(self.observer_dateinfo, self.observer_info),
            moonset=self._safe_moon_event(self.observer_dateinfo, self.observer_info, "set"),
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
        self._mark_astro_dirty(False)

    # Prayer Time Calculations
    def calculate_prayer_times(self) -> None:
        """Compute prayer times from current observer, astronomy, and method state.

        Raises
        ------
        ValueError
            Raised when astronomy is stale while ``auto_calculate`` is disabled.
        """
        can_calculate = self.auto_calculate or not self.datetime_modified
        if not can_calculate:
            raise ValueError(
                "Since `auto_calculate` is False, prayer times cannot be calculated until astronomical parameters are updated. "
                "Call `calculate_astro()` first."
            )

        self.times_of_prayer: PrayerTimes = pt.calculate_prayer_times(self.observer_dateinfo, self.observer_info, self.sun_info, self.method)

        # Prayers have been calculated so the flag is set to False
        self._mark_prayers_dirty(False)

    # Set the method of calculating prayer times among the available default options
    # The default option (from creation) is the Jaʿfarī method.
    def set_prayer_method(self, method_key: str = 'JAFARI', asr_type: int = 0) -> None:
        """Select one of the predefined prayer profiles.

        Parameters
        ----------
        method_key : str, default='JAFARI'
            Method alias defined in ``prayer_times.DEFAULT_PRAYER_METHODS``.
        asr_type : int, default=0
            ʿAṣr selector (``0`` standard, ``1`` Ḥanafī).

        Raises
        ------
        ValueError
            Raised for unknown method aliases or invalid ``asr_type`` values.
        """

        method_key = method_key.strip().upper()
    
        for method in pt.DEFAULT_PRAYER_METHODS:
            method_keys = method.keys or ()
            if method_key in (key.upper() for key in method_keys):
                if asr_type not in (0, 1):
                    raise ValueError(f"'asr_type' must be either 0 or 1. Invalid value: {asr_type}")

                if asr_type == 0:
                    self.method = method
                else:
                    self.method = replace(method, asr_type=asr_type)

                self._refresh_prayers()
                
                return

        valid_options = [method.keys[0] for method in pt.DEFAULT_PRAYER_METHODS if method.keys]
        raise ValueError(
            f"Invalid prayer method '{method_key}'. Valid options are: {', '.join(valid_options)}")
    
    # Helper function to validate and set angle
    def __validate_and_set(self, attribute_name: str, value: Optional[float]) -> None:
        """Validate an optional angle input and update method configuration."""
        if value is None:
            return

        if not isinstance(value, (int, float)) or isinstance(value, bool):
            raise ValueError(f"{attribute_name} must be a number. Invalid value: {value}")

        if value <= 0:
            raise ValueError(f"{attribute_name} must be greater than 0. Invalid value: {value}")

        self.method = replace(self.method, **{attribute_name: Angle(float(value))})

    # Alows user to set their own solar hour angles for prayer time calculations.
    def set_custom_prayer_angles(
        self,
        fajr_angle: Optional[float] = None,
        maghrib_angle: Optional[float] = None,
        isha_angle: Optional[float] = None,
    ) -> None:
        """Override default solar angles used by prayer-time computations.

        Parameters
        ----------
        fajr_angle : float, optional
            Fajr angle in degrees below the horizon.
        maghrib_angle : float, optional
            Maghrib angle in degrees below the horizon.
        isha_angle : float, optional
            ʿIshāʾ angle in degrees below the horizon.

        Raises
        ------
        ValueError
            Raised when a provided value is non-numeric or not strictly positive.
        """

        # Validate and set each angle
        self.__validate_and_set('fajr_angle', fajr_angle)
        self.__validate_and_set('maghrib_angle', maghrib_angle)
        self.__validate_and_set('isha_angle', isha_angle)

        self.method = replace(self.method, name="Custom")

        self._refresh_prayers()
    
    # Separated from angles since it is defined by shadow ratio
    def set_asr_type(self, asr_type: int = 0) -> None:
        """Set the ʿAṣr shadow-ratio rule.

        Parameters
        ----------
        asr_type : int, default=0
            ``0`` for the majority rule, ``1`` for the Ḥanafī rule.

        Raises
        ------
        ValueError
            Raised when ``asr_type`` is not ``0`` or ``1``.
        """

        if asr_type in (0, 1):
            self.method = replace(self.method, asr_type=asr_type)
        else:
            raise ValueError(f"'asr_type' must be either 0 or 1. Check documentation to understand each type. Invalid value: {asr_type}")

        # Method is now custom
        self.method = replace(self.method, name="Custom")

        self._refresh_prayers()

    # Set to either 0 (sunset to sunrise; the majority method) or 1 (sunset to fajr, the 'Jaʿfarī' method)
    def set_midnight_type(self, midnight_type: int = 0) -> None:
        """Set the Islamic-midnight definition.

        Parameters
        ----------
        midnight_type : int, default=0
            ``0`` for midpoint between sunset and sunrise, ``1`` for midpoint
            between sunset and Fajr.

        Raises
        ------
        ValueError
            Raised when ``midnight_type`` is not ``0`` or ``1``.
        """

        if midnight_type in (0, 1):
           self.method = replace(self.method, midnight_type=midnight_type)
        else:
            raise ValueError(f"'midnight_type' must be either 0 or 1. Check documentation to understand each type. Invalid value: {midnight_type}")

        # Method is now custom
        self.method = replace(self.method, name="Custom")

        self._refresh_prayers()
        
    def set_extreme_latitude_rule(self, rule: str = 'ANGLEBASED') -> None:
        """Set the fallback strategy for extreme-latitude prayer calculations.

        Parameters
        ----------
        rule : str, default='ANGLEBASED'
            One of ``'NONE'``, ``'NEARESTLAT'``, ``'MIDDLENIGHT'``,
            ``'ONESEVENTH'``, or ``'ANGLEBASED'``.

        Raises
        ------
        TypeError
            Raised when ``rule`` is not a string.
        ValueError
            Raised when ``rule`` is not a supported option.
        """
        EXTREME_LATITUDE_RULES: List[str] = ['NONE', 'NEARESTLAT', 'MIDDLENIGHT', 'ONESEVENTH', 'ANGLEBASED']
        
        if not isinstance(rule, str):
            raise TypeError("'rule' must be a string type.")
        normalized_rule = rule.strip().upper()
        if normalized_rule not in EXTREME_LATITUDE_RULES:
            raise ValueError(f'The value of {rule} for \'rule\' is not a valid rule. Please select from among the following: {EXTREME_LATITUDE_RULES}')
        
        self.method = replace(self.method, extreme_lats=normalized_rule)

        self._refresh_prayers()

    # Return Observer Parameters
    def observer(self) -> ObserverInfo:
        """Return observer location and environmental configuration."""

        return self.observer_info

    # Return date and time information
    def dates_times(self) -> DateTimeInfo:
        """Return datetime, Hijri, Julian-day, and delta-T context.

        Raises
        ------
        ValueError
            Raised when astronomy context is stale and ``auto_calculate`` is off.
        """

        can_print = self.auto_calculate or not self.datetime_modified

        if not can_print:
            raise ValueError("Cannot print dates and times without calculating the astronomical parameters. Set `auto_calculate` to `True` or call `calculate_astro()`.")

        return self.observer_dateinfo

    # Return prayer times
    def prayer_times(self) -> PrayerTimes:
        """Return cached prayer-time output dataclass.

        Raises
        ------
        ValueError
            Raised when prayer results are stale and ``auto_calculate`` is off.
        """

        can_print: bool = self.auto_calculate or not self.prayers_modified

        if not can_print:
            raise ValueError("`auto_calculate` has been set to False and a change to the prayer methods and/or the date & time may have been made. `calculate_prayer_times()` must be called first.")
        
        return self.times_of_prayer
    
    # Return Mecca information
    def mecca(self) -> MeccaInfo:
        """Return distance and bearing from observer to Mecca."""

        mecca_distance, mecca_direction = ce.haversine(self.observer_info.latitude.decimal, self.observer_info.longitude.decimal, te.MECCA_LAT, te.MECCA_LONG)
        mecca_direction %= 360
        
        return MeccaInfo(
            distance=Distance(mecca_distance, DistanceUnits.KILOMETRE),
            angle=Angle(mecca_direction),
            cardinal=ce.get_cardinal_direction(mecca_direction)
        )
    
    # Return sun properties and position values
    def sun(self) -> SunInfo:
        """Return cached solar event and coordinate summary.

        Raises
        ------
        ValueError
            Raised when astronomy context is stale and ``auto_calculate`` is off.
        """

        can_print = self.auto_calculate or not self.datetime_modified

        if not can_print:
            raise ValueError("Cannot print dates and times without calculating the astronomical parameters. First call `calculate_astro()`.")
        
        return self.sun_info
    
    # Return moon properties and position values
    def moon(self) -> MoonInfo:
        """Return cached lunar event and coordinate summary.

        Raises
        ------
        ValueError
            Raised when astronomy context is stale and ``auto_calculate`` is off.
        """

        can_print = self.auto_calculate or not self.datetime_modified

        if not can_print:
            raise ValueError("Cannot print dates and times without calculating the astronomical parameters. First call `calculate_astro()`.")

        return self.moon_info
    
    def moonphases(self) -> List[Tuple[str, datetime]]:
        """Return upcoming lunar phases as ``(name, datetime)`` tuples."""

        # Find Next New Moon (and the rest of the phases)
        phases: Tuple[datetime, datetime, datetime, datetime] = me.next_phases_of_moon_utc(self.observer_dateinfo.date)
        phase_names = ["New Moon", "First Quarter", "Full Moon", "Last Quarter"]

        return [(phase_names[i], phase + timedelta(hours=self.utc_offset)) for i, phase in enumerate(phases)]

    # Calculate Next New Moon Visibilities
    def visibilities(self, days: int = 3, criterion: int = 1) -> Visibilities:
        """Compute new-moon visibility predictions for the configured observer.

        Parameters
        ----------
        days : int, default=3
            Number of consecutive days to evaluate from the nearest new moon.
        criterion : int, default=1
            Visibility classifier: ``0`` (Odeh) or ``1`` (Yallop/HMNAO TN 69).

        Returns
        -------
        Visibilities
            Dataclass containing criterion name, event datetimes, q-values,
            and textual classifications.
        """

        if not isinstance(days, int) or isinstance(days, bool):
            raise TypeError(f"'days' must be of type `int`, but got `{type(days).__name__}`.")
        
        if days < 1:
            raise ValueError(f"'days' must be greater than 0. Invalid value: {days}.")

        if not isinstance(criterion, int) or isinstance(criterion, bool):
            raise TypeError(f"'criterion' must be of type `int`, but got `{type(criterion).__name__}`.")

        if criterion not in (0, 1):
            raise ValueError(f"'criterion' must be either 0 or 1. Invalid value: {criterion}.")

        visibilities: Visibilities = fast_astro.compute_visibilities(self.observer_dateinfo.date, self.observer_dateinfo.utc_offset, 
                                              self.observer_info.latitude.decimal, self.observer_info.longitude.decimal, 
                                              self.observer_info.elevation.in_unit(DistanceUnits.METRE), 
                                              self.observer_info.temperature, self.observer_info.pressure, 
                                              days, criterion)

        return visibilities
