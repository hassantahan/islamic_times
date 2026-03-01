"""
Module for Moon astronomical calculations.

This module provides functions to:
  - Calculate the Moon's position (declination, right ascension, altitude, azimuth, etc.).
  - Compute the Moon's illumination percentage.
  - Estimate the moonset time for a given observer and date.
  - Predict the visibility of the new crescent moon and classify it as visible or invisible.
  - Calculate the Moon's phase at a given date.

This module should not be used directly unless the user is familiar with the underlying calculations.

References:
  - Meeus, Jean. *Astronomical Algorithms*, 2nd ed. VA: Willmann-Bell, 1998.
  - Odeh, Mohammad Sh. "New criterion for lunar crescent visibility." Experimental astronomy 18, no. 1 (2004): 39-64.
  - Yallop, Bernard D. "A method for predicting the first sighting of the new Crescent Moon." RGO NAO Technical Note 69 (1997).
"""

from islamic_times import sun_equations as se
from islamic_times import time_equations as te
from islamic_times import calculation_equations as ce
from islamic_times.it_dataclasses import Angle, DateTimeInfo, Distance, DistanceUnits, ObserverInfo, RightAscension
from datetime import datetime, timedelta
from dataclasses import dataclass, replace
from typing import List, Sequence, Tuple
from warnings import warn
import numpy as np
import math

__MOON_NUTATION_ARGUMENTS_LR = [
#	D		 M		 M'		 F
	0,		 0,		 1,		 0,
	2,		 0,		-1,		 0,
	2,		 0,		 0,		 0,
	0,		 0,		 2,		 0,
	0,		 1,		 0,		 0,
	0,		 0,		 0,		 2,
	2,		 0,		-2,		 0,
	2,		-1,		-1,		 0,
	2,		 0,		 1,		 0,
	2,		-1,		 0,		 0,
	0,		 1,		-1,		 0,
	1,		 0,		 0,		 0,
	0,		 1,		 1,		 0,
	2,		 0,		 0,		-2,
	0,		 0,		 1,		 2,
	0,		 0,		 1,		-2,
	4,		 0,		-1,		 0,
	0,		 0,		 3,		 0,
	4,		 0,		-2,		 0,
	2,		 1,		-1,		 0,
	2,		 1,		 0,		 0,
	1,		 0,		-1,		 0,
	1,		 1,		 0,		 0,
	2,		-1,		 1,		 0,
	2,		 0,		 2,		 0,
	4,		 0,		 0,		 0,
	2,		 0,		-3,		 0,
	0,		 1,		-2,		 0,
	2,		 0,		-1,		 2,
	2,		-1,		-2,		 0,
	1,		 0,		 1,		 0,
	2,		-2,		 0,		 0,
	0,		 1,		 2,		 0,
	0,		 2,		 0,		 0,
	2,		-2,		-1,		 0,
	2,		 0,		 1,		-2,
	2,		 0,		 0,		 2,
	4,		-1,		-1,		 0,
	0,		 0,		 2,		 2,
	3,		 0,		-1,		 0,
	2,		 1,		 1,		 0,
	4,		-1,		-2,		 0,
	0,		 2,		-1,		 0,
	2,		 2,		-1,		 0,
	2,		 1,		-2,		 0,
	2,		-1,		 0,		-2,
	4,		 0,		 1,		 0,
	0,		 0,		 4,		 0,
	4,		-1,		 0,		 0,
	1,		 0,		-2,		 0,
	2,		 1,		 0,		-2,
	0,		 0,		 2,		-2,
	1,		 1,		 1,		 0,
	3,		 0,		-2,		 0,
	4,		 0,		-3,		 0,
	2,		-1,		 2,		 0,
	0,		 2,		 1,		 0,
	1,		 1,		-1,		 0,
	2,		 0,		 3,		 0,
	2,		 0,		-1,		-2
]

__MOON_NUTATION_ARGUMENTS_B = [
#	D		 M		 M'		 F
	0,		 0,		 0,		 1,
	0,		 0,		 1,		 1,
	0,		 0,		 1,		-1,
	2,		 0,		 0,		-1,
	2,		 0,		-1,		 1,
	2,		 0,		-1,		-1,
	2,		 0,		 0,		 1,
	0,		 0,		 2,		 1,
	2,		 0,		 1,		-1,
	0,		 0,		 2,		-1,
	2,		-1,		 0,		-1,
	2,		 0,		-2,		-1,
	2,		 0,		 1,		 1,
	2,		 1,		 0,		-1,
	2,		-1,		-1,		 1,
	2,		-1,		 0,		 1,
	2,		-1,		-1,		-1,
	0,		 1,		-1,		-1,
	4,		 0,		-1,		-1,
	0,		 1,		 0,		 1,
	0,		 0,		 0,		 3,
	0,		 1,		-1,		 1,
	1,		 0,		 0,		 1,
	0,		 1,		 1,		 1,
	0,		 1,		 1,		-1,
	0,		 1,		 0,		-1,
	1,		 0,		 0,		-1,
	0,		 0,		 3,		 1,
	4,		 0,		 0,		-1,
	4,		 0,		-1,		 1,
	0,		 0,		 1,		-3,
	4,		 0,		-2,		 1,
	2,		 0,		 0,		-3,
	2,		 0,		 2,		-1,
	2,		-1,		 1,		-1,
	2,		 0,		-2,		 1,
	0,		 0,		 3,		-1,
	2,		 0,		 2,		 1,
	2,		 0,		-3,		-1,
	2,		 1,		-1,		 1,
	2,		 1,		 0,		 1,
	4,		 0,		 0,		 1,
	2,		-1,		 1,		 1,
	2,		-2,		 0,		-1,
	0,		 0,		 1,		 3,
	2,		 1,		 1,		-1,
	1,		 1,		 0,		-1,
	1,		 1,		 0,		 1,
	0,		 1,		-2,		-1,
	2,		 1,		-1,		-1,
	1,		 0,		 1,		 1,
	2,		-1,		-2,		-1,
	0,		 1,		 2,		 1,
	4,		 0,		-2,		-1,
	4,		-1,		-1,		-1,
	1,		 0,		 1,		-1,
	4,		 0,		 1,		-1,
	1,		 0,		-1,		-1,
	4,		-1,		 0,		-1,
	2,		-2,		 0,		 1
]

__MOON_NUTATION_COEFF_LR = [
	#	l		 	r
	6288774,	-20905355,
	1274027,	-3699111,
	658314,		-2955968,
	213618,		-569925,
	-185116,	48888,
	-114332,	-3149,
	58793,		246158,
	57066,		-152138,
	53322,		-170733,
	45758,		-204586,
	-40923,		-129620,
	-34720,		108743,
	-30383,		104755,
	15327,		10321,
	-12528,		0,
	10980,		79661,
	10675,		-34782,
	10034,		-23210,
	8548,		-21636,
	-7888,		24208,
	-6766,		30824,
	-5163,		-8379,
	4987,		-16675,
	4036,		-12831,
	3994,		-10445,
	3861,		-11650,
	3665,		14403,
	-2689,		-7003,
	-2602,		0,
	2390,		10056,
	-2348,		6322,
	2236,		-9884,
	-2120,		5751,
	-2069,		0,
	2048,		-4950,
	-1773,		4130,
	-1595,		0,
	1215,		-3958,
	-1110,		0,
	-892,		3258,
	-810,		2616,
	759,		-1897,
	-713,		-2117,
	-700,		2354,
	691,		0,
	596,		0,
	549,		-1423,
	537,		-1117,
	520,		-1571,
	-487,		-1739,
	-399,		0,
	-381,		-4421,
	351,		0,
	-340,		0,
	330,		0,
	327,		0,
	-323,		1165,
	299,		0,
	294,		0,
	0,			8752
]

__MOON_NUTATION_COEFF_B = [
	5128122,
	280602,
	277693,
	173237,
	55413,
	46271,
	32573,
	17198,
	9266,
	8822,
	8216,
	4324,
	4200,
	-3359,
	2463,
	2211,
	2065,
	-1870,
	1828,
	-1794,
	-1749,
	-1565,
	-1491,
	-1475,
	-1410,
	-1344,
	-1335,
	1107,
	1021,
	833,
	777,
	671,
	607,
	596,
	491,
	-451,
	439,
	422,
	421,
	-366,
	-351,
	331,
	315,
	302,
	-283,
	-229,
	223,
	223,
	-220,
	-220,
	-185,
	181,
	-177,
	176,
	166,
	-164,
	132,
	-119,
	115,
	107
]

@dataclass(frozen=True, slots=True)
class Moon:
	"""
	A class to compute the position of the Moon in the sky based on given astronomical parameters.
    """

	# Orbital elements
	true_longitude: Angle
	true_latitude: Angle
	geocentric_distance: Distance

	# Nutation and obliquity
	lunar_nutation: Tuple[List[float], float, float, float]
	omega: Angle
	apparent_longitude: Angle
	deltaPsi: Angle
	true_obliquity: Angle

	# Apparent coordinates
	right_ascension: RightAscension
	declination: Angle

	# Hour angles
	greenwich_hour_angle: Angle
	local_hour_angle: Angle

	# Topocentric quantities
	eh_parallax: Angle
	topocentric_ascension: RightAscension
	top_declination: Angle
	topocentric_local_hour_angle: Angle

	# Horizontal coordinates
	true_altitude: Angle
	true_azimuth: Angle
	apparent_altitude: Angle

# Chapter 47
def moon_nutation(jde: float) -> Tuple[np.ndarray, float, float, float]:
	"""
	Calculate the Moon's nutation in for a given Julian Ephemeris Day. See Chapter 47 of *Astronomical Algorithms* for more information.

	Parameters:
		jde (float): The Julian Ephemeris Day.

	Returns:
		List[float]: The nutation in longitude and obliquity (in degrees).
	"""
	warn('This particular function will no longer be supported in python. The proper function is in the C extension but cannot be called from python. Use "islamic_times.astro_core.compute_moon()" or "islamic_times.sun_equations.compute_moon()" instead.', DeprecationWarning)
	from islamic_times._legacy_py_impl import moon_equations as legacy_me
	return legacy_me.moon_nutation(jde)

def moonpos(observer_date: DateTimeInfo, observer: ObserverInfo, deltaPsi: Angle, ecliptic: Angle) -> Moon:
	"""
	Calculate the various solar positional parameters for a given Julian Ephemeris Day, ΔT, and observer coordinates. See Chapter 47 of *Astronomical Algorthims* for more information.

	Parameters:
		observer_date (DateTimeInfo): The date and time of the observer.
		observer (ObserverInfo): The observer's coordinates and elevation.

	Returns:
		Moon: A `Moon` object that contains various attributes that describe its position and parameters. 

	Notes: 
		- The temperature and pressure are used for atmospheric refraction calculations. Currently, this feature is enabled.
	"""

	import islamic_times.astro_core as fast_astro
	the_moon: Moon = fast_astro.compute_moon(observer_date.jde, observer_date.deltaT, 
						 observer.latitude.decimal, observer.longitude.decimal, 
						 observer.elevation.in_unit(DistanceUnits.METRE), 0, 0, 
						 deltaPsi.decimal, ecliptic.decimal)

	return the_moon

# Known limitation:
# this wrapper returns the four phases following the nearest new moon event
# from the native implementation.
# Chapter 49
def next_phases_of_moon_utc(date: datetime) -> Tuple[datetime, datetime, datetime, datetime]:
	"""
	Calculate the next four phases of the Moon (New Moon, First Quarter, Full Moon, Last Quarter) after a given date. See Chapter 49 of *Astronomical Algorithms* for more information.

	Parameters:
		date (datetime): The date to calculate the next phases from.

	Returns:
		Tuple ([datetime, datetime, datetime, datetime]): The next four phases of the Moon.
	"""
	import islamic_times.astro_core as fast_astro
	return fast_astro.next_phases_of_moon_utc(date)

# Refer to Chapter 48 of AA
def moon_illumination(sun_dec: Angle, sun_ra: RightAscension, moon_dec: Angle, moon_ra: RightAscension, sun_earth_distance: Distance, moon_earth_distance: Distance) -> float:
	"""
	Calculate the fraction of the Moon illuminated for a given date. See Chapter 48 of *Astronomical Algorithms* for more information.

	Parameters:
		sun_dec (Angle): The Sun's declination.
		sun_ra (RightAscension): The Sun's Right Ascension.
		moon_dec (Angle): The Moon's declination.
		moon_ra (RightAscension): The Moon's Right Ascension.
		sun_earth_distance (Distance): The Sun-Earth distance.
		moon_earth_distance (Distance): The Moon-Earth distance.

	Returns:
		float: The fraction of the Moon illuminated.
	"""
	# Eqs. 48.2
	cos_psi: float = math.sin(sun_dec.radians) * math.sin(moon_dec.radians) + math.cos(sun_dec.radians) * math.cos(moon_dec.radians) * math.cos(sun_ra.radians - moon_ra.radians)
	psi: Angle = Angle(math.degrees(math.acos(cos_psi)))
	sin_psi: float = math.sin(psi.radians)
	
	# Eq. 48.3
	phase_angle = Angle(math.degrees(math.atan2(sun_earth_distance.to(DistanceUnits.AU).value * sin_psi, moon_earth_distance.to(DistanceUnits.AU).value - sun_earth_distance.to(DistanceUnits.AU).value * cos_psi)))
	
	# Eq. 48.1
	fraction_illuminated = (1 + math.cos(phase_angle.radians)) / 2
	return fraction_illuminated

def _split_solar_nutation_triplets(
    sun_nutation: Sequence[float] | float | None,
) -> tuple[List[float] | None, List[float] | None]:
	"""Split optional sun nutation buffers into two 3-value lists.

	Passing ``None`` (or legacy ``np.inf``) requests native auto-resolution.
	"""
	if sun_nutation is None:
		return None, None

	if isinstance(sun_nutation, (int, float)) and math.isinf(float(sun_nutation)):
		return None, None

	values = list(sun_nutation)
	if len(values) < 6:
		raise ValueError("`sun_nutation` must contain at least 6 values: first 3 deltaPsi and last 3 true_obliquity.")

	return values[:3], values[-3:]


def find_moon_transit(
	observer_date: DateTimeInfo,
	observer: ObserverInfo,
	sun_nutation: Sequence[float] | float | None = None,
) -> datetime:
	"""
	Calculate transit of the moon (specifically culmination) for a given date and observer coordinates. See Chapter 15 of *Astronomical Algorithms* for more information.

	Parameters:
		observer_date (DateTimeInfo): The date and time of the observer.
		observer (ObserverInfo): The observer's coordinates and elevation.
		sun_nutation (Sequence[float] | float | None): Optional precomputed
			solar nutation buffers. Pass ``None`` (default) to auto-resolve in C.

	Returns:
		datetime: The time of moonset.
	"""
	import islamic_times.astro_core as fast_astro

	delPsi, true_obliquity = _split_solar_nutation_triplets(sun_nutation)

	moon_transit: datetime = fast_astro.find_moon_transit(observer_date.jd, observer_date.deltaT, 
									observer.latitude.decimal, observer.longitude.decimal, 
									observer.elevation.in_unit(DistanceUnits.METRE), observer.temperature, observer.pressure, 
									observer_date.utc_offset, delPsi, true_obliquity)
	
	return moon_transit.replace(tzinfo=observer_date.date.tzinfo)

# Refer to Chapter 15 of AA
def moonrise_or_moonset(observer_date: DateTimeInfo, observer: ObserverInfo, rise_or_set: str = 'set') -> datetime | float:
	"""
	Calculate the time of moonset for a given date and observer coordinates. See Chapter 15 of *Astronomical Algorithms* for more information.

	Parameters:
		observer_date (DateTimeInfo): The date and time of the observer.
		observer (ObserverInfo): The observer's coordinates and elevation.
		rise_or_set (str): Specify whether to find the moonrise or moonset time. Default is 'set'.

	Returns:
		datetime: The time of moonset.
	"""
	warn('This particular function will no longer be supported in python. The proper function is in the C extension but cannot be called from python. Use "islamic_times.astro_core.find_proper_moontime()" or "islamic_times.sun_equations.find_proper_moontime()" instead.', DeprecationWarning)
	from islamic_times._legacy_py_impl import moon_equations as legacy_me
	return legacy_me.moonrise_or_moonset(observer_date, observer, rise_or_set)

# Normalize moonrise/moonset output so it matches the observer-local calendar day.
def find_proper_moontime(
	observer_date: DateTimeInfo,
	observer: ObserverInfo,
	rise_or_set: str = 'set',
	sun_nutation: Sequence[float] | float | None = None,
) -> datetime:
	"""
	Determines the proper local time for a setting or rising moon. It finds the time that corresponds to the reference date given.

	Parameters:
		observer_date (DateTimeInfo): The date and time of the observer.
		observer (ObserverInfo): The observer's coordinates and elevation.
		rise_or_set (str): Find either the setting or rising option. Default is set to 'set'.
		sun_nutation (Sequence[float] | float | None): Optional precomputed
			nutation/obliquity buffers. Pass ``None`` (default) to auto-resolve in C.

	Returns:
		datetime: The date and time of the moon event.

	Raises:
		ValueError: If `rise_or_set` is not set correctly to either 'rise' or 'set'.
		ArithmeticError: If the moon event does not exist for the given location at the given date & time.
    """
	import islamic_times.astro_core as fast_astro

	if rise_or_set not in ['rise', 'set', 'moonrise', 'moonset']:
		raise ValueError("Invalid value for rise_or_set. Please use 'rise' or 'set'.")
	
	if rise_or_set in ["set", "moonset"]:
		event = 's'
	else:
		event = 'r'

	delPsi, true_obliquity = _split_solar_nutation_triplets(sun_nutation)

	try:
		moontime: datetime = fast_astro.find_proper_moontime(
			observer_date.jd,
			observer_date.deltaT,
			observer.latitude.decimal,
			observer.longitude.decimal,
			observer.elevation.in_unit(DistanceUnits.METRE),
			observer.temperature,
			observer.pressure,
			observer_date.utc_offset,
			delPsi,
			true_obliquity,
			event,
		)
	except ValueError as exc:
		raise ArithmeticError("Moon event does not exist for the given location at the given date & time.") from exc
	
	return moontime.replace(tzinfo=observer_date.date.tzinfo)

# Visibility calculations either:
# Criterion 0: Odeh, 2006
# Criterion 1: HMNAO TN No. 69, a.k.a. Yallop, 1997
# Criterion 2: Shaukat, n.d. (Yallop q-values with Shaukat classification thresholds)
def calculate_visibility(sun_az: Angle, sun_alt: Angle, moon_az: Angle, moon_alt: Angle, moon_pi: Angle, criterion: int = 0) -> float:
	"""
	Calculate the visibility of the Moon's crescent for a given date and observer coordinates.

	Parameters:
		sun_az (Angle): The Sun's azimuth.
		sun_alt (Angle): The Sun's altitude.
		moon_az (Angle): The Moon's azimuth.
		moon_alt (Angle): The Moon's altitude.
		moon_pi (Angle): The Moon's parallax.
		criterion (int): The criterion of visibility calculation to use (0, 1, or 2).
			Default is 0 which uses Odeh, 2006. When set to 1 or 2, it uses
			the HMNAO TN No. 69 (Yallop, 1997) and Shaukat, n.d. q-value kernel respectively.
	
	Returns:
		float: The visibility of the Moon's crescent.
	"""
	warn('This particular function will no longer be supported in python. The proper function is in the C extension but cannot be called from python. Use "islamic_times.islamic_times.ITLocation.visibilities()" instead.', DeprecationWarning)
	from islamic_times._legacy_py_impl import moon_equations as legacy_me
	return legacy_me.calculate_visibility(sun_az, sun_alt, moon_az, moon_alt, moon_pi, criterion)

# Shaukat criteria mapping adapted from:
# https://moonsighting.com/faq_ms.html#Criteria
def calculate_visibility_shaukat(sun_az: Angle, sun_long: Angle, moon_long: Angle, moon_lat: Angle, moon_az: Angle, moon_pi: Angle, moon_illumin: float):
	"""Compute crescent-visibility q-value using a Shaukat-style criterion.

	Parameters
	----------
	sun_az : Angle
		Sun azimuth.
	sun_long : Angle
		Sun ecliptic longitude.
	moon_long : Angle
		Moon ecliptic longitude.
	moon_lat : Angle
		Moon ecliptic latitude.
	moon_az : Angle
		Moon azimuth.
	moon_pi : Angle
		Moon equatorial horizontal parallax.
	moon_illumin : float
		Fractional lunar illumination in the range [0, 1].

	Returns
	-------
	float
		Dimensionless visibility score. Larger values indicate better visibility.

	Notes
	-----
	This Python helper is deprecated in favor of native visibility routines used
	by ``ITLocation.visibilities()``.
	"""
	warn('This particular function will no longer be supported in python. The proper function is in the C extension but cannot be called from python. Use "islamic_times.islamic_times.ITLocation.visibilities()" instead.', DeprecationWarning)
	from islamic_times._legacy_py_impl import moon_equations as legacy_me
	return legacy_me.calculate_visibility_shaukat(sun_az, sun_long, moon_long, moon_lat, moon_az, moon_pi, moon_illumin)

# Classification according to Odeh, Yallop, or Shaukat thresholds.
def classify_visibility(q: float, criterion: int = 1) -> str:
	"""
	Classify the visibility of the Moon's crescent based on the given q value.

	Parameters:
		q (float): The q value.
		criterion (int): The criterion of visibility classification to use (0, 1, or 2).

	Returns:
		str: The classification of the Moon's crescent visibility.
	"""
	warn('This particular function will no longer be supported in python. The proper function is in the C extension but cannot be called from python. Use "islamic_times.islamic_times.ITLocation.visibilities()" instead.', DeprecationWarning)
	from islamic_times._legacy_py_impl import moon_equations as legacy_me
	return legacy_me.classify_visibility(q, criterion)
