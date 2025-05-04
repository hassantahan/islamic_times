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
from islamic_times.it_dataclasses import *
from datetime import datetime, timedelta
from dataclasses import dataclass, replace
from typing import List, Tuple
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
def moon_nutation(jde: float) -> Tuple[float, float, float]:
	"""
	Calculate the Moon's nutation in for a given Julian Ephemeris Day. See Chapter 47 of *Astronomical Algorithms* for more information.

	Parameters:
		jde (float): The Julian Ephemeris Day.

	Returns:
		List[float]: The nutation in longitude and obliquity (in degrees).
	"""
	warn('This particular function will no longer be supported in python. The proper function is in the C extension but cannot be called from python. Use "islamic_times.astro_core.compute_moon()" or "islamic_times.sun_equations.compute_moon()" instead.', DeprecationWarning)

	t = (jde - te.J2000) / te.JULIAN_CENTURY
	t2 = t ** 2
	t3 = t ** 3
	t4 = t ** 4

	fundamental_arguments = np.array([
		297.8501921 + 445267.1114034 * t - 0.0018819 * t2 + t3 / 545868 - t4 / 113065000, 	# D
		357.5291092 + 35999.0502909 * t - 0.0001536 * t2 + t3 / 24490000, 					# M
		134.9633964 + 477198.8675055 * t + 0.0087414 * t2 + t3 / 69699 - t4 / 14712000, 	# M'						
		93.2720950 + 483202.0175233 * t - 0.0036539 * t2 - t3 / 3526000 + t4 / 863310000, 	# F
		218.3164477 + 481267.88123421 * t - 0.0015786 * t2 + t3 / 538841 - t4 / 65194000 	# L'
	]) % 360
	
	a = np.array([
		119.75 + 131.849 * t,
		53.09 + 479264.290 * t,
		313.45 + 481266.484 * t
	]) % 360

	sum_l = 0
	sum_r = 0
	sum_b = 0

	eccentricity = 1 - 0.002516 * t - 0.0000074 * t2

	args_lr = np.array(__MOON_NUTATION_ARGUMENTS_LR).reshape(-1, 4)
	coeff_lr = np.array(__MOON_NUTATION_COEFF_LR).reshape(-1, 2)

	args_b = np.array(__MOON_NUTATION_ARGUMENTS_B).reshape(-1, 4)
	coeff_b = np.array(__MOON_NUTATION_COEFF_B)

	temp_lr = np.dot(args_lr, fundamental_arguments[:4])
	temp_b = np.dot(args_b, fundamental_arguments[:4])

	eccentricity_exponent = np.abs(args_lr[:, 1])
	eccentricity_comp_lr = np.where(eccentricity_exponent == 0, 1, eccentricity ** eccentricity_exponent)

	eccentricity_exponent_b = np.abs(args_b[:, 1])
	eccentricity_comp_b = np.where(eccentricity_exponent_b == 0, 1, eccentricity ** eccentricity_exponent_b)

	sum_l = np.sum(eccentricity_comp_lr * coeff_lr[:, 0] * ce.sin(temp_lr))
	sum_r = np.sum(eccentricity_comp_lr * coeff_lr[:, 1] * np.cos(np.radians(temp_lr)))
	sum_b = np.sum(eccentricity_comp_b * coeff_b * ce.sin(temp_b))

	sum_l += 3958 * ce.sin(a[0]) 													+ \
			 1962 * ce.sin(fundamental_arguments[4] - fundamental_arguments[3]) 	+ \
			 318 * ce.sin(a[1])

	sum_b += -2235 * ce.sin(fundamental_arguments[4])						 		+ \
			 382 * ce.sin(a[2])						 								+ \
			 175 * ce.sin(a[0] - fundamental_arguments[3])					 		+ \
			 175 * ce.sin(a[0] + fundamental_arguments[3])					 		+ \
			 127 * ce.sin(fundamental_arguments[4] - fundamental_arguments[2]) 		- \
			 115 * ce.sin(fundamental_arguments[4] + fundamental_arguments[2])

	return (fundamental_arguments, sum_l, sum_b, sum_r)

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

# TODO: Properly get the next phases instead of all the phases after the next new moon ==> Using round (temp fix)
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
def moon_illumination(sun_dec: Angle, sun_ra: Angle, moon_dec: Angle, moon_ra: Angle, sun_earth_distance: Distance, moon_earth_distance: Distance) -> float:
	"""
	Calculate the fraction of the Moon illuminated for a given date. See Chapter 48 of *Astronomical Algorithms* for more information.

	Parameters:
		sun_dec (Angle): The Sun's declination.
		sun_ra (Angle): The Sun's Right Ascension.
		moon_dec (Angle): The Moon's declination.
		moon_ra (Angle): The Moon's Right Ascension.
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

def find_moon_transit(observer_date: DateTimeInfo, observer: ObserverInfo, sun_nutation = np.inf) -> datetime:
	"""
	Calculate transit of the moon (specifically culmination) for a given date and observer coordinates. See Chapter 15 of *Astronomical Algorithms* for more information.

	Parameters:
		observer_date (DateTimeInfo): The date and time of the observer.
		observer (ObserverInfo): The observer's coordinates and elevation.

	Returns:
		datetime: The time of moonset.
	"""
	import islamic_times.astro_core as fast_astro

	# To be later explained
	if (sun_nutation == np.inf):
		delPsi = [-123456.0, -123456.0, -123456.0]
		true_obliquity = [-123456.0, -123456.0, -123456.0]
	else:
		delPsi = list(sun_nutation[:3])
		true_obliquity = list(sun_nutation[-3:])

	moon_transit: datetime = fast_astro.find_moon_transit(observer_date.jd, observer_date.deltaT, 
									observer.latitude.decimal, observer.longitude.decimal, 
									observer.elevation.in_unit(DistanceUnits.METRE), observer.temperature, observer.pressure, 
									observer_date.utc_offset, delPsi, true_obliquity)
	
	return moon_transit.replace(tzinfo=observer_date.date.tzinfo)

# Refer to Chapter 15 of AA
def moonrise_or_moonset(observer_date: DateTimeInfo, observer: ObserverInfo, rise_or_set: str = 'set') -> datetime:
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

	# Validate the input
	if rise_or_set not in ['rise', 'set', 'moonrise', 'moonset']:
		raise ValueError("Invalid value for rise_or_set. Please use 'rise' or 'set'.")

	# First find the Year Month Day at UT 0h from JDE
	ymd = datetime(observer_date.date.year, observer_date.date.month, observer_date.date.day)
	new_jd = te.gregorian_to_jd(observer_date.date) - te.fraction_of_day(observer_date.date)
	new_deltaT = te.delta_t_approx(ymd.year, ymd.month)

	# Calculate new moon params with the new_jd
	moon_params: List[Moon] = []
	for i in range(3):
		ymd_temp = te.jd_to_gregorian(new_jd + i - 1, observer_date.utc_offset)
		delT_temp = te.delta_t_approx(ymd_temp.year, ymd_temp.month)
		sun_params = se.sunpos(replace(observer_date, date=ymd_temp, jd=(new_jd + i - 1), deltaT=delT_temp), observer)
		delPsi = sun_params.nutation[0]
		moon_params.append(
						moonpos(
							replace(observer_date, date=ymd_temp, jd=(new_jd + i - 1), deltaT=delT_temp),
							observer, delPsi, sun_params.true_obliquity
						)
					)

	# Compute H0: the hour angle corresponding to the Moon's apparent altitude for rise/set
	h_zero: Angle = Angle(0.7275 * moon_params[1].eh_parallax.decimal - 0.566667)

	cosH_zero: float = (math.sin(h_zero.radians) - math.sin(observer.latitude.radians) * math.sin(moon_params[1].declination.radians)) / (
		math.cos(observer.latitude.radians) * math.cos(moon_params[1].declination.radians))

	# Moon never rises or sets
	if abs(cosH_zero) < 1:
		H_zero = Angle(math.degrees(math.acos(cosH_zero)))
	else:
		return math.inf

	if math.isnan(H_zero.decimal):
		return math.inf

	# Greenwich Mean Sidereal Time (GMST)
	sidereal_time: Angle = te.greenwich_mean_sidereal_time(new_jd)

	# Compute transit
	m0: float = (moon_params[1].right_ascension.decimal_degrees.decimal - observer.longitude.decimal - sidereal_time.decimal) / 360

	# Choose rise or set
	if rise_or_set in ['rise', 'moonrise']:
		m_event = m0 - H_zero.decimal / 360
	else:
		m_event = m0 + H_zero.decimal / 360

	# Iteratively refine the result
	for _ in range(3):
		theta_event: Angle = Angle((sidereal_time.decimal + 360.985647 * m_event) % 360)
		n_event: float = m_event + new_deltaT / 86400

		# Interpolated declination and right ascension
		interp_dec: Angle = Angle(ce.interpolation(n_event,
										moon_params[0].declination.decimal,
										moon_params[1].declination.decimal,
										moon_params[2].declination.decimal)
									)

		interp_ra = RightAscension(ce.interpolation(n_event,
										moon_params[0].right_ascension.decimal_degrees.decimal,
										moon_params[1].right_ascension.decimal_degrees.decimal,
										moon_params[2].right_ascension.decimal_degrees.decimal) / 15
									)

		# Local hour angle
		local_hour_angle = Angle((theta_event.decimal - (-observer.longitude.decimal) - interp_ra.decimal_degrees.decimal) % 360)

		# Altitude of the Moon
		moon_alt = Angle(math.degrees(math.asin(
			math.sin(observer.latitude.radians) * math.sin(interp_dec.radians) +
			math.cos(observer.latitude.radians) * math.cos(interp_dec.radians) *
			math.cos(local_hour_angle.radians)
		)))

		# Correction to m_event
		deltaM = (moon_alt.decimal - h_zero.decimal) / (
			360 * math.cos(interp_dec.radians) * math.cos(observer.latitude.radians) * math.sin(local_hour_angle.radians))

		m_event += deltaM

	# Compute final moonrise or moonset time
	event_dt = datetime(ymd.year, ymd.month, ymd.day) + timedelta(days=m_event) - timedelta(hours=observer_date.utc_offset)

	return event_dt

# This is necessary because UTC offsets for coords not near UTC, but also not using local TZ.
def find_proper_moontime(observer_date: DateTimeInfo, observer: ObserverInfo, rise_or_set: str = 'set', sun_nutation = np.inf) -> datetime:
	"""
	Determines the proper local time for a setting or rising moon. It finds the time that corresponds to the reference date given.

	Parameters:
		observer_date (DateTimeInfo): The date and time of the observer.
		observer (ObserverInfo): The observer's coordinates and elevation.
		rise_or_set (str): Find either the setting or rising option. Default is set to 'set'.
		sun_nutation (List[float]): The nutation in longitude and obliquity (in degrees). Default is `np.inf`. This should only be used if user understands the C code.

	Returns:
		datetime: The date and time of the moon event.

	Raises:
		ValueError: If `rise_or_set` is not set correctly to either 'rise' or 'set'.
		ArithmeticError: If the moon event does not exist for the given location at the given date & time.
    """
	import islamic_times.astro_core as fast_astro

	if rise_or_set not in ['rise', 'set', 'moonrise', 'moonset']:
		raise ValueError("Invalid value for rise_or_set. Please use 'rise' or 'set'.")
	
	if rise_or_set in ["set", "sunset"]:
		event = 's'
	else:
		event = 'r'

	# To be later explained
	if (sun_nutation == np.inf):
		delPsi = [-123456.0, -123456.0, -123456.0]
		true_obliquity = [-123456.0, -123456.0, -123456.0]
	else:
		delPsi = list(sun_nutation[:3])
		true_obliquity = list(sun_nutation[-3:])

	try:
		moontime: datetime = fast_astro.find_proper_moontime(observer_date.jd, observer_date.deltaT, 
										observer.latitude.decimal, observer.longitude.decimal, 
										observer.elevation.in_unit(DistanceUnits.METRE), observer.temperature, observer.pressure, 
										observer_date.utc_offset, delPsi, true_obliquity, event)
	except:
		raise ArithmeticError("Moon event does not exist for the given location at the given date & time.")
	
	return moontime.replace(tzinfo=observer_date.date.tzinfo)

# Visibility calculations either:
# Criterion 0: Odeh, 2006
# Criterion 1: HMNAO TN No. 69, a.k.a. Yallop, 1997
def calculate_visibility(sun_az: Angle, sun_alt: Angle, moon_az: Angle, moon_alt: Angle, moon_pi: Angle, criterion: int = 0) -> float:
	"""
	Calculate the visibility of the Moon's crescent for a given date and observer coordinates.

	Parameters:
		sun_az (Angle): The Sun's azimuth.
		sun_alt (Angle): The Sun's altitude.
		moon_az (Angle): The Moon's azimuth.
		moon_alt (Angle): The Moon's altitude.
		moon_pi (Angle): The Moon's parallax.
		criterion (int): The criterion of visibility calculation to use (0 or 1). Default is 0 which uses Odeh, 2006. When set to 1, it uses HMNAO TN No. 69 (a.k.a. Yallop, 1997).
	
	Returns:
		float: The visibility of the Moon's crescent.
	"""
	warn('This particular function will no longer be supported in python. The proper function is in the C extension but cannot be called from python. Use "islamic_times.islamic_times.ITLocation.visibilities()" instead.', DeprecationWarning)

	daz = Angle(abs(sun_az.decimal - moon_az.decimal))
	arcv = Angle(abs(sun_alt.decimal - moon_alt.decimal))
	arcl = Angle(math.degrees(math.acos(math.cos(daz.radians) * math.cos(arcv.radians))))

	semi_diameter = 0.27245 * (moon_pi.decimal * 60)
	semi_diameter_prime = semi_diameter * (1 + math.sin(moon_alt.radians) * math.sin(moon_pi.radians))

	w_prime = semi_diameter_prime * (1 - math.cos(arcl.radians))

	if criterion == 0:
		q_value = (arcv.decimal - (-0.1018 * w_prime ** 3 + 0.7319 * w_prime ** 2 - 6.3226 * w_prime + 7.1651))
	else:
		q_value = (arcv.decimal - (11.8371 - 6.3226 * w_prime + 0.7319 * w_prime ** 2 - 0.1018 * w_prime ** 3)) / 10

	return q_value

# According to Shaukat (yet to find his paper; reference: https://moonsighting.com/faq_ms.html#Criteria)
def calculate_visibility_shaukat(sun_az: Angle, sun_long: Angle, moon_long: Angle, moon_lat: Angle, moon_az: Angle, moon_pi: Angle, moon_illumin: float):
	warn('This particular function will no longer be supported in python. The proper function is in the C extension but cannot be called from python. Use "islamic_times.islamic_times.ITLocation.visibilities()" instead.', DeprecationWarning)

	semi_diameter: float = 0.27245 * (moon_pi.decimal * 60)
	width: float = 2 * semi_diameter * moon_illumin

	arcv: Angle = Angle(math.degrees(math.acos(
		(math.cos(moon_long.radians - sun_long.radians) * math.cos(moon_lat.radians)) / \
		math.cos(sun_az.radians - moon_az.radians)
	)))

	q_value: float = (arcv.decimal - (11.8371 - 6.3226 * width + 0.7319 * width ** 2 - 0.1018 * width ** 3)) / 10

	return q_value

# Classification according to Odeh, 2006 or HMNAO TN No.69
def classify_visibility(q: float, criterion: int = 1) -> str:
	"""
	Classify the visibility of the Moon's crescent based on the given q value.

	Parameters:
		q (float): The q value.
		criterion (int): The criterion of visibility classification to use (0 or 1). Default is 1. When set to 1, it uses HMNAO TN No. 69 (a.k.a. Yallop, 1997).

	Returns:
		str: The classification of the Moon's crescent visibility.
	"""
	warn('This particular function will no longer be supported in python. The proper function is in the C extension but cannot be called from python. Use "islamic_times.islamic_times.ITLocation.visibilities()" instead.', DeprecationWarning)
	
	if q == -999: 
		return "Moonset before the new moon."
	elif q == -998: 
		return "Moonset before sunset."
	# Only for extreme latitudes
	elif q == -997:
		return "Moonset & Sunset don't exist."
	elif q == -996:
		return "Sunset doesn't exist."
	elif q == -995:
		return "Moonset doesn't exist."

	if criterion == 0:
		if q >= 5.65:
			return "A: Crescent is visible by naked eyes."
		elif 5.65 > q >= 2:
			return "B: Crescent is visible by optical aid, and it could be seen by naked eyes."
		elif 2 > q >= -0.96:
			return "C: Crescent is visible by optical aid only."
		elif -0.96 > q:
			return "D: Crescent is not visible even by optical aid."
		else:
			raise ValueError("Invalid q value. Must be a float.")
	else:
		if q > 0.216:
			return "A: Easily visible."
		elif 0.216 >= q > -0.014:
			return "B: Visible under perfect conditions."
		elif -0.014 >= q > -0.160:
			return "C: May need optical aid to find crescent."
		elif -0.160 >= q > -0.232:
			return "D: Will need optical aid to find crescent."
		elif -0.232 >= q > -0.293:
			return "E: Not visible with a [conventional] telescope."
		elif -0.293 >= q:
			return "F: Not visible; below the Danjon limit."
		else:
			raise ValueError("Invalid q value. Must be a float.")