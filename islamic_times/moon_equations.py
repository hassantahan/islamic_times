'''
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
'''

from islamic_times import sun_equations as se
from islamic_times import time_equations as te
from islamic_times import calculation_equations as ce
from datetime import datetime, timedelta
from typing import List, Tuple
import numpy as np
import math

__moon_nutation_arguments_lr = [
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

__moon_nutation_arguments_b = [
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

__moon_nutation_coeff_lr = [
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

__moon_nutation_coeff_b = [
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

__moon_phase_corrections_coeff = [
	-0.40720,	-0.40614,	-0.62801,
	 0.17241,	 0.17302,	 0.17172,
	 0.01608,	 0.01614,	-0.01183,
	 0.01039,	 0.01043,	 0.00862,
	 0.00739,	 0.00734,	 0.00804,
	-0.00514,	-0.00515,	 0.00454,
	 0.00208,	 0.00209,	 0.00204,
	-0.00111,	-0.00111,	-0.00180,
	-0.00057,	-0.00057,	-0.00070,
	 0.00056,	 0.00056,	-0.00040,
	-0.00042,	-0.00042,	-0.00034,
	 0.00042,	 0.00042,	 0.00032,
	 0.00038,	 0.00038,	 0.00032,
	-0.00024,	-0.00024,	-0.00028,
	-0.00017,	-0.00017,	 0.00027,
	-0.00007,	-0.00007,	-0.00017,
	 0.00004,	 0.00004,	-0.00005,
	 0.00004,	 0.00004,	 0.00004,
	 0.00003,	 0.00003,	-0.00004,
	 0.00003,	 0.00003,	 0.00004,
	-0.00003,	-0.00003,	 0.00003,
	 0.00003,	 0.00003,	 0.00003,
	-0.00002,	-0.00002,	 0.00002,
	-0.00002,	-0.00002,	 0.00002,
	 0.00002,	 0.00002,	-0.00002
]

__moon_phase_corrections_arg = [
	#New & Full Moon
	[
	#	E	 M	 M'	 F
	#	Note: The number for E is it's power, not coeff
		0,	 0,	  1,	 0,
		1,	 1,	  0,	 0,
		0,	 0,	  2,	 0,
		0,	 0,	  0,	 2,
		1,	-1,	  1,	 0,
		1,	 1,	  1,	 0,
		2,	 2,	  0,	 0,
		0,	 0,	  1,	-2,
		0,	 0,	  1,	 2,
		1,	 1,	  2,	 0,
		0,	 0,	  3,	 0,
		1,	 1,	  0,	 2,
		1,	 1,	  0,	-2,
		1,	-1,	  2,	 0,
		9,	 9,	  9,	 9,		#OMEGA
		0,	 2,	  1,	 0,
		0,	 0,	  2,	-2,
		0,	 3,	  0,	 0,
		0,	 1,	  1,	-2,
		0,	 0,	  2,	 2,
		0,	 1,	  1,	 2,
		0,	-1,	  1,	 2,
		0,	-1,	  1,	-2,
		0,	 1,	  3,	 0,
		0,	 0,	  4,	 0
	],

	#First & Last Quarter Moon
	[
	#	E	 M	 M'	 F
	#	Note: The number for E is it's power, not coeff
		0,	 0,	  1,	 0,
		1,	 1,	  0,	 0,
		1,	 1,	  1,	 0,
		0,	 0,	  2,	 0,
		0,	 0,	  0,	 2,
		1,	-1,	  1,	 0,
		2,	 2,	  0,	 0,
		0,	 0,	  1,	-2,
		0,	 0,	  1,	 2,
		0,	 0,	  3,	 0,
		1,	-1,	  2,	 0,
		1,	 1,	  0,	 2,
		1,	 1,	  0,	-2,
		2,	 2,	  1,	 0,
		1,	 1,	  2,	 0,
		9,	 9,	  9,	 9,		#OMEGA
		0,	-1,	  1,	-2,
		0,	 0,	  2,	 2,
		0,	 1,	  1,	 2,
		0,	-2,	  1,	 0,
		0,	 1,	  1,	-2,
		0,	 3,	  0,	 0,
		0,	 0,	  2,	-2,
		0,	-1,	  1,	 2,
		0,	 1,	  3,	 0
	]
]

__a_sin_term_phases_coeff = [
	[299.77,	0.1074080,	-0.009173],
	[251.88,	0.0163210],
	[251.83,	26.651886],
	[349.42,	36.412478],
	[84.660,	18.206239],
	[141.74,	53.303771],
	[207.14,	2.4537320],
	[154.84,	7.3068600],
	[34.520,	27.261239],
	[207.19,	0.1218240],
	[291.34,	1.8443790],
	[161.72,	24.198154],
	[239.56,	25.513099],
	[331.55,	3.5925180]
]

__a_coeffs = [
	0.000325,	0.000165,	0.000164,
	0.000126,	0.000110,	0.000062,
	0.000060,	0.000056,	0.000047,
	0.000042,	0.000040,	0.000037,
	0.000035,	0.000023
]

# Chapter 47
def moon_nutation(jde: float) -> Tuple[float, float, float]:
	'''
	Calculate the Moon's nutation in for a given Julian Ephemeris Day. See Chapter 47 of *Astronomical Algorithms* for more information.

	Parameters:
		jde (float): The Julian Ephemeris Day.

	Returns:
		List[float]: The nutation in longitude and obliquity (in degrees).
	'''

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

	args_lr = np.array(__moon_nutation_arguments_lr).reshape(-1, 4)
	coeff_lr = np.array(__moon_nutation_coeff_lr).reshape(-1, 2)

	args_b = np.array(__moon_nutation_arguments_b).reshape(-1, 4)
	coeff_b = np.array(__moon_nutation_coeff_b)

	temp_lr = np.dot(args_lr, fundamental_arguments[:4])
	temp_b = np.dot(args_b, fundamental_arguments[:4])

	eccentricity_exponent = np.abs(args_lr[:, 1])
	eccentricity_comp_lr = np.where(eccentricity_exponent == 0, 1, eccentricity ** eccentricity_exponent)

	eccentricity_exponent_b = np.abs(args_b[:, 1])
	eccentricity_comp_b = np.where(eccentricity_exponent_b == 0, 1, eccentricity ** eccentricity_exponent_b)

	sum_l = np.sum(eccentricity_comp_lr * coeff_lr[:, 0] * np.sin(np.radians(temp_lr)))
	sum_r = np.sum(eccentricity_comp_lr * coeff_lr[:, 1] * np.cos(np.radians(temp_lr)))
	sum_b = np.sum(eccentricity_comp_b * coeff_b * np.sin(np.radians(temp_b)))

	sum_l += 3958 * ce.sin(np.radians(a[0])) + \
			 1962 * ce.sin(np.radians(fundamental_arguments[4] - fundamental_arguments[3])) + \
			 318 * ce.sin(np.radians(a[1]))

	sum_b += -2235 * ce.sin(np.radians(fundamental_arguments[4])) + \
			 382 * ce.sin(np.radians(a[2])) + \
			 175 * ce.sin(np.radians(a[0] - fundamental_arguments[3])) + \
			 175 * ce.sin(np.radians(a[0] + fundamental_arguments[3])) + \
			 127 * ce.sin(np.radians(fundamental_arguments[4] - fundamental_arguments[2])) - \
			 115 * ce.sin(np.radians(fundamental_arguments[4] + fundamental_arguments[2]))

	return (fundamental_arguments, sum_l, sum_b, sum_r)

def moonpos(jde: float, deltaT: float, local_latitude: float, local_longitude: float, deltaPsi: float, ecliptic: float, elev: float) -> List[float]:
	'''
	Calculate the various solar positional parameters for a given Julian Ephemeris Day, ΔT, and observer coordinates. See Chapter 47 of the Astronomical Algorthims for more information.

	Parameters:
		jde (float): The Julian Ephemeris Day.
		deltaT (float): The difference between Terrestrial Time and Universal Time (ΔT).
		local_latitude (float): The observer's latitude (°).
		local_longitude (float): The observer's longitude (°).
		deltaPsi (float): The nutation in longitude (°).
		ecliptic (float): The observer's ecliptic (°).
		elev (float): The observer's elevation above sea level (m).

	Returns:
		list (List[float]): The Moon's position in the sky. The list needs to be restructured.
	Notes: 
		- The temperature and pressure are used for atmospheric refraction calculations. Currently, this feature is enabled.
	'''

	t = (jde - te.J2000) / te.JULIAN_CENTURY

	# Calculation nutations
	nut = moon_nutation(jde)

	# Rect Coordinates + Distance
	longitude = nut[0][4] + nut[1] / 10 ** 6 + deltaPsi
	latitude = nut[2] / 10 ** 6
	distance = 385000.56 + nut[3] / 10 ** 3

	# Place in the sky (pg. 93)
	ascension = (np.rad2deg(math.atan2((ce.sin(longitude) * ce.cos(ecliptic) - ce.tan(latitude) * ce.sin(ecliptic)), ce.cos(longitude)))) % 360
	declination = np.rad2deg(math.asin(ce.sin(latitude) * ce.cos(ecliptic) + ce.cos(latitude) * ce.sin(ecliptic) * ce.sin(longitude)))
	eh_parallax = np.rad2deg(np.arcsin(te.EARTH_RADIUS_KM / distance))

	# Local Hour Angle Calculations (pg. 88 & 92)
	mean_greenwich_sidereal_time = te.greenwich_mean_sidereal_time(jde - deltaT / 86400)
	st_correction = ce.decimal_to_dms(deltaPsi)[2] * ce.cos(ecliptic) / 15
	app_greenwich_sidereal_time = mean_greenwich_sidereal_time + (st_correction / 3600)
	local_hour_angle = (app_greenwich_sidereal_time - -1 * local_longitude - ascension) % 360

	# Modify RA and Declination to their topocentric equivalents
	top_ascension, top_declination = correct_ra_dec(ascension, declination, local_hour_angle, eh_parallax, local_latitude, elev / 1000)

	# Altitude & Azimuth calculations
	altitude = np.rad2deg(np.asin(ce.sin(local_latitude) * ce.sin(declination) + ce.cos(local_latitude) * ce.cos(declination)* ce.cos(local_hour_angle)))
	azimuth = np.rad2deg(np.arctan2(-1 * ce.cos(declination) * ce.sin(local_hour_angle), ce.sin(declination) * ce.cos(local_latitude) - ce.cos(declination) * ce.sin(local_latitude) * ce.cos(local_hour_angle))) % 360

	# Correct for atmospheric refraction (taken from https://en.wikipedia.org/wiki/Atmospheric_refraction)
	refraction = 1.02 / (ce.tan(altitude + 10.3 / (altitude + 5.11))) #+ 0.0019279 - 0.000034 * elev
	altitude += refraction / 60

	# Correct for parallax
	parallax_correction = -np.rad2deg(np.deg2rad(eh_parallax) * ce.cos(altitude))
	altitude += parallax_correction
	
	return [
		longitude, 						# 0:  deg. decimal (lambda)
		latitude, 						# 1:  deg. decimal (beta)
		distance,						# 2:  km
		ecliptic,						# 3:  deg. decimal
		ascension,						# 4:  deg. decimal
		declination,					# 5:  deg. decimal
		top_ascension,					# 6:  deg. decimal
		top_declination,				# 7:  deg. decimal
		eh_parallax,					# 8:  deg. decimal
		altitude,						# 9:  deg. decimal (apparent)
		azimuth,						# 10: deg. decimal (apparent)
		mean_greenwich_sidereal_time,	# 11: deg. decimal
		local_hour_angle				# 12: deg. decimal
	]

# Fixing RA and Dec for apparency pg. 279
def correct_ra_dec(ra: float, dec: float, lha: float, parallax: float, lat: float, elev: float, dist: float = te.EARTH_RADIUS_KM) -> Tuple[float, float]:
	'''
	Correct the Moon's Right Ascension and Declination for apparent position. Currently not used for much.

	Parameters:
		ra (float): The Moon's Right Ascension (°).
		dec (float): The Moon's Declination (°).
		lha (float): The Local Hour Angle (°).
		parallax (float): The Moon's parallax (°).
		lat (float): The observer's latitude (°).
		elev (float): The observer's elevation above sea level (m).
		dist (float): The observer's distance from the Moon (km).

	Returns:
		tuple (Tuple[float, float]): The topocentric Right Ascension and Declination (°).
	'''
	a = dist
	f = 1 / 298.257
	b = a * (1 - f)

	u = np.rad2deg(np.arctan2(b * ce.tan(lat), a))
	p_sin_psi_prime = b / a * ce.sin(u) + elev / dist * ce.sin(lat)
	p_cos_psi_prime = ce.cos(u) + elev / dist * ce.cos(lat)

	temp_num = -1 * p_cos_psi_prime * ce.sin(parallax) * ce.sin(lha)
	temp_denom = ce.cos(dec) - p_cos_psi_prime * ce.sin(parallax) * ce.sin(lha)
	deltaA = np.rad2deg(np.arctan2(temp_num, temp_denom))

	temp_num = (ce.sin(dec) - p_sin_psi_prime * ce.sin(parallax)) * ce.cos(deltaA)

	ascension_prime = ra + deltaA
	declination_prime = np.rad2deg(np.arctan2(temp_num, temp_denom))

	return (ascension_prime, declination_prime)

# TODO: Properly get the next phases instead of all the phases after the next new moon ==> Using round (temp fix)
# Chapter 49
def next_phases_of_moon_utc(date: datetime) -> List[datetime]:
	'''
	Calculate the next four phases of the Moon (New Moon, First Quarter, Full Moon, Last Quarter) after a given date. See Chapter 49 of *Astronomical Algorithms* for more information.

	Parameters:
		date (datetime): The date to calculate the next phases from.

	Returns:
		list (List[datetime]): The next four phases of the Moon.
	'''

	# Find the day of the year.
	day_of_year = date.timetuple().tm_yday

	# Preserve the sign relative to year 2000.
	p_sign = np.sign(date.year - 2000)

	# Calculate k based on eq. 49.2.
	k_temp = ((abs(date.year - 2000) + p_sign * day_of_year / te.TROPICAL_YEAR) * 12.3685)
	# Determine k for each phase:
	k_array = np.array([
		p_sign * np.round(k_temp),
		p_sign * (np.round(k_temp - 0.25) + 0.25),
		p_sign * (np.round(k_temp - 0.5) + 0.5),
		p_sign * (np.round(k_temp - 0.75) + 0.75)
	])

	# Prepare output array.
	moon_phases = [None] * 4 # np.empty(4)

	# Pre-convert A-term coefficients to NumPy arrays.
	a_sin_term_coeffs = [np.array(item) for item in __a_sin_term_phases_coeff]
	a_coeffs = np.array(__a_coeffs)

	# Loop over each phase (New Moon, First Quarter, Full Moon, Last Quarter)
	for p, k in enumerate(k_array):
		# Julian centuries (eq. 49.3)
		t = k / 1236.85
		t2 = t * t
		t3 = t * t2
		t4 = t * t3

		# Mean Julian Ephemeris Date.
		jde = 2451550.09766 + 29.530588861 * k + 0.00015437 * t2 - 0.000000150 * t3 + 0.00000000073 * t4
		phase = jde  # start with JDE

		# Fundamental arguments (in degrees)
		fundamental_arguments = np.array([
			1 - 0.002516 * t - 0.0000074 * t2,                                                   # E
			2.5534 + 29.10535670 * k - 0.0000014 * t2 - 0.00000011 * t3,                         # M
			201.5643 + 385.81693528 * k + 0.0107582 * t2 + 0.00001238 * t3 - 0.000000058 * t4,   # M'
			160.7108 + 390.67050284 * k + 0.0016118 * t2 + 0.00000227 * t3 - 0.000000011 * t4,   # F
			124.7746 - 1.5637558 * k + 0.0020672 * t2 + 0.00000215 * t3                          # Omega
		])
		# Normalize each to [0, 360)
		fundamental_arguments %= 360

		# === A-Term Corrections ===
		# For each set in __a_sin_term_phases_coeff, use multipliers: [1, k, t2, ...] (if more than 2 numbers).
		a_sin_args = []
		for coeff_array in a_sin_term_coeffs:
			L = coeff_array.size
			if L == 2:
				multipliers = np.array([1, k])
			else:
				multipliers = np.concatenate(([1, k], np.full(L - 2, t2)))
			# Sum the product and reduce modulo 360.
			a_sin_val = np.dot(coeff_array, multipliers) % 360
			a_sin_args.append(a_sin_val)

		# === W Correction (for quarter phases) ===
		# (Note: here we convert degrees to radians for the trigonometric functions.)
		w = (0.00306 
				- 0.00038 * fundamental_arguments[0] * np.cos(np.radians(fundamental_arguments[1]))
				+ 0.00026 * np.cos(np.radians(fundamental_arguments[2]))
				- 0.00002 * np.cos(np.radians(fundamental_arguments[2] - fundamental_arguments[1]))
				+ 0.00002 * np.cos(np.radians(fundamental_arguments[2] + fundamental_arguments[1]))
				+ 0.00002 * np.cos(np.radians(2 * fundamental_arguments[3])))

		# === Periodic Term Corrections ===
		# Reshape the correction coefficients and arguments for vectorized operations.
		coeff_array = np.array(__moon_phase_corrections_coeff).reshape(-1, 3)
		arg_array = np.array(__moon_phase_corrections_arg[p % 2]).reshape(-1, 4)

		# Determine the shift "s" for the coefficients based on phase.
		# p == 0 (New Moon)  -> s = 0; p == 1 or 3 (First/Last Quarter) -> s = 2; p == 2 (Full Moon) -> s = 1.
		if p in (1, 3):
			s = 2
		elif p == 2:
			s = 1
		else:
			s = 0

		# For non-Omega terms (where the first number in a row is not 9):
		mask = arg_array[:, 0] != 9
		periodic_corr = 0.0
		if np.any(mask):
			# For each non-Omega term, compute:
			#   sin_argument = dot(fundamental_arguments[1:4], arg_row[1:4])
			sin_arguments = np.dot(arg_array[mask, 1:4], fundamental_arguments[1:4])
			sin_vals = np.sin(np.radians(sin_arguments))
			# Multiply by the coefficient and the eccentricity factor (fundamental_arguments[0] raised to arg_row[0])
			periodic_corr += np.sum(
				coeff_array[mask, s] * (fundamental_arguments[0] ** arg_array[mask, 0]) * sin_vals
			)
		# For Omega terms (first number equals 9), the sine argument is simply fundamental_arguments[4].
		if np.any(~mask):
			sin_omega = np.sin(np.radians(fundamental_arguments[4]))
			periodic_corr += np.sum(
				coeff_array[~mask, s] * (fundamental_arguments[0] ** 0) * sin_omega
			)
		phase += periodic_corr

		# === Additional A-Term Corrections ===
		# Use vectorized dot product (after converting the a_sin_args list to an array).
		a_term_corr = np.dot(a_coeffs, np.sin(np.radians(a_sin_args)))
		phase += a_term_corr

		# === Apply W Factor for First and Last Quarters ===
		if p == 1:
			phase += w
		elif p == 3:
			phase -= w

		# === Adjust from Terrestrial Dynamical Time (TD) to UT.
		phase -= te.delta_t_approx(date.year, date.month) / 86400

		# === Convert from JD to Gregorian.
		moon_phases[p] = te.jd_to_gregorian(phase)

	return moon_phases

# Refer to Chapter 48 of AA
def moon_illumination(sun_dec: float, sun_ra: float, moon_dec: float, moon_ra: float, sun_earth_distance: float, moon_earth_distance: float) -> float:
	'''
	Calculate the fraction of the Moon illuminated for a given date. See Chapter 48 of *Astronomical Algorithms* for more information.

	Parameters:
		sun_dec (float): The Sun's declination (°).
		sun_ra (float): The Sun's Right Ascension (°).
		moon_dec (float): The Moon's declination (°).
		moon_ra (float): The Moon's Right Ascension (°).
		sun_earth_distance (float): The Sun-Earth distance (AU).
		moon_earth_distance (float): The Moon-Earth distance (AU).

	Returns:
		float: The fraction of the Moon illuminated.
	'''
	# Eqs. 48.2
	cos_psi = ce.sin(sun_dec) * ce.sin(moon_dec) + ce.cos(sun_dec) * ce.cos(moon_dec) * ce.cos(sun_ra - moon_ra)
	psi = np.rad2deg(np.arccos(cos_psi))
	sin_psi = ce.sin(psi)
	
	# Eq. 48.3
	phase_angle = np.rad2deg(np.arctan2(sun_earth_distance * sin_psi, moon_earth_distance - sun_earth_distance * cos_psi))
	
	# Eq. 48.1
	fraction_illuminated = (1 + ce.cos(phase_angle)) / 2
	return fraction_illuminated

# Refer to Chapter 15 of AA
def calculate_moonset(date: datetime, lat: float, long: float, elev: float, utc_diff: float) -> datetime:
	'''
	Calculate the time of moonset for a given date and observer coordinates. See Chapter 15 of *Astronomical Algorithms* for more information.

	Parameters:
		date (datetime): The date to calculate the moonset for.
		lat (float): The observer's latitude (°).
		long (float): The observer's longitude (°).
		elev (float): The observer's elevation above sea level (m).
		utc_diff (float): The observer's difference from UTC (hours).

	Returns:
		datetime: The time of moonset.
	'''

	# First find the Year Month Day at UT 0h from JDE
	ymd = datetime(date.year, date.month, date.day)
	new_jd = te.gregorian_to_jd(date) - te.fraction_of_day(date)
	new_deltaT = te.delta_t_approx(ymd.year, ymd.month)
	new_jde = new_jd + new_deltaT / 86400

	# Calculate new sun and moon params with the new_jd
	moon_params = []
	for i in range(3):
		ymd_temp = te.jd_to_gregorian(new_jd + i - 1, utc_diff)
		delT_temp = te.delta_t_approx(ymd_temp.year, ymd_temp.month)
		sun_params = se.sunpos(new_jde + i - 1, delT_temp, lat, long)
		delPsi, delEps = se.sun_nutation(new_jde + i - 1)
		moon_params.append(moonpos(new_jde + i - 1, delT_temp, lat, long, delPsi, sun_params[13], elev))

	# Find H0
	h_zero = 0.7275 * moon_params[1][8] - 0.566667
	cosH_zero = (ce.sin(h_zero) - ce.sin(lat) * ce.sin(moon_params[1][5])) / (ce.cos(lat) * ce.cos(moon_params[1][5]))
	H_zero = np.rad2deg(np.arccos(cosH_zero))

	# No moonset/rise
	if math.isnan(H_zero):
		return np.inf

	# GMST
	sidereal_time = te.greenwich_mean_sidereal_time(new_jd)

	# Transit
	m0 = (moon_params[1][4] + -1 * long - sidereal_time) / 360
	if m0 < 0: 
		m0 += 1
	elif m0 > 1: 
		m0 -= 1

	# Setting
	m2 = m0 + H_zero / 360
	if m2 < 0: 
		m2 += 1
	elif m2 > 1: 
		m2 -= 1

	# Minor corrective steps
	for _ in range(3):
		little_theta_zero = (sidereal_time + 360.985647 * m2) % 360

		n = m2 + new_deltaT / 86400
		interpolated_moon_dec = ce.interpolation(n, moon_params[0][5], moon_params[1][5], moon_params[2][5])
		interpolated_moon_ra = ce.interpolation(n, moon_params[0][4], moon_params[1][4], moon_params[2][4])

		lunar_local_hour_angle = (little_theta_zero - -1 * long - interpolated_moon_ra) % 360
		moon_alt = np.rad2deg(np.arcsin(ce.sin(lat) * ce.sin(interpolated_moon_dec) + ce.cos(lat) * ce.cos(interpolated_moon_dec) * ce.cos(lunar_local_hour_angle)))
		deltaM = (moon_alt - h_zero) / (360 * ce.cos(interpolated_moon_dec) * ce.cos(lat) * ce.sin(lunar_local_hour_angle))

		m2 += deltaM

	# Offset days
	local_hours = m2 * 24 - utc_diff

	# Combine
	moonset_dt = datetime(ymd.year, ymd.month, ymd.day) + timedelta(hours=(local_hours))

	return moonset_dt

# Visibility calculations either:
# Type 0: Odeh, 2006
# Type 1: HMNAO TN No. 69, a.k.a. Yallop, 1997
def calculate_visibility(sun_az: float, sun_alt: float, moon_az: float, moon_alt: float, moon_pi: float, type: int = 0) -> float:
	'''
	Calculate the visibility of the Moon's crescent for a given date and observer coordinates.

	Parameters:
		sun_az (float): The Sun's azimuth (°).
		sun_alt (float): The Sun's altitude (°).
		moon_az (float): The Moon's azimuth (°).
		moon_alt (float): The Moon's altitude (°).
		moon_pi (float): The Moon's parallax (°).
		type (int): The type of visibility calculation to use (0 or 1). Default is 0 which uses Odeh, 2006. When set to 1, it uses HMNAO TN No. 69 (a.k.a. Yallop, 1997).
	
	Returns:
		float: The visibility of the Moon's crescent.
	'''

	arcl = ce.calculate_angle_diff(sun_az, sun_alt, moon_az, moon_alt)
	arcv = np.abs(sun_alt - moon_alt)
	#daz = sun_az - moon_az
	moon_pi *= 206265 / 60

	semi_diameter = 0.27245 * moon_pi
	semi_diameter_prime = semi_diameter * (1 + ce.sin(moon_alt) * ce.sin(moon_pi / 60))

	w_prime = semi_diameter_prime * (1 - ce.cos(arcl))

	if type == 0:
		q_value = (arcv - (-0.1018 * w_prime ** 3 + 0.7319 * w_prime ** 2 - 6.3226 * w_prime + 7.1651))
	else:
		q_value = (arcv - (11.8371 - 6.3226 * w_prime + 0.7319 * w_prime ** 2 - 0.1018 * w_prime ** 3)) / 10

	return q_value

# Classification according to Odeh, 2006 or HMNAO TN No.69
def classify_visibility(q: float, type: int = 0) -> str:
	'''
	Classify the visibility of the Moon's crescent based on the given q value.

	Parameters:
		q (float): The q value.
		type (int): The type of visibility classification to use (0 or 1). Default is 0 which uses Odeh, 2006. When set to 1, it uses HMNAO TN No. 69 (a.k.a. Yallop, 1997).

	Returns:
		str: The classification of the Moon's crescent visibility.
	'''
	
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

	if type == 0:
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