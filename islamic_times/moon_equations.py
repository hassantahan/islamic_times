'''
Module for Moon astronomical calculations.

This module provides functions to:
  - Calculate the Moon's position (declination, right ascension, altitude, azimuth, etc.).
  - Compute the Moon's illumination percentage.
  - Estimate the moonset time for a given observer and date.

This module should not be used directly unless the user is familiar with the underlying calculations.
'''

from islamic_times import sun_equations as se
from islamic_times import time_equations as te
from islamic_times import calculation_equations as ce
from datetime import datetime, timedelta
from typing import List
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
def moon_nutation(jde: float) -> List[float]:
	'''
	Calculate the Moon's nutation in for a given Julian Ephemeris Day. See Chapter 22 of the Astronomical Almanac for more information.
	'''

	t = (jde - te.J2000) / te.JULIAN_CENTURY
	t2 = t ** 2
	t3 = t ** 3
	t4 = t ** 4

	fundamental_arguments = [
		297.8501921 + 445267.1114034 * t - 0.0018819 * t2 + t3 / 545868 - t4 / 113065000, 	# D
		357.5291092 + 35999.0502909 * t - 0.0001536 * t2 + t3 / 24490000, 					# M
		134.9633964 + 477198.8675055 * t + 0.0087414 * t2 + t3 / 69699 - t4 / 14712000, 	# M'
		93.2720950 + 483202.0175233 * t - 0.0036539 * t2 - t3 / 3526000 + t4 / 863310000, 	# F
		218.3164477 + 481267.88123421 * t - 0.0015786 * t2 + t3 / 538841 - t4 / 65194000 	# L'
	] 
	
	a = [
		119.75 + 131.849 * t,
		53.09 + 479264.290 * t,
		313.45 + 481266.484 * t
	]

	for i in range(len(fundamental_arguments)):
		fundamental_arguments[i] %= 360
	
	for i in range(len(a)):
		a[i] %= 360

	sum_l = 0
	sum_r = 0
	sum_b = 0

	eccentricity = 1 - 0.002516 * t - 0.0000074 * t2

	for i in np.arange(0, np.size(__moon_nutation_arguments_lr), 4):
		temp =  __moon_nutation_arguments_lr[i] 		* 	fundamental_arguments[0] + \
				__moon_nutation_arguments_lr[i + 1] 	* 	fundamental_arguments[1] + \
				__moon_nutation_arguments_lr[i + 2] 	* 	fundamental_arguments[2] + \
				__moon_nutation_arguments_lr[i + 3] 	* 	fundamental_arguments[3]

		if (__moon_nutation_arguments_lr[i + 1] == 0):
			eccentricity_compensation = 1
		else:
			pow = np.abs(__moon_nutation_arguments_lr[i + 1])
			eccentricity_compensation = eccentricity ** pow
		
		sum_l += eccentricity_compensation * __moon_nutation_coeff_lr[int(i / 2)] * ce.sin(temp)
		sum_r += eccentricity_compensation * __moon_nutation_coeff_lr[int(i / 2) + 1] * ce.cos(temp)
	
	for i in np.arange(0, np.size(__moon_nutation_arguments_b), 4):
		temp =  __moon_nutation_arguments_b[i] * fundamental_arguments[0] + \
				__moon_nutation_arguments_b[i + 1] * fundamental_arguments[1] + \
				__moon_nutation_arguments_b[i + 2] * fundamental_arguments[2] + \
				__moon_nutation_arguments_b[i + 3] * fundamental_arguments[3]

		if (__moon_nutation_arguments_b[i + 1] == 0):
			eccentricity_compensation = 1
		else:
			pow = np.abs(__moon_nutation_arguments_b[i + 1])
			eccentricity_compensation = eccentricity ** pow

		sum_b += eccentricity_compensation * __moon_nutation_coeff_b[int(i / 4)] * ce.sin(temp)


	sum_l += 3958 * ce.sin(a[0]) + 1962 * ce.sin(fundamental_arguments[4] - fundamental_arguments[3]) + 318 * ce.sin(a[1])

	sum_b += -2235 * ce.sin(fundamental_arguments[4]) + 382 * ce.sin(a[2]) + 175 * ce.sin(a[0] - fundamental_arguments[3]) \
			+ 175 * ce.sin(a[0] + fundamental_arguments[3]) + 127 * ce.sin(fundamental_arguments[4] - fundamental_arguments[2]) \
			- 115 * ce.sin(fundamental_arguments[4] + fundamental_arguments[2])

	return [fundamental_arguments, sum_l, sum_b, sum_r]

def moonpos(jde: float, deltaT: float, local_latitude: float, local_longitude: float, deltaPsi: float, ecliptic: float, elev: float) -> List[float]:
	'''
    Calculate the various solar positional parameters for a given Julian Ephemeris Day, ΔT, and observer coordinates. See Chapter 25 of the Astronomical Algorthims for more information.

    Note: The temperature and pressure are used for atmospheric refraction calculations. Currently, this feature is disabled.
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
def correct_ra_dec(ra: float, dec: float, lha: float, parallax: float, lat: float, elev: float, dist: float = te.EARTH_RADIUS_KM) -> float | float:
	#
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

	return ascension_prime, declination_prime

# TODO: Properly get the next phases instead of all the phases after the next new moon ==> Using round (temp fix)
# Chapter 49
def next_phases_of_moon_utc(date: datetime) -> List[datetime]:
	# Find the day of the year
	day_of_year = date.timetuple().tm_yday

	# Preserve the sign
	p_sign = np.sign(date.year - 2000)

	# Calculate k based on eq. 49.2
	# k = p_sign * np.ceil((abs(date.year - 2000) + p_sign * day_of_year / TROPICAL_YEAR) * 12.3685)
	k_temp = ((abs(date.year - 2000) + p_sign * day_of_year / te.TROPICAL_YEAR) * 12.3685)

	k_array = [
		p_sign * (np.round(k_temp)), 
		p_sign * (np.round(k_temp - 0.25) + 0.25), 
		p_sign * (np.round(k_temp - 0.5) + 0.5), 
		p_sign * (np.round(k_temp - 0.75) + 0.75)
	]

	# Construct phase array
	moon_phases = [0, 0, 0, 0]

	for p, k in enumerate(k_array): 
		# Julian centuries based on eq. 49.3
		t = k / 1236.85
		t2 = t ** 2
		t3 = t ** 3
		t4 = t ** 4

		# JDE calculation
		jde = 2451550.09766 + 29.530588861 * k + 0.00015437 * t2 - 0.000000150 * t3 + 0.00000000073 * t4

		moon_phases[p] += jde

		# Construct the sine arguments
		fundamental_arguments = [
			1 - 0.002516 * t - 0.0000074 * t2,													# E
			2.5534 + 29.10535670 * k - 0.0000014 * t2 - 0.00000011 * t3,						# M
			201.5643 + 385.81693528 * k + 0.0107582 * t2 + 0.00001238 * t3 - 0.000000058 * t4,	# M'
			160.7108 + 390.67050284 * k + 0.0016118 * t2 + 0.00000227 * t3 - 0.000000011 * t4,	# F
			124.7746 - 1.5637558 * k + 0.0020672 * t2 + 0.00000215 * t3							# Omega
		]

		# Bound the fundamental arguments
		for i in range(5):
			fundamental_arguments[i] %= 360
		
		# Construct the A terms
		a_sin_args = []
		for a_sin_coeff_arrays in __a_sin_term_phases_coeff:
			temp = 0
			j = 0
			for a_sin_coeff in a_sin_coeff_arrays:
				if j == 0:
					temp += a_sin_coeff
				elif j == 1:
					temp += a_sin_coeff * k
				else:
					temp += a_sin_coeff * t2
				j += 1
			temp %= 360
			a_sin_args.append(temp)

		# Calculate W for First & Last Quarter phases pg. 352
		w = 0.00306 - 0.00038 * fundamental_arguments[0] * ce.cos(fundamental_arguments[1]) \
			+ 0.00026 * ce.cos(fundamental_arguments[2]) - 0.00002 * ce.cos(fundamental_arguments[2] - fundamental_arguments[1]) \
			+ 0.00002 * ce.cos(fundamental_arguments[2] + fundamental_arguments[1]) + 0.00002 * ce.cos(2 * fundamental_arguments[3])

		# Calculate the correction for the first group of periodic terms
		temp = 0
		i = 0
		j = 0
		# Iterate over the first group by 3 and the other group by 4
		# TODO: change the __moon_phase_corrections_arg[0] to account for the different phases (p) ==> done
		while i < np.size(__moon_phase_corrections_coeff) and j < np.size(__moon_phase_corrections_arg):
			sin_coeff = 0
			
			# Shift for the __moon_phase_corrections_coeff array
			s = 0
			if p == 1 or p == 3:
				s = 2
			elif p == 2:
				s = 1

			# Note: __moon_phase_corrections_arg has its first [] set to p % 2 since 0 & 2
			# are the new and full moons which use the first array set in __moon_phase_corrections_arg
			# Therefore, 1 & 3 are for the First and Last Quarters

			# Check if Omega term
			if __moon_phase_corrections_arg[p % 2][j] != 9:
				# Go through the sine argument terms and add them up
				sin_argument = 0
				for l in np.arange(0, 3):
					sin_argument += fundamental_arguments[l + 1] * __moon_phase_corrections_arg[p % 2][j + l + 1]
				# Take the sine of the sum of the arguments
				sin_coeff = ce.sin(sin_argument)
				# Add to the correction term the product of the phase coefficient, the eccentricity factor, and the sine coefficient
				temp += __moon_phase_corrections_coeff[i + s] * pow(fundamental_arguments[0], __moon_phase_corrections_arg[p % 2][j]) * sin_coeff
			else:
				# Omega term
				temp += __moon_phase_corrections_coeff[i + s] * pow(fundamental_arguments[0], 0) * ce.sin(fundamental_arguments[4])

			i += 3
			j += 4

		moon_phases[p] += temp

		# Calculate the 14 additional corrections from the a terms
		temp = 0
		for i, a_sin_arg in enumerate(a_sin_args):
			temp += __a_coeffs[i] * ce.sin(a_sin_arg)

		moon_phases[p] += temp

		# TODO: implement the w factor correction for first (+w) and last (-w) quarters
		if p == 1:
			moon_phases[p] += w
		elif p == 3:
			moon_phases[p] -= w

		# Convert from TD to UT in this approximation.
		moon_phases[p] -= te.delta_t_approx(date.year, date.month) / 86400
		
		# Convert from JD to Gregorian
		moon_phases[p] = te.jd_to_gregorian(moon_phases[p])

	return moon_phases

# Refer to Chapter 48 of AA
def moon_illumination(sun_dec: float, sun_ra: float, moon_dec: float, moon_ra: float, sun_earth_distance: float, moon_earth_distance: float) -> float:
	
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
			return "Invalid Input"
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
			return "Invalid Input"
		