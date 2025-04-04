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
from islamic_times.dataclasses import *
from datetime import datetime, timedelta
from dataclasses import dataclass, replace, field
from typing import List, Tuple
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

__MOON_PHASE_CORRECTIONS_COEFF = [
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

__MOON_PHASE_CORRECTIONS_ARG = [
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

__A_SIN_TERM_PHASES_COEFF = [
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

__A_COEFFS = [
	0.000325,	0.000165,	0.000164,
	0.000126,	0.000110,	0.000062,
	0.000060,	0.000056,	0.000047,
	0.000042,	0.000040,	0.000037,
	0.000035,	0.000023
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

	# jde: float
	# deltaT: float
	# local_latitude: Angle
	# local_longitude: Angle
	# elevation: Distance
	# deltaPsi: Angle
	# ecliptic: Angle

	# lunar_nutation: Tuple = field(init=False)
	# true_longitude: Angle = field(init=False)
	# true_latitude: Angle = field(init=False)
	# geocentric_distance: Distance = field(init=False)
	# apparent_longitude: Angle = field(init=False)
	# right_ascension: RightAscension = field(init=False)
	# declination: Angle = field(init=False)
	# eh_parallax: Angle = field(init=False)
	# greenwich_hour_angle: Angle = field(init=False)
	# local_hour_angle: Angle = field(init=False)
	# topocentric_ascension: RightAscension = field(init=False)
	# top_declination: Angle = field(init=False)
	# top_local_hour_angle: Angle = field(init=False)
	# true_altitude: Angle = field(init=False)
	# true_azimuth: Angle = field(init=False)
	# apparent_altitude: Angle = field(init=False)

	# def __post_init__(self):
	# 	self._compute_nutation()
	# 	self._compute_position()
	# 	self._compute_equatorial()
	# 	self._compute_local_hour_angle()
	# 	self._compute_topocentric()
	# 	self._compute_horizontal_coordinates()

	# def _compute_nutation(self):
	# 	object.__setattr__(self, 'lunar_nutation', moon_nutation(self.jde))

	# def _compute_position(self):
	# 	lon = self.lunar_nutation[0][4] + self.lunar_nutation[1] / 1e6
	# 	lat = self.lunar_nutation[2] / 1e6
	# 	dist = 385000.56 + self.lunar_nutation[3] / 1e3

	# 	object.__setattr__(self, 'true_longitude', Angle(lon))
	# 	object.__setattr__(self, 'true_latitude', Angle(lat))
	# 	object.__setattr__(self, 'geocentric_distance', Distance(dist, DistanceUnits.KILOMETRE))
	# 	object.__setattr__(self, 'apparent_longitude', Angle(lon + self.deltaPsi.decimal))

	# def _compute_equatorial(self):
	# 	ra = math.atan2(
	# 		math.sin(math.radians(self.apparent_longitude.decimal)) * math.cos(self.ecliptic.radians) -
	# 		math.tan(math.radians(self.true_latitude.decimal)) * math.sin(self.ecliptic.radians),
	# 		math.cos(math.radians(self.apparent_longitude.decimal))
	# 	)

	# 	dec = math.asin(
	# 		math.sin(math.radians(self.true_latitude.decimal)) * math.cos(self.ecliptic.radians) +
	# 		math.cos(math.radians(self.true_latitude.decimal)) * math.sin(self.ecliptic.radians) *
	# 		math.sin(math.radians(self.apparent_longitude.decimal))
	# 	)

	# 	parallax = math.degrees(math.asin(
	# 		math.sin(math.radians(8.794 / 3600)) / (self.geocentric_distance.in_unit(DistanceUnits.AU))
	# 	))

	# 	object.__setattr__(self, 'right_ascension', RightAscension(math.degrees(ra) % 360 / 15))
	# 	object.__setattr__(self, 'declination', Angle(math.degrees(dec)))
	# 	object.__setattr__(self, 'eh_parallax', Angle(parallax))

	# def _compute_local_hour_angle(self):
	# 	gst = te.greenwich_mean_sidereal_time(self.jde - self.deltaT / 86400).decimal
	# 	arcsec = self.deltaPsi.dms[2]  # arcseconds directly from DMS
	# 	st_corr = arcsec * math.cos(self.ecliptic.radians) / 15
	# 	app_gst = gst + (st_corr / 3600)
	# 	lha = (app_gst + self.local_longitude.decimal - self.right_ascension.decimal_degrees.decimal) % 360

	# 	object.__setattr__(self, 'greenwich_hour_angle', Angle(app_gst))
	# 	object.__setattr__(self, 'local_hour_angle', Angle(lha))

	# def _compute_topocentric(self):
	# 	top_ra, top_dec = ce.correct_ra_dec(
	# 		self.right_ascension,
	# 		self.declination,
	# 		self.local_hour_angle,
	# 		self.eh_parallax,
	# 		self.local_latitude,
	# 		self.elevation
	# 	)

	# 	top_lha = (te.greenwich_mean_sidereal_time(self.jde - self.deltaT / 86400).decimal +
	# 				self.local_longitude.decimal - top_ra.decimal_degrees.decimal) % 360

	# 	object.__setattr__(self, 'topocentric_ascension', top_ra)
	# 	object.__setattr__(self, 'top_declination', top_dec)
	# 	object.__setattr__(self, 'top_local_hour_angle', Angle(top_lha))

	# def _compute_horizontal_coordinates(self):
	# 	alt = math.asin(
	# 		math.sin(self.local_latitude.radians) * math.sin(self.top_declination.radians) +
	# 		math.cos(self.local_latitude.radians) * math.cos(self.top_declination.radians) *
	# 		math.cos(self.top_local_hour_angle.radians)
	# 	)

	# 	az = math.atan2(
	# 		-math.cos(self.top_declination.radians) * math.sin(self.top_local_hour_angle.radians),
	# 		math.sin(self.top_declination.radians) * math.cos(self.local_latitude.radians) -
	# 		math.cos(self.top_declination.radians) * math.sin(self.local_latitude.radians) *
	# 		math.cos(self.top_local_hour_angle.radians)
	# 	)

	# 	alt_deg = math.degrees(alt)
	# 	az_deg = math.degrees(az) % 360

	# 	refraction = 1.02 / math.tan(math.radians(alt_deg + 10.3 / (alt_deg + 5.11))) + 0.0019279 - 0.000034 * self.elevation.value

	# 	object.__setattr__(self, 'true_altitude', Angle(alt_deg))
	# 	object.__setattr__(self, 'true_azimuth', Angle(az_deg))
	# 	object.__setattr__(self, 'apparent_altitude', Angle(alt_deg + refraction / 60))

# Chapter 47
def moon_nutation(jde: float) -> Tuple[float, float, float]:
	"""
	Calculate the Moon's nutation in for a given Julian Ephemeris Day. See Chapter 47 of *Astronomical Algorithms* for more information.

	Parameters:
		jde (float): The Julian Ephemeris Day.

	Returns:
		List[float]: The nutation in longitude and obliquity (in degrees).
	"""

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
		Moon (obj): A `Moon` object that contains various attributes that describe its position. 

	Notes: 
		- The temperature and pressure are used for atmospheric refraction calculations. Currently, this feature is enabled.
	"""

	# the_moon = Moon(
	# 	jde=observer_date.jde,
	# 	deltaT=observer_date.deltaT,
	# 	local_latitude=observer.latitude,
	# 	local_longitude=observer.longitude,
	# 	elevation=observer.elevation,
	# 	deltaPsi=deltaPsi,
	# 	ecliptic=ecliptic
	# )

	# print(the_moon)
	import islamic_times.astro_core as fast_astro
	the_moon: Moon = fast_astro.compute_moon(observer_date.jde, observer_date.deltaT, 
						 observer.latitude.decimal, observer.longitude.decimal, 
						 observer.elevation.in_unit(DistanceUnits.METRE), 0, 0, 
						 deltaPsi.decimal, ecliptic.decimal)

	return the_moon

# TODO: Properly get the next phases instead of all the phases after the next new moon ==> Using round (temp fix)
# Chapter 49
def next_phases_of_moon_utc(date: datetime) -> List[datetime]:
	"""
	Calculate the next four phases of the Moon (New Moon, First Quarter, Full Moon, Last Quarter) after a given date. See Chapter 49 of *Astronomical Algorithms* for more information.

	Parameters:
		date (datetime): The date to calculate the next phases from.

	Returns:
		list (List[datetime]): The next four phases of the Moon.
	"""

	# Find the day of the year.
	day_of_year = date.timetuple().tm_yday

	# Preserve the sign relative to year 2000.
	p_sign = np.sign(date.year - 2000)
	if not (p_sign in [1, -1]):
		p_sign = 1

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
	moon_phases = [None] * 4

	# Pre-convert A-term coefficients to NumPy arrays.
	a_sin_term_coeffs = [np.array(item) for item in __A_SIN_TERM_PHASES_COEFF]
	a_coeffs = np.array(__A_COEFFS)

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
		# For each set in __A_SIN_TERM_PHASES_COEFF, use multipliers: [1, k, t2, ...] (if more than 2 numbers).
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
				- 0.00038 * fundamental_arguments[0] * ce.sin(fundamental_arguments[1]))	\
				+ 0.00026 * ce.sin(fundamental_arguments[2])								\
				- 0.00002 * ce.sin(fundamental_arguments[2] - fundamental_arguments[1])		\
				+ 0.00002 * ce.sin(fundamental_arguments[2] + fundamental_arguments[1])		\
				+ 0.00002 * ce.sin(2 * fundamental_arguments[3])

		# === Periodic Term Corrections ===
		# Reshape the correction coefficients and arguments for vectorized operations.
		coeff_array = np.array(__MOON_PHASE_CORRECTIONS_COEFF).reshape(-1, 3)
		arg_array = np.array(__MOON_PHASE_CORRECTIONS_ARG[p % 2]).reshape(-1, 4)

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
			sin_vals = ce.sin(sin_arguments)
			# Multiply by the coefficient and the eccentricity factor (fundamental_arguments[0] raised to arg_row[0])
			periodic_corr += np.sum(
				coeff_array[mask, s] * (fundamental_arguments[0] ** arg_array[mask, 0]) * sin_vals
			)
		# For Omega terms (first number equals 9), the sine argument is simply fundamental_arguments[4].
		if np.any(~mask):
			sin_omega = ce.sin(fundamental_arguments[4])
			periodic_corr += np.sum(
				coeff_array[~mask, s] * (fundamental_arguments[0] ** 0) * sin_omega
			)
		phase += periodic_corr

		# === Additional A-Term Corrections ===
		# Use vectorized dot product (after converting the a_sin_args list to an array).
		a_term_corr = np.dot(a_coeffs, ce.sin(a_sin_args))
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

def find_moon_transit(observer_date: DateTimeInfo, observer: ObserverInfo) -> datetime:
	"""
	Calculate transit of the moon (specifically culmination) for a given date and observer coordinates. See Chapter 15 of *Astronomical Algorithms* for more information.

	Parameters:
		observer_date (DateTimeInfo): The date and time of the observer.
		observer (ObserverInfo): The observer's coordinates and elevation.

	Returns:
		datetime: The time of moonset.
	"""

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
		delPsi = sun_params.delta_obliquity
		moon_params.append(
							moonpos(
								replace(observer_date, date=ymd_temp, jd=(new_jd + i - 1), deltaT=delT_temp),
								observer, delPsi, sun_params.true_obliquity
							)
						)

	# Greenwich Mean Sidereal Time (GMST)
	sidereal_time: Angle = te.greenwich_mean_sidereal_time(new_jd)

	# Transit approximation (m0)
	m0 = (moon_params[1].right_ascension.decimal_degrees.decimal - observer.longitude.decimal - sidereal_time.decimal) / 360

	# Iteratively refine transit time
	for _ in range(3):
		little_theta_zero = Angle((sidereal_time.decimal + 360.985647 * m0) % 360)
		n = m0 + new_deltaT / 86400

		interpolated_moon_ra = RightAscension(ce.interpolation(n,
												moon_params[0].right_ascension.decimal_degrees.decimal,
												moon_params[1].right_ascension.decimal_degrees.decimal,
												moon_params[2].right_ascension.decimal_degrees.decimal) / 15
											)

		lunar_local_hour_angle = Angle((little_theta_zero.decimal - (-observer.longitude.decimal) - interpolated_moon_ra.decimal_degrees.decimal) % 360)

		m0 -= lunar_local_hour_angle.decimal / 360

	# Normalize m0 to [0,1]
	m0 %= 1

	# Compute final transit datetime
	moon_transit_dt = datetime(ymd.year, ymd.month, ymd.day) + timedelta(days=m0) - timedelta(hours=observer_date.utc_offset)

	return moon_transit_dt.replace(tzinfo=observer_date.date.tzinfo)

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
		delPsi = sun_params.delta_obliquity
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
def find_proper_moontime(observer_date: DateTimeInfo, observer: ObserverInfo, rise_or_set: str = 'set') -> datetime:
	"""
    Determines the proper local time for a setting or rising moon. It finds the time that corresponds to the reference date given.

    Parameters:
        observer_date (DateTimeInfo): The date and time of the observer.
		observer (ObserverInfo): The observer's coordinates and elevation.
        rise_or_set (str): Find either the setting or rising option. Default is set to 'set'.

    Returns:
        datetime: The date and time of the moon event. If the moon event is not found (does not set or rise), returns `math.inf`.

    Raises:
        ValueError: If `rise_or_set` is not set correctly to either 'rise' or 'set'.
    """

	if rise_or_set not in ['rise', 'set', 'moonrise', 'moonset']:
		raise ValueError("Invalid value for rise_or_set. Please use 'rise' or 'set'.")

	if observer_date.utc_offset == 0:
			temp_utc_offset = math.floor(observer.longitude.decimal / 15) - 1
	else:
		temp_utc_offset = observer_date.utc_offset * -1

	temp_moonset: datetime = moonrise_or_moonset(observer_date, observer, rise_or_set)
	date_doy = observer_date.date.timetuple().tm_yday

	i = 1
	while(True):
		if temp_moonset == math.inf:
			return datetime.min
		
		temp_moonset_doy = (temp_moonset + timedelta(hours=temp_utc_offset)).timetuple().tm_yday
		if (temp_moonset_doy < date_doy and temp_moonset.year == observer_date.date.year) or ((temp_moonset + timedelta(hours=temp_utc_offset)).year < observer_date.date.year):
			temp_moonset = moonrise_or_moonset(
										replace(observer_date, date=observer_date.date + timedelta(days=i)), 
										observer, rise_or_set
									)
			i += 1
		else: 
			return temp_moonset.replace(tzinfo=observer_date.date.tzinfo)

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