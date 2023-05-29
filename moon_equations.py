import numpy as np
import math
from time_equations import *
from sun_equations import *

moon_nutation_arguments_lr = [
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

moon_nutation_arguments_b = [
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

moon_nutation_coeff_lr = [
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

moon_nutation_coeff_b = [
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

def moon_nutation(julian_day):
	t = (julian_day - J2000) / JULIAN_CENTURY
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

	for i in range(5):
		fundamental_arguments[i] = bound_angle_deg(fundamental_arguments[i])
	
	for i in range(3):
		a[i] = bound_angle_deg(a[i])

	sum_l = 0
	sum_r = 0
	sum_b = 0

	eccentricity = 1 - 0.002516 * t - 0.0000074 * t2

	for i in np.arange(0, np.size(moon_nutation_arguments_lr), 4):
		temp =  moon_nutation_arguments_lr[i] * fundamental_arguments[0] + \
				moon_nutation_arguments_lr[i + 1] * fundamental_arguments[1] + \
				moon_nutation_arguments_lr[i + 2] * fundamental_arguments[2] + \
				moon_nutation_arguments_lr[i + 3] * fundamental_arguments[3]

		if (moon_nutation_arguments_lr[i + 1] == 0):
			eccentricity_compensation = 1
		else:
			pow = np.abs(moon_nutation_arguments_lr[i + 1])
			eccentricity_compensation = eccentricity ** pow
		
		sum_l += eccentricity_compensation * moon_nutation_coeff_lr[int(i / 2)] * sin(temp)
		sum_r += eccentricity_compensation * moon_nutation_coeff_lr[int(i / 2) + 1] * cos(temp)
	
	for i in np.arange(0, np.size(moon_nutation_arguments_b), 4):
		temp =  moon_nutation_arguments_b[i] * fundamental_arguments[0] + \
				moon_nutation_arguments_b[i + 1] * fundamental_arguments[1] + \
				moon_nutation_arguments_b[i + 2] * fundamental_arguments[2] + \
				moon_nutation_arguments_b[i + 3] * fundamental_arguments[3]

		if (moon_nutation_arguments_b[i + 1] == 0):
			eccentricity_compensation = 1
		else:
			pow = np.abs(moon_nutation_arguments_b[i + 1])
			eccentricity_compensation = eccentricity ** pow

		sum_b += eccentricity_compensation * moon_nutation_coeff_b[int(i / 4)] * sin(temp)


	sum_l += 3958 * sin(a[0]) + 1962 * sin(fundamental_arguments[4] - fundamental_arguments[3]) + 318 * sin(a[1])

	sum_b += -2235 * sin(fundamental_arguments[4]) + 382 * sin(a[2]) + 175 * sin(a[0] - fundamental_arguments[3]) \
			+ 175 * sin(a[0] + fundamental_arguments[3]) + 127 * sin(fundamental_arguments[4] - fundamental_arguments[2]) \
			- 115 * sin(fundamental_arguments[4] + fundamental_arguments[2])

	return [fundamental_arguments, sum_l, sum_b, sum_r]

def moonpos(julian_day, local_latitude, local_longitude):
	t = (julian_day - J2000) / JULIAN_CENTURY
	t2 = t ** 2
	t3 = t ** 3
	t4 = t ** 4

	nut = moon_nutation(julian_day)
	longitude = nut[0][4] + nut[1] / 1000000
	latitude = nut[2] / 1000000
	distance = 385000.56 + nut[3] / 1000

	sun_factors = sunpos(julian_day, local_latitude, local_longitude)
	ecliptic = sun_factors[13]

	ascension = bound_angle_deg(np.rad2deg(math.atan2((sin(longitude) * cos(ecliptic) - tan(latitude) * sin(ecliptic)), cos(longitude))))
	declination = np.rad2deg(math.asin(sin(latitude) * cos(ecliptic) + cos(latitude) * sin(ecliptic) * sin(longitude)))

	eh_parallax = np.rad2deg(np.arcsin(6378.14 / distance))

	greenwich_hour_angle = bound_angle_deg(siderial_time(julian_day))

	delta_psi = sun_nutation(julian_day)[0]
	st_correction = decimal_to_dms(delta_psi)[2] * cos(ecliptic) / 15
	greenwich_hour_angle += (st_correction / 240)
	
	local_hour_angle = greenwich_hour_angle + local_longitude - ascension #- 15 * ((fraction_of_day(today) + 5/24) % 1)

	altitude = np.rad2deg(math.asin(sin(local_latitude) * sin(declination) + cos(local_latitude) * cos(declination)* cos(local_hour_angle)))

	azimuth = np.rad2deg(np.arccos((sin(declination) * cos(local_latitude) - cos(declination) * sin(local_latitude) * cos(local_hour_angle)) / cos(altitude)))
	if local_hour_angle >= 0:
		azimuth = 360 - azimuth

	return [
		longitude, 			# 0: lambda
		latitude, 			# 1: beta
		distance,			# 2:
		ecliptic,			# 3:
		ascension,			# 4:
		declination,		# 5:
		eh_parallax,		# 6:
		altitude,			# 7:
		azimuth				# 8:
	]