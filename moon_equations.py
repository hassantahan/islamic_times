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

moon_phase_corrections_coeff = [
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

moon_phase_corrections_arg = [
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

a_sin_term_phases_coeff = [
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

a_coeffs = [
	0.000325,	0.000165,	0.000164,
	0.000126,	0.000110,	0.000062,
	0.000060,	0.000056,	0.000047,
	0.000042,	0.000040,	0.000037,
	0.000035,	0.000023
]

# Chapter 47
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

def moonpos(julian_day, local_latitude, local_longitude, deltaPsi, ecliptic, elev):
	t = (julian_day - J2000) / JULIAN_CENTURY

	# Calculation nutations
	nut = moon_nutation(julian_day)

	# Rect Coordinates + Distance
	longitude = nut[0][4] + nut[1] / 10 ** 6 + deltaPsi
	latitude = nut[2] / 10 ** 6
	distance = 385000.56 + nut[3] / 10 ** 3

	# Place in the sky
	ascension = bound_angle_deg(np.rad2deg(math.atan2((sin(longitude) * cos(ecliptic) - tan(latitude) * sin(ecliptic)), cos(longitude))))
	declination = np.rad2deg(math.asin(sin(latitude) * cos(ecliptic) + cos(latitude) * sin(ecliptic) * sin(longitude)))
	eh_parallax = np.rad2deg(np.arcsin(EARTH_RADIUS_KM / distance))

	# Local Hour Angle Calculations (see sun_equations.py for more info)
	greenwich_hour_angle = bound_angle_deg(siderial_time(julian_day))
	delta_psi = sun_nutation(julian_day)[0]
	st_correction = decimal_to_dms(delta_psi)[2] * cos(ecliptic) / 15
	greenwich_hour_angle += (st_correction / 240)
	local_hour_angle = greenwich_hour_angle + local_longitude - ascension

	# Modify RA and Declination to their apparent equivalents
	app_ascension, app_declination = correct_ra_dec(ascension, declination, local_hour_angle, eh_parallax, local_latitude, elev / 1000)

	# Final calculations
	altitude = np.rad2deg(math.asin(sin(local_latitude) * sin(declination) + cos(local_latitude) * cos(declination)* cos(local_hour_angle)))
	azimuth = np.rad2deg(np.arccos((sin(declination) * cos(local_latitude) - cos(declination) * sin(local_latitude) * cos(local_hour_angle)) / cos(altitude)))

	# Possibly useless
	if local_hour_angle >= 0:
		azimuth = 360 - azimuth

	return [
		longitude, 			# 0: lambda
		latitude, 			# 1: beta
		distance,			# 2:
		ecliptic,			# 3:
		ascension,			# 4:
		declination,		# 5:
		app_ascension,		# 6
		app_declination,	# 7
		eh_parallax,		# 8:
		altitude,			# 9:
		azimuth				# 10:
	]

# Fixing RA and Dec for apparency pg. 279
def correct_ra_dec(ra, dec, lha, parallax, lat, elev, dist = EARTH_RADIUS_KM):
	a = dist
	f = 1 / 298.257
	b = a * (1 - f)

	u = np.rad2deg(np.arctan2(b * tan(lat), a))
	p_sin_psi_prime = b / a * sin(u) + elev / dist * sin(lat)
	p_cos_psi_prime = cos(u) + elev / dist * cos(lat)


	temp_num = -1 * p_cos_psi_prime * sin(parallax) * sin(lha)
	temp_denom = cos(dec) - p_cos_psi_prime * sin(parallax) * sin(lha)
	deltaA = np.rad2deg(np.arctan2(temp_num, temp_denom))

	temp_num = (sin(dec) - p_sin_psi_prime * sin(parallax)) * cos(deltaA)

	ascension_prime = ra + deltaA
	declination_prime = np.rad2deg(np.arctan2(temp_num, temp_denom))

	return ascension_prime, declination_prime

# TODO: Properly get the next phases instead of all the phases after the next new moon ==> Using round (temp fix)
# Chapter 49
def next_phases_of_moon_utc(date):
	# Find the day of the year
	day_of_year = date.timetuple().tm_yday

	# Preserve the sign
	p_sign = np.sign(date.year - 2000)

	# Calculate k based on eq. 49.2
	# k = p_sign * np.ceil((abs(date.year - 2000) + p_sign * day_of_year / TROPICAL_YEAR) * 12.3685)
	k_temp = ((abs(date.year - 2000) + p_sign * day_of_year / TROPICAL_YEAR) * 12.3685)

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
			fundamental_arguments[i] = bound_angle_deg(fundamental_arguments[i])
		
		# Construct the A terms
		a_sin_args = []
		for a_sin_coeff_arrays in a_sin_term_phases_coeff:
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
			temp = bound_angle_deg(temp)
			a_sin_args.append(temp)

		# Calculate W for First & Last Quarter phases pg. 352
		w = 0.00306 - 0.00038 * fundamental_arguments[0] * cos(fundamental_arguments[1]) \
			+ 0.00026 * cos(fundamental_arguments[2]) - 0.00002 * cos(fundamental_arguments[2] - fundamental_arguments[1]) \
			+ 0.00002 * cos(fundamental_arguments[2] + fundamental_arguments[1]) + 0.00002 * cos(2 * fundamental_arguments[3])

		# Calculate the correction for the first group of periodic terms
		temp = 0
		i = 0
		j = 0
		# Iterate over the first group by 3 and the other group by 4
		# TODO: change the moon_phase_corrections_arg[0] to account for the different phases (p) ==> done
		while i < np.size(moon_phase_corrections_coeff) and j < np.size(moon_phase_corrections_arg):
			# print('Term', i // 3 + 1,':')
			sin_coeff = 0
			
			# Shift for the moon_phase_corrections_coeff array
			s = 0
			if p == 1 or p == 3:
				s = 2
			elif p == 2:
				s = 1

			# Note: moon_phase_corrections_arg has its first [] set to p % 2 since 0 & 2
			# are the new and full moons which use the first array set in moon_phase_corrections_arg
			# Therefore, 1 & 3 are for the First and Last Quarters

			# Check if Omega term
			if moon_phase_corrections_arg[p % 2][j] != 9:
				# Go through the sine argument terms and add them up
				sin_argument = 0
				for l in np.arange(0, 3):
					# print(f'{fundamental_arguments[l + 1]:.4f}', moon_phase_corrections_arg[p % 2][j + l + 1])
					sin_argument += fundamental_arguments[l + 1] * moon_phase_corrections_arg[p % 2][j + l + 1]
				# Take the sine of the sum of the arguments
				sin_coeff = sin(sin_argument)
				# print(moon_phase_corrections_coeff[i + s], f'{fundamental_arguments[0]:.7f}', moon_phase_corrections_arg[p % 2][j])
				# Add to the correction term the product of the phase coefficient, the eccentricity factor, and the sine coefficient
				temp += moon_phase_corrections_coeff[i + s] * pow(fundamental_arguments[0], moon_phase_corrections_arg[p % 2][j]) * sin_coeff
			else:
				# Omega term
				# print(moon_phase_corrections_coeff[i + s], f'{fundamental_arguments[0]:.7f}', moon_phase_corrections_arg[p % 2][j], f'{fundamental_arguments[4]:.4f}')
				temp += moon_phase_corrections_coeff[i + s] * pow(fundamental_arguments[0], 0) * sin(fundamental_arguments[4])
			# print(f'Temp: {temp:.5f}')
			i += 3
			j += 4

		moon_phases[p] += temp

		# Calculate the 14 additional corrections from the a terms
		temp = 0
		for i, a_sin_arg in enumerate(a_sin_args):
			temp += a_coeffs[i] * sin(a_sin_arg)

		moon_phases[p] += temp

		# TODO: implement the w factor correction for first (+w) and last (-w) quarters
		if p == 1:
			moon_phases[p] += w
		elif p == 3:
			moon_phases[p] -= w

		# Convert from TD to UT in this approximation.
		moon_phases[p] -= delta_t_approx(date.year) / 86400
		
		# Convert from JD to Gregorian
		moon_phases[p] = jd_to_gregorian(moon_phases[p])

	return moon_phases

# Chapter 48
def moon_illumination(sun_dec, sun_ra, sun_long, moon_dec, moon_ra, moon_lat, moon_long, sun_earth_distance, moon_earth_distance):
	
	# Eqs. 48.2
	cos_psi = sin(sun_dec) * sin(moon_dec) + cos(sun_dec) * cos(moon_dec) * cos(sun_ra - moon_ra)
	psi = np.rad2deg(np.arccos(cos_psi))
	sin_psi = sin(psi)
	
	# Eq. 48.3
	phase_angle = np.rad2deg(np.arctan2(sun_earth_distance * sin_psi, moon_earth_distance - sun_earth_distance * cos_psi))
	# print(phase_angle)
	
	# Eq. 48.1
	fraction_illuminated = (1 + cos(phase_angle)) / 2
	return fraction_illuminated

# Visibility calculations from HMNAO TN No. 69
def calculate_visibility(sun_az, sun_alt, moon_az, moon_alt, moon_pi):
	
	# print(sun_az, sun_alt, moon_az, moon_alt, moon_pi)

	arcl = calculate_angle_diff(sun_az, sun_alt, moon_az, moon_alt)
	arcv = np.abs(sun_alt - moon_alt)
	daz = sun_az - moon_az
	moon_pi *= 206265 / 60

	semi_diameter = 0.27245 * moon_pi
	semi_diameter_prime = semi_diameter * (1 + sin(moon_alt) * sin(moon_pi / 60))

	w_prime = semi_diameter_prime * (1 - cos(arcl))

	q_value = (arcv - (11.8371 - 6.3226 * w_prime + 0.7319 * w_prime ** 2 - 0.1018 * w_prime ** 3)) / 10

	# print(arcl, arcv, daz, cos(arcl) - cos(arcv) * cos(daz))
	# print(moon_pi, w_prime)

	return q_value

# Classification according to HMNAO TN No.69
def classify_visibility(q):
    if q > 0.216:
        return "Easily visible"
    elif 0.216 >= q > -0.014:
        return "Visible under perfect conditions"
    elif -0.014 >= q > -0.160:
        return "May need optical aid"
    elif -0.160 >= q > -0.232:
        return "Will need optical aid"
    elif -0.232 >= q:
        return "Not visible"
    else:
        return "Invalid Input"