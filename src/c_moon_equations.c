#define PY_SSIZE_T_CLEAN
#include <math.h>
#include <string.h>
#include "c_moon_equations.h"

/* ================================
   Definitions & Helper Constants
   ================================ */

#define MOON_NUTATION_ARGUMENTS_LR_SIZE 60
#define ASTRONOMICAL_UNIT_KM 149597870.7
#define EARTH_RADIUS_KM 6378.14


/* ================================
   Nutation, Obliquity, and Moon Phases array terms
   ================================ */

static const int __MOON_NUTATION_ARGUMENTS_LR[] = {
//	D		 M		 M'		 F
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
};
    
static const int __MOON_NUTATION_ARGUMENTS_B[] = {
//  D		 M		 M'		 F
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
};
    
static const int __MOON_NUTATION_COEFF_LR[] = {
//  	   l		 	 r
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
};

static const int __MOON_NUTATION_COEFF_B[] = {
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
};
    
static const double __MOON_PHASE_CORRECTIONS_COEFF[] = {
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
};
    
static const int __MOON_PHASE_CORRECTIONS_ARG[2][100] = {
        //New & Full Moon
        {
        //	E	 M	 M'	 F
        //	Note: The number for E is it's power, not coeff
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
            9,	 9,	  9,	 9,		//OMEGA
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
        },
    
        //First & Last Quarter Moon
        {
        //	E	 M	 M'	 F
        //	Note: The number for E is it's power, not coeff
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
            9,	 9,	  9,	 9,		//OMEGA
            0,	-1,	  1,	-2,
            0,	 0,	  2,	 2,
            0,	 1,	  1,	 2,
            0,	-2,	  1,	 0,
            0,	 1,	  1,	-2,
            0,	 3,	  0,	 0,
            0,	 0,	  2,	-2,
            0,	-1,	  1,	 2,
            0,	 1,	  3,	 0
        }
};
    
static const double __A_SIN_TERM_PHASES_COEFF[][14] = {
        {299.77,	0.1074080,	-0.009173},
        {251.88,	0.0163210},
        {251.83,	26.651886},
        {349.42,	36.412478},
        {84.660,	18.206239},
        {141.74,	53.303771},
        {207.14,	2.4537320},
        {154.84,	7.3068600},
        {34.520,	27.261239},
        {207.19,	0.1218240},
        {291.34,	1.8443790},
        {161.72,	24.198154},
        {239.56,	25.513099},
        {331.55,	3.5925180}
};
    
static const double __A_COEFFS[] = {
        0.000325,	0.000165,	0.000164,
        0.000126,	0.000110,	0.000062,
        0.000060,	0.000056,	0.000047,
        0.000042,	0.000040,	0.000037,
        0.000035,	0.000023
};

/* ================================
   core method: moon nuation
   ================================ */

MoonNutationResult moon_nutation(double jde) {
    double t = (jde - J2000) / JULIAN_CENTURY;
    double t2 = t * t;
    double t3 = t2 * t;
    double t4 = t3 * t;

    double fundamental_arguments[5] = {
        normalize_angle(297.8501921 + 445267.1114034  * t - 0.0018819 * t2 + t3 / 545868  - t4 / 113065000),  // D
        normalize_angle(357.5291092 + 35999.0502909   * t - 0.0001536 * t2 + t3 / 24490000),                  // M
        normalize_angle(134.9633964 + 477198.8675055  * t + 0.0087414 * t2 + t3 / 69699   - t4 / 14712000),   // M'
        normalize_angle(93.2720950  + 483202.0175233  * t - 0.0036539 * t2 - t3 / 3526000 + t4 / 863310000),  // F
        normalize_angle(218.3164477 + 481267.88123421 * t - 0.0015786 * t2 + t3 / 538841  - t4 / 65194000)    // L'
    };

    double a[3] = {
        normalize_angle(119.75 + 131.849    * t),
        normalize_angle(53.09  + 479264.290 * t),
        normalize_angle(313.45 + 481266.484 * t)
    };

    double sum_l = 0.0;
    double sum_r = 0.0;
    double sum_b = 0.0;

    double eccentricity = 1 - 0.002516 * t - 0.0000074 * t2;

    for (int i = 0; i < MOON_NUTATION_ARGUMENTS_LR_SIZE; i++) {
        int d = __MOON_NUTATION_ARGUMENTS_LR[i * 4];
        int m = __MOON_NUTATION_ARGUMENTS_LR[i * 4 + 1];
        int mp = __MOON_NUTATION_ARGUMENTS_LR[i * 4 + 2];
        int f = __MOON_NUTATION_ARGUMENTS_LR[i * 4 + 3];

        double temp_lr = d * fundamental_arguments[0] + m * fundamental_arguments[1] +
                         mp * fundamental_arguments[2] + f * fundamental_arguments[3];
        temp_lr = normalize_angle(temp_lr);

        double eccentricity_comp = (m == 0) ? 1.0 : pow(eccentricity, abs(m));

        sum_l += eccentricity_comp * __MOON_NUTATION_COEFF_LR[i * 2] * sin(temp_lr * M_PI / 180.0);
        sum_r += eccentricity_comp * __MOON_NUTATION_COEFF_LR[i * 2 + 1] * cos(temp_lr * M_PI / 180.0);
    }

    for (int i = 0; i < sizeof(__MOON_NUTATION_ARGUMENTS_B) / (4 * sizeof(int)); i++) {
        int d = __MOON_NUTATION_ARGUMENTS_B[i * 4];
        int m = __MOON_NUTATION_ARGUMENTS_B[i * 4 + 1];
        int mp = __MOON_NUTATION_ARGUMENTS_B[i * 4 + 2];
        int f = __MOON_NUTATION_ARGUMENTS_B[i * 4 + 3];

        double temp_b = d * fundamental_arguments[0] + m * fundamental_arguments[1] +
                        mp * fundamental_arguments[2] + f * fundamental_arguments[3];
        temp_b = normalize_angle(temp_b);

        double eccentricity_comp = (m == 0) ? 1.0 : pow(eccentricity, abs(m));

        sum_b += eccentricity_comp * __MOON_NUTATION_COEFF_B[i] * sin(temp_b * M_PI / 180.0);
    }

    sum_l += 3958 * sin(a[0] * M_PI / 180.0) +
             1962 * sin((fundamental_arguments[4] - fundamental_arguments[3]) * M_PI / 180.0) +
              318 * sin(a[1] * M_PI / 180.0);

    sum_b += -2235 * sin(fundamental_arguments[4] * M_PI / 180.0) +
               382 * sin(a[2] * M_PI / 180.0) +
               175 * sin((a[0] - fundamental_arguments[3]) * M_PI / 180.0) +
               175 * sin((a[0] + fundamental_arguments[3]) * M_PI / 180.0) +
               127 * sin((fundamental_arguments[4] - fundamental_arguments[2]) * M_PI / 180.0) -
               115 * sin((fundamental_arguments[4] + fundamental_arguments[2]) * M_PI / 180.0);

    MoonNutationResult result;
    for (int i = 0; i < 5; i++) {
        result.fundamental_arguments[i] = fundamental_arguments[i];
    }
    result.sum_l = sum_l;
    result.sum_b = sum_b;
    result.sum_r = sum_r;

    return result;
}


/* ================================
   core method: moon result
   ================================ */

void compute_moon_result(double jde, double deltaT, double local_latitude, double local_longitude,
                        double elevation, double temperature, double pressure, 
                        double deltaPsi, double ecliptic, 
                        MoonResult* result) {
    // Time variables
    double t = (jde - J2000) / JULIAN_CENTURY;
    double t2 = t * t;
    double t3 = t2 * t;
    double t4 = t3 * t;

    // Compute nutation
    MoonNutationResult nutation_result = moon_nutation(jde);

    // Orbital elements
    result->true_longitude = normalize_angle(nutation_result.fundamental_arguments[4] + nutation_result.sum_l / 1000000.0);
    result->true_latitude = nutation_result.sum_b / 1000000.0;
    result->geocentric_distance = 385000.56 + nutation_result.sum_r / 1000.0;

    // Nutation and obliquity
    result->lunar_nutation = nutation_result;
    result->omega = normalize_angle(125.04452 - 1934.136261 * t + 0.0020708 * t2 + t3 / 450000.0);
    result->apparent_longitude = result->true_longitude + deltaPsi;
    result->deltaPsi = deltaPsi;
    result->true_obliquity = ecliptic;

    // Apparent coordinates
    double ra_deg, dec_deg;
    compute_equitorial_coordinates(result->apparent_longitude, result->true_obliquity, 
                                    result->true_latitude, &ra_deg, &dec_deg);
    result->right_ascension = ra_deg;
    result->declination = dec_deg;

    // Hour angles
    double jd = jde - deltaT / 86400.0;
    double gmst = greenwich_mean_sidereal_time(jd);
    double gst, gha_deg, lha_deg;
    compute_gha_lha(jd, result->true_obliquity, deltaPsi, gmst, local_longitude, ra_deg,
                    &gst, &gha_deg, &lha_deg);

    result->greenwich_hour_angle = gha_deg;
    result->local_hour_angle = lha_deg;

    // Parallax
    result->eh_parallax = DEGREES(asin(sin(RADIANS(8.794 / 3600.0)) / (result->geocentric_distance / ASTRONOMICAL_UNIT_KM)));

    // Topocentric corrections
    double top_ra = result->right_ascension;
    double top_dec = result->declination;
    correct_ra_dec(&top_ra, &top_dec, result->local_hour_angle, result->eh_parallax, local_latitude, elevation / 1000.0, EARTH_RADIUS_KM);

    result->topocentric_ascension = top_ra;
    result->top_declination = top_dec;
    result->topocentric_local_hour_angle = normalize_angle(gst + local_longitude - top_ra);

    // Horizontal coordinates
    double true_alt, true_az;
    compute_horizontal_coordinates(top_ra, top_dec, result->topocentric_local_hour_angle,
                                    local_latitude, &true_az, &true_alt);
    result->true_altitude = true_alt;
    result->true_azimuth = true_az;

    // Apparent altitude (refraction correction)
    double refraction_correction = 0.0167 / tan((result->true_altitude + 10.3 / (result->true_altitude + 5.11)) * M_PI / 180.0);
    result->apparent_altitude = result->true_altitude + refraction_correction;
}

/* Python wrapper */
PyObject* py_compute_moon(PyObject* self, PyObject* const* args, Py_ssize_t nargs) {
    if (nargs != 9) {
        char error_message[100];
        snprintf(error_message, sizeof(error_message), "Expected 9 arguments but received %zu", nargs);
        PyErr_SetString(PyExc_TypeError, error_message);
        return NULL;
    }
    
    // double jde, deltaT, latitude, longitude, elevation, temperature, pressure;
    // if (!PyArg_ParseTuple(args, "dd", &jde, &deltaT, &latitude, &longitude, 
    //                     &elevation, &temperature, &pressure)) return NULL;

    double jde = PyFloat_AsDouble(args[0]);
    double deltaT = PyFloat_AsDouble(args[1]);
    double latitude = PyFloat_AsDouble(args[2]);
    double longitude = PyFloat_AsDouble(args[3]);
    double elevation = PyFloat_AsDouble(args[4]);
    double temperature = PyFloat_AsDouble(args[5]);
    double pressure = PyFloat_AsDouble(args[6]);
    double deltaPsi = PyFloat_AsDouble(args[7]);
    double ecliptic = PyFloat_AsDouble(args[8]);

    if (PyErr_Occurred()) return NULL;
    
    MoonResult result;
    compute_moon_result(jde, deltaT, latitude, longitude, 
                        elevation, temperature, pressure, 
                        deltaPsi, ecliptic,
                        &result);

    if (!MoonType || !PyType_Check(MoonType)) {
        PyErr_SetString(PyExc_RuntimeError, "MoonType is not a valid Python type");
        return NULL;
    }

    PyObject *arg_list = PyList_New(5);
    for (int i = 0; i < 5; ++i) {
        PyList_SET_ITEM(arg_list, i, PyFloat_FromDouble(result.lunar_nutation.fundamental_arguments[i]));
    }

    PyObject *lunar_nutation = PyTuple_Pack(4,
        arg_list,
        PyFloat_FromDouble(result.lunar_nutation.sum_l),
        PyFloat_FromDouble(result.lunar_nutation.sum_b),
        PyFloat_FromDouble(result.lunar_nutation.sum_r)
    );

    PyObject* moon = PyObject_CallFunctionObjArgs((PyObject *)MoonType,
        ANGLE(result.true_longitude),
        ANGLE(result.true_latitude),
        DIST(result.geocentric_distance, UNIT("KILOMETRE")),
        lunar_nutation,
        ANGLE(result.omega),
        ANGLE(result.apparent_longitude),
        ANGLE(result.deltaPsi),
        ANGLE(result.true_obliquity),
        RA(result.right_ascension),
        ANGLE(result.declination),
        ANGLE(result.greenwich_hour_angle),
        ANGLE(result.local_hour_angle),
        ANGLE(result.eh_parallax),
        RA(result.topocentric_ascension),
        ANGLE(result.top_declination),
        ANGLE(result.topocentric_local_hour_angle),
        ANGLE(result.true_altitude),
        ANGLE(result.true_azimuth),
        ANGLE(result.apparent_altitude),
        NULL
    );

    if (!moon) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create Moon object");
        return NULL;
    }

    return moon;
}