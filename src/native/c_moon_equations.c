#define PY_SSIZE_T_CLEAN
#include "c_moon_equations.h"
#include "c_event_solver.h"

/* ================================
   Definitions & Helper Constants
   ================================ */

#define MOON_NUTATION_ARGUMENTS_LR_SIZE 60
#define NUM_CORRECTION_TERMS 25
#define NUM_A_SIN_TERMS 14
#define NUM_A_COEFFS 14
#define MAX_EVENT_SEARCH_DAYS 370

/* ================================
   Nutation, Obliquity, and Moon Phases array terms
   ================================ */

static const int MOON_NUTATION_ARGUMENTS_LR[] = {
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
    
static const int MOON_NUTATION_ARGUMENTS_B[] = {
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
    
static const int MOON_NUTATION_COEFF_LR[] = {
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

static const int MOON_NUTATION_COEFF_B[] = {
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

/* MOON_PHASE_CORRECTIONS_COEFF: 25 rows x 3 columns */
static const double MOON_PHASE_CORRECTIONS_COEFF[NUM_CORRECTION_TERMS][3] = {
    {-0.40720, -0.40614, -0.62801},
    { 0.17241,  0.17302,  0.17172},
    { 0.01608,  0.01614, -0.01183},
    { 0.01039,  0.01043,  0.00862},
    { 0.00739,  0.00734,  0.00804},
    {-0.00514, -0.00515,  0.00454},
    { 0.00208,  0.00209,  0.00204},
    {-0.00111, -0.00111, -0.00180},
    {-0.00057, -0.00057, -0.00070},
    { 0.00056,  0.00056, -0.00040},
    {-0.00042, -0.00042, -0.00034},
    { 0.00042,  0.00042,  0.00032},
    { 0.00038,  0.00038,  0.00032},
    {-0.00024, -0.00024, -0.00028},
    {-0.00017, -0.00017,  0.00027},
    {-0.00007, -0.00007, -0.00017},
    { 0.00004,  0.00004, -0.00005},
    { 0.00004,  0.00004,  0.00004},
    { 0.00003,  0.00003, -0.00004},
    { 0.00003,  0.00003,  0.00004},
    {-0.00003, -0.00003,  0.00003},
    { 0.00003,  0.00003,  0.00003},
    {-0.00002, -0.00002,  0.00002},
    {-0.00002, -0.00002,  0.00002},
    { 0.00002,  0.00002, -0.00002}
};

/* MOON_PHASE_CORRECTIONS_ARG: 2 sets (for different phases) of 25 rows x 4 columns */
static const int MOON_PHASE_CORRECTIONS_ARG[2][NUM_CORRECTION_TERMS][4] = {
    {   // New & Full Moon
        {0,  0,  1,  0},
        {1,  1,  0,  0},
        {0,  0,  2,  0},
        {0,  0,  0,  2},
        {1, -1,  1,  0},
        {1,  1,  1,  0},
        {2,  2,  0,  0},
        {0,  0,  1, -2},
        {0,  0,  1,  2},
        {1,  1,  2,  0},
        {0,  0,  3,  0},
        {1,  1,  0,  2},
        {1,  1,  0, -2},
        {1, -1,  2,  0},
        {9,  9,  9,  9},  // OMEGA
        {0,  2,  1,  0},
        {0,  0,  2, -2},
        {0,  3,  0,  0},
        {0,  1,  1, -2},
        {0,  0,  2,  2},
        {0,  1,  1,  2},
        {0, -1,  1,  2},
        {0, -1,  1, -2},
        {0,  1,  3,  0},
        {0,  0,  4,  0}
    },
    {   // First & Last Quarter Moon
        {0,  0,  1,  0},
        {1,  1,  0,  0},
        {1,  1,  1,  0},
        {0,  0,  2,  0},
        {0,  0,  0,  2},
        {1, -1,  1,  0},
        {2,  2,  0,  0},
        {0,  0,  1, -2},
        {0,  0,  1,  2},
        {0,  0,  3,  0},
        {1, -1,  2,  0},
        {1,  1,  0,  2},
        {1,  1,  0, -2},
        {2,  2,  1,  0},
        {1,  1,  2,  0},
        {9,  9,  9,  9},  // OMEGA
        {0, -1,  1, -2},
        {0,  0,  2,  2},
        {0,  1,  1,  2},
        {0, -2,  1,  0},
        {0,  1,  1, -2},
        {0,  3,  0,  0},
        {0,  0,  2, -2},
        {0, -1,  1,  2},
        {0,  1,  3,  0}
    }
};

/* A_SIN_TERM_PHASES_COEFF: variable-length rows (using up to 3 coefficients here) */
static const double A_SIN_TERM_PHASES_COEFF[NUM_A_SIN_TERMS][14] = {
    {299.77,  0.1074080, -0.009173},
    {251.88,  0.0163210},
    {251.83, 26.651886},
    {349.42, 36.412478},
    {84.660,  18.206239},
    {141.74, 53.303771},
    {207.14,  2.4537320},
    {154.84,  7.3068600},
    {34.520,  27.261239},
    {207.19,  0.1218240},
    {291.34,  1.8443790},
    {161.72, 24.198154},
    {239.56, 25.513099},
    {331.55,  3.5925180}
};

/* Lengths of each row in A_SIN_TERM_PHASES_COEFF */
static const int A_SIN_TERM_LENGTH[NUM_A_SIN_TERMS] = {
    3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2
};

/* A_COEFFS array */
static const double A_COEFFS[NUM_A_COEFFS] = {
    0.000325, 0.000165, 0.000164,
    0.000126, 0.000110, 0.000062,
    0.000060, 0.000056, 0.000047,
    0.000042, 0.000040, 0.000037,
    0.000035, 0.000023
};

/* ================================
   Core calculation: lunar nutation
   ================================ */

/*
 * Compute lunar nutation series terms for one JDE instant.
 * Input:
 *   - jde: Julian Ephemeris Day
 * Output:
 *   - MoonNutationResult with fundamental arguments and summed terms.
 */
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
        int d = MOON_NUTATION_ARGUMENTS_LR[i * 4];
        int m = MOON_NUTATION_ARGUMENTS_LR[i * 4 + 1];
        int mp = MOON_NUTATION_ARGUMENTS_LR[i * 4 + 2];
        int f = MOON_NUTATION_ARGUMENTS_LR[i * 4 + 3];

        double temp_lr = d * fundamental_arguments[0] + m * fundamental_arguments[1] +
                         mp * fundamental_arguments[2] + f * fundamental_arguments[3];
        temp_lr = normalize_angle(temp_lr);
        double temp_lr_rad = RADIANS(temp_lr);

        double eccentricity_comp = (m == 0) ? 1.0 : pow(eccentricity, abs(m));

        sum_l += eccentricity_comp * MOON_NUTATION_COEFF_LR[i * 2] * sin(temp_lr_rad);
        sum_r += eccentricity_comp * MOON_NUTATION_COEFF_LR[i * 2 + 1] * cos(temp_lr_rad);
    }

    const size_t n_args = sizeof(MOON_NUTATION_ARGUMENTS_B) / (4 * sizeof(int));
    for (size_t i = 0; i < n_args; i++) {
        int d = MOON_NUTATION_ARGUMENTS_B[i * 4];
        int m = MOON_NUTATION_ARGUMENTS_B[i * 4 + 1];
        int mp = MOON_NUTATION_ARGUMENTS_B[i * 4 + 2];
        int f = MOON_NUTATION_ARGUMENTS_B[i * 4 + 3];

        double temp_b = d * fundamental_arguments[0] + m * fundamental_arguments[1] +
                        mp * fundamental_arguments[2] + f * fundamental_arguments[3];
        temp_b = normalize_angle(temp_b);

        double eccentricity_comp = (m == 0) ? 1.0 : pow(eccentricity, abs(m));

        sum_b += eccentricity_comp * MOON_NUTATION_COEFF_B[i] * sin(RADIANS(temp_b));
    }

    sum_l += 3958 * sin(RADIANS(a[0])) +
             1962 * sin(RADIANS(fundamental_arguments[4] - fundamental_arguments[3])) +
              318 * sin(RADIANS(a[1]));

    sum_b += -2235 * sin(RADIANS(fundamental_arguments[4])) +
               382 * sin(RADIANS(a[2])) +
               175 * sin(RADIANS(a[0] - fundamental_arguments[3])) +
               175 * sin(RADIANS(a[0] + fundamental_arguments[3])) +
               127 * sin(RADIANS(fundamental_arguments[4] - fundamental_arguments[2])) -
               115 * sin(RADIANS(fundamental_arguments[4] + fundamental_arguments[2]));

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
   Core calculation: full lunar result
   ================================ */

/*
 * Populate MoonResult for one observer/time.
 * Units:
 *   - jde: Julian Ephemeris Day
 *   - deltaT: seconds
 *   - latitude/longitude: degrees
 *   - elevation: metres
 *   - temperature: Celsius
 *   - pressure: kPa
 *   - deltaPsi/ecliptic: degrees
 */
void compute_moon_result(double jde, double deltaT, double local_latitude, double local_longitude,
                        double elevation, double temperature, double pressure, 
                        double deltaPsi, double ecliptic, 
                        MoonResult* result) {
    // Time variables
    double t = (jde - J2000) / JULIAN_CENTURY;
    double t2 = t * t;
    double t3 = t2 * t;

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
    compute_gha_lha(result->true_obliquity, deltaPsi, gmst, local_longitude, ra_deg,
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

/* Python wrapper helpers for constructing dataclass-compatible objects. */
static PyObject* create_angle_obj_moon(AstroCoreState* state, double value) {
    PyObject* py_value = PyFloat_FromDouble(value);
    PyObject* angle_obj;
    if (!py_value) {
        return NULL;
    }

    angle_obj = PyObject_CallFunctionObjArgs(state->angle_type, py_value, NULL);
    Py_DECREF(py_value);
    return angle_obj;
}

static PyObject* create_ra_obj_from_degrees_moon(AstroCoreState* state, double value_deg) {
    PyObject* py_hours = PyFloat_FromDouble(value_deg / 15.0);
    PyObject* ra_obj;
    if (!py_hours) {
        return NULL;
    }

    ra_obj = PyObject_CallFunctionObjArgs(state->right_ascension_type, py_hours, NULL);
    Py_DECREF(py_hours);
    return ra_obj;
}

static PyObject* create_distance_obj_moon(AstroCoreState* state, double value, const char* unit_name) {
    PyObject* py_value = NULL;
    PyObject* unit_obj = NULL;
    PyObject* dist_obj = NULL;

    py_value = PyFloat_FromDouble(value);
    if (!py_value) {
        goto error;
    }

    unit_obj = PyObject_GetAttrString(state->distance_units_type, unit_name);
    if (!unit_obj) {
        goto error;
    }

    dist_obj = PyObject_CallFunctionObjArgs(state->distance_type, py_value, unit_obj, NULL);

error:
    Py_XDECREF(py_value);
    Py_XDECREF(unit_obj);
    return dist_obj;
}

PyObject* py_compute_moon(PyObject* self, PyObject* const* args, Py_ssize_t nargs) {
    #define SET_TUPLE_ITEM_OR_FAIL(tuple_obj, index, value_expr)                 \
        do {                                                                      \
            PyObject* _tmp = (value_expr);                                        \
            if (!_tmp) {                                                          \
                goto error;                                                       \
            }                                                                     \
            PyTuple_SET_ITEM((tuple_obj), (index), _tmp);                        \
        } while (0)

    if (nargs != 9) {
        char error_message[100];
        snprintf(error_message, sizeof(error_message), "Expected 9 arguments but received %zd", nargs);
        PyErr_SetString(PyExc_TypeError, error_message);
        return NULL;
    }
    
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

    AstroCoreState* state = (AstroCoreState*)PyModule_GetState(self);
    if (!state) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to access astro_core module state.");
        return NULL;
    }
    if (!state->moon_type || !PyType_Check(state->moon_type)) {
        PyErr_SetString(PyExc_RuntimeError, "MoonType is not a valid Python type");
        return NULL;
    }

    PyObject* arg_list = PyList_New(5);
    PyObject* sum_l_obj = NULL;
    PyObject* sum_b_obj = NULL;
    PyObject* sum_r_obj = NULL;
    PyObject* lunar_nutation = NULL;
    PyObject* args_tuple = NULL;
    PyObject* moon = NULL;
    if (!arg_list) {
        return NULL;
    }

    for (int i = 0; i < 5; ++i) {
        PyObject* fundamental_argument = PyFloat_FromDouble(result.lunar_nutation.fundamental_arguments[i]);
        if (!fundamental_argument) {
            Py_DECREF(arg_list);
            return NULL;
        }
        PyList_SET_ITEM(arg_list, i, fundamental_argument);
    }

    sum_l_obj = PyFloat_FromDouble(result.lunar_nutation.sum_l);
    sum_b_obj = PyFloat_FromDouble(result.lunar_nutation.sum_b);
    sum_r_obj = PyFloat_FromDouble(result.lunar_nutation.sum_r);
    if (!sum_l_obj || !sum_b_obj || !sum_r_obj) {
        goto error;
    }

    lunar_nutation = PyTuple_New(4);
    if (!lunar_nutation) {
        goto error;
    }

    PyTuple_SET_ITEM(lunar_nutation, 0, arg_list);
    arg_list = NULL;
    PyTuple_SET_ITEM(lunar_nutation, 1, sum_l_obj);
    sum_l_obj = NULL;
    PyTuple_SET_ITEM(lunar_nutation, 2, sum_b_obj);
    sum_b_obj = NULL;
    PyTuple_SET_ITEM(lunar_nutation, 3, sum_r_obj);
    sum_r_obj = NULL;

    args_tuple = PyTuple_New(19);
    if (!args_tuple) {
        goto error;
    }
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 0, create_angle_obj_moon(state, result.true_longitude));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 1, create_angle_obj_moon(state, result.true_latitude));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 2, create_distance_obj_moon(state, result.geocentric_distance, "KILOMETRE"));
    PyTuple_SET_ITEM(args_tuple, 3, lunar_nutation);
    lunar_nutation = NULL;
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 4, create_angle_obj_moon(state, result.omega));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 5, create_angle_obj_moon(state, result.apparent_longitude));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 6, create_angle_obj_moon(state, result.deltaPsi));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 7, create_angle_obj_moon(state, result.true_obliquity));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 8, create_ra_obj_from_degrees_moon(state, result.right_ascension));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 9, create_angle_obj_moon(state, result.declination));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 10, create_angle_obj_moon(state, result.greenwich_hour_angle));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 11, create_angle_obj_moon(state, result.local_hour_angle));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 12, create_angle_obj_moon(state, result.eh_parallax));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 13, create_ra_obj_from_degrees_moon(state, result.topocentric_ascension));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 14, create_angle_obj_moon(state, result.top_declination));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 15, create_angle_obj_moon(state, result.topocentric_local_hour_angle));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 16, create_angle_obj_moon(state, result.true_altitude));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 17, create_angle_obj_moon(state, result.true_azimuth));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 18, create_angle_obj_moon(state, result.apparent_altitude));

    moon = PyObject_CallObject(state->moon_type, args_tuple);
    Py_DECREF(args_tuple);
    #undef SET_TUPLE_ITEM_OR_FAIL

    if (!moon) {
        return NULL;
    }

    return moon;

error:
    Py_XDECREF(arg_list);
    Py_XDECREF(sum_l_obj);
    Py_XDECREF(sum_b_obj);
    Py_XDECREF(sum_r_obj);
    Py_XDECREF(lunar_nutation);
    Py_XDECREF(args_tuple);
    #undef SET_TUPLE_ITEM_OR_FAIL
    return NULL;
}


/* ================================
   Moon transit / culmination solver
   ================================ */

/*
 * Solve moon transit time for the reference civil date.
 * Returns 0 on success; non-zero values indicate failure.
 */
int find_moon_transit(datetime date, double utc_offset, double local_latitude, double local_longitude,
                    double elevation, double temperature, double pressure, 
                    double deltaPsi[3], double true_obliquity[3],
                    datetime* moon_event) {

    double new_jd = gregorian_to_jd(date, 0) - fraction_of_day_datetime(date);
    double new_deltaT = delta_t_approx(date.year, date.month);

    SunResult temp_sun_param;
    MoonResult moon_params[3];
    for (int i = 0; i < 3; i++) {
        datetime temp_date;
        jd_to_gregorian(new_jd + i - 1, utc_offset, &temp_date);
        double t_deltaT = delta_t_approx(temp_date.year, temp_date.month);
        double t_jde = (new_jd + i - 1) + t_deltaT / SECONDS_IN_DAY;
        if (deltaPsi[i] == CALCULATE_SUN_PARAMS_FOR_MOON_TIME) {   
            compute_sun_result(t_jde, t_deltaT, 
                local_latitude, local_longitude, elevation, temperature, pressure, 
                &temp_sun_param);
            deltaPsi[i] = temp_sun_param.nutation_longitude;
            true_obliquity[i] = temp_sun_param.true_obliquity;
        }
        compute_moon_result(t_jde, t_deltaT, 
            local_latitude, local_longitude, elevation, temperature, pressure,
            deltaPsi[i], true_obliquity[i], &moon_params[i]);
    }
    
    double sidereal_time = greenwich_mean_sidereal_time(new_jd);
    double right_ascension_deg[3] = {
        moon_params[0].right_ascension,
        moon_params[1].right_ascension,
        moon_params[2].right_ascension
    };
    double m0 = (moon_params[1].right_ascension - local_longitude - sidereal_time) / 360.0;
    m0 = refine_transit_fraction(m0, sidereal_time, local_longitude, new_deltaT, right_ascension_deg);

    // Construct transit datetime
    moon_event->year = date.year;
    moon_event->month = date.month;
    moon_event->day = date.day;
    moon_event->hour = 0; moon_event->minute = 0; moon_event->second = 0; moon_event->microsecond = 0; 
    *moon_event = add_days(*moon_event, m0);

    return 0;
}

/* Python wrapper for find_moon_transit.
 * The deltaT argument is accepted for API compatibility but recomputed natively.
 */
PyObject* py_find_moon_transit(PyObject* self, PyObject* args) {
    PyObject *deltaPsi_obj, *true_obliquity_obj;
    double jd, unused_deltaT, latitude, longitude, 
           elevation, temperature, pressure, 
           utc_offset;
    double deltaPsi[3], true_obliquity[3];
    if (!PyArg_ParseTuple(args, "ddddddddOO", &jd, &unused_deltaT, &latitude, &longitude, 
                                              &elevation, &temperature, &pressure, 
                                              &utc_offset, &deltaPsi_obj, &true_obliquity_obj))
        return NULL;
    (void)unused_deltaT;
    
    if (!PySequence_Check(deltaPsi_obj) || PySequence_Size(deltaPsi_obj) != 3 || !PySequence_Check(true_obliquity_obj) || PySequence_Size(true_obliquity_obj) != 3) {
        PyErr_SetString(PyExc_TypeError, "Expected two sequences of 3 numeric values.");
        return NULL;
    }

    for (int i = 0; i < 3; ++i) {
        PyObject* item1 = PySequence_GetItem(deltaPsi_obj, i);
        PyObject* item2 = PySequence_GetItem(true_obliquity_obj, i);
        if (!item1 || !item2) {
            Py_XDECREF(item1);
            Py_XDECREF(item2);
            return NULL;
        }

        deltaPsi[i] = PyFloat_AsDouble(item1);
        true_obliquity[i] = PyFloat_AsDouble(item2);
        if (PyErr_Occurred()) {
            Py_DECREF(item1);
            Py_DECREF(item2);
            return NULL;
        }

        Py_DECREF(item1);
        Py_DECREF(item2);
    }

    datetime reference_dt;
    jd_to_gregorian(jd, utc_offset, &reference_dt);
    datetime transit_dt;
    int status = find_moon_transit(reference_dt, utc_offset, latitude, longitude,
                                    elevation, temperature, pressure, 
                                    deltaPsi, true_obliquity,
                                    &transit_dt);

    if (status == 0)
        return datetime_to_pydatetime(add_days(transit_dt, utc_offset / 24.0));
    else
        return NULL;
}


/* ================================
   Moonrise and moonset solver
   ================================ */

/*
 * Solve moonrise/moonset time for a reference civil date.
 * Return codes:
 *   0  -> success
 *  -1  -> event does not occur on this date/location
 *  -2  -> invalid event type
 */
int moonrise_or_moonset(datetime date, double utc_offset, double local_latitude, double local_longitude,
                        double elevation, double temperature, double pressure, char event_type,
                        double deltaPsi[3], double true_obliquity[3], 
                        datetime* moon_event) {
    
    double lat_rad = RADIANS(local_latitude);
    double new_jd = gregorian_to_jd(date, 0) - fraction_of_day_datetime(date);
    double new_deltaT = delta_t_approx(date.year, date.month);

    SunResult temp_sun_param;
    MoonResult moon_params[3];
    for (int i = 0; i < 3; i++) {
        datetime temp_date;
        jd_to_gregorian(new_jd + i - 1, utc_offset, &temp_date);
        double t_deltaT = delta_t_approx(temp_date.year, temp_date.month);
        double t_jde = (new_jd + i - 1) + t_deltaT / SECONDS_IN_DAY;
        if (deltaPsi[i] == CALCULATE_SUN_PARAMS_FOR_MOON_TIME) {   
            compute_sun_result(t_jde, t_deltaT, 
                local_latitude, local_longitude, elevation, temperature, pressure, 
                &temp_sun_param);
            deltaPsi[i] = temp_sun_param.nutation_longitude;
            true_obliquity[i] = temp_sun_param.true_obliquity;
        }
        compute_moon_result(t_jde, t_deltaT, 
            local_latitude, local_longitude, elevation, temperature, pressure,
            deltaPsi[i], true_obliquity[i], &moon_params[i]);
    }

    double declination_deg[3] = {
        moon_params[0].declination,
        moon_params[1].declination,
        moon_params[2].declination
    };
    double right_ascension_deg[3] = {
        moon_params[0].right_ascension,
        moon_params[1].right_ascension,
        moon_params[2].right_ascension
    };

    double horizon_altitude_deg = 0.7275 * moon_params[1].eh_parallax - 0.566667;
    double H_zero_deg;
    if (compute_event_hour_angle(lat_rad, declination_deg[1], horizon_altitude_deg, &H_zero_deg) != 0) {
        return -1;
    }
    
    double sidereal_time = greenwich_mean_sidereal_time(new_jd);
    double m0 = (moon_params[1].right_ascension - local_longitude - sidereal_time) / 360.0;

    double m_event;
    if (event_type == 'r')
        m_event = m0 - H_zero_deg / 360;
    else if (event_type == 's')
        m_event = m0 + H_zero_deg / 360;
    else
        return -2;
    
    m_event = refine_altitude_event_fraction(
        m_event,
        sidereal_time,
        local_longitude,
        new_deltaT,
        lat_rad,
        horizon_altitude_deg,
        declination_deg,
        right_ascension_deg
    );

    moon_event->year = date.year;
    moon_event->month = date.month;
    moon_event->day = date.day;
    moon_event->hour = 0; moon_event->minute = 0; moon_event->second = 0; moon_event->microsecond = 0; 
    *moon_event = add_days(*moon_event, m_event);

    return 0;
}

typedef struct {
    double latitude;
    double longitude;
    double elevation;
    double temperature;
    double pressure;
    char event;
    double* deltaPsi;
    double* true_obliquity;
} MoonEventSearchContext;

static int solve_moon_event_for_day(datetime date, double utc_offset, void* ctx, datetime* out_event) {
    MoonEventSearchContext* context = (MoonEventSearchContext*)ctx;
    return moonrise_or_moonset(
        date,
        utc_offset,
        context->latitude,
        context->longitude,
        context->elevation,
        context->temperature,
        context->pressure,
        context->event,
        context->deltaPsi,
        context->true_obliquity,
        out_event
    );
}

/*
 * Find the first moonrise/moonset event that belongs to the reference civil day.
 * Returns INVALID_DATETIME when no valid event is found.
 */
datetime find_proper_moontime(double jd, double utc_offset, double latitude, double longitude, double elevation,
    double temperature, double pressure, char event, double deltaPsi[3], double true_obliquity[3]) {
    datetime reference_dt;
    jd_to_gregorian(jd, utc_offset, &reference_dt);

    MoonEventSearchContext context = {
        .latitude = latitude,
        .longitude = longitude,
        .elevation = elevation,
        .temperature = temperature,
        .pressure = pressure,
        .event = event,
        .deltaPsi = deltaPsi,
        .true_obliquity = true_obliquity
    };

    return find_event_on_reference_day(
        reference_dt,
        utc_offset,
        longitude,
        MAX_EVENT_SEARCH_DAYS,
        solve_moon_event_for_day,
        &context
    );
}

/* Python wrapper for find_proper_moontime. */
PyObject* py_find_proper_moontime(PyObject* self, PyObject* args) {
    PyObject *deltaPsi_obj, *true_obliquity_obj;
    double jd, unused_deltaT, latitude, longitude, 
           elevation, temperature, pressure, 
           utc_offset;
    double deltaPsi[3], true_obliquity[3];
    int event_code;
    char event;
    if (!PyArg_ParseTuple(args, "ddddddddOOC", &jd, &unused_deltaT, &latitude, &longitude, 
                                               &elevation, &temperature, &pressure, 
                                               &utc_offset, &deltaPsi_obj, &true_obliquity_obj, 
                                               &event_code))
        return NULL;
    (void)unused_deltaT;

    if (!PySequence_Check(deltaPsi_obj) || PySequence_Size(deltaPsi_obj) != 3 ||
        !PySequence_Check(true_obliquity_obj) || PySequence_Size(true_obliquity_obj) != 3) {

        PyErr_SetString(PyExc_TypeError, "Expected two sequences of 3 numeric values.");
        return NULL;
    }

    event = (char)event_code;

    for (int i = 0; i < 3; ++i) {
        PyObject* item1 = PySequence_GetItem(deltaPsi_obj, i);
        PyObject* item2 = PySequence_GetItem(true_obliquity_obj, i);
        if (!item1 || !item2) {
            Py_XDECREF(item1);
            Py_XDECREF(item2);
            return NULL;
        }

        deltaPsi[i] = PyFloat_AsDouble(item1);
        true_obliquity[i] = PyFloat_AsDouble(item2);
        if (PyErr_Occurred()) {
            Py_DECREF(item1);
            Py_DECREF(item2);
            return NULL;
        }

        Py_DECREF(item1);
        Py_DECREF(item2);
    }
    
    return datetime_to_pydatetime(
        find_proper_moontime(jd, utc_offset, latitude, longitude, elevation, temperature, pressure, event, deltaPsi, true_obliquity));
}


/* ================================
   Calculate next moon phases
   ================================ */

void next_phases_of_moon_utc(datetime date, datetime phases[4]) {
    int day = day_of_year(date.year, date.month, date.day);
    int p_sign = (date.year - 2000 > 0) ? 1 : ((date.year - 2000 < 0) ? -1 : 1);

    double k_temp = ((fabs((double)(date.year - 2000)) +
                        p_sign * ((double)day) / TROPICAL_YEAR) * 12.3685);
    double k_array[4];
    k_array[0] = p_sign * round(k_temp);
    k_array[1] = p_sign * (round(k_temp - 0.25) + 0.25);
    k_array[2] = p_sign * (round(k_temp - 0.50) + 0.50);
    k_array[3] = p_sign * (round(k_temp - 0.75) + 0.75);

    // Loop over each phase (New Moon, First Quarter, Full Moon, Last Quarter)
    for (int p = 0; p < 4; p++) {
        double k = k_array[p];
        double t  = k / 1236.85;
        double t2 = t * t;
        double t3 = t2 * t;
        double t4 = t3 * t;

        /* Mean Julian Ephemeris Date */
        double phase_jde = 2451550.09766 + 29.530588861 * k +
                            0.00015437 * t2 - 0.000000150 * t3 +
                            0.00000000073 * t4;

        /* Fundamental arguments:
            F[0]=E, F[1]=M, F[2]=M', F[3]=F, F[4]=Omega */
        double F[5];
        F[0] = 1 - 0.002516 * t - 0.0000074 * t2;
        F[1] = 2.5534 + 29.10535670 * k - 0.0000014 * t2 - 0.00000011 * t3;
        F[2] = 201.5643 + 385.81693528 * k + 0.0107582 * t2 +
                0.00001238 * t3 - 0.000000058 * t4;
        F[3] = 160.7108 + 390.67050284 * k + 0.0016118 * t2 +
                0.00000227 * t3 - 0.000000011 * t4;
        F[4] = 124.7746 - 1.5637558 * k + 0.0020672 * t2 + 0.00000215 * t3;
        for (int i = 0; i < 5; i++) {
            F[i] = normalize_angle(F[i]);
        }

        /* A-Term Corrections */
        double a_sin_args[NUM_A_SIN_TERMS];
        for (int i = 0; i < NUM_A_SIN_TERMS; i++) {
            int L = A_SIN_TERM_LENGTH[i];
            double sum = 0.0;
            for (int j = 0; j < L; j++) {
                double multiplier = (j == 0) ? 1.0 : (j == 1 ? k : t2);
                sum += A_SIN_TERM_PHASES_COEFF[i][j] * multiplier;
            }
            a_sin_args[i] = normalize_angle(sum);
        }

        /* W Correction (for quarter phases) */
        double w = (0.00306 
                    - 0.00038 * F[0] * sin(RADIANS(F[1])))
                    + 0.00026 * sin(RADIANS(F[2]))
                    - 0.00002 * sin(RADIANS(F[2] - F[1]))
                    + 0.00002 * sin(RADIANS(F[2] + F[1]))
                    + 0.00002 * sin(RADIANS(2 * F[3]));

        /* Periodic Term Corrections */
        int s = (p == 1 || p == 3) ? 2 : (p == 2 ? 1 : 0);
        double periodic_corr = 0.0;
        for (int i = 0; i < NUM_CORRECTION_TERMS; i++) {
            double term = 0.0;
            if (MOON_PHASE_CORRECTIONS_ARG[p % 2][i][0] != 9) {
                double dot = MOON_PHASE_CORRECTIONS_ARG[p % 2][i][1] * F[1] +
                                MOON_PHASE_CORRECTIONS_ARG[p % 2][i][2] * F[2] +
                                MOON_PHASE_CORRECTIONS_ARG[p % 2][i][3] * F[3];
                double sin_val = sin(RADIANS(dot));
                term = MOON_PHASE_CORRECTIONS_COEFF[i][s] *
                        pow(F[0], MOON_PHASE_CORRECTIONS_ARG[p % 2][i][0]) *
                        sin_val;
            } else {
                double sin_omega = sin(RADIANS(F[4]));
                term = MOON_PHASE_CORRECTIONS_COEFF[i][s] * sin_omega;
            }
            periodic_corr += term;
        }
        phase_jde += periodic_corr;

        /* Additional A-Term Corrections */
        double a_term_corr = 0.0;
        for (int i = 0; i < NUM_A_COEFFS; i++) {
            a_term_corr += A_COEFFS[i] * sin(RADIANS(a_sin_args[i]));
        }
        phase_jde += a_term_corr;

        /* Apply W factor for quarter phases */
        if (p == 1)
            phase_jde += w;
        else if (p == 3)
            phase_jde -= w;

        /* Adjust from Terrestrial Dynamical Time to UT */
        phase_jde -= delta_t_approx(date.year, date.month) / 86400.0;

        jd_to_gregorian(phase_jde, 0, &phases[p]);
    }
}

PyObject* py_next_phases_of_moon_utc(PyObject* self, PyObject* args) {
    ENSURE_PYDATETIME();
    PyObject *input_datetime;

    if (!PyArg_ParseTuple(args, "O!", PyDateTimeAPI->DateTimeType, &input_datetime)) {
        return NULL;
    }

    datetime date;
    fill_in_datetime_values(&date, input_datetime);

    // Call method
    datetime moon_phases[4];
    next_phases_of_moon_utc(date, moon_phases);

    // Convert to Python datetime
    PyObject* py_moon_phases[4];
    for (int i = 0; i < 4; i++) {
        py_moon_phases[i] = datetime_to_pydatetime(moon_phases[i]);
        if (!py_moon_phases[i]) {
            for (int j = 0; j < i; j++) {
                Py_DECREF(py_moon_phases[j]);
            }
            return NULL;
        }
    }

    // Create Python tuple
    PyObject* result = PyTuple_New(4);
    if (!result) return NULL;

    for (int i = 0; i < 4; i++) {
        PyTuple_SET_ITEM(result, i, py_moon_phases[i]);  // Steals reference
    }

    return result;
}
