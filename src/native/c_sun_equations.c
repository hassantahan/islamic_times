#define PY_SSIZE_T_CLEAN
#include "c_sun_equations.h"

/* ================================
   Definitions & Helper Constants
   ================================ */

#define OBLIQUITY_TERMS_SIZE 10
#define SUN_NUTATION_COEFFICIENTS_LOOP 63
#define SUN_NUTATION_ARGUMENTS_LOOP 5
/* Keep loop constants aligned with the tabulated term dimensions above. */


/* ================================
   Nutation and Obliquity array terms
   ================================ */

static const int SUN_NUTATION_ARGUMENTS[] = {
    0,  0,  0,  0,  1,
   -2,  0,  0,  2,  2,
    0,  0,  0,  2,  2,
    0,  0,  0,  0,  2,
    0,  1,  0,  0,  0,
    0,  0,  1,  0,  0,
   -2,  1,  0,  2,  2,
    0,  0,  0,  2,  1,
    0,  0,  1,  2,  2,
   -2, -1,  0,  2,  2,
   -2,  0,  1,  0,  0,
   -2,  0,  0,  2,  1,
    0,  0, -1,  2,  2,
    2,  0,  0,  0,  0,
    0,  0,  1,  0,  1,
    2,  0, -1,  2,  2,
    0,  0, -1,  0,  1,
    0,  0,  1,  2,  1,
   -2,  0,  2,  0,  0,
    0,  0, -2,  2,  1,
    2,  0,  0,  2,  2,
    0,  0,  2,  2,  2,
    0,  0,  2,  0,  0,
   -2,  0,  1,  2,  2,
    0,  0,  0,  2,  0,
   -2,  0,  0,  2,  0,
    0,  0, -1,  2,  1,
    0,  2,  0,  0,  0,
    2,  0, -1,  0,  1,
   -2,  2,  0,  2,  2,
    0,  1,  0,  0,  1,
   -2,  0,  1,  0,  1,
    0, -1,  0,  0,  1,
    0,  0,  2, -2,  0,
    2,  0, -1,  2,  1,
    2,  0,  1,  2,  2,
    0,  1,  0,  2,  2,
   -2,  1,  1,  0,  0,
    0, -1,  0,  2,  2,
    2,  0,  0,  2,  1,
    2,  0,  1,  0,  0,
   -2,  0,  2,  2,  2,
   -2,  0,  1,  2,  1,
    2,  0, -2,  0,  1,
    2,  0,  0,  0,  1,
    0, -1,  1,  0,  0,
   -2, -1,  0,  2,  1,
   -2,  0,  0,  0,  1,
    0,  0,  2,  2,  1,
   -2,  0,  2,  0,  1,
   -2,  1,  0,  2,  1,
    0,  0,  1, -2,  0,
   -1,  0,  1,  0,  0,
   -2,  1,  0,  0,  0,
    1,  0,  0,  0,  0,
    0,  0,  1,  2,  0,
    0,  0, -2,  2,  2,
   -1, -1,  1,  0,  0,
    0,  1,  1,  0,  0,
    0, -1,  1,  2,  2,
    2, -1, -1,  2,  2,
    0,  0,  3,  2,  2,
    2, -1,  0,  2,  2
};

static const double SUN_NUTATION_COEFFICIENTS[] = {
   -171996,  -174.2,   92025,     8.9,          //  0,  0,  0,  0,  1 
    -13187,    -1.6,    5736,    -3.1,          // -2,  0,  0,  2,  2 
     -2274,     -.2,     977,     -.5,          //  0,  0,  0,  2,  2 
      2062,      .2,    -895,      .5,          //  0,  0,  0,  0,  2 
      1426,    -3.4,      54,     -.1,          //  0,  1,  0,  0,  0 
       712,      .1,      -7,       0,          //  0,  0,  1,  0,  0 
      -517,     1.2,     224,     -.6,          // -2,  1,  0,  2,  2 
      -386,     -.4,     200,       0,          //  0,  0,  0,  2,  1 
      -301,      .0,     129,     -.1,          //  0,  0,  1,  2,  2 
       217,     -.5,     -95,      .3,          // -2, -1,  0,  2,  2 
      -158,       0,       0,       0,          // -2,  0,  1,  0,  0 
       129,      .1,     -70,       0,          // -2,  0,  0,  2,  1 
       123,       0,     -53,       0,          //  0,  0, -1,  2,  2 
        63,       0,       0,       0,          //  2,  0,  0,  0,  0 
        63,      .1,     -33,       0,          //  0,  0,  1,  0,  1 
       -59,       0,      26,       0,          //  2,  0, -1,  2,  2 
       -58,     -.1,      32,       0,          //  0,  0, -1,  0,  1 
       -51,       0,      27,       0,          //  0,  0,  1,  2,  1 
        48,       0,       0,       0,          // -2,  0,  2,  0,  0 
        46,       0,     -24,       0,          //  0,  0, -2,  2,  1 
       -38,       0,      16,       0,          //  2,  0,  0,  2,  2 
       -31,       0,      13,       0,          //  0,  0,  2,  2,  2 
        29,       0,       0,       0,          //  0,  0,  2,  0,  0 
        29,       0,     -12,       0,          // -2,  0,  1,  2,  2 
        26,       0,       0,       0,          //  0,  0,  0,  2,  0 
       -22,       0,       0,       0,          // -2,  0,  0,  2,  0 
        21,       0,     -10,       0,          //  0,  0, -1,  2,  1 
        17,     -.1,       0,       0,          //  0,  2,  0,  0,  0 
        16,       0,      -8,       0,          //  2,  0, -1,  0,  1 
       -16,      .1,       7,       0,          // -2,  2,  0,  2,  2 
       -15,       0,       9,       0,          //  0,  1,  0,  0,  1 
       -13,       0,       7,       0,          // -2,  0,  1,  0,  1 
       -12,       0,       6,       0,          //  0, -1,  0,  0,  1 
        11,       0,       0,       0,          //  0,  0,  2, -2,  0 
       -10,       0,       5,       0,          //  2,  0, -1,  2,  1 
        -8,       0,       3,       0,          //  2,  0,  1,  2,  2 
         7,       0,      -3,       0,          //  0,  1,  0,  2,  2 
        -7,       0,       0,       0,          // -2,  1,  1,  0,  0 
        -7,       0,       3,       0,          //  0, -1,  0,  2,  2 
        -7,       0,       3,       0,          //  2,  0,  0,  2,  1 
         6,       0,       0,       0,          //  2,  0,  1,  0,  0 
         6,       0,      -3,       0,          // -2,  0,  2,  2,  2 
         6,       0,      -3,       0,          // -2,  0,  1,  2,  1 
        -6,       0,       3,       0,          //  2,  0, -2,  0,  1 
        -6,       0,       3,       0,          //  2,  0,  0,  0,  1 
         5,       0,       0,       0,          //  0, -1,  1,  0,  0 
        -5,       0,       3,       0,          // -2, -1,  0,  2,  1 
        -5,       0,       3,       0,          // -2,  0,  0,  0,  1 
        -5,       0,       3,       0,          //  0,  0,  2,  2,  1 
         4,       0,       0,       0,          // -2,  0,  2,  0,  1 
         4,       0,       0,       0,          // -2,  1,  0,  2,  1 
         4,       0,       0,       0,          //  0,  0,  1, -2,  0 
        -4,       0,       0,       0,          // -1,  0,  1,  0,  0 
        -4,       0,       0,       0,          // -2,  1,  0,  0,  0 
        -4,       0,       0,       0,          //  1,  0,  0,  0,  0 
         3,       0,       0,       0,          //  0,  0,  1,  2,  0 
        -3,       0,       0,       0,          // -1, -1,  1,  0,  0 
        -3,       0,       0,       0,          //  0,  1,  1,  0,  0 
        -3,       0,       0,       0,          //  0, -1,  1,  2,  2 
        -3,       0,       0,       0,          //  2, -1, -1,  2,  2 
        -3,       0,       0,       0,          //  0,  0, -2,  2,  2 
        -3,       0,       0,       0,          //  0,  0,  3,  2,  2 
        -3,       0,       0,       0           //  2, -1,  0,  2,  2 
};

static const double OBLIQUITY_TERMS[] = {
    -4680.93,
    -1.55,
    1999.25,
    -51.38,
    -249.67,
    -39.05,
    7.12,
    27.87,
    5.79,
    2.45
};


/* ================================
   Core calculation: solar nutation
   ================================ */

/*
 * Compute nutation in longitude and obliquity for a JDE instant.
 * Inputs:
 *   - jde: Julian Ephemeris Day
 * Outputs (degrees via out-pointers):
 *   - deltaPsi
 *   - deltaEpsilon
 */
void sun_nutation(double jde, double* deltaPsi, double* deltaEpsilon) {
    double t = (jde - J2000) / JULIAN_CENTURY;
    double t2 = t * t;
    double t3 = t * t2;

    // Mean arguments of lunar-solar motion
    double ta[5];
    ta[0] = normalize_angle(297.850363 + 445267.11148 * t - 0.0019142 * t2 + t3 / 189474.0);
    ta[1] = normalize_angle(357.52772 + 35999.05034 * t - 0.0001603 * t2 - t3 / 300000.0);
    ta[2] = normalize_angle(134.96298 + 477198.867398 * t + 0.0086972 * t2 + t3 / 56250.0);
    ta[3] = normalize_angle(93.27191 + 483202.017538 * t - 0.0036825 * t2 + t3 / 327270.0);
    ta[4] = normalize_angle(125.04452 - 1934.136261 * t + 0.0020708 * t2 + t3 / 450000.0);

    for (int i = 0; i < 5; ++i) {
        ta[i] *= M_PI / 180.0; // convert to radians
    }

    // Nutation sums
    double dp = 0.0;
    double de = 0.0;

    for (int i = 0; i < SUN_NUTATION_COEFFICIENTS_LOOP; ++i) {
        double arg = 0.0;
        for (int j = 0; j < SUN_NUTATION_ARGUMENTS_LOOP; ++j) {
            int local_arg = SUN_NUTATION_ARGUMENTS[i * 5 + j];
            if (local_arg != 0) {
                arg += local_arg * ta[j];
            }
        }
        
        dp += (SUN_NUTATION_COEFFICIENTS[i * 4 + 0] + SUN_NUTATION_COEFFICIENTS[i * 4 + 1] * t) * sin(arg);
        de += (SUN_NUTATION_COEFFICIENTS[i * 4 + 2] + SUN_NUTATION_COEFFICIENTS[i * 4 + 3] * t) * cos(arg);
    }

    *deltaPsi = dp / (3600.0 * 10000.0);
    *deltaEpsilon = de / (3600.0 * 10000.0);
}


/* ================================
   Core calculation: mean obliquity
   ================================ */

/*
 * Compute mean obliquity of the ecliptic in degrees.
 * Input t is Julian millennia from J2000.
 */
double oblique_eq(double t) {
    // See Jean Meeus, Astronomical Algorithms, Chapter 22.
    double u = t / 10.0;
    double eps = 23.0 + 26.0 / 60.0 + (21.448 / 3600.0);

    // Compute powers of u
    double power = u;
    for (int i = 0; i < OBLIQUITY_TERMS_SIZE; ++i) {
        eps += OBLIQUITY_TERMS[i] * power / 3600.0;
        power *= u;
    }

    return eps;
}


/* ================================
   Core calculation: full solar result
   ================================ */

/*
 * Populate the SunResult struct for one observer/time.
 * Units:
 *   - jde: Julian Ephemeris Day
 *   - deltaT: seconds
 *   - latitude/longitude: degrees
 *   - elevation: metres
 *   - temperature: Celsius
 *   - pressure: kPa
 */
void compute_sun_result(double jde, double deltaT, double local_latitude, double local_longitude,
                        double elevation, double temperature, double pressure,
                        SunResult* result) {
    // Compute time variable
    double jd = jde - deltaT / SECONDS_IN_DAY;
    result->t = (jde - J2000) / JULIAN_MILLENNIUM;
    double t = result->t;
    double t2 = t * t;
    double t3 = t2 * t;
    double t4 = t3 * t;
    double t5 = t4 * t;
    
    // Orbital Elements
    // Mean longitude
    double mean_long = normalize_angle(280.4664567 + 360007.6982779 * t + 0.03032028 * t2 + t3 / 49931.0 - t4 / 15300.0 - t5 / 2000000.0);
    result->mean_longitude = mean_long;
    
    // Mean anomaly
    double mean_anom = normalize_angle(357.52911 + 359990.50340 * t - 0.001603 * t2 - t3 / 30000.0);
    result->mean_anomaly = mean_anom;
    double mean_anom_rad = RADIANS(mean_anom);
    
    // Eccentricity of Earth's orbit
    double ecc = 0.016708634 - 0.00042037 * t - 0.000001267 * t2;
    result->earth_orbit_eccentricity = ecc;
    
    // Sun centre (the equation of the centre)
    double sun_center_val = (1.914602 -  0.04817 * t - 0.000014 * t2) * sin(1 * mean_anom_rad) +
                            (0.019993 - 0.000101 * t)                 * sin(2 * mean_anom_rad) +
                             0.000289                                 * sin(3 * mean_anom_rad);
    result->sun_centre = sun_center_val;
    
    // True longitude and anomaly
    double true_long = mean_long + sun_center_val;
    result->true_longitude = true_long;
    double true_anom = mean_anom + sun_center_val;
    result->true_anomaly = true_anom;

    // Geocentric distance (in AU)
    double geo_dist = (1.000001018 * (1 - ecc * ecc)) / (1 + ecc * cos(RADIANS(true_anom)));
    result->geocentric_distance = geo_dist;

    // Nutation and Obliquity
    double dp, de;
    sun_nutation(jde, &dp, &de);
    result->nutation_longitude = dp;
    result->nutation_obliquity = de;
    
    // Mean obliquity
    result->mean_obliquity = oblique_eq(result->t);
    
    // True obliquity
    result->true_obliquity = result->mean_obliquity + de;
    
    // Compute omega
    double omega = 125.04452 - 19341.36261 * t + 0.020708 * t2 + t3 / 45000.0;
    result->omega = omega;
    double omega_rad = RADIANS(omega);
    
    // Apparent longitude
    double app_long = true_long - 0.00569 - 0.00478 * sin(omega_rad);
    result->apparent_longitude = app_long;
    

    // True equatorial coordinates
    double ra, dec;
    compute_equitorial_coordinates(app_long, result->true_obliquity, 0, &ra, &dec);
    result->true_right_ascension = ra;
    result->true_declination = dec;
    
    // Apparent equatorial coordinates
    double epsilon_corr = result->true_obliquity + 0.00256 * cos(omega_rad); // obliquity correction
    double app_ra, app_dec;
    compute_equitorial_coordinates(app_long, epsilon_corr, 0, &app_ra, &app_dec);

    result->apparent_right_ascension = app_ra;
    result->apparent_declination = app_dec;
    

    // Hour Angles
    double gst, gha, lha;
    double gmst = greenwich_mean_sidereal_time(jd);
    compute_gha_lha(result->true_obliquity, dp, gmst, local_longitude, app_ra, &gst, &gha, &lha);
    result->greenwich_hour_angle = gha;
    result->local_hour_angle = lha;
    

    // Parallax
    result->eh_parallax = DEGREES(asin(sin(RADIANS(8.794 / 3600)) / result->geocentric_distance));


    // Topocentric Corrections
    double top_ra = result->apparent_right_ascension; 
    double top_dec = result->apparent_declination;
    correct_ra_dec(&top_ra, &top_dec, result->local_hour_angle, result->eh_parallax, local_latitude, elevation / 1000, EARTH_RADIUS_KM);

    result->topocentric_ascension = top_ra;
    result->topocentric_declination = top_dec;
    result->topocentric_local_hour_angle = normalize_angle(gst + local_longitude - top_ra);
    

    // Horizontal Coordinates
    // True
    double true_alt, true_az;
    compute_horizontal_coordinates(top_ra, top_dec, result->topocentric_local_hour_angle, 
                                    local_latitude, &true_az, &true_alt);
    result->true_altitude = true_alt;
    result->true_azimuth = true_az;

    // Apparent (disabled)
    // double refraction_factor = 1.02 / tan((result->true_altitude + 10.3 / (result->true_altitude + 5.11)) * M_PI / 180) * pressure / 101.325 *  283 / (273 + temperature);
    result->apparent_altitude = result->true_altitude;//+ refraction_factor / 60;
}


/* Python wrapper helpers for constructing dataclass-compatible objects. */
static PyObject* create_angle_obj(double value) {
    PyObject* py_value = PyFloat_FromDouble(value);
    PyObject* angle_obj;
    if (!py_value) {
        return NULL;
    }

    angle_obj = PyObject_CallFunctionObjArgs(AngleType, py_value, NULL);
    Py_DECREF(py_value);
    return angle_obj;
}

static PyObject* create_ra_obj_from_degrees(double value_deg) {
    PyObject* py_hours = PyFloat_FromDouble(value_deg / 15.0);
    PyObject* ra_obj;
    if (!py_hours) {
        return NULL;
    }

    ra_obj = PyObject_CallFunctionObjArgs(RightAscensionType, py_hours, NULL);
    Py_DECREF(py_hours);
    return ra_obj;
}

static PyObject* create_distance_obj(double value, const char* unit_name) {
    PyObject* py_value = NULL;
    PyObject* unit_obj = NULL;
    PyObject* dist_obj = NULL;

    py_value = PyFloat_FromDouble(value);
    if (!py_value) {
        goto error;
    }

    unit_obj = PyObject_GetAttrString(DistanceUnitsType, unit_name);
    if (!unit_obj) {
        goto error;
    }

    dist_obj = PyObject_CallFunctionObjArgs(DistanceType, py_value, unit_obj, NULL);

error:
    Py_XDECREF(py_value);
    Py_XDECREF(unit_obj);
    return dist_obj;
}

PyObject* py_compute_sun(PyObject* self, PyObject* const* args, Py_ssize_t nargs) {
    #define SET_TUPLE_ITEM_OR_FAIL(tuple_obj, index, value_expr)                 \
        do {                                                                      \
            PyObject* _tmp = (value_expr);                                        \
            if (!_tmp) {                                                          \
                goto error;                                                       \
            }                                                                     \
            PyTuple_SET_ITEM((tuple_obj), (index), _tmp);                        \
        } while (0)

    if (nargs != 7) {
        PyErr_SetString(PyExc_TypeError, "Expected 7 arguments");
        return NULL;
    }

    double jde = PyFloat_AsDouble(args[0]);
    double deltaT = PyFloat_AsDouble(args[1]);
    double latitude = PyFloat_AsDouble(args[2]);
    double longitude = PyFloat_AsDouble(args[3]);
    double elevation = PyFloat_AsDouble(args[4]);
    double temperature = PyFloat_AsDouble(args[5]);
    double pressure = PyFloat_AsDouble(args[6]);

    if (PyErr_Occurred()) return NULL;

    SunResult result;
    compute_sun_result(jde, deltaT, latitude, longitude, elevation, temperature, pressure, &result);

    if (!SunType || !PyType_Check(SunType)) {
        PyErr_SetString(PyExc_RuntimeError, "SunType is not a valid Python type");
        return NULL;
    }

    PyObject* args_tuple = PyTuple_New(26);
    PyObject* nutation_tuple = NULL;
    PyObject* sun = NULL;
    if (!args_tuple) {
        return NULL;
    }

    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 0, create_angle_obj(result.mean_longitude));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 1, create_angle_obj(result.mean_anomaly));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 2, PyFloat_FromDouble(result.earth_orbit_eccentricity));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 3, create_angle_obj(result.sun_centre));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 4, create_angle_obj(result.true_longitude));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 5, create_angle_obj(result.true_anomaly));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 6, create_distance_obj(result.geocentric_distance, "AU"));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 7, create_angle_obj(result.omega));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 8, create_angle_obj(result.apparent_longitude));

    nutation_tuple = PyTuple_New(2);
    if (!nutation_tuple) {
        goto error;
    }
    SET_TUPLE_ITEM_OR_FAIL(nutation_tuple, 0, create_angle_obj(result.nutation_longitude));
    SET_TUPLE_ITEM_OR_FAIL(nutation_tuple, 1, create_angle_obj(result.nutation_obliquity));
    PyTuple_SET_ITEM(args_tuple, 9, nutation_tuple);
    nutation_tuple = NULL;

    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 10, create_angle_obj(result.nutation_obliquity));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 11, create_angle_obj(result.mean_obliquity));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 12, create_angle_obj(result.true_obliquity));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 13, create_ra_obj_from_degrees(result.true_right_ascension));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 14, create_angle_obj(result.true_declination));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 15, create_ra_obj_from_degrees(result.apparent_right_ascension));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 16, create_angle_obj(result.apparent_declination));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 17, create_angle_obj(result.greenwich_hour_angle));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 18, create_angle_obj(result.local_hour_angle));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 19, create_angle_obj(result.eh_parallax));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 20, create_ra_obj_from_degrees(result.topocentric_ascension));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 21, create_angle_obj(result.topocentric_declination));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 22, create_angle_obj(result.topocentric_local_hour_angle));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 23, create_angle_obj(result.true_altitude));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 24, create_angle_obj(result.true_azimuth));
    SET_TUPLE_ITEM_OR_FAIL(args_tuple, 25, create_angle_obj(result.apparent_altitude));

    sun = PyObject_CallObject(SunType, args_tuple);
    Py_DECREF(args_tuple);
    #undef SET_TUPLE_ITEM_OR_FAIL
    return sun;

error:
    Py_XDECREF(nutation_tuple);
    Py_DECREF(args_tuple);
    #undef SET_TUPLE_ITEM_OR_FAIL
    return NULL;
}


/* ================================
   Sun transit / culmination solver
   ================================ */

/*
 * Solve transit time for the reference civil date.
 * Returns 0 on success; non-zero values indicate failure.
 */
int find_sun_transit(datetime date, double utc_offset, double local_latitude, double local_longitude,
                    double elevation, double temperature, double pressure, datetime* sun_event) {

    double new_jd = gregorian_to_jd(date, 0) - fraction_of_day_datetime(date);
    double new_deltaT = delta_t_approx(date.year, date.month);

    SunResult sun_params[3];
    for (int i = 0; i < 3; i++) {
        datetime temp_date;
        jd_to_gregorian(new_jd + i - 1, utc_offset, &temp_date);
        double t_deltaT = delta_t_approx(temp_date.year, temp_date.month);
        double t_jde = (new_jd + i - 1) + t_deltaT / SECONDS_IN_DAY;
        compute_sun_result(t_jde, t_deltaT, 
            local_latitude, local_longitude, elevation, temperature, pressure, 
            &sun_params[i]
        );
    }
    
    double sidereal_time = greenwich_mean_sidereal_time(new_jd);
    double m0 = (sun_params[1].apparent_right_ascension - local_longitude - sidereal_time) / 360.0;
    
    for (int i = 0; i < 3; i++) {
        double theta_event_deg = normalize_angle(sidereal_time + 360.985647 * m0);
        double n_event = m0 + new_deltaT / SECONDS_IN_DAY;
        double interp_ra_event_deg = angle_interpolation(n_event,
                                    sun_params[0].apparent_right_ascension,
                                    sun_params[1].apparent_right_ascension,
                                    sun_params[2].apparent_right_ascension);

        double lha_event_deg = normalize_angle(theta_event_deg - (-local_longitude) - interp_ra_event_deg);
        m0 -= lha_event_deg / 360;
    }

    // Canonical modulo into [0, 1)
    m0 = fmod(m0, 1.0);
    if (m0 < 0.0) {
        m0 += 1.0;
    }

    // Construct transit datetime
    sun_event->year = date.year;
    sun_event->month = date.month;
    sun_event->day = date.day;
    sun_event->hour = 0; sun_event->minute = 0; sun_event->second = 0; sun_event->microsecond = 0; 
    *sun_event = add_days(*sun_event, m0);

    return 0;
}

/* Python wrapper for find_sun_transit.
 * The deltaT argument from Python is accepted for API compatibility but
 * recomputed internally per-day inside native logic.
 */
PyObject* py_find_sun_transit(PyObject* self, PyObject* args) {
    double jd, unused_deltaT, latitude, longitude, elevation, temperature, pressure, utc_offset;
    if (!PyArg_ParseTuple(args, "dddddddd", &jd, &unused_deltaT, &latitude, &longitude, &elevation, &temperature, &pressure, &utc_offset))
        return NULL;
    (void)unused_deltaT;
    
    datetime reference_dt;
    jd_to_gregorian(jd, utc_offset, &reference_dt);
    datetime transit_dt;
    int status = find_sun_transit(reference_dt, utc_offset, latitude, longitude,
                                    elevation, temperature, pressure, &transit_dt);

    if (status == 0)
        return datetime_to_pydatetime(add_days(transit_dt, utc_offset / 24.0));
    else
        return NULL;
}


/* ================================
   Sunrise and sunset solver
   ================================ */

/*
 * Solve sunrise/sunset time for a reference civil date.
 * Return codes:
 *   0  -> success
 *  -1  -> event does not occur on this date/location
 *  -2  -> invalid event type
 */
int sunrise_or_sunset(datetime date, double utc_offset, double local_latitude, double local_longitude,
                        double elevation, double temperature, double pressure, char event_type, double angle_deg, datetime* sun_event) {
    
    double lat_rad = RADIANS(local_latitude);
    double new_jd = gregorian_to_jd(date, 0) - fraction_of_day_datetime(date);
    double new_deltaT = delta_t_approx(date.year, date.month);

    SunResult sun_params[3];
    for (int i = 0; i < 3; i++) {
        datetime temp_date;
        jd_to_gregorian(new_jd + i - 1, utc_offset, &temp_date);
        double t_deltaT = delta_t_approx(temp_date.year, temp_date.month);
        double t_jde = (new_jd + i - 1) + t_deltaT / SECONDS_IN_DAY;
        compute_sun_result(t_jde, t_deltaT, 
            local_latitude, local_longitude, elevation, temperature, pressure, 
            &sun_params[i]
        );
    }

    double h_zero_rad = RADIANS(-angle_deg);
    double cosH_zero = (sin(h_zero_rad) - sin(lat_rad) * sin(RADIANS(sun_params[1].apparent_declination))) / \
                        (cos(lat_rad) * cos(RADIANS(sun_params[1].apparent_declination)));

    double H_zero_rad, H_zero_deg;
    if (cosH_zero < 1.0 && cosH_zero > -1.0) {
        H_zero_rad = acos(cosH_zero);
        H_zero_deg = H_zero_rad * 180 / M_PI;
    }
    else
        return -1;
    
    double sidereal_time = greenwich_mean_sidereal_time(new_jd);
    double m0 = (sun_params[1].apparent_right_ascension - local_longitude - sidereal_time) / 360.0;

    double m_event;
    if (event_type == 'r')
        m_event = m0 - H_zero_deg / 360;
    else if (event_type == 's')
        m_event = m0 + H_zero_deg / 360;
    else
        return -2;
    
    for (int i = 0; i < 3; i++) {
        double theta_event_deg = normalize_angle(sidereal_time + 360.985647 * m_event);
        double n_event = m_event + new_deltaT / SECONDS_IN_DAY;
        double interp_dec_event_rad = angle_interpolation(n_event,
                                    sun_params[0].apparent_declination,
                                    sun_params[1].apparent_declination,
                                    sun_params[2].apparent_declination) * M_PI / 180;
        double interp_ra_event_deg = angle_interpolation(n_event,
                                    sun_params[0].apparent_right_ascension,
                                    sun_params[1].apparent_right_ascension,
                                    sun_params[2].apparent_right_ascension);

        double lha_event_rad = normalize_angle(theta_event_deg - (-local_longitude) - interp_ra_event_deg) * M_PI / 180;
        double sun_alt_deg = asin(sin(lat_rad) * sin(interp_dec_event_rad) + 
                                cos(lat_rad) * cos(interp_dec_event_rad) * cos(lha_event_rad)) * 180 / M_PI;

        double deltaM = (sun_alt_deg - h_zero_rad * 180 / M_PI) / (360 * cos(interp_dec_event_rad) * cos(lat_rad) * sin(lha_event_rad));
        m_event += deltaM;
    }

    sun_event->year = date.year;
    sun_event->month = date.month;
    sun_event->day = date.day;
    sun_event->hour = 0; sun_event->minute = 0; sun_event->second = 0; sun_event->microsecond = 0; 
    *sun_event = add_days(*sun_event, m_event);

    return 0;
}

/*
 * Find the first sunrise/sunset event that belongs to the reference civil day.
 * Returns INVALID_DATETIME when no valid event is found.
 */
datetime find_proper_suntime(double jd, double utc_offset, double latitude, double longitude, double elevation, 
    double temperature, double pressure, double angle_deg, char event) {
    // Get gregorian datetime from JD
    datetime reference_dt;
    jd_to_gregorian(jd, utc_offset, &reference_dt);

    // Calculate UTC Offset estimate if not given
    double temp_utc_offset = utc_offset;
    if (utc_offset == 0)
        temp_utc_offset = floor(longitude / 15) - 1;
    
    // Set the reference day of year
    int reference_doy = day_of_year(reference_dt.year, reference_dt.month, reference_dt.day);

    int status = 0;
    int i = 0;
    while(1) {
        // Shift reference datetime
        datetime new_datetime;
        new_datetime = add_days(reference_dt, i);

        // Set temp_suntime by sending in the shifted reference datetime
        datetime temp_suntime;
        status = sunrise_or_sunset(new_datetime, temp_utc_offset, latitude, longitude, elevation, 
                                temperature, pressure, 
                                event, angle_deg, &temp_suntime
        );

        if (status != 0)
            return INVALID_DATETIME;

        datetime temp_suntime_with_estimate_offset = add_days(temp_suntime, (double)i + temp_utc_offset / 24.0);

        int temp_suntime_doy = day_of_year(temp_suntime_with_estimate_offset.year, 
                                temp_suntime_with_estimate_offset.month, temp_suntime_with_estimate_offset.day);

        if ((temp_suntime_doy < reference_doy && temp_suntime.year == reference_dt.year) || 
                                (temp_suntime_with_estimate_offset.year < reference_dt.year)) {
            i++;
        }
        else {
            return add_days(temp_suntime, utc_offset / 24.0);
        }
    }
}

/* Python wrapper for find_proper_suntime. */
PyObject* py_find_proper_suntime(PyObject* self, PyObject* args) {
    double jd, latitude, longitude, elevation, temperature, pressure, utc_offset, angle_deg;
    int event_code;
    char event;
    if (!PyArg_ParseTuple(args, "ddddddddC", &jd, &latitude, &longitude, 
                                              &elevation, &temperature, &pressure, 
                                              &utc_offset, &angle_deg, &event_code))
        return NULL;

    event = (char)event_code;
    
    return datetime_to_pydatetime(
        find_proper_suntime(jd, utc_offset, latitude, longitude, elevation, temperature, pressure, angle_deg, event));
}

/* Specialized rise/set solver used by visibility calculations.
 * Also returns per-day deltaPsi and obliquity values.
 */

int sunrise_or_sunset_w_nutation(datetime date, double utc_offset, double local_latitude, double local_longitude,
    double elevation, double temperature, double pressure, char event_type, double angle_deg, 
    double* deltaPsi, double* true_obliquity, datetime* sun_event) {

    double lat_rad = RADIANS(local_latitude);
    double new_jd = gregorian_to_jd(date, 0) - fraction_of_day_datetime(date);
    double new_deltaT = delta_t_approx(date.year, date.month);

    SunResult sun_params[3];
    for (int i = 0; i < 3; i++) {
        datetime temp_date;
        jd_to_gregorian(new_jd + i - 1, utc_offset, &temp_date);
        double t_deltaT = delta_t_approx(temp_date.year, temp_date.month);
        double t_jde = (new_jd + i - 1) + t_deltaT / SECONDS_IN_DAY;
        compute_sun_result(t_jde, t_deltaT, 
            local_latitude, local_longitude, elevation, temperature, pressure, 
            &sun_params[i]
        );
        deltaPsi[i] = sun_params[i].nutation_longitude;
        true_obliquity[i] = sun_params[i].true_obliquity;
    }

    double h_zero_rad = RADIANS(-angle_deg);
    double cosH_zero = (sin(h_zero_rad) - sin(lat_rad) * sin(RADIANS(sun_params[1].apparent_declination))) / \
                        (cos(lat_rad) * cos(RADIANS(sun_params[1].apparent_declination)));

    double H_zero_rad, H_zero_deg;
    if (cosH_zero < 1.0 && cosH_zero > -1.0) {
        H_zero_rad = acos(cosH_zero);
        H_zero_deg = H_zero_rad * 180 / M_PI;
    }
    else
        return -1;
    
    double sidereal_time = greenwich_mean_sidereal_time(new_jd);
    double m0 = (sun_params[1].apparent_right_ascension - local_longitude - sidereal_time) / 360.0;

    double m_event;
    if (event_type == 'r')
        m_event = m0 - H_zero_deg / 360;
    else if (event_type == 's')
        m_event = m0 + H_zero_deg / 360;
    else
        return -2;
    
    for (int i = 0; i < 3; i++) {
        double theta_event_deg = normalize_angle(sidereal_time + 360.985647 * m_event);
        double n_event = m_event + new_deltaT / SECONDS_IN_DAY;
        double interp_dec_event_rad = angle_interpolation(n_event,
                                    sun_params[0].apparent_declination,
                                    sun_params[1].apparent_declination,
                                    sun_params[2].apparent_declination) * M_PI / 180;
        double interp_ra_event_deg = angle_interpolation(n_event,
                                    sun_params[0].apparent_right_ascension,
                                    sun_params[1].apparent_right_ascension,
                                    sun_params[2].apparent_right_ascension);

        double lha_event_rad = normalize_angle(theta_event_deg - (-local_longitude) - interp_ra_event_deg) * M_PI / 180;
        double sun_alt_deg = asin(sin(lat_rad) * sin(interp_dec_event_rad) + 
                                cos(lat_rad) * cos(interp_dec_event_rad) * cos(lha_event_rad)) * 180 / M_PI;

        double deltaM = (sun_alt_deg - h_zero_rad * 180 / M_PI) / (360 * cos(interp_dec_event_rad) * cos(lat_rad) * sin(lha_event_rad));
        m_event += deltaM;
    }

    sun_event->year = date.year;
    sun_event->month = date.month;
    sun_event->day = date.day;
    sun_event->hour = 0; sun_event->minute = 0; sun_event->second = 0; sun_event->microsecond = 0; 
    *sun_event = add_days(*sun_event, m_event);

    return 0;
}

/*
 * Companion helper for visibility workflows that also exposes nutation terms.
 * Returns INVALID_DATETIME when no valid event is found.
 */
datetime find_proper_suntime_w_nutation(double jd, double utc_offset, double latitude, double longitude, double elevation, 
    double temperature, double pressure, double angle_deg, char event, double* deltaPsi, double* true_obliquity) {
    // Get gregorian datetime from JD
    datetime reference_dt;
    jd_to_gregorian(jd, utc_offset, &reference_dt);

    // Calculate UTC Offset estimate if not given
    double temp_utc_offset = utc_offset;
    if (utc_offset == 0)
        temp_utc_offset = floor(longitude / 15) - 1;
    
    // Set the reference day of year
    int reference_doy = day_of_year(reference_dt.year, reference_dt.month, reference_dt.day);

    int status = 0;
    int i = 0;
    while(1) {
        // Shift reference datetime
        datetime new_datetime;
        new_datetime = add_days(reference_dt, i);

        // Set temp_suntime by sending in the shifted reference datetime
        datetime temp_suntime;
        status = sunrise_or_sunset_w_nutation(new_datetime, temp_utc_offset, latitude, longitude, elevation, 
                                temperature, pressure, event, angle_deg, deltaPsi, true_obliquity, &temp_suntime);

        if (status != 0)
            return INVALID_DATETIME;

        datetime temp_suntime_with_estimate_offset = add_days(temp_suntime, (double)i + temp_utc_offset / 24.0);

        int temp_suntime_doy = day_of_year(temp_suntime_with_estimate_offset.year, 
                                temp_suntime_with_estimate_offset.month, temp_suntime_with_estimate_offset.day);

        if ((temp_suntime_doy < reference_doy && temp_suntime.year == reference_dt.year) || 
                                (temp_suntime_with_estimate_offset.year < reference_dt.year)) {
            i++;
        }
        else {
            return add_days(temp_suntime, utc_offset / 24.0);
        }
    }
}
