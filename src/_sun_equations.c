#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <math.h>
#include <string.h>
#include "astro_core.h"
#include "_time_equations.h"
#include "_calculation_equations.h"

/* ================================
   Definitions & Helper Constants
   ================================ */

#define M_PI 3.14159265358979
#define DEG2RAD M_PI / 180
#define J2000 2451545.0
#define JULIAN_CENTURY 36525.0  // as used in your _compute_time()
#define JULIAN_MILLENNIUM 365250.0  // as used in your _compute_time()
#define OBLIQUITY_TERMS_SIZE 10
#define SUN_NUTATION_COEFFICIENTS_LOOP 63
#define SUN_NUTATION_ARGUMENTS_LOOP 5
#define ASTRONOMICAL_UNIT_KM 149597870.7
#define EARTH_RADIUS_KM 6378.14


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
   SunResult Structure
   ================================ */

typedef struct {
    double t;
    double mean_longitude;
    double mean_anomaly;
    double earth_orbit_eccentricity;
    double sun_centre;
    double true_longitude;
    double true_anomaly;
    double geocentric_distance;
    double omega;
    double apparent_longitude;
    double nutation_longitude;
    double nutation_obliquity;
    double mean_obliquity;
    double true_obliquity;
    double true_right_ascension;
    double true_declination;
    double apparent_right_ascension;
    double apparent_declination;
    double greenwich_hour_angle;
    double local_hour_angle;
    double eh_parallax;
    double topocentric_ascension;
    double topocentric_declination;
    double topocentric_local_hour_angle;
    double true_altitude;
    double true_azimuth;
    double apparent_altitude;
} SunResult;


/* ================================
   Core Calculation: sun_nutation
   ================================ */

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

    for (int i = 0; i < SUN_NUTATION_COEFFICIENTS_LOOP; ++i) {  // Replace 63 with the actual count of terms
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
   Core Calculation: oblique_eq
   ================================ */

double oblique_eq(double t) {
    // Calculate the obliquity of the ecliptic for a given Julian Ephemeris Day. See Chapter 22 of *Astronomical Algorthims* for more information.
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
   Core Calculation: compute_sun_result
   ================================ */

void compute_sun_result(double jde, double deltaT, double local_latitude, double local_longitude,
                        double elevation, double temperature, double pressure,
                        SunResult* result) {
    // Compute time variable
    result->t = (jde - J2000) / JULIAN_MILLENNIUM;
    double t = result->t;
    double t2 = t * t;
    double t3 = t2 * t;
    double t4 = t3 * t;
    double t5 = t4 * t;
    
    // --- Orbital Elements ---
    // Mean longitude
    double mean_long = normalize_angle(280.4664567 + 360007.6982779 * t + 0.03032028 * t2 + t3 / 49931.0 - t4 / 15300.0 - t5 / 2000000.0);
    result->mean_longitude = mean_long;
    
    // Mean anomaly
    double mean_anom = normalize_angle(357.52911 + 359990.50340 * t - 0.001603 * t2 - t3 / 30000.0);
    result->mean_anomaly = mean_anom;
    
    // Eccentricity of Earth's orbit
    double ecc = 0.016708634 - 0.00042037 * t - 0.000001267 * t2;
    result->earth_orbit_eccentricity = ecc;
    
    // Sun centre (the equation of the centre)
    double sun_center_val = (1.914602 - 0.04817 * t - 0.000014 * t2) * sin(mean_anom * M_PI / 180.0) +
                            (0.019993 - 0.000101 * t) * sin(2 * mean_anom * M_PI / 180.0) +
                            0.000289 * sin(3 * mean_anom * M_PI / 180.0);
    result->sun_centre = sun_center_val;
    
    // True longitude and anomaly
    double true_long = mean_long + sun_center_val;
    result->true_longitude = true_long;
    double true_anom = mean_anom + sun_center_val;
    result->true_anomaly = true_anom;
    
    // Geocentric distance (in AU)
    double geo_dist = (1.000001018 * (1 - ecc * ecc)) / (1 + ecc * cos(true_anom * M_PI / 180.0));
    result->geocentric_distance = geo_dist;
    
    // --- Nutation and Obliquity ---
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
    
    // Apparent longitude
    double app_long = true_long - 0.00569 - 0.00478 * sin(omega * M_PI / 180.0);
    result->apparent_longitude = app_long;
    
    // --- Apparent Coordinates (Right Ascension & Declination) ---
    // Using the simplified formulas:
    double ra = atan2(cos(result->mean_obliquity * M_PI / 180.0) * sin(true_long * M_PI / 180.0),
                      cos(true_long * M_PI / 180.0));
    double dec = asin(sin(result->mean_obliquity * M_PI / 180.0) * sin(true_long * M_PI / 180.0));
    
    // Convert RA (in radians) to degrees; note that in your Python code RA is later divided by 15 to get hours.
    result->true_right_ascension = normalize_angle(ra * 180.0 / M_PI);
    result->true_declination = dec * 180.0 / M_PI;
    
    // Apparent RA/Dec
    double epsilon_corr = result->true_obliquity + 0.00256 * cos(omega * M_PI / 180.0); // obliquity correction
    double app_long_rad = app_long * M_PI / 180.0; // apparent longitude in radians
    double epsilon_corr_rad = epsilon_corr * M_PI / 180.0; // corrected obliquity in radians

    double app_ra = atan2(cos(epsilon_corr_rad) * sin(app_long_rad), cos(app_long_rad));
    double app_dec = asin(sin(epsilon_corr_rad) * sin(app_long_rad));

    result->apparent_right_ascension = normalize_angle(app_ra * 180.0 / M_PI);
    result->apparent_declination = app_dec * 180.0 / M_PI;
    
    // --- Hour Angles ---
    double sidereal_time = greenwich_mean_sidereal_time(jde - deltaT / 86400);
    double st_correction = dp * 3600 * cos(result->true_obliquity * M_PI / 180) / 15;
    double gha = normalize_angle(sidereal_time + st_correction / 240);
    result->greenwich_hour_angle = gha;

    result->local_hour_angle = normalize_angle(gha + local_longitude - result->apparent_right_ascension);
    
    // --- Topocentric Corrections ---
    result->eh_parallax = asin(sin(8.794 / 3600 * M_PI / 180) / result->geocentric_distance) * 180 / M_PI;
    double top_ra = result->apparent_right_ascension; 
    double top_dec = result->apparent_declination;

    correct_ra_dec(&top_ra, &top_dec, result->local_hour_angle, result->eh_parallax, local_latitude, elevation / 1000, EARTH_RADIUS_KM);

    result->topocentric_ascension = top_ra;
    result->topocentric_declination = top_dec;
    result->topocentric_local_hour_angle = normalize_angle(gha + local_longitude - result->topocentric_ascension);
    
    // --- Horizontal Coordinates ---
    // True
    double top_dec_rad = top_dec * M_PI / 180;
    double lat_rad = local_latitude * M_PI / 180;
    double lha_rad = result->local_hour_angle * M_PI / 180;
    
    double true_alt_rad = asin(
        sin(lat_rad) * sin(top_dec_rad) + 
        cos(lat_rad) * cos(top_dec_rad) * cos(lha_rad)
    );

    double true_az_ra = atan2(
        -cos(top_dec_rad) * sin(lha_rad),
        sin(top_dec_rad) * cos(lat_rad) - 
        cos(top_dec_rad) * sin(lat_rad) * cos(lha_rad)
    );

    result->true_altitude = true_alt_rad * 180 / M_PI;
    result->true_azimuth = normalize_angle(true_az_ra * 180 / M_PI);

    // Apparent
    // double refraction_factor = 1.02 / tan((result->true_altitude + 10.3 / (result->true_altitude + 5.11)) * M_PI / 180) * pressure / 101.325 *  283 / (273 + temperature);
    result->apparent_altitude = result->true_altitude ;//+ refraction_factor / 60;
}

/* ================================
   Python Wrapper Functions
   ================================ */

// compute_sun(jde, deltaT, latitude, longitude, elevation, temperature, pressure)
PyObject* py_compute_sun(PyObject* self, PyObject* const* args, Py_ssize_t nargs) {
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

    PyObject* sun = PyObject_CallFunctionObjArgs((PyObject *)SunType,
        // FLOAT(jde),
        // FLOAT(deltaT),
        // ANGLE(latitude),
        // ANGLE(longitude),
        // DIST(elevation, UNIT("METRE")),
        // FLOAT(temperature),
        // FLOAT(pressure),
        // FLOAT(result.t),
        ANGLE(result.mean_longitude),
        ANGLE(result.mean_anomaly),
        FLOAT(result.earth_orbit_eccentricity),
        ANGLE(result.sun_centre),
        ANGLE(result.true_longitude),
        ANGLE(result.true_anomaly),
        DIST(result.geocentric_distance, UNIT("AU")),
        ANGLE(result.omega),
        ANGLE(result.apparent_longitude),
        Py_BuildValue("(OO)", ANGLE(result.nutation_longitude), ANGLE(result.nutation_obliquity)),
        ANGLE(result.nutation_obliquity),
        ANGLE(result.mean_obliquity),
        ANGLE(result.true_obliquity),
        RA(result.true_right_ascension),
        ANGLE(result.true_declination),
        RA(result.apparent_right_ascension),
        ANGLE(result.apparent_declination),
        ANGLE(result.greenwich_hour_angle),
        ANGLE(result.local_hour_angle),
        ANGLE(result.eh_parallax),
        RA(result.topocentric_ascension),
        ANGLE(result.topocentric_declination),
        ANGLE(result.topocentric_local_hour_angle),
        ANGLE(result.true_altitude),
        ANGLE(result.true_azimuth),
        ANGLE(result.apparent_altitude),
        NULL
    );

    return sun;
}

// Stub for find_sun_transit: 
// For demonstration, this dummy implementation returns jd + 0.5 days.
PyObject* py_find_sun_transit(PyObject* self, PyObject* args) {
    double jd, deltaT, longitude;
    if (!PyArg_ParseTuple(args, "ddd", &jd, &deltaT, &longitude))
        return NULL;
    
    double transit_jd = jd + 0.5;  // Dummy calculation
    return PyFloat_FromDouble(transit_jd);
}

// Stub for find_proper_suntime:
// For demonstration, returns jd + 0.25 for sunrise and jd + 0.75 for sunset.
PyObject* py_find_proper_suntime(PyObject* self, PyObject* args) {
    double jd, deltaT, latitude, longitude;
    const char* event;
    if (!PyArg_ParseTuple(args, "dddds", &jd, &deltaT, &latitude, &longitude, &event))
        return NULL;
    
    double suntime = 0.0;
    if (strcmp(event, "rise") == 0 || strcmp(event, "sunrise") == 0) {
        suntime = jd + 0.25;
    } else if (strcmp(event, "set") == 0 || strcmp(event, "sunset") == 0) {
        suntime = jd + 0.75;
    } else {
        PyErr_SetString(PyExc_ValueError, "Invalid event");
        return NULL;
    }
    return PyFloat_FromDouble(suntime);
}