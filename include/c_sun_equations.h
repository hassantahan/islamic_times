#ifndef SUN_EQUATIONS_H
#define SUN_EQUATIONS_H

#include "astro_core.h"
#include "c_time_equations.h"
#include "c_calculation_equations.h"

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

void compute_sun_result(double jde, double deltaT, double local_latitude, double local_longitude,
    double elevation, double temperature, double pressure,
    SunResult* result);

int sunrise_or_sunset_w_nutation(datetime date, double utc_offset, double local_latitude, double local_longitude,
    double elevation, double temperature, double pressure, char event_type, double angle_deg, double* deltaPsi, double* true_obliquity, datetime* sun_event);

PyObject* py_compute_sun(PyObject* self, PyObject* const* args, Py_ssize_t nargs);
PyObject* py_find_sun_transit(PyObject* self, PyObject* args);
PyObject* py_find_proper_suntime(PyObject* self, PyObject* args);

#endif // SUN_EQUATIONS_H