#ifndef MOON_EQUATIONS_H
#define MOON_EQUATIONS_H

// Includes
#include "c_sun_equations.h"

#define CALCULATE_SUN_PARAMS_FOR_MOON_TIME -123456.0

/* ================================
   Helper structs for calculation
   ================================ */

typedef struct {
    double fundamental_arguments[5];
    double sum_l;
    double sum_b;
    double sum_r;
} MoonNutationResult;

   
typedef struct {
    // Orbital elements
	double true_longitude;
	double true_latitude;
	double geocentric_distance;

	// Nutation and obliquity
	MoonNutationResult lunar_nutation;
	double omega;
	double apparent_longitude;
	double deltaPsi;
	double true_obliquity;

	// Apparent coordinates
	double right_ascension;
	double declination;

	// Hour angles
	double greenwich_hour_angle;
	double local_hour_angle;

	// Topocentric quantities
	double eh_parallax;
	double topocentric_ascension;
	double top_declination;
	double topocentric_local_hour_angle;

	// Horizontal coordinates
	double true_altitude;
	double true_azimuth;
	double apparent_altitude;
} MoonResult;


MoonNutationResult moon_nutation(double jde);
void compute_moon_result(double jde, double deltaT, double local_latitude, double local_longitude,
    double elevation, double temperature, double pressure, 
    double deltaPsi, double ecliptic, 
    MoonResult* result);
datetime find_proper_moontime(double jd, double utc_offset, double latitude, double longitude, double elevation, 
		double temperature, double pressure, char event, double deltaPsi[3], double true_obliquity[3]);
void next_phases_of_moon_utc(datetime date, datetime phases[4]);

/* Python wrappers */
PyObject* py_compute_moon(PyObject* self, PyObject* const* args, Py_ssize_t nargs);
PyObject* py_find_moon_transit(PyObject* self, PyObject* args);
PyObject* py_find_proper_moontime(PyObject* self, PyObject* args); 
PyObject* py_next_phases_of_moon_utc(PyObject* self, PyObject* args);

#endif