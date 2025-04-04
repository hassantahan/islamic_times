#ifndef MOON_EQUATIONS_H
#define MOON_EQUATIONS_H

#include "astro_core.h"
#include "c_time_equations.h"
#include "c_calculation_equations.h"

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

/* Python wrappers */
PyObject* py_compute_moon(PyObject* self, PyObject* const* args, Py_ssize_t nargs);

#endif