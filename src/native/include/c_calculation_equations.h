#ifndef ISLAMIC_TIMES_NATIVE_CALCULATION_EQUATIONS_H
#define ISLAMIC_TIMES_NATIVE_CALCULATION_EQUATIONS_H

#include "astro_core.h"

/* Normalize angle to [0, 360). */
double normalize_angle(double angle_deg);
/* Interpolate three angles while respecting wrap-around discontinuities. */
double angle_interpolation(double n, double ang1_deg, double ang2_deg, double ang3_deg);
/* Convert ecliptic coordinates to geocentric equatorial right ascension/declination.
 * The exported name keeps legacy "equitorial" spelling for ABI/API compatibility.
 */
void compute_equitorial_coordinates(double apparent_longitude_deg, double true_obliquity_deg, double true_latitude_deg, double* ra_deg, double* dec_deg);
/* Apply parallax/topocentric correction to right ascension and declination in-place. */
void correct_ra_dec(double* ra_deg, double* dec_deg, double lha_deg, double parallax_deg, double lat_deg, double elev_km, double dist_km);
/* Compute Greenwich sidereal time, Greenwich hour angle, and local hour angle. */
void compute_gha_lha(double true_obliquity_deg, double nut_lon_deg, double gmst_deg, double observer_longitude_deg, double ra_deg, double* gst_deg, double* gha_deg, double* lha_deg);
/* Convert equatorial coordinates to topocentric azimuth/altitude. */
void compute_horizontal_coordinates(double ra_deg, double dec_deg, double lha_deg, double lat_deg, double* az_deg, double* alt_deg);
/* Convert declination/hour-angle to geocentric horizontal azimuth/altitude. */
void geocentric_horizontal_coordinates(double lat_deg, double body_dec_deg, double body_lha_deg, double* body_geo_alt_deg, double* body_geo_az_deg);

#endif  // ISLAMIC_TIMES_NATIVE_CALCULATION_EQUATIONS_H
