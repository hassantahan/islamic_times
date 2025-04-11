#ifndef CALCULATION_EQUATIONS_H
#define CALCULATION_EQUATIONS_H

#include "astro_core.h"

double normalize_angle(double julian_day);
double angle_interpolation(double n, double ang1_deg, double ang2_deg, double ang3_deg);
void compute_equitorial_coordinates(double apparent_longitude_deg, double true_obliquity_deg, double true_latitude_deg, double* ra_deg, double* dec_deg);
void correct_ra_dec(double* ra_deg, double* dec_deg, double lha_deg, double parallax_deg, double lat_deg, double elev_km, double dist_km);
void compute_gha_lha(double true_obliquity_deg, double nut_lon_deg, double gmst_deg, double observer_longitude_deg, double ra_deg, double* gst_deg, double* gha_deg, double* lha_deg);
void compute_horizontal_coordinates(double ra_deg, double dec_deg, double lha_deg, double lat_deg, double* az_deg, double* alt_deg);
void geocentric_horizontal_coordinates(double lat_deg,double body_dec_deg, double body_lha_deg, double* body_geo_alt_deg, double* body_geo_az_deg);

#endif // CALCULATION_EQUATIONS_H