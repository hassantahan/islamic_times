#ifndef CALCULATION_EQUATIONS_H
#define CALCULATION_EQUATIONS_H

#include "astro_core.h"
#include "datetime.h"

double normalize_angle(double julian_day);
double angle_interpolation(double n, double ang1_deg, double ang2_deg, double ang3_deg);
void correct_ra_dec(double* ra_deg, double* dec_deg, double lha_deg, double parallax_deg, double lat_deg, double elev_km, double dist_km);

#endif // CALCULATION_EQUATIONS_H