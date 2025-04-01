#ifndef CALCULATION_EQUATIONS_H
#define CALCULATION_EQUATIONS_H

double normalize_angle(double julian_day);
void correct_ra_dec(double* ra_deg, double* dec_deg, double lha_deg, double parallax_deg, double lat_deg, double elev_km, double dist_km);

#endif // CALCULATION_EQUATIONS_H