#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <math.h>
#include <string.h>

/* ================================
   Definitions & Helper Constants
   ================================ */

#define M_PI 3.14159265358979
#define EARTH_FLATTENING_FACTOR 0.003352810665

/*  ================================
    Proper modding of angles
    ================================ */

double normalize_angle(double angle_deg) {
    double res = fmod(angle_deg, 360.0);
    if (res < 0)
        res += 360.0;
    return res;
}

/* ================================
   Properly interpret angles
   ================================ */

double angle_diff(double ang1_deg, double ang2_deg) {
    double diff = normalize_angle(ang2_deg - ang1_deg);
    if (diff > 180.0)
        diff -= 360.0;
    
    return diff;
}

double angle_interpolation(double n, double ang1_deg, double ang2_deg, double ang3_deg) {
    double a = angle_diff(ang1_deg, ang2_deg);
    double b = angle_diff(ang2_deg, ang3_deg);

    double c = b - a;

    double val = ang2_deg + n / 2 * (a + b + n * c);

    return normalize_angle(val + 360);
}


/*  ================================
    correct_ra_dec (turn ra/dec into topocentric counterparts)
    ================================ */

void correct_ra_dec(double* ra_deg, double* dec_deg, double lha_deg, double parallax_deg, double lat_deg, double elev_km, double dist_km) {
    // Correct the Moon's Right Ascension and Declination for apparent position. See Chapter 40 of *Astronomical Algorthims* for more information.
    double a = dist_km;
    double b = a * (1 - EARTH_FLATTENING_FACTOR);

    double ra_rad = *ra_deg * M_PI / 180;
    double dec_rad = *dec_deg * M_PI / 180;
    double lha_rad = lha_deg * M_PI / 180;
    double lat_rad = lat_deg * M_PI / 180;
    double parallax_rad = parallax_deg * M_PI / 180;

    double u_rad = atan2(b / a * tan(lat_rad), 1);
    double p_sin_phi_prime = b / a * sin(u_rad) + elev_km / dist_km * sin(lat_rad);
    double p_cos_phi_prime = cos(u_rad) + elev_km / dist_km * cos(lat_rad);

    double numerator_delA = -1 * p_cos_phi_prime * sin(parallax_rad) * sin(lha_rad);
    double denominator = cos(dec_rad) - p_cos_phi_prime * sin(parallax_rad) * cos(lha_rad);
    double delA = atan2(numerator_delA, denominator) * 180 / M_PI;

    // Right ascension corrected
    *ra_deg += delA;

    double numerator_dec = (sin(dec_rad) - p_sin_phi_prime * sin(parallax_rad)) * cos(delA * M_PI / 180);

    // Declination corrected
    *dec_deg = atan2(numerator_dec, denominator) * 180 / M_PI;

}