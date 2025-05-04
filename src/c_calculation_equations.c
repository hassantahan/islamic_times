#define PY_SSIZE_T_CLEAN
#include "c_calculation_equations.h"

/* ================================
   Definitions & Helper Constants
   ================================ */

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


/* ================================
    Calculation of Geocentric Equitorial Coordinates 
    from Ecliptic Coordinates
   ================================ */

void compute_equitorial_coordinates(double apparent_longitude_deg, double true_obliquity_deg, double true_latitude_deg, double* ra_deg, double* dec_deg) {
    // See Chapter 13 of *Astronomical Algorithms* (pg. 93)
    double apparent_longitude_rad = RADIANS(apparent_longitude_deg);
    double true_obliquity_rad = RADIANS(true_obliquity_deg);
    double true_latitude_rad = RADIANS(true_latitude_deg);

    double ra_rad = atan2(sin(apparent_longitude_rad) * cos(true_obliquity_rad) - tan(true_latitude_rad) * sin(true_obliquity_rad),
                    cos(apparent_longitude_rad));

    double dec_rad = asin(sin(true_latitude_rad) * cos(true_obliquity_rad) + 
                        cos(true_latitude_rad) * sin(true_obliquity_rad) * sin(apparent_longitude_rad));

    *ra_deg = normalize_angle(DEGREES(ra_rad));
    *dec_deg = DEGREES(dec_rad);
}

/*  ================================
    correct_ra_dec 
    (turn ra/dec into topocentric counterparts)
    ================================ */

void correct_ra_dec(double* ra_deg, double* dec_deg, double lha_deg, double parallax_deg, double lat_deg, double elev_km, double dist_km) {
    // Correct a celestial body's Right Ascension and Declination for apparent position. See Chapter 40 of *Astronomical Algorthims* for more information.
    double a = dist_km;
    double b = a * (1 - EARTH_FLATTENING_FACTOR);

    double dec_rad = RADIANS(*dec_deg);
    double lha_rad = RADIANS(lha_deg);
    double lat_rad = RADIANS(lat_deg);
    double parallax_rad = RADIANS(parallax_deg);

    double u_rad = atan2(b / a * tan(lat_rad), 1);
    double p_sin_phi_prime = b / a * sin(u_rad) + elev_km / dist_km * sin(lat_rad);
    double p_cos_phi_prime = cos(u_rad) + elev_km / dist_km * cos(lat_rad);

    double numerator_delA = -1 * p_cos_phi_prime * sin(parallax_rad) * sin(lha_rad);
    double denominator = cos(dec_rad) - p_cos_phi_prime * sin(parallax_rad) * cos(lha_rad);
    double delA = DEGREES(atan2(numerator_delA, denominator));
    
    // Right ascension corrected
    *ra_deg = normalize_angle(*ra_deg + delA);
    
    // Declination corrected
    double numerator_dec = (sin(dec_rad) - p_sin_phi_prime * sin(parallax_rad)) * cos(RADIANS(delA));
    *dec_deg = DEGREES(atan2(numerator_dec, denominator));
}

/*  ================================
    GHA and LHA computation
    ================================ */

void compute_gha_lha(double true_obliquity_deg, double nut_lon_deg, double gmst_deg, double observer_longitude_deg, double ra_deg,
    double* gst_deg, double* gha_deg, double* lha_deg) {
    double st_correction = nut_lon_deg * 3600.0 * cos(RADIANS(true_obliquity_deg)) / 15.0;
    *gst_deg = normalize_angle(gmst_deg + st_correction / 240.0);

    *gha_deg = normalize_angle(*gst_deg - ra_deg);
    *lha_deg = normalize_angle(*gst_deg + observer_longitude_deg - ra_deg);
}

/* ================================
    Calculation of Horizontal Coordinates
    from Geocentric Equatorial Coordinates
   ================================ */

void compute_horizontal_coordinates(double ra_deg, double dec_deg, double lha_deg, double lat_deg, double* az_deg, double* alt_deg) {
    // See Chapter 13 of *Astronomical Algorithms* (pg. 93)

    double dec_rad = RADIANS(dec_deg);
    double lha_rad = RADIANS(lha_deg);
    double lat_rad = RADIANS(lat_deg);

    double az_rad = atan2(sin(lha_rad),
                          cos(lha_rad) * sin(lat_rad) - tan(dec_rad) * cos(lat_rad)) + M_PI;

    double alt_rad = asin(sin(lat_rad) * sin(dec_rad) + 
                          cos(lat_rad) * cos(dec_rad) * cos(lha_rad));
    
    *az_deg = normalize_angle(DEGREES(az_rad));
    *alt_deg = DEGREES(alt_rad);
}

void geocentric_horizontal_coordinates(double lat_deg, double body_dec_deg, double body_lha_deg, double* body_geo_alt_deg, double* body_geo_az_deg) {
    double lat_rad = RADIANS(lat_deg);
    double body_dec_rad = RADIANS(body_dec_deg);
    double body_lha_rad = RADIANS(body_lha_deg);
    
    double tan_phi = pow(1 - EARTH_FLATTENING_FACTOR, 2) * tan(lat_rad);
    double geo_lat_rad = atan2(tan_phi, 1);

    // Geocentric Altitude
    double body_geo_alt_rad = asin(sin(geo_lat_rad) * sin(body_dec_rad) + 
                                   cos(geo_lat_rad) * cos(body_dec_rad) * cos(body_lha_rad));
    
    *body_geo_alt_deg = DEGREES(body_geo_alt_rad);

    // Geocentric Azimuth
    double cos_az = (sin(body_dec_rad) - sin(body_geo_alt_rad) * sin(geo_lat_rad)) / \
                    (cos(body_geo_alt_rad) * cos(geo_lat_rad));

    double sin_az = (-cos(body_dec_rad) * sin(body_lha_rad)) / (cos(body_geo_alt_rad));
    double body_geo_az_rad = atan2(sin_az, cos_az);

    *body_geo_az_deg = normalize_angle(DEGREES(body_geo_az_rad));
}