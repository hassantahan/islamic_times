#define PY_SSIZE_T_CLEAN
#include "c_event_solver.h"

double normalize_event_fraction(double value) {
    double normalized = fmod(value, 1.0);
    if (normalized < 0.0) {
        normalized += 1.0;
    }
    return normalized;
}

double refine_transit_fraction(
    double initial_fraction,
    double sidereal_time_deg,
    double local_longitude_deg,
    double delta_t_seconds,
    const double right_ascension_deg[3]
) {
    double event_fraction = initial_fraction;

    for (int i = 0; i < 3; ++i) {
        double theta_event_deg = normalize_angle(sidereal_time_deg + 360.985647 * event_fraction);
        double n_event = event_fraction + delta_t_seconds / SECONDS_IN_DAY;
        double interp_ra_event_deg = angle_interpolation(
            n_event,
            right_ascension_deg[0],
            right_ascension_deg[1],
            right_ascension_deg[2]
        );

        double lha_event_deg = normalize_angle(theta_event_deg - (-local_longitude_deg) - interp_ra_event_deg);
        event_fraction -= lha_event_deg / 360.0;
    }

    return normalize_event_fraction(event_fraction);
}

int compute_event_hour_angle(
    double latitude_rad,
    double declination_deg,
    double horizon_altitude_deg,
    double* hour_angle_deg
) {
    double declination_rad = RADIANS(declination_deg);
    double horizon_altitude_rad = RADIANS(horizon_altitude_deg);
    double cos_h0 = (sin(horizon_altitude_rad) - sin(latitude_rad) * sin(declination_rad)) /
                    (cos(latitude_rad) * cos(declination_rad));

    if (!(cos_h0 < 1.0 && cos_h0 > -1.0)) {
        return -1;
    }

    *hour_angle_deg = DEGREES(acos(cos_h0));
    return 0;
}

double refine_altitude_event_fraction(
    double initial_fraction,
    double sidereal_time_deg,
    double local_longitude_deg,
    double delta_t_seconds,
    double latitude_rad,
    double horizon_altitude_deg,
    const double declination_deg[3],
    const double right_ascension_deg[3]
) {
    double event_fraction = initial_fraction;
    double horizon_altitude_rad = RADIANS(horizon_altitude_deg);

    for (int i = 0; i < 3; ++i) {
        double theta_event_deg = normalize_angle(sidereal_time_deg + 360.985647 * event_fraction);
        double n_event = event_fraction + delta_t_seconds / SECONDS_IN_DAY;
        double interp_dec_event_rad = RADIANS(
            angle_interpolation(
                n_event,
                declination_deg[0],
                declination_deg[1],
                declination_deg[2]
            )
        );
        double interp_ra_event_deg = angle_interpolation(
            n_event,
            right_ascension_deg[0],
            right_ascension_deg[1],
            right_ascension_deg[2]
        );

        double lha_event_rad = RADIANS(
            normalize_angle(theta_event_deg - (-local_longitude_deg) - interp_ra_event_deg)
        );
        double event_altitude_deg = DEGREES(
            asin(
                sin(latitude_rad) * sin(interp_dec_event_rad) +
                cos(latitude_rad) * cos(interp_dec_event_rad) * cos(lha_event_rad)
            )
        );

        double delta_m = (event_altitude_deg - DEGREES(horizon_altitude_rad)) /
                         (360.0 * cos(interp_dec_event_rad) * cos(latitude_rad) * sin(lha_event_rad));
        event_fraction += delta_m;
    }

    return event_fraction;
}

datetime find_event_on_reference_day(
    datetime reference_dt,
    double utc_offset,
    double longitude_deg,
    int max_search_days,
    DailyEventSolverFn solver_fn,
    void* solver_ctx
) {
    double temp_utc_offset = utc_offset;
    if (utc_offset == 0.0) {
        temp_utc_offset = floor(longitude_deg / 15.0) - 1.0;
    }

    int reference_doy = day_of_year(reference_dt.year, reference_dt.month, reference_dt.day);

    for (int i = 0; i <= max_search_days; ++i) {
        datetime shifted_dt = add_days(reference_dt, i);
        datetime raw_event_dt;

        int status = solver_fn(shifted_dt, temp_utc_offset, solver_ctx, &raw_event_dt);
        if (status != 0) {
            return INVALID_DATETIME;
        }

        datetime normalized_event_dt = add_days(raw_event_dt, (double)i + temp_utc_offset / 24.0);
        int normalized_doy = day_of_year(
            normalized_event_dt.year,
            normalized_event_dt.month,
            normalized_event_dt.day
        );

        if ((normalized_doy < reference_doy && raw_event_dt.year == reference_dt.year) ||
            (normalized_event_dt.year < reference_dt.year)) {
            continue;
        }

        return add_days(raw_event_dt, utc_offset / 24.0);
    }

    return INVALID_DATETIME;
}
