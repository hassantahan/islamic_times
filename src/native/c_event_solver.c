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
    double lat_sin = sin(latitude_rad);
    double lat_cos = cos(latitude_rad);
    double dec_sin = sin(declination_rad);
    double dec_cos = cos(declination_rad);
    double cos_h0 = (sin(horizon_altitude_rad) - lat_sin * dec_sin) / (lat_cos * dec_cos);

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
    double target_altitude_rad = RADIANS(horizon_altitude_deg);
    double target_altitude_deg = horizon_altitude_deg;
    double lat_sin = sin(latitude_rad);
    double lat_cos = cos(latitude_rad);

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
        double dec_sin = sin(interp_dec_event_rad);
        double dec_cos = cos(interp_dec_event_rad);
        double lha_sin = sin(lha_event_rad);
        double event_altitude_deg = DEGREES(
            asin(lat_sin * dec_sin + lat_cos * dec_cos * cos(lha_event_rad))
        );

        double delta_m = (event_altitude_deg - target_altitude_deg) /
                         (360.0 * dec_cos * lat_cos * lha_sin);
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
    datetime reference_day = reference_dt;
    reference_day.hour = 0;
    reference_day.minute = 0;
    reference_day.second = 0;
    reference_day.microsecond = 0;

    if (utc_offset == 0.0) {
        temp_utc_offset = floor(longitude_deg / 15.0) - 1.0;
    }

    for (int i = 0; i <= max_search_days; ++i) {
        datetime shifted_dt = add_days(reference_dt, i);
        datetime raw_event_dt;

        int status = solver_fn(shifted_dt, temp_utc_offset, solver_ctx, &raw_event_dt);
        if (status != 0) {
            return INVALID_DATETIME;
        }

        datetime normalized_event_dt = add_days(raw_event_dt, (double)i + temp_utc_offset / 24.0);
        datetime normalized_event_day = normalized_event_dt;
        normalized_event_day.hour = 0;
        normalized_event_day.minute = 0;
        normalized_event_day.second = 0;
        normalized_event_day.microsecond = 0;
        if (compare_datetime(&normalized_event_day, &reference_day) < 0) {
            continue;
        }

        return add_days(raw_event_dt, utc_offset / 24.0);
    }

    return INVALID_DATETIME;
}
