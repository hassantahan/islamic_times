#ifndef ISLAMIC_TIMES_NATIVE_EVENT_SOLVER_H
#define ISLAMIC_TIMES_NATIVE_EVENT_SOLVER_H

#include "c_time_equations.h"

typedef int (*DailyEventSolverFn)(datetime date, double utc_offset, void* ctx, datetime* out_event);

/* Normalize an event fraction into [0, 1). */
double normalize_event_fraction(double value);

/* Shared Newton-style refinement used by solar/lunar transit solvers. */
double refine_transit_fraction(
    double initial_fraction,
    double sidereal_time_deg,
    double local_longitude_deg,
    double delta_t_seconds,
    const double right_ascension_deg[3]
);

/* Compute event hour angle in degrees for a target horizon altitude. */
int compute_event_hour_angle(
    double latitude_rad,
    double declination_deg,
    double horizon_altitude_deg,
    double* hour_angle_deg
);

/* Shared Newton-style refinement used by rise/set solvers. */
double refine_altitude_event_fraction(
    double initial_fraction,
    double sidereal_time_deg,
    double local_longitude_deg,
    double delta_t_seconds,
    double latitude_rad,
    double horizon_altitude_deg,
    const double declination_deg[3],
    const double right_ascension_deg[3]
);

/* Shared "find first event on reference civil day" search loop. */
datetime find_event_on_reference_day(
    datetime reference_dt,
    double utc_offset,
    double longitude_deg,
    int max_search_days,
    DailyEventSolverFn solver_fn,
    void* solver_ctx
);

#endif  // ISLAMIC_TIMES_NATIVE_EVENT_SOLVER_H
