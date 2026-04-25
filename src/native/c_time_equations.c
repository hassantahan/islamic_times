#define PY_SSIZE_T_CLEAN
#include "c_time_equations.h"
#include "c_time_scale_data.h"

/* ================================
   Definitions & Helper Constants
   ================================ */

/* ================================
   GMST calculation
   Input: Julian Day (UT-based)
   Output: GMST angle in degrees [0, 360)
   ================================ */

double greenwich_mean_sidereal_time(double julian_day){
    double t = (julian_day - J2000) / JULIAN_CENTURY;
    double t2 = t * t;
    double t3 = t2 * t;

    double theta_zero = 280.46061837 + 360.98564736629 * (julian_day - J2000) 
                                     + 0.000387933 * t2 
                                     - t3 / 38710000.0;

    return normalize_angle(theta_zero);
}

/* Python wrapper: greenwich_mean_sidereal_time(jd) -> float */
PyObject* py_greenwich_mean_sidereal_time(PyObject *self, PyObject *args) {
    double jd;
    if (!PyArg_ParseTuple(args, "d", &jd)) return NULL;
    double result = greenwich_mean_sidereal_time(jd);
    return Py_BuildValue("d", result);
}

/* ================================
   Gregorian -> Julian Day
   UTC offset is provided in hours.
   ================================ */

double gregorian_to_jd(datetime date, double utc_offset) {
    double y = date.year;
    double m = date.month;
    double frac = fraction_of_day_datetime(date);
    double d = (double)date.day + frac;

    if (m <= 2) {
        y--;
        m += 12;
    }

    int a = (int)(floor(y / 100));
    int b = 2 - a + (int)(floor(a / 4));

    double jd = (int)(floor(365.25 * (y + 4716)) + (int)(floor(30.6001 * (m + 1)))) + d + b - 1524.5 - utc_offset / 24.0;

    return jd;
}

/* Python wrapper: gregorian_to_jd(datetime, utc_offset) -> float */
PyObject* py_gregorian_to_jd(PyObject *self, PyObject *args) {
    ENSURE_PYDATETIME()
    PyObject *input_datetime;
    double utc_offset;

    // Parse argument
    if (!PyArg_ParseTuple(args, "O!d", PyDateTimeAPI->DateTimeType, &input_datetime, &utc_offset)) {
        return NULL;
    }

    datetime date;
    fill_in_datetime_values(&date, input_datetime);

    double result = gregorian_to_jd(date, utc_offset);

    return Py_BuildValue("d", result);
}


/* ================================
   Julian Day -> Gregorian datetime
   Returns status code: 0 success, 1 invalid input.
   ================================ */

int jd_to_gregorian(double jd, double utc_offset, datetime* date) {
    if (jd > 0) {
        jd += 0.5 + utc_offset / 24.0;
        
        int z = (int)jd;
        double f = jd - z;

        int a, alpha;
        if (z < 2299161) {
            a = (int)z;
        }
        else {
            alpha = (int)floor((z - 1867216.25) / 36524.25);
            a = z + 1 + alpha - (int)floor(alpha / 4);
        }

        int b = a + 1524;
        int c = (int)floor((b - 122.1) / 365.25);
        int d = (int)floor(365.25 * c);
        int e = (int)floor((b - d) / 30.6001);

        double d_day = b - d - floor(30.6001 * e) + f;
        
        if (e < 14) 
            date->month = e - 1;
        else
            date->month = e - 13;

        if (date->month > 2)
            date->year = c - 4716;
        else
            date->year = c - 4715;

        date->day = (int)floor(d_day);
        double f_day = d_day - date->day;
        
        double d_hour = f_day * 24;
        date->hour = (int)d_hour;
        
        double d_minute = (d_hour - date->hour) * 60;
        date->minute = (int)d_minute;
        
        double d_second = (d_minute - date->minute) * 60;
        date->second = (int)d_second;
        
        double d_microsecond = (d_second - date->second) * 1000000;
        date->microsecond = (int)d_microsecond;

        return 0;
    }
    else {
        return 1;
    }
}

/* Python wrapper: jd_to_gregorian(jd, utc_offset) -> datetime */
PyObject* py_jd_to_gregorian(PyObject *self, PyObject *args) {
    ENSURE_PYDATETIME()
    double jd;
    double utc_offset; 

    if (!PyArg_ParseTuple(args, "dd", &jd, &utc_offset)) return NULL;

    datetime date;
    int status = jd_to_gregorian(jd, utc_offset, &date);
    if (!status) {
        return datetime_to_pydatetime(date);
    }
    else {
        PyErr_SetString(PyExc_ValueError, "'jd' must be a positive Julian Day value.");
        return NULL;
    }
}

/* ================================
   Fraction-of-day helpers
   ================================ */

double fraction_of_day_jd(double jd, double utc_offset) {
    return fmod((jd + 0.5 + utc_offset / 24), 1);
}

/* ================================
   Delta-T approximation (seconds)
   Modern dates use bundled official USNO data. Historical/out-of-range dates
   fall back to the NASA eclipse polynomial set.
   ================================ */

#define FUTURE_DELTA_T_SMOOTHING_DECAY_YEARS 20.0
#define FUTURE_DELTA_T_CONSTRAINT_RAMP_YEARS 3.0
#define HUBER_DELTA_T_Q_MS2_PER_YEAR 0.058
#define HUBER_DELTA_T_M_YEARS 2500.0
#define HUBER_DELTA_T_BASE_YEAR 2005.0
#define USNO_DELTA_T_TAIL_SLOPE_POINTS 4

static double decimal_year_from_datetime(datetime date) {
    int doy = day_of_year(date.year, date.month, date.day);
    int days_this_year = is_leap_year(date.year) ? 366 : 365;
    if (doy < 1 || days_this_year <= 0) {
        return (double)date.year;
    }
    return (double)date.year + (((double)(doy - 1)) + fraction_of_day_datetime(date)) / (double)days_this_year;
}

static double decimal_year_from_jd(double jd_utc) {
    datetime date;
    if (jd_to_gregorian(jd_utc, 0.0, &date) != 0) {
        return 2000.0;
    }
    return decimal_year_from_datetime(date);
}

static double delta_t_polynomial_for_decimal_year(double y) {
    int year = (int)floor(y);
    double u, t, dt;

    if (year < -500) {
        u = (year - 1820) / 100.0;
        dt = -20.0 + 32.0 * u * u;
    } else if (year < 500) {
        u = y / 100.0;
        dt = 10583.6 - 1014.41 * u + 33.78311 * pow(u, 2)
           - 5.952053 * pow(u, 3) - 0.1798452 * pow(u, 4)
           + 0.022174192 * pow(u, 5) + 0.0090316521 * pow(u, 6);
    } else if (year < 1600) {
        u = (y - 1000.0) / 100.0;
        dt = 1574.2 - 556.01 * u + 71.23472 * pow(u, 2)
           + 0.319781 * pow(u, 3) - 0.8503463 * pow(u, 4)
           - 0.005050998 * pow(u, 5) + 0.0083572073 * pow(u, 6);
    } else if (year < 1700) {
        t = y - 1600.0;
        dt = 120.0 - 0.9808 * t - 0.01532 * pow(t, 2)
           + pow(t, 3) / 7129.0;
    } else if (year < 1800) {
        t = y - 1700.0;
        dt = 8.83 + 0.1603 * t - 0.0059285 * pow(t, 2)
           + 0.00013336 * pow(t, 3) - pow(t, 4) / 1174000.0;
    } else if (year < 1860) {
        t = y - 1800.0;
        dt = 13.72 - 0.332447 * t + 0.0068612 * pow(t, 2)
           + 0.0041116 * pow(t, 3) - 0.00037436 * pow(t, 4)
           + 0.0000121272 * pow(t, 5) - 0.0000001699 * pow(t, 6)
           + 0.000000000875 * pow(t, 7);
    } else if (year < 1900) {
        t = y - 1860.0;
        dt = 7.62 + 0.5737 * t - 0.251754 * pow(t, 2)
           + 0.01680668 * pow(t, 3) - 0.0004473624 * pow(t, 4)
           + pow(t, 5) / 233174.0;
    } else if (year < 1920) {
        t = y - 1900.0;
        dt = -2.79 + 1.494119 * t - 0.0598939 * pow(t, 2)
           + 0.0061966 * pow(t, 3) - 0.000197 * pow(t, 4);
    } else if (year < 1941) {
        t = y - 1920.0;
        dt = 21.20 + 0.84493 * t - 0.076100 * pow(t, 2)
           + 0.0020936 * pow(t, 3);
    } else if (year < 1961) {
        t = y - 1950.0;
        dt = 29.07 + 0.407 * t - pow(t, 2) / 233.0
           + pow(t, 3) / 2547.0;
    } else if (year < 1986) {
        t = y - 1975.0;
        dt = 45.45 + 1.067 * t - pow(t, 2) / 260.0
           - pow(t, 3) / 718.0;
    } else if (year < 2005) {
        t = y - 2000.0;
        dt = 63.86 + 0.3345 * t - 0.060374 * pow(t, 2)
           + 0.0017275 * pow(t, 3) + 0.000651814 * pow(t, 4)
           + 0.00002373599 * pow(t, 5);
    } else if (year < 2050) {
        t = y - 2000.0;
        dt = 62.92 + 0.32217 * t + 0.005589 * pow(t, 2);
    } else if (year < 2150) {
        dt = -20.0 + 32.0 * pow((y - 1820.0) / 100.0, 2)
           - 0.5628 * (2150.0 - y);
    } else {
        u = (year - 1820) / 100.0;
        dt = -20.0 + 32.0 * u * u;
    }

    return dt;
}

static double delta_t_polynomial_approx(int year, int month) {
    double y = year + (month - 0.5) / 12.0;
    return delta_t_polynomial_for_decimal_year(y);
}

static double delta_t_polynomial_slope_seconds_per_year(double decimal_year) {
    const double half_step_years = 1.0 / 365.25;
    return (
        delta_t_polynomial_for_decimal_year(decimal_year + half_step_years)
        - delta_t_polynomial_for_decimal_year(decimal_year - half_step_years)
    ) / (2.0 * half_step_years);
}

static double huber_delta_t_standard_error_seconds(double decimal_year) {
    double n = decimal_year - HUBER_DELTA_T_BASE_YEAR;
    if (n <= 0.0) {
        return 0.0;
    }

    double variance_term = (n * HUBER_DELTA_T_Q_MS2_PER_YEAR / 3.0) * (1.0 + n / HUBER_DELTA_T_M_YEARS);
    if (variance_term <= 0.0) {
        return 0.0;
    }

    return 365.25 * n * sqrt(variance_term) / 1000.0;
}

static double usno_delta_t_boundary_jd_utc(void) {
    return USNO_DELTA_T_TABLE[USNO_DELTA_T_TABLE_LEN - 1].mjd + 2400000.5;
}

static double estimate_usno_delta_t_tail_slope_seconds_per_year(void) {
    int point_count = USNO_DELTA_T_TAIL_SLOPE_POINTS;
    if (point_count > USNO_DELTA_T_TABLE_LEN) {
        point_count = USNO_DELTA_T_TABLE_LEN;
    }
    if (point_count < 2) {
        return 0.0;
    }

    int start_index = USNO_DELTA_T_TABLE_LEN - point_count;
    double sum_x = 0.0;
    double sum_y = 0.0;
    double sum_xx = 0.0;
    double sum_xy = 0.0;

    for (int i = start_index; i < USNO_DELTA_T_TABLE_LEN; ++i) {
        double jd_utc = USNO_DELTA_T_TABLE[i].mjd + 2400000.5;
        double x = decimal_year_from_jd(jd_utc);
        double y = USNO_DELTA_T_TABLE[i].delta_t_seconds;

        sum_x += x;
        sum_y += y;
        sum_xx += x * x;
        sum_xy += x * y;
    }

    double count = (double)point_count;
    double denominator = count * sum_xx - sum_x * sum_x;
    if (fabs(denominator) < 1e-12) {
        return 0.0;
    }
    return (count * sum_xy - sum_x * sum_y) / denominator;
}

static double delta_t_future_smoothed_for_jd(double jd_utc) {
    double boundary_jd = usno_delta_t_boundary_jd_utc();
    if (jd_utc <= boundary_jd) {
        return USNO_DELTA_T_TABLE[USNO_DELTA_T_TABLE_LEN - 1].delta_t_seconds;
    }

    double decimal_year = decimal_year_from_jd(jd_utc);
    double boundary_decimal_year = decimal_year_from_jd(boundary_jd);
    double delta_years = decimal_year - boundary_decimal_year;

    double baseline = delta_t_polynomial_for_decimal_year(decimal_year);
    double baseline_boundary = delta_t_polynomial_for_decimal_year(boundary_decimal_year);
    double boundary_delta_t = USNO_DELTA_T_TABLE[USNO_DELTA_T_TABLE_LEN - 1].delta_t_seconds;
    double boundary_offset = boundary_delta_t - baseline_boundary;

    double baseline_slope = delta_t_polynomial_slope_seconds_per_year(boundary_decimal_year);
    double table_slope = estimate_usno_delta_t_tail_slope_seconds_per_year();
    double tau = FUTURE_DELTA_T_SMOOTHING_DECAY_YEARS;

    double raw_correction = (
        boundary_offset + (table_slope - baseline_slope + boundary_offset / tau) * delta_years
    ) * exp(-delta_years / tau);

    double sigma = huber_delta_t_standard_error_seconds(decimal_year);
    double sigma_boundary = huber_delta_t_standard_error_seconds(boundary_decimal_year);
    double envelope = sigma;
    if (fabs(boundary_offset) > sigma_boundary) {
        envelope += (fabs(boundary_offset) - sigma_boundary) * exp(-delta_years / tau);
    }

    if (fabs(raw_correction) > envelope) {
        double limited_correction = copysign(envelope, raw_correction);
        double ramp_weight = 1.0 - exp(-delta_years / FUTURE_DELTA_T_CONSTRAINT_RAMP_YEARS);
        raw_correction += ramp_weight * (limited_correction - raw_correction);
    }

    return baseline + raw_correction;
}

static int interpolate_usno_delta_t(double mjd_utc, double* delta_t_seconds) {
    if (USNO_DELTA_T_TABLE_LEN <= 0) {
        return 0;
    }

    if (mjd_utc < USNO_DELTA_T_TABLE[0].mjd || mjd_utc > USNO_DELTA_T_TABLE[USNO_DELTA_T_TABLE_LEN - 1].mjd) {
        return 0;
    }

    if (mjd_utc == USNO_DELTA_T_TABLE[USNO_DELTA_T_TABLE_LEN - 1].mjd) {
        *delta_t_seconds = USNO_DELTA_T_TABLE[USNO_DELTA_T_TABLE_LEN - 1].delta_t_seconds;
        return 1;
    }

    for (int i = 0; i < USNO_DELTA_T_TABLE_LEN - 1; ++i) {
        const DeltaTTableEntry left = USNO_DELTA_T_TABLE[i];
        const DeltaTTableEntry right = USNO_DELTA_T_TABLE[i + 1];
        if (mjd_utc < left.mjd || mjd_utc > right.mjd) {
            continue;
        }
        if (right.mjd == left.mjd) {
            *delta_t_seconds = right.delta_t_seconds;
            return 1;
        }
        double alpha = (mjd_utc - left.mjd) / (right.mjd - left.mjd);
        *delta_t_seconds = left.delta_t_seconds + alpha * (right.delta_t_seconds - left.delta_t_seconds);
        return 1;
    }

    *delta_t_seconds = USNO_DELTA_T_TABLE[USNO_DELTA_T_TABLE_LEN - 1].delta_t_seconds;
    return 1;
}

static double tai_minus_utc_for_jd(double jd_utc) {
    const TaiUtcTableEntry* current = &USNO_TAI_UTC_TABLE[0];

    for (int i = 0; i < USNO_TAI_UTC_TABLE_LEN; ++i) {
        if (jd_utc >= USNO_TAI_UTC_TABLE[i].jd_start) {
            current = &USNO_TAI_UTC_TABLE[i];
        } else {
            break;
        }
    }

    double mjd_utc = jd_utc - 2400000.5;
    return current->base_seconds + (mjd_utc - current->mjd_reference) * current->slope_seconds_per_day;
}

double tt_minus_utc_for_jd(double jd_utc) {
    return 32.184 + tai_minus_utc_for_jd(jd_utc);
}

double delta_t_for_jd(double jd_utc) {
    double mjd_utc = jd_utc - 2400000.5;
    double delta_t_seconds;
    if (interpolate_usno_delta_t(mjd_utc, &delta_t_seconds)) {
        return delta_t_seconds;
    }

    if (USNO_DELTA_T_TABLE_LEN > 0 && jd_utc > usno_delta_t_boundary_jd_utc()) {
        return delta_t_future_smoothed_for_jd(jd_utc);
    }

    datetime date;
    if (jd_to_gregorian(jd_utc, 0.0, &date) != 0) {
        return delta_t_polynomial_approx(2000, 1);
    }
    return delta_t_polynomial_for_decimal_year(decimal_year_from_datetime(date));
}

double ut1_minus_utc_for_jd(double jd_utc) {
    return tt_minus_utc_for_jd(jd_utc) - delta_t_for_jd(jd_utc);
}

void resolve_time_scales_for_jd(double jd_utc, double* tt_minus_utc, double* ut1_minus_utc, double* delta_t) {
    double resolved_tt_utc = tt_minus_utc_for_jd(jd_utc);
    double resolved_delta_t = delta_t_for_jd(jd_utc);
    double resolved_ut1_utc = resolved_tt_utc - resolved_delta_t;

    if (tt_minus_utc != NULL) {
        *tt_minus_utc = resolved_tt_utc;
    }
    if (ut1_minus_utc != NULL) {
        *ut1_minus_utc = resolved_ut1_utc;
    }
    if (delta_t != NULL) {
        *delta_t = resolved_delta_t;
    }
}

double delta_t_approx(int year, int month) {
    datetime month_start = {
        .year = year,
        .month = month,
        .day = 1,
        .hour = 0,
        .minute = 0,
        .second = 0,
        .microsecond = 0,
    };
    int month_length = days_in_month(year, month);
    double midpoint_jd = gregorian_to_jd(month_start, 0.0) + ((double)month_length / 2.0);
    return delta_t_for_jd(midpoint_jd);
}

/* Python wrapper: delta_t_approx(year, month) -> float */
PyObject* py_delta_t_approx(PyObject *self, PyObject *args) {
    int year, month; 

    if (!PyArg_ParseTuple(args, "ii", &year, &month)) return NULL;

    double result = delta_t_approx(year, month);

    return Py_BuildValue("d", result);
}

PyObject* py_resolve_time_scales(PyObject *self, PyObject *args) {
    double jd_utc;
    if (!PyArg_ParseTuple(args, "d", &jd_utc)) {
        return NULL;
    }

    double tt_minus_utc;
    double ut1_minus_utc;
    double delta_t;
    resolve_time_scales_for_jd(jd_utc, &tt_minus_utc, &ut1_minus_utc, &delta_t);

    return Py_BuildValue("ddd", tt_minus_utc, ut1_minus_utc, delta_t);
}
