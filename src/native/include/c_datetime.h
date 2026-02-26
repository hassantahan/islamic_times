#ifndef ISLAMIC_TIMES_NATIVE_DATETIME_H
#define ISLAMIC_TIMES_NATIVE_DATETIME_H

#include "astro_core.h"

/* Sentinel used by event solvers when no valid datetime exists. */
#define INVALID_DATETIME (datetime){.year = -9999}

/* Naive Gregorian timestamp container (timezone is handled separately). */
typedef struct {
    int year;
    int month;
    int day;
    int hour;
    int minute;
    int second;
    int microsecond;
} datetime;

/* Gregorian leap-year rule. */
int is_leap_year(int year);
/* Day-of-year in [1, 366], or -1 for obviously invalid month/day inputs.
 * Performs coarse validation only (month/day combinations are not fully checked).
 */
int day_of_year(int year, int month, int day);
/* Month length for Gregorian calendar. */
int days_in_month(int year, int month);
/* Carry/borrow normalize datetime fields in-place. */
void normalize_datetime(datetime *dt);
/* Add fractional days and return normalized datetime. */
datetime add_days(datetime dt, double days_to_add);
/* Fractional day represented by hour/minute/second/microsecond. */
double fraction_of_day_datetime(datetime date);
/* Convert Python datetime object fields into native datetime struct. */
void fill_in_datetime_values(datetime* date, PyObject* input_datetime);
/* Convert native datetime struct into Python datetime object. */
PyObject* datetime_to_pydatetime(datetime dt);
/* Compare datetimes: -1 earlier, 0 equal, 1 later. */
int compare_datetime(const datetime* a, const datetime* b);

#endif  // ISLAMIC_TIMES_NATIVE_DATETIME_H
