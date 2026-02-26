#include "c_datetime.h"

// Cumulative days up to the first day of each month (0-indexed, non-leap year).
static const int month_days[] = { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 };
static const int days_per_month[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

/* Fractional UTC day represented by the datetime fields. */
double fraction_of_day_datetime(datetime date) {
    return (((double)date.hour / 24.0) +
            ((double)date.minute / 1440.0) +
            ((double)date.second / 86400.0) +
            ((double)date.microsecond / 86400000000.0));
}

/* Gregorian leap-year rule. */
int is_leap_year(int year) {
    return (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0));
}

/* Day-of-year in [1, 366], or -1 for obviously invalid month/day inputs.
 * This helper performs coarse validation only (it does not reject month/day
 * combinations like April 31).
 */
int day_of_year(int year, int month, int day) {
    if (month < 1 || month > 12 || day < 1 || day > 31) return -1;
    int doy = month_days[month - 1] + day;
    if (month > 2 && is_leap_year(year)) doy += 1;
    return doy;
}

/* Month length for Gregorian calendar. */
int days_in_month(int year, int month) {
    return (month == 2 && is_leap_year(year)) ? 29 : days_per_month[month - 1];
}

/* Carry/borrow normalize datetime fields into canonical ranges. */
void normalize_datetime(datetime *dt) {
    // Normalize microseconds to seconds
    dt->second += dt->microsecond / 1000000;
    dt->microsecond %= 1000000;
    if (dt->microsecond < 0) {
        dt->second -= 1;
        dt->microsecond += 1000000;
    }

    // Normalize seconds to minutes
    dt->minute += dt->second / 60;
    dt->second %= 60;
    if (dt->second < 0) {
        dt->minute -= 1;
        dt->second += 60;
    }

    // Normalize minutes to hours
    dt->hour += dt->minute / 60;
    dt->minute %= 60;
    if (dt->minute < 0) {
        dt->hour -= 1;
        dt->minute += 60;
    }

    // Normalize hours to days
    dt->day += dt->hour / 24;
    dt->hour %= 24;
    if (dt->hour < 0) {
        dt->day -= 1;
        dt->hour += 24;
    }

    // Normalize day, month, year
    while (dt->day > days_in_month(dt->year, dt->month)) {
        dt->day -= days_in_month(dt->year, dt->month);
        dt->month++;
        if (dt->month > 12) {
            dt->month = 1;
            dt->year++;
        }
    }

    while (dt->day <= 0) {
        dt->month--;
        if (dt->month <= 0) {
            dt->month = 12;
            dt->year--;
        }
        dt->day += days_in_month(dt->year, dt->month);
    }
}

/* Add a fractional number of days and return the normalized result. */
datetime add_days(datetime dt, double days_to_add) {
    int whole_days = (int)days_to_add;
    double fractional_day = days_to_add - (double)whole_days;

    dt.day += whole_days;

    // Convert fractional day to h:m:s:us and let normalization handle overflow/underflow.
    double total_hours = fractional_day * 24.0;
    int whole_hours = (int)floor(total_hours);
    
    double total_minutes = (total_hours - whole_hours) * 60.0;
    int whole_minutes = (int)floor(total_minutes);

    double total_seconds = (total_minutes - whole_minutes) * 60.0;
    int whole_seconds = (int)floor(total_seconds);

    double total_micro = (total_seconds - whole_seconds) * 1e6;
    int whole_micro = (int)floor(total_micro);

    dt.hour += whole_hours;
    dt.minute += whole_minutes;
    dt.second += whole_seconds;
    dt.microsecond += whole_micro;

    normalize_datetime(&dt);
    return dt;
}

/* Extract a Python datetime object into the native datetime struct. */
void fill_in_datetime_values(datetime* date, PyObject* input_datetime) {
    int year        = PyDateTime_GET_YEAR(input_datetime);
    int month       = PyDateTime_GET_MONTH(input_datetime);
    int day         = PyDateTime_GET_DAY(input_datetime);
    int hour        = PyDateTime_DATE_GET_HOUR(input_datetime);
    int minute      = PyDateTime_DATE_GET_MINUTE(input_datetime);
    int second      = PyDateTime_DATE_GET_SECOND(input_datetime);
    int microsecond = PyDateTime_DATE_GET_MICROSECOND(input_datetime);

    date->year = year;
    date->month = month;
    date->day = day;
    date->hour = hour;
    date->minute = minute;
    date->second = second;
    date->microsecond = microsecond;
}

/* Build a Python datetime object from native datetime fields. */
PyObject* datetime_to_pydatetime(datetime dt) {
    ENSURE_PYDATETIME()
    // Build the Python datetime object from the struct fields
    PyObject *py_dt = PyDateTime_FromDateAndTime(
        dt.year,
        dt.month,
        dt.day,
        dt.hour,
        dt.minute,
        dt.second,
        dt.microsecond
    );

    if (!py_dt)
        return NULL;

    return py_dt;
}

/*
 * Compare two datetimes.
 * Returns:
 *   1  if *a is later than *b
 *  -1  if *a is earlier than *b
 *   0  if equal
 */
int compare_datetime(const datetime* a, const datetime* b) {
    if (a->year != b->year) return (a->year < b->year) ? -1 : 1;
    if (a->month != b->month) return (a->month < b->month) ? -1 : 1;
    if (a->day != b->day) return (a->day < b->day) ? -1 : 1;
    if (a->hour != b->hour) return (a->hour < b->hour) ? -1 : 1;
    if (a->minute != b->minute) return (a->minute < b->minute) ? -1 : 1;
    if (a->second != b->second) return (a->second < b->second) ? -1 : 1;
    if (a->microsecond != b->microsecond) return (a->microsecond < b->microsecond) ? -1 : 1;
    return 0;
}
