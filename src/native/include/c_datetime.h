#ifndef C_DATETIME_H
#define C_DATETIME_H

#include "astro_core.h"

#define INVALID_DATETIME (datetime){.year = -9999}

typedef struct {
    int year;
    int month;
    int day;
    int hour;
    int minute;
    int second;
    int microsecond;
} datetime;

int is_leap_year(int year);
int day_of_year(int year, int month, int day);
int days_in_month(int year, int month);
void normalize_datetime(datetime *dt);
datetime add_days(datetime dt, double days_to_add);
double fraction_of_day_datetime(datetime date);
void fill_in_datetime_values(datetime* date, PyObject* input_datetime);
PyObject* datetime_to_pydatetime(datetime dt);
int compare_datetime(const datetime* a, const datetime* b);

#endif