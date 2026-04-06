#ifndef ISLAMIC_TIMES_NATIVE_TIME_EQUATIONS_H
#define ISLAMIC_TIMES_NATIVE_TIME_EQUATIONS_H

// Includes
#include "c_datetime.h"
#include "c_calculation_equations.h"

/* Core time-conversion and sidereal helpers (all angles in degrees). */
double greenwich_mean_sidereal_time(double julian_day);
/* Convert Gregorian datetime + UTC offset (hours) to Julian Day. */
double gregorian_to_jd(datetime date, double utc_offset);
/* Convert Julian Day + UTC offset (hours) to datetime; 0 on success. */
int jd_to_gregorian(double jd, double utc_offset, datetime* date);
/* Fractional day from Julian Day in local offset context. */
double fraction_of_day_jd(double jd, double utc_offset);
/* Fractional day from datetime fields. */
double fraction_of_day_datetime(datetime date);
/* Approximate Delta-T (TT-UT) in seconds. */
double delta_t_approx(int year, int month);
/* TT-UTC in seconds for a UTC Julian Day. */
double tt_minus_utc_for_jd(double jd_utc);
/* UT1-UTC in seconds for a UTC Julian Day. */
double ut1_minus_utc_for_jd(double jd_utc);
/* Delta-T (TT-UT1) in seconds for a UTC Julian Day. */
double delta_t_for_jd(double jd_utc);
/* Resolve TT-UTC, UT1-UTC, and Delta-T for a UTC Julian Day. */
void resolve_time_scales_for_jd(double jd_utc, double* tt_minus_utc, double* ut1_minus_utc, double* delta_t);

/* Python wrappers exposed by astro_core. */
PyObject* py_greenwich_mean_sidereal_time(PyObject* self, PyObject* args);
PyObject* py_gregorian_to_jd(PyObject *self, PyObject *args);
PyObject* py_jd_to_gregorian(PyObject *self, PyObject *args);
PyObject* py_delta_t_approx(PyObject *self, PyObject *args);
PyObject* py_resolve_time_scales(PyObject *self, PyObject *args);

#endif  // ISLAMIC_TIMES_NATIVE_TIME_EQUATIONS_H
