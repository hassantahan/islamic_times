#ifndef TIME_EQUATIONS_H
#define TIME_EQUATIONS_H

// Includes
#include "c_datetime.h"
#include "c_calculation_equations.h"


double greenwich_mean_sidereal_time(double julian_day);
double gregorian_to_jd(datetime date, double utc_offset);
int jd_to_gregorian(double jd, double utc_offset, datetime* date);
double fraction_of_day_jd(double jd, double utc_offset);
double fraction_of_day_datetime(datetime date);
double delta_t_approx(int year, int month);

PyObject* py_greenwich_mean_sidereal_time(PyObject* self, PyObject* args);
PyObject* py_gregorian_to_jd(PyObject *self, PyObject *args);
PyObject* py_jd_to_gregorian(PyObject *self, PyObject *args);
PyObject* py_delta_t_approx(PyObject *self, PyObject *args);

#endif // TIME_EQUATIONS_H