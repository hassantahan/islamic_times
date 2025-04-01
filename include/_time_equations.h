#ifndef TIME_EQUATIONS_H
#define TIME_EQUATIONS_H

double greenwich_mean_sidereal_time(double julian_day);
double gregorian_to_jd(int year, int month, int day, double fraction_of_day, double zone);
int jd_to_gregorian(double jd, double zone, int* year, int* month, int* day, int* hour, int* minute, int* second, int* microsecond);

PyObject* py_greenwich_mean_sidereal_time(PyObject* self, PyObject* args);
PyObject* py_gregorian_to_jd(PyObject *self, PyObject *args);
PyObject* py_jd_to_gregorian(PyObject *self, PyObject *args);

#endif // TIME_EQUATIONS_H