#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <math.h>
#include <string.h>
#include "astro_core.h"

/* ================================
   Definitions & Helper Constants
   ================================ */

#define M_PI 3.14159265358979
#define J2000 2451545.0
#define JULIAN_CENTURY 36525.0
#define JULIAN_MILLENNIUM 365250.0

/* ================================
   GMST Calculation
   ================================ */

double greenwich_mean_sidereal_time(double julian_day){
    double t = (julian_day - J2000) / JULIAN_CENTURY;
    double t2 = t * t;
    double t3 = t2 * t;

    double theta_zero = 280.46061837 + 360.98564736629 * (julian_day - J2000) + 0.000387933 * t2 - t3 / 38710000.0;

    return theta_zero;
}

/* Python wrapper */
PyObject* py_greenwich_mean_sidereal_time(PyObject *self, PyObject *args) {
    double jd;
    if (!PyArg_ParseTuple(args, "d", &jd)) return NULL;
    double result = greenwich_mean_sidereal_time(jd);
    return Py_BuildValue("d", result);
}

/* ================================
   Gregorian 2 JD Calculation
   ================================ */

   double gregorian_to_jd(int year, int month, int day, double fraction_of_day, double zone){
    double y = year;
    double m = month;
    double d = (double)day + fraction_of_day;
    
    if (m <= 2) {
        y--;
        m += 12;
    }

    int a = (int)(floor(y / 100));
    int b = 2 - a + (int)(floor(a / 4));

    double jd = (int)(floor(365.25 * (y + 4716)) + (int)(floor(30.6001 * (m + 1)))) + d + b - 1524.5 - zone / 24;

    return jd;
}

/* Python wrapper */
PyObject* py_gregorian_to_jd(PyObject *self, PyObject *args) {
    PyObject *input_datetime;
    double zone;

    // Parse argument as a PyObject
    if (!PyArg_ParseTuple(args, "O!d", datetime_datetime, &input_datetime, &zone)) {
        return NULL;
    }

    // Extract datetime fields
    int year = PyLong_AsLong(PyObject_GetAttrString(input_datetime, "year"));
    int month = PyLong_AsLong(PyObject_GetAttrString(input_datetime, "month"));
    int day = PyLong_AsLong(PyObject_GetAttrString(input_datetime, "day"));
    int hour = PyLong_AsLong(PyObject_GetAttrString(input_datetime, "hour"));
    int minute = PyLong_AsLong(PyObject_GetAttrString(input_datetime, "minute"));
    int second = PyLong_AsLong(PyObject_GetAttrString(input_datetime, "second"));
    int microsecond = PyLong_AsLong(PyObject_GetAttrString(input_datetime, "microsecond"));

    double fraction_of_day = ((double)hour + (double)minute / 60 + (double)second / 3600 + (double)microsecond / (1000000 * 3600));
    double result = gregorian_to_jd(year, month, day, fraction_of_day, zone);

    return Py_BuildValue("d", result);
}


/* ================================
   JD 2 Gregorian Calculation
   ================================ */

int jd_to_gregorian(double jd, double zone, int* year, int* month, int* day, int* hour, int* minute, int* second, int* microsecond) {
    if (jd > 0) {
        jd += 0.5 - zone / 24;
        
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
            *month = e - 1;
        else
            *month = e - 13;

        if (*month > 2)
            *year = c - 4716;
        else
            *year = c - 4715;

        *day = (int)floor(d_day);
        double f_day = d_day - *day;
        
        double d_hour = f_day * 24;
        *hour = (int)d_hour;
        
        double d_minute = (d_hour - *hour) * 60;
        *minute = (int)d_minute;
        
        double d_second = (d_minute - *minute) * 60;
        *second = (int)d_second;
        
        double d_microsecond = (d_second - *second) * 1000000;
        *microsecond = (int)d_microsecond;

        return 0;
    }
    else {
        return 1;
    }
}

/* Python wrapper */
PyObject* py_jd_to_gregorian(PyObject *self, PyObject *args) {
    double jd;
    double zone; 

    if (!PyArg_ParseTuple(args, "dd", &jd, &zone)) return NULL;

    int year, month, day, hour, minute, second, microsecond;
    int status = jd_to_gregorian(jd, zone, &year, &month, &day, &hour, &minute, &second, &microsecond);
    if (!status) {
        return PyObject_CallFunction(datetime_datetime, "iiiiiii",
                year, month, day, hour, minute, second, microsecond);
    }
    else {
        return PyObject_CallFunction(datetime_datetime, "iiiiiii",
            1, 1, 1, 0, 0, 0, 0);
    } 
}