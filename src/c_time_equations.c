#define PY_SSIZE_T_CLEAN
#include "c_time_equations.h"

/* ================================
   Definitions & Helper Constants
   ================================ */

/* ================================
   GMST Calculation
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

/* Python wrapper */
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
   JD 2 Gregorian Calculation
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

/* Python wrapper */
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
        return NULL;
    }
}

/* ================================
   Functions for geting the fraction of day as float 
   ================================ */

double fraction_of_day_jd(double jd, double utc_offset) {
    return fmod((jd + 0.5 + utc_offset / 24), 1);
}

/* ================================
   DeltaT approximator; from: https://eclipse.gsfc.nasa.gov/LEcat5/deltatpoly.html
   ================================ */

double delta_t_approx(int year, int month) {
    double y = year + (month - 0.5) / 12.0;
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

/* Python wrapper */
PyObject* py_delta_t_approx(PyObject *self, PyObject *args) {
    int year, month; 

    if (!PyArg_ParseTuple(args, "ii", &year, &month)) return NULL;

    double result = delta_t_approx(year, month);

    return Py_BuildValue("d", result);
}