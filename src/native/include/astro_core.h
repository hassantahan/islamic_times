// Shared native extension types/constants for islamic_times.astro_core.
#ifndef ISLAMIC_TIMES_NATIVE_ASTRO_CORE_H
#define ISLAMIC_TIMES_NATIVE_ASTRO_CORE_H

#include <Python.h>
#include <datetime.h>
#include <numpy/arrayobject.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define J2000 2451545.0
#define JULIAN_CENTURY 36525.0
#define JULIAN_MILLENNIUM 365250.0
#define SECONDS_IN_DAY 86400.0
#define TROPICAL_YEAR 365.24219878

#define STANDARD_SET_AND_RISE 0.8333333

#define ASTRONOMICAL_UNIT_KM 149597870.7
#define EARTH_RADIUS_KM 6378.14

#define RADIANS(value) ((value) * (M_PI) / 180.0)
#define DEGREES(value) ((value) * (180.0) / M_PI)

#define FLOAT(value) PyFloat_FromDouble(value)

// Safely import a type/class object from a Python module.
// PyObject_GetAttrString already returns a new reference.
#define IMPORT_TYPE(var, mod, name)                                       \
    do {                                                                  \
        var = PyObject_GetAttrString((mod), (name));                      \
        if (!var) return NULL;                                            \
    } while (0)

#define ENSURE_PYDATETIME() if (!PyDateTimeAPI) { PyDateTime_IMPORT; }
#define ENSURE_NUMPY() do { \
    if (PyArray_API == NULL) { \
        import_array(); \
    } \
} while(0)

typedef struct {
    PyObject* sun_type;
    PyObject* moon_type;
    PyObject* visibilities_type;
    PyObject* angle_type;
    PyObject* distance_type;
    PyObject* distance_units_type;
    PyObject* right_ascension_type;
} AstroCoreState;

#endif  // ISLAMIC_TIMES_NATIVE_ASTRO_CORE_H
