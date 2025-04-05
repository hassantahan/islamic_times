// astro_types.h
#ifndef ASTRO_TYPES_H
#define ASTRO_TYPES_H

#include <Python.h>
#include <datetime.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define J2000 2451545.0
#define JULIAN_CENTURY 36525.0
#define JULIAN_MILLENNIUM 365250.0
#define SECONDS_IN_DAY 86400.0
#define TROPICAL_YEAR 365.24219878

#define RADIANS(value) ((value) * (M_PI) / 180.0)
#define DEGREES(value) ((value) * (180.0) / M_PI)

#define ANGLE(value) PyObject_CallFunctionObjArgs(AngleType, PyFloat_FromDouble(value), NULL)
#define DIST(value, unit) PyObject_CallFunctionObjArgs(DistanceType, PyFloat_FromDouble(value), unit, NULL)
#define RA(value) PyObject_CallFunctionObjArgs(RightAscensionType, PyFloat_FromDouble(value / 15.0), NULL)
#define FLOAT(value) PyFloat_FromDouble(value)
#define UNIT(name) PyObject_GetAttrString(DistanceUnitsType, (name))

// Macro for safe type import + incref
#define IMPORT_TYPE(var, mod, name)                                       \
    do {                                                                  \
        var = PyObject_GetAttrString((mod), (name));                      \
        if (!var) return NULL;                                            \
        Py_INCREF(var);                                                   \
    } while (0)

#define ENSURE_PYDATETIME() if (!PyDateTimeAPI) { PyDateTime_IMPORT; }


extern PyObject *datetime_datetime;
extern PyObject *SunType;
extern PyObject *MoonType;
extern PyObject *AngleType;
extern PyObject *DistanceType;
extern PyObject *DistanceUnitsType;
extern PyObject *RightAscensionType;


/* ================================
   
   ================================ */

#endif
