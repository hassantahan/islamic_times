// astro_types.h

#ifndef ASTRO_TYPES_H
#define ASTRO_TYPES_H

#include <Python.h>

#define ANGLE(value) PyObject_CallFunctionObjArgs(AngleType, PyFloat_FromDouble(value), NULL)
#define DIST(value, unit) PyObject_CallFunctionObjArgs(DistanceType, PyFloat_FromDouble(value), unit, NULL)
#define RA(value) PyObject_CallFunctionObjArgs(RightAscensionType, PyFloat_FromDouble(value / 15.0), NULL)
#define FLOAT(value) PyFloat_FromDouble(value)
#define UNIT(name) PyObject_GetAttrString(DistanceUnitsType, (name))

extern PyObject *datetime_datetime;
extern PyObject *SunType;
extern PyObject *AngleType;
extern PyObject *DistanceType;
extern PyObject *DistanceUnitsType;
extern PyObject *RightAscensionType;

// Macro for safe type import + incref
#define IMPORT_TYPE(var, mod, name)                                      \
    do {                                                                 \
        var = PyObject_GetAttrString((mod), (name));                      \
        if (!var) return NULL;                                            \
        Py_INCREF(var);                                                   \
    } while (0)

#endif
