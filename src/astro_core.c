#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "_time_equations.h"
#include "_calculation_equations.h"
#include "_sun_equations.h"
#include "astro_core.h"

PyObject *datetime_datetime = NULL;
PyObject *SunType = NULL;
PyObject *AngleType = NULL;
PyObject *DistanceType = NULL;
PyObject *DistanceUnitsType = NULL;
PyObject *RightAscensionType = NULL;

/* ============================
   Cleaning up types
   ============================ */

static void cleanup_types(void) {
    // Optional cleanup at exit
    Py_XDECREF(SunType);
    Py_XDECREF(AngleType);
    Py_XDECREF(DistanceType);
    Py_XDECREF(DistanceUnitsType);
    Py_XDECREF(RightAscensionType);
}

/* ============================
   Python-exposed Methods Table
   ============================ */

static PyMethodDef AstroCoreMethods[] = {
    {"greenwich_mean_sidereal_time", py_greenwich_mean_sidereal_time, METH_VARARGS, "Compute GMST give a Julian Day."},
    {"jd_to_gregorian", py_jd_to_gregorian, METH_VARARGS, "Go from a Julian day to a Gregorian date and time."},
    {"gregorian_to_jd", py_gregorian_to_jd, METH_VARARGS, "Go from a Gregorian date and time to Julian Day."},
    {"compute_sun", (PyCFunction)(void(*)(void))py_compute_sun, METH_FASTCALL, "Compute sun parameters"},
    {NULL, NULL, 0, NULL}  /* Sentinel */
};

/* ============================
   Module Definition
   ============================ */

static struct PyModuleDef astro_core_module = {
    PyModuleDef_HEAD_INIT,
    "astro_core",                              // Module name
    "Core astronomical calculations module.",  // Module docstring
    -1,                                        // Size of per-interpreter state or -1
    AstroCoreMethods
};

/* ============================
   Module Initialization
   ============================ */

PyMODINIT_FUNC PyInit_astro_core(void) {
    PyObject *m = PyModule_Create(&astro_core_module);
    if (m == NULL)
        return NULL;

    PyObject *datetime_module = PyImport_ImportModule("datetime");
    if (!datetime_module)
        return NULL;

    datetime_datetime = PyObject_GetAttrString(datetime_module, "datetime");
    Py_DECREF(datetime_module);
    if (!datetime_datetime)
        return NULL;

    PyObject *mod_sun = PyImport_ImportModule("islamic_times.sun_equations");
    PyObject *mod_dc  = PyImport_ImportModule("islamic_times.dataclasses");

    if (!mod_sun || !mod_dc) return NULL;

    // Clean and safe import
    IMPORT_TYPE(SunType, mod_sun, "Sun");
    IMPORT_TYPE(AngleType, mod_dc, "Angle");
    IMPORT_TYPE(DistanceType, mod_dc, "Distance");
    IMPORT_TYPE(DistanceUnitsType, mod_dc, "DistanceUnits");
    IMPORT_TYPE(RightAscensionType, mod_dc, "RightAscension");

    Py_DECREF(mod_sun);
    Py_DECREF(mod_dc);

    // Register cleanup when interpreter shuts down
    if (Py_AtExit(cleanup_types) != 0) {
        PyErr_WarnEx(PyExc_RuntimeWarning, "Failed to register cleanup", 1);
    }

    return m;
}