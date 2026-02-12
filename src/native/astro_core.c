#define PY_SSIZE_T_CLEAN
#include "astro_core.h"
#include "c_time_equations.h"
#include "c_calculation_equations.h"
#include "c_sun_equations.h"
#include "c_moon_equations.h"
#include "c_visibilities.h"

PyObject *SunType = NULL;
PyObject *MoonType = NULL;
PyObject *VisibilitiesType = NULL;
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
    SunType = NULL;
    Py_XDECREF(MoonType);
    MoonType = NULL;
    Py_XDECREF(VisibilitiesType);
    VisibilitiesType = NULL;
    Py_XDECREF(AngleType);
    AngleType = NULL;
    Py_XDECREF(DistanceType);
    DistanceType = NULL;
    Py_XDECREF(DistanceUnitsType);
    DistanceUnitsType = NULL;
    Py_XDECREF(RightAscensionType);
    RightAscensionType = NULL;
}

/* ============================
   Python-exposed Methods Table
   ============================ */

static PyMethodDef AstroCoreMethods[] = {
    {"greenwich_mean_sidereal_time", py_greenwich_mean_sidereal_time, METH_VARARGS, "Compute GMST give a Julian Day."},
    {"jd_to_gregorian", py_jd_to_gregorian, METH_VARARGS, "Find the Gregorian date and time given a Julian day."},
    {"gregorian_to_jd", py_gregorian_to_jd, METH_VARARGS, "Find the Julian Day given a Gregorian date and time."},
    {"delta_t_approx", py_delta_t_approx, METH_VARARGS, "Find the approximate ΔT given a year and month."},
    {"compute_sun", (PyCFunction)(PyObject *(*)(PyObject *, PyObject *const *, Py_ssize_t))py_compute_sun, METH_FASTCALL, "Compute the sun's position and parameters."},
    {"find_sun_transit", py_find_sun_transit, METH_VARARGS, "Compute the time of the transit of the sun for the given date."},
    {"find_proper_suntime", py_find_proper_suntime, METH_VARARGS, "Compute sunrise or sunset for the given date."},
    {"compute_moon", (PyCFunction)(PyObject *(*)(PyObject *, PyObject *const *, Py_ssize_t))py_compute_moon, METH_FASTCALL, "Compute the moon's position and parameters."},
    {"find_moon_transit", py_find_moon_transit, METH_VARARGS, "Compute the time of the transit of the moon for the given date."},
    {"find_proper_moontime", py_find_proper_moontime, METH_VARARGS, "Compute moonrise or moonset for the given date."},
    {"next_phases_of_moon_utc", py_next_phases_of_moon_utc, METH_VARARGS, "Compute the nearest phases of the moon."},
    {"compute_visibilities", py_compute_visibilities, METH_VARARGS, "Compute the new moon crescent visibility for a given number of days according to a specified critierion."},
    {"compute_visibilities_batch", compute_visibilities_batch_py, METH_VARARGS, "Batch computation of new moon visibilities."},
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
    PyObject *mod_sun = NULL;
    PyObject *mod_moon = NULL;
    PyObject *mod_dc = NULL;

    if (m == NULL) {
        return NULL;
    }

    // Import and initialize datetime C-API
    PyDateTime_IMPORT;
    if (PyDateTimeAPI == NULL) {
        Py_DECREF(m);
        return NULL;
    }

    import_array();
    if (PyArray_API == NULL) {
        Py_DECREF(m);
        return NULL;
    }

    // islamic_times module imports
    mod_sun = PyImport_ImportModule("islamic_times.sun_equations");
    mod_moon = PyImport_ImportModule("islamic_times.moon_equations");
    mod_dc = PyImport_ImportModule("islamic_times.it_dataclasses");
    if (!mod_sun || !mod_moon || !mod_dc) {
        goto error;
    }

    // Class/type imports
    SunType = PyObject_GetAttrString(mod_sun, "Sun");
    MoonType = PyObject_GetAttrString(mod_moon, "Moon");
    VisibilitiesType = PyObject_GetAttrString(mod_dc, "Visibilities");
    AngleType = PyObject_GetAttrString(mod_dc, "Angle");
    DistanceType = PyObject_GetAttrString(mod_dc, "Distance");
    DistanceUnitsType = PyObject_GetAttrString(mod_dc, "DistanceUnits");
    RightAscensionType = PyObject_GetAttrString(mod_dc, "RightAscension");
    if (!SunType || !MoonType || !VisibilitiesType || !AngleType || !DistanceType || !DistanceUnitsType || !RightAscensionType) {
        goto error;
    }

    Py_DECREF(mod_sun);
    Py_DECREF(mod_moon);
    Py_DECREF(mod_dc);

    // Register cleanup when interpreter shuts down
    if (Py_AtExit(cleanup_types) != 0) {
        PyErr_WarnEx(PyExc_RuntimeWarning, "Failed to register cleanup", 1);
    }

    return m;

error:
    Py_XDECREF(mod_sun);
    Py_XDECREF(mod_moon);
    Py_XDECREF(mod_dc);
    cleanup_types();
    Py_DECREF(m);
    return NULL;
}
