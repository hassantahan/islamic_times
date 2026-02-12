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
    {"greenwich_mean_sidereal_time", py_greenwich_mean_sidereal_time, METH_VARARGS,
     "greenwich_mean_sidereal_time(jd) -> float\nReturn GMST in degrees for Julian Day jd."},
    {"jd_to_gregorian", py_jd_to_gregorian, METH_VARARGS,
     "jd_to_gregorian(jd, utc_offset) -> datetime\nConvert Julian Day to Gregorian datetime."},
    {"gregorian_to_jd", py_gregorian_to_jd, METH_VARARGS,
     "gregorian_to_jd(dt, utc_offset) -> float\nConvert Gregorian datetime to Julian Day."},
    {"delta_t_approx", py_delta_t_approx, METH_VARARGS,
     "delta_t_approx(year, month) -> float\nApproximate Delta-T (TT-UT) in seconds."},
    {"compute_sun", (PyCFunction)(PyObject *(*)(PyObject *, PyObject *const *, Py_ssize_t))py_compute_sun, METH_FASTCALL,
     "compute_sun(...) -> Sun\nCompute solar coordinates and derived parameters."},
    {"find_sun_transit", py_find_sun_transit, METH_VARARGS,
     "find_sun_transit(...) -> datetime\nCompute local solar transit time."},
    {"find_proper_suntime", py_find_proper_suntime, METH_VARARGS,
     "find_proper_suntime(...) -> datetime\nCompute local sunrise/sunset for a target angle."},
    {"compute_moon", (PyCFunction)(PyObject *(*)(PyObject *, PyObject *const *, Py_ssize_t))py_compute_moon, METH_FASTCALL,
     "compute_moon(...) -> Moon\nCompute lunar coordinates and derived parameters."},
    {"find_moon_transit", py_find_moon_transit, METH_VARARGS,
     "find_moon_transit(...) -> datetime\nCompute local lunar transit time."},
    {"find_proper_moontime", py_find_proper_moontime, METH_VARARGS,
     "find_proper_moontime(...) -> datetime\nCompute local moonrise/moonset."},
    {"next_phases_of_moon_utc", py_next_phases_of_moon_utc, METH_VARARGS,
     "next_phases_of_moon_utc(dt) -> tuple[datetime, datetime, datetime, datetime]\nReturn next New/First/Full/Last phases in UTC."},
    {"compute_visibilities", py_compute_visibilities, METH_VARARGS,
     "compute_visibilities(...) -> Visibilities\nCompute crescent visibility q-values and classes."},
    {"compute_visibilities_batch", compute_visibilities_batch_py, METH_VARARGS,
     "compute_visibilities_batch(...) -> numpy.ndarray\nBatch crescent visibility computation across coordinate arrays."},
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
