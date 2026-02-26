#define PY_SSIZE_T_CLEAN
#include "astro_core.h"
#include "c_time_equations.h"
#include "c_calculation_equations.h"
#include "c_sun_equations.h"
#include "c_moon_equations.h"
#include "c_visibilities.h"

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
    {NULL, NULL, 0, NULL}
};

static int astro_core_traverse(PyObject* module, visitproc visit, void* arg) {
    AstroCoreState* state = (AstroCoreState*)PyModule_GetState(module);
    if (!state) {
        return 0;
    }

    Py_VISIT(state->sun_type);
    Py_VISIT(state->moon_type);
    Py_VISIT(state->visibilities_type);
    Py_VISIT(state->angle_type);
    Py_VISIT(state->distance_type);
    Py_VISIT(state->distance_units_type);
    Py_VISIT(state->right_ascension_type);
    return 0;
}

static int astro_core_clear(PyObject* module) {
    AstroCoreState* state = (AstroCoreState*)PyModule_GetState(module);
    if (!state) {
        return 0;
    }

    Py_CLEAR(state->sun_type);
    Py_CLEAR(state->moon_type);
    Py_CLEAR(state->visibilities_type);
    Py_CLEAR(state->angle_type);
    Py_CLEAR(state->distance_type);
    Py_CLEAR(state->distance_units_type);
    Py_CLEAR(state->right_ascension_type);
    return 0;
}

static void astro_core_free(void* module) {
    astro_core_clear((PyObject*)module);
}

static int import_cached_type(PyObject* source_module, const char* type_name, PyObject** out_type) {
    PyObject* type_obj = PyObject_GetAttrString(source_module, type_name);
    if (!type_obj) {
        return -1;
    }
    *out_type = type_obj;
    return 0;
}

static int astro_core_load_runtime_types(PyObject* module) {
    PyObject* mod_sun = NULL;
    PyObject* mod_moon = NULL;
    PyObject* mod_dc = NULL;
    AstroCoreState* state = (AstroCoreState*)PyModule_GetState(module);

    if (!state) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to allocate astro_core module state.");
        return -1;
    }

    mod_sun = PyImport_ImportModule("islamic_times.sun_equations");
    mod_moon = PyImport_ImportModule("islamic_times.moon_equations");
    mod_dc = PyImport_ImportModule("islamic_times.it_dataclasses");
    if (!mod_sun || !mod_moon || !mod_dc) {
        goto error;
    }

    if (import_cached_type(mod_sun, "Sun", &state->sun_type) < 0 ||
        import_cached_type(mod_moon, "Moon", &state->moon_type) < 0 ||
        import_cached_type(mod_dc, "Visibilities", &state->visibilities_type) < 0 ||
        import_cached_type(mod_dc, "Angle", &state->angle_type) < 0 ||
        import_cached_type(mod_dc, "Distance", &state->distance_type) < 0 ||
        import_cached_type(mod_dc, "DistanceUnits", &state->distance_units_type) < 0 ||
        import_cached_type(mod_dc, "RightAscension", &state->right_ascension_type) < 0) {
        goto error;
    }

    Py_DECREF(mod_sun);
    Py_DECREF(mod_moon);
    Py_DECREF(mod_dc);
    return 0;

error:
    Py_XDECREF(mod_sun);
    Py_XDECREF(mod_moon);
    Py_XDECREF(mod_dc);
    astro_core_clear(module);
    return -1;
}

/* ============================
   Module Definition
   ============================ */

static struct PyModuleDef astro_core_module = {
    PyModuleDef_HEAD_INIT,
    "astro_core",
    "Core astronomical calculations module.",
    sizeof(AstroCoreState),
    AstroCoreMethods,
    NULL,
    astro_core_traverse,
    astro_core_clear,
    astro_core_free
};

/* ============================
   Module Initialization
   ============================ */

PyMODINIT_FUNC PyInit_astro_core(void) {
    PyObject* module = PyModule_Create(&astro_core_module);
    if (module == NULL) {
        return NULL;
    }

    PyDateTime_IMPORT;
    if (PyDateTimeAPI == NULL) {
        Py_DECREF(module);
        return NULL;
    }

    import_array();
    if (PyArray_API == NULL) {
        Py_DECREF(module);
        return NULL;
    }

    if (astro_core_load_runtime_types(module) < 0) {
        Py_DECREF(module);
        return NULL;
    }

    return module;
}
