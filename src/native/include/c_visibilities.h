#ifndef ISLAMIC_TIMES_NATIVE_VISIBILITIES_H
#define ISLAMIC_TIMES_NATIVE_VISIBILITIES_H

// Includes
#include "c_moon_equations.h"

typedef struct {
    /* Criterion score (or sentinel for no-event conditions). */
    double q_value;
    /* Best-evaluation datetime in observer-local context. */
    datetime best_dt;
    /* Human-readable category label corresponding to q_value. */
    const char* classification;
} VisibilityResult;

/* Python wrappers exposed by astro_core. */
PyObject* py_compute_visibilities(PyObject* self, PyObject* args);
PyObject* compute_visibilities_batch_py(PyObject* self, PyObject* args);

#endif  // ISLAMIC_TIMES_NATIVE_VISIBILITIES_H
