#ifndef C_VISIBILITIES_H
#define C_VISIBILITIES_H

// Includes
#include "c_moon_equations.h"

typedef struct {
    double q_value;
    datetime best_dt;
    const char* classification;
} VisibilityResult;

PyObject* py_compute_visibilities(PyObject* self, PyObject* args);
PyObject* compute_visibilities_batch_py(PyObject* self, PyObject* args);

#endif