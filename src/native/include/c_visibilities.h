#ifndef ISLAMIC_TIMES_NATIVE_VISIBILITIES_H
#define ISLAMIC_TIMES_NATIVE_VISIBILITIES_H

#include "astro_core.h"

/* Python wrappers exposed by astro_core. */
PyObject* py_compute_visibilities(PyObject* self, PyObject* args);
PyObject* py_compute_crescent_evening(PyObject* self, PyObject* args);
PyObject* compute_visibilities_batch_py(PyObject* self, PyObject* args);
PyObject* compute_visibilities_batch_codes_py(PyObject* self, PyObject* args);

#endif  // ISLAMIC_TIMES_NATIVE_VISIBILITIES_H
