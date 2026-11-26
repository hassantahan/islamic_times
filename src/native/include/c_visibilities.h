#ifndef ISLAMIC_TIMES_NATIVE_VISIBILITIES_H
#define ISLAMIC_TIMES_NATIVE_VISIBILITIES_H

#include "astro_core.h"

/* Legacy reference-date forecast, anchored to the nearest new moon. */
PyObject* py_compute_visibilities(PyObject* self, PyObject* args);
/* Explicit UTC conjunction and local reference day. Returns sunset, moonset,
 * and best time in the solver's fixed local offset, using None for missing times.
 * The final pair contains the criterion score or no-event code and its label. */
PyObject* py_compute_crescent_evening(PyObject* self, PyObject* args);
/* Batch results are flat and location-major, with consecutive days per point.
 * Longitudes and latitudes are in degrees and UTC offset is in hours. */
PyObject* compute_visibilities_batch_py(PyObject* self, PyObject* args);
PyObject* compute_visibilities_batch_codes_py(PyObject* self, PyObject* args);

#endif  // ISLAMIC_TIMES_NATIVE_VISIBILITIES_H
