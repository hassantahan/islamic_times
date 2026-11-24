#define PY_SSIZE_T_CLEAN
#include <limits.h>
#include <stdint.h>
#include <numpy/arrayobject.h>
#include "c_moon_equations.h"
#include "c_visibilities.h"

#define MAX_VISIBILITY_DAYS 366

typedef struct {
    double q_value;    /* Criterion score, or a legacy no-event code. */
    datetime best_dt;  /* Observer-local solver time. */
    datetime sunset;
    datetime moonset;
} VisibilityResult;

/* Shared visibility-score kernels used by Odeh/Yallop/Shaukat classification paths. */

double yallop_odeh(double sun_az, double sun_alt, double moon_az, double moon_alt, double moon_pi, int criterion) {
    double daz_deg = fabs(sun_az - moon_az);
    double arcv_deg = fabs(sun_alt - moon_alt);
    double arcl_deg = DEGREES(acos(cos(RADIANS(arcv_deg)) * cos(RADIANS(daz_deg))));

    double semi_diameter = 0.27245 * moon_pi * 60.0;
    double semi_diameter_prime = semi_diameter * (1 + sin(RADIANS(moon_alt)) * sin(RADIANS(moon_pi)));

    double w_prime = semi_diameter_prime * (1 - cos(RADIANS(arcl_deg)));

    if (criterion)
        return (arcv_deg - (11.8371 - 6.3226 * w_prime + 0.7319 * pow(w_prime, 2) - 0.1018 * pow(w_prime, 3))) / 10.0;
    
    return arcv_deg - (-0.1018 * pow(w_prime, 3) + 0.7319 * pow(w_prime, 2) - 6.3226 * w_prime + 7.1651);
}

const char* classify_visibility(double q, int criterion) {
    // Special moon/sun conditions
    if (q == -999.0) 
        return "Moonset before the new moon.";
    else if (q == -998.0) 
        return "Moonset before sunset.";
    else if (q == -997.0) 
        return "Moonset & Sunset don't exist.";
    else if (q == -996.0) 
        return "Sunset doesn't exist.";
    else if (q == -995.0) 
        return "Moonset doesn't exist.";
    else if (q == -994.0)
        return "The nominal best time is before conjunction.";

    // Visibility classification
    if (criterion == 0) {
        if (q >= 5.65)
            return "A: Crescent is visible by naked eyes.";
        else if (q >= 2.0)
            return "B: Crescent is visible by optical aid, and it could be seen by naked eyes.";
        else if (q >= -0.96)
            return "C: Crescent is visible by optical aid only.";
        else
            return "D: Crescent is not visible even by optical aid.";
    } else if (criterion == 1) {
        if (q > 0.216)
            return "A: Easily visible.";
        else if (q > -0.014)
            return "B: Visible under perfect conditions.";
        else if (q > -0.160)
            return "C: May need optical aid to find crescent.";
        else if (q > -0.232)
            return "D: Will need optical aid to find crescent.";
        else if (q > -0.293)
            return "E: Not visible with a [conventional] telescope.";
        else
            return "F: Not visible; below the Danjon limit.";
    }

    if (q > 0.27)
        return "A: Easily visible";
    else if (q > -0.024)
        return "B: Visible under perfect conditions.";
    else if (q > -0.212)
        return "C: Optical aid needed to find the moon.";
    else if (q > -0.48)
        return "D: Visible with optical aid only.";
    else
        return "F: Not visible.";
}

/* ================================
   Crescent visibility solver
   Computes q-values and classes for
   consecutive days after new moon.
   ================================ */

/* Evaluate one local evening. The legacy mapper gives the conjunction-day
 * moonset checks precedence; the detailed API also records the solar event. */
static void evaluate_visibility(
    datetime conjunction, double jd, double offset, double lat, double lon, double elev, double temp, double press,
    int criterion, int conjunction_day, int detailed, VisibilityResult* result
) {
    *result = (VisibilityResult){
        .best_dt = conjunction, .sunset = INVALID_DATETIME, .moonset = INVALID_DATETIME
    };
    if (conjunction_day) {
        result->moonset = find_proper_moontime(jd, offset, lat, lon, elev, temp, press, 's', NULL, NULL);
        if (!detailed && result->moonset.year == INVALID_DATETIME.year) {
            result->q_value = -995.0;
            return;
        }
        if (!detailed && compare_datetime(&result->moonset, &conjunction) < 0) {
            result->best_dt = result->moonset;
            result->q_value = -999.0;
            return;
        }
    }

    double nutation[3], obliquity[3];
    result->sunset = find_proper_suntime_w_nutation(
        jd, offset, lat, lon, elev, temp, press, STANDARD_SET_AND_RISE, 's', nutation, obliquity
    );
    int has_sunset = result->sunset.year != INVALID_DATETIME.year;
    if (!conjunction_day) {
        result->moonset = find_proper_moontime(
            jd, offset, lat, lon, elev, temp, press, 's',
            has_sunset ? nutation : NULL, has_sunset ? obliquity : NULL
        );
    }

    int has_moonset = result->moonset.year != INVALID_DATETIME.year;
    if (!has_sunset && !has_moonset) {
        result->q_value = -997.0;
        return;
    }
    if (!has_sunset) {
        result->best_dt = result->moonset;
        result->q_value = -996.0;
        return;
    }
    if (!has_moonset) {
        result->best_dt = result->sunset;
        result->q_value = -995.0;
        return;
    }
    if (conjunction_day && compare_datetime(&result->moonset, &conjunction) < 0) {
        result->best_dt = result->moonset;
        result->q_value = -999.0;
        return;
    }
    if (compare_datetime(&result->moonset, &result->sunset) < 0) {
        result->best_dt = result->moonset;
        result->q_value = -998.0;
        return;
    }

    /* Evaluate at Yallop's best time, four ninths of the sunset-to-moonset lag. */
    double sunset_jd = gregorian_to_jd(result->sunset, offset);
    double moonset_jd = gregorian_to_jd(result->moonset, offset);
    double best_jd = sunset_jd + 4.0 / 9.0 * (moonset_jd - sunset_jd);
    if (detailed && conjunction_day && best_jd < gregorian_to_jd(conjunction, offset)) {
        result->q_value = -994.0;
        return;
    }
    double tt_utc, ut1_utc;
    resolve_time_scales_for_jd(best_jd, &tt_utc, &ut1_utc, NULL);
    double jde = best_jd + tt_utc / SECONDS_IN_DAY;
    double jd_ut1 = best_jd + ut1_utc / SECONDS_IN_DAY;
    SunResult sun;
    MoonResult moon;
    compute_sun_result(jde, jd_ut1, lat, lon, elev, temp, press, &sun);
    compute_moon_result(jde, jd_ut1, lat, lon, elev, temp, press, sun.nutation_longitude, sun.true_obliquity, &moon);

    if (criterion == 0) {
        /* Odeh uses apparent topocentric horizontal coordinates. */
        result->q_value = yallop_odeh(
            sun.true_azimuth, sun.apparent_altitude, moon.true_azimuth, moon.apparent_altitude, moon.eh_parallax, 0
        );
    } else {
        /* Yallop and Shaukat share the same geocentric score. */
        double sun_alt, sun_az, moon_alt, moon_az;
        geocentric_horizontal_coordinates(lat, sun.apparent_declination, sun.local_hour_angle, &sun_alt, &sun_az);
        geocentric_horizontal_coordinates(lat, moon.declination, moon.local_hour_angle, &moon_alt, &moon_az);
        result->q_value = yallop_odeh(sun_az, sun_alt, moon_az, moon_alt, moon.eh_parallax, 1);
    }
    jd_to_gregorian(best_jd, offset, &result->best_dt);
}

static void compute_visibilities(
    datetime conjunction, double jd, double offset, double lat, double lon, double elev, double temp, double press,
    VisibilityResult* results, int days, int criterion
) {
    for (int i = 0; i < days; i++) {
        evaluate_visibility(conjunction, jd + i, offset, lat, lon, elev, temp, press, criterion, i == 0, 0, &results[i]);
    }
}

/* Compact category-code mapper used by batch code API.
 * Codes follow the ordered palette used in mapper:
 *  0:-999, 1:-998, 2:-997, 3:-996, 4:-995,
 *  criterion 0: 5:D, 6:C, 7:B, 8:A
 *  criterion 1: 5:F, 6:E, 7:D, 8:C, 9:B, 10:A
 *  criterion 2: 5:F, 6:D, 7:C, 8:B, 9:A
 */
static uint8_t classify_visibility_code(double q, int criterion) {
    if (q == -999.0) return 0;
    if (q == -998.0) return 1;
    if (q == -997.0) return 2;
    if (q == -996.0) return 3;
    if (q == -995.0) return 4;

    if (criterion == 0) {
        if (q >= 5.65) return 8;
        if (q >= 2.0) return 7;
        if (q >= -0.96) return 6;
        return 5;
    }

    if (criterion == 1) {
        if (q > 0.216) return 10;
        if (q > -0.014) return 9;
        if (q > -0.160) return 8;
        if (q > -0.232) return 7;
        if (q > -0.293) return 6;
        return 5;
    }

    if (q > 0.27) return 9;
    if (q > -0.024) return 8;
    if (q > -0.212) return 7;
    if (q > -0.48) return 6;
    return 5;
}

/* Python wrapper: returns Visibilities dataclass instance. */
PyObject* py_compute_visibilities(PyObject* self, PyObject* args) {
    ENSURE_PYDATETIME();
    PyObject *input_datetime;
    double utc_offset, lat, lon, elev, temp, press;
    int days, criterion;

    if (!PyArg_ParseTuple(args, "O!ddddddii", PyDateTimeAPI->DateTimeType, &input_datetime, &utc_offset, &lat, &lon,
                                                                        &elev, &temp, &press, &days, &criterion))
        return NULL;

    if (days <= 0 || days > MAX_VISIBILITY_DAYS) {
        PyErr_Format(PyExc_ValueError, "'days' must be in the range [1, %d].", MAX_VISIBILITY_DAYS);
        return NULL;
    }
    if (criterion < 0 || criterion > 2) {
        PyErr_SetString(PyExc_ValueError, "'criterion' must be 0 (Odeh), 1 (Yallop), or 2 (Shaukat).");
        return NULL;
    }

    // Parse observer datetime.
    datetime date;
    fill_in_datetime_values(&date, input_datetime);

    // Use nearest new moon as the visibility baseline.
    datetime moon_phases[4];
    next_phases_of_moon_utc(date, moon_phases);

    // New moon is phase index 0.
    datetime conjunction = add_days(moon_phases[0], utc_offset / 24.0);
    double jd = gregorian_to_jd(conjunction, utc_offset);
    

    // Allocate native result buffer and run computation.
    VisibilityResult* c_results = malloc(sizeof(VisibilityResult) * (size_t)days);
    if (!c_results) {
        return PyErr_NoMemory();
    }
    compute_visibilities(conjunction, jd, utc_offset, lat, lon, elev, temp, press,
                        c_results, days, criterion);

    // Prepare Python sequences.
    PyObject* py_dates = PyTuple_New(days);
    PyObject* py_q_values = PyTuple_New(days);
    PyObject* py_classifications = PyTuple_New(days);
    if (!py_dates || !py_q_values || !py_classifications) {
        Py_XDECREF(py_dates);
        Py_XDECREF(py_q_values);
        Py_XDECREF(py_classifications);
        free(c_results);
        return NULL;
    }

    for (int i = 0; i < days; i++) {
        PyObject* py_dt = datetime_to_pydatetime(c_results[i].best_dt);
        PyObject* py_q = PyFloat_FromDouble(c_results[i].q_value);
        PyObject* py_class = PyUnicode_FromString(classify_visibility(c_results[i].q_value, criterion));

        if (!py_dt || !py_q || !py_class) {
            Py_XDECREF(py_dt);
            Py_XDECREF(py_q);
            Py_XDECREF(py_class);
            Py_DECREF(py_dates);
            Py_DECREF(py_q_values);
            Py_DECREF(py_classifications);
            free(c_results);
            return NULL;
        }

        PyTuple_SET_ITEM(py_dates, i, py_dt);
        PyTuple_SET_ITEM(py_q_values, i, py_q);
        PyTuple_SET_ITEM(py_classifications, i, py_class);
    }

    free(c_results);

    // Build Visibilities dataclass instance.
    const char* criterion_name = "Shaukat";
    if (criterion == 0) {
        criterion_name = "Odeh";
    } else if (criterion == 1) {
        criterion_name = "Yallop";
    }
    PyObject* py_criterion = PyUnicode_FromString(criterion_name);
    if (!py_criterion) {
        Py_DECREF(py_dates);
        Py_DECREF(py_q_values);
        Py_DECREF(py_classifications);
        return NULL;
    }

    PyObject* args_tuple = Py_BuildValue("OOOO", py_criterion, py_dates, py_q_values, py_classifications);
    Py_DECREF(py_criterion);  // Py_BuildValue increases refcount
    Py_DECREF(py_dates);
    Py_DECREF(py_q_values);
    Py_DECREF(py_classifications);
    if (!args_tuple) return NULL;

    AstroCoreState* state = (AstroCoreState*)PyModule_GetState(self);
    if (!state) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to access astro_core module state.");
        Py_DECREF(args_tuple);
        return NULL;
    }
    if (!state->visibilities_type || !PyType_Check(state->visibilities_type)) {
        Py_DECREF(args_tuple);
        PyErr_SetString(PyExc_RuntimeError, "VisibilitiesType is not a valid Python type");
        return NULL;
    }

    PyObject* result = PyObject_CallObject(state->visibilities_type, args_tuple);
    Py_DECREF(args_tuple);

    return result;
}

static PyObject* visibility_time_or_none(datetime value) {
    if (value.year == INVALID_DATETIME.year) {
        Py_RETURN_NONE;
    }
    return datetime_to_pydatetime(value);
}

/* Compute metadata for one explicit observer-local date and UTC conjunction. */
PyObject* py_compute_crescent_evening(PyObject* self, PyObject* args) {
    ENSURE_PYDATETIME();
    PyObject* conjunction_obj;
    PyObject* reference_obj;
    double offset, lat, lon, elev, temp, press;
    int criterion, conjunction_day;
    if (!PyArg_ParseTuple(
        args, "O!O!ddddddip", PyDateTimeAPI->DateTimeType, &conjunction_obj,
        PyDateTimeAPI->DateTimeType, &reference_obj, &offset, &lat, &lon, &elev, &temp, &press,
        &criterion, &conjunction_day
    )) {
        return NULL;
    }
    if (criterion < 0 || criterion > 2) {
        PyErr_SetString(PyExc_ValueError, "'criterion' must be 0 (Odeh), 1 (Yallop), or 2 (Shaukat).");
        return NULL;
    }

    datetime conjunction_utc, reference;
    fill_in_datetime_values(&conjunction_utc, conjunction_obj);
    fill_in_datetime_values(&reference, reference_obj);
    datetime conjunction = add_days(conjunction_utc, offset / 24.0);
    double jd = gregorian_to_jd(reference, offset);
    VisibilityResult result;
    evaluate_visibility(
        conjunction, jd, offset, lat, lon, elev, temp, press, criterion, conjunction_day, 1, &result
    );

    PyObject* sunset = visibility_time_or_none(result.sunset);
    PyObject* moonset = visibility_time_or_none(result.moonset);
    PyObject* best_time = visibility_time_or_none(result.q_value <= -994.0 ? INVALID_DATETIME : result.best_dt);
    if (!sunset || !moonset || !best_time) {
        Py_XDECREF(sunset);
        Py_XDECREF(moonset);
        Py_XDECREF(best_time);
        return NULL;
    }
    return Py_BuildValue(
        "NNNds", sunset, moonset, best_time, result.q_value, classify_visibility(result.q_value, criterion)
    );
}

static PyArrayObject* new_visibility_array(npy_intp size, char type) {
    npy_intp dims[] = {size};
    if (type == 'r' || type == 'i') {
        return (PyArrayObject*)PyArray_SimpleNew(1, dims, type == 'r' ? NPY_DOUBLE : NPY_UINT8);
    }
    PyObject* text = PyUnicode_FromString("U100");
    if (!text) return NULL;

    PyArray_Descr* dtype = NULL;
    int converted = PyArray_DescrConverter(text, &dtype);
    Py_DECREF(text);
    if (!converted) return NULL;
    return (PyArrayObject*)PyArray_NewFromDescr(&PyArray_Type, dtype, 1, dims, NULL, NULL, 0, NULL);
}

static PyObject* compute_visibilities_batch_impl(
    PyObject* lats_obj, PyObject* lons_obj, PyObject* new_moon_obj,
    int days, int criterion, double offset, double elev, double temp, double press, char type
) {
    if (days < 1 || days > MAX_VISIBILITY_DAYS) {
        PyErr_Format(PyExc_ValueError, "'days' must be in the range [1, %d].", MAX_VISIBILITY_DAYS);
        return NULL;
    }
    if (criterion < 0 || criterion > 2) {
        PyErr_SetString(PyExc_ValueError, "'criterion' must be 0 (Odeh), 1 (Yallop), or 2 (Shaukat).");
        return NULL;
    }
    if (type != 'r' && type != 'c' && type != 'i') {
        PyErr_SetString(PyExc_ValueError, "Unknown type flag (expected 'r', 'c', or 'i')");
        return NULL;
    }
    if (!PyArray_Check(lats_obj) || !PyArray_Check(lons_obj)) {
        PyErr_SetString(PyExc_TypeError, "lat/lon must be NumPy arrays");
        return NULL;
    }
    PyArrayObject* lats_arr = (PyArrayObject*)lats_obj;
    PyArrayObject* lons_arr = (PyArrayObject*)lons_obj;
    if (PyArray_TYPE(lats_arr) != NPY_DOUBLE || !PyArray_ISCARRAY_RO(lats_arr) ||
        PyArray_TYPE(lons_arr) != NPY_DOUBLE || !PyArray_ISCARRAY_RO(lons_arr)) {
        PyErr_SetString(PyExc_ValueError, "lat/lon arrays must be float64 and C-contiguous");
        return NULL;
    }
    npy_intp count = PyArray_SIZE(lats_arr);
    if (count != PyArray_SIZE(lons_arr)) {
        PyErr_SetString(PyExc_ValueError, "Latitude and longitude arrays must have same size");
        return NULL;
    }
    if (count > INT_MAX / days) {
        PyErr_SetString(PyExc_OverflowError, "Requested batch size exceeds native loop limits.");
        return NULL;
    }

    npy_intp size = count * days;
    PyArrayObject* output = new_visibility_array(size, type);
    if (!output || size == 0) return (PyObject*)output;

    uint8_t* codes = type == 'i' ? (uint8_t*)PyArray_DATA(output) : NULL;
    double* scores = NULL;
    if (type == 'r') {
        scores = (double*)PyArray_DATA(output);
    } else if (type == 'c') {
        scores = PyMem_New(double, size);
        if (!scores) {
            Py_DECREF(output);
            return PyErr_NoMemory();
        }
    }

    datetime conjunction_utc;
    fill_in_datetime_values(&conjunction_utc, new_moon_obj);
    datetime conjunction = add_days(conjunction_utc, offset / 24.0);
    double jd = gregorian_to_jd(conjunction, offset);
    const double* lats = (const double*)PyArray_DATA(lats_arr);
    const double* lons = (const double*)PyArray_DATA(lons_arr);
    VisibilityResult location[MAX_VISIBILITY_DAYS];

    /* Stream one location into the NumPy output at a time. */
    Py_BEGIN_ALLOW_THREADS
    for (npy_intp i = 0; i < count; i++) {
        compute_visibilities(conjunction, jd, offset, lats[i], lons[i], elev, temp, press, location, days, criterion);
        for (int day = 0; day < days; day++) {
            npy_intp index = i * days + day;
            if (codes) {
                codes[index] = classify_visibility_code(location[day].q_value, criterion);
            } else {
                scores[index] = location[day].q_value;
            }
        }
    }
    Py_END_ALLOW_THREADS

    if (type == 'c') {
        for (npy_intp i = 0; i < size; i++) {
            PyObject* label = PyUnicode_FromString(classify_visibility(scores[i], criterion));
            int status = label ? PyArray_SETITEM(output, PyArray_GETPTR1(output, i), label) : -1;
            Py_XDECREF(label);
            if (status < 0) {
                PyMem_Free(scores);
                Py_DECREF(output);
                return NULL;
            }
        }
        PyMem_Free(scores);
    }
    return (PyObject*)output;
}

/* Python wrapper: returns either raw q-values or category labels as NumPy array. */
PyObject* compute_visibilities_batch_py(PyObject* self, PyObject* args) {
    ENSURE_PYDATETIME();
    ENSURE_NUMPY();

    PyObject *lats_obj, *lons_obj, *new_moon_obj;
    int days, criterion, type_code;
    double utc_offset, elev, temp, press;
    if (!PyArg_ParseTuple(args,
            "OOO!iiddddC",
            &lats_obj, &lons_obj, PyDateTimeAPI->DateTimeType, &new_moon_obj,
            &days, &criterion,
            &utc_offset, &elev, &temp, &press,
            &type_code)) {
        return NULL;
    }

    return compute_visibilities_batch_impl(
        lats_obj,
        lons_obj,
        new_moon_obj,
        days,
        criterion,
        utc_offset,
        elev,
        temp,
        press,
        (char)type_code
    );
}

/* Python wrapper: returns category codes as uint8 NumPy array. */
PyObject* compute_visibilities_batch_codes_py(PyObject* self, PyObject* args) {
    ENSURE_PYDATETIME();
    ENSURE_NUMPY();

    PyObject *lats_obj, *lons_obj, *new_moon_obj;
    int days, criterion;
    double utc_offset, elev, temp, press;
    if (!PyArg_ParseTuple(args,
            "OOO!iidddd",
            &lats_obj, &lons_obj, PyDateTimeAPI->DateTimeType, &new_moon_obj,
            &days, &criterion,
            &utc_offset, &elev, &temp, &press)) {
        return NULL;
    }

    return compute_visibilities_batch_impl(
        lats_obj,
        lons_obj,
        new_moon_obj,
        days,
        criterion,
        utc_offset,
        elev,
        temp,
        press,
        'i'
    );
}
