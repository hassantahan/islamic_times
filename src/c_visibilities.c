#define PY_SSIZE_T_CLEAN
#include "c_visibilities.h"

/* Different Visibility Criteria */

double yallop_odeh(double sun_az, double sun_alt, double moon_az, double moon_alt, double moon_pi, int criterion) {
    double daz_deg = fabs(sun_az - moon_az);
    double arcv_deg = fabs(sun_alt - moon_alt);
    double arcl_deg = DEGREES(acos(cos(RADIANS(arcv_deg)) * cos(RADIANS(daz_deg))));

    double semi_diameter = 0.27245 * moon_pi * 60.0;
    double semi_diameter_prime = semi_diameter * (1 + sin(RADIANS(moon_alt)) * sin(RADIANS(moon_pi)));

    double w_prime = semi_diameter_prime * (1 - cos(RADIANS(arcl_deg)));

    if (criterion)
        return (arcv_deg - (11.8371 - 6.3226 * w_prime + 0.7319 * pow(w_prime, 2) - 0.1018 * pow(w_prime, 3))) / 10;
    
    return arcv_deg - (-0.1018 * pow(w_prime, 3) + 0.7319 * pow(w_prime, 2) - 6.3226 * w_prime + 7.1651);
}

double shaukat() {
    return 0.0;
}

/* ================================
    Calculate New Moon Crescent 
    Visibilities According to Different
    Criteria 
    ================================ */

void compute_visibilities(datetime date, double utc_offset, double lat, double lon, double elev, double temp, double press,
    double* q_values, datetime* best_dts, int days, int criterion) {

    // Get New Moon Date from moon_phases list
    datetime moon_phases[4];
    next_phases_of_moon_utc(date, moon_phases);

    // New Moon is the first element
    datetime new_moon_dt = moon_phases[0];

    // Get the JD and ΔT for the new moon
    double jd_new_moon = gregorian_to_jd(new_moon_dt, utc_offset);
    double deltaT_new_moon = delta_t_approx(new_moon_dt.year, new_moon_dt.month);

    // Initialize best_jds array
    double* best_jds = malloc(sizeof(double) * days);
    
    // Find the local moonset on that day
    double temp1[3] = {CALCULATE_SUN_PARAMS_FOR_MOON_TIME, CALCULATE_SUN_PARAMS_FOR_MOON_TIME, CALCULATE_SUN_PARAMS_FOR_MOON_TIME}; 
    double temp2[3] = {CALCULATE_SUN_PARAMS_FOR_MOON_TIME, CALCULATE_SUN_PARAMS_FOR_MOON_TIME, CALCULATE_SUN_PARAMS_FOR_MOON_TIME}; 
    datetime nm_moonset = find_proper_moontime(jd_new_moon, utc_offset, lat, lon, elev, temp, press, 's', temp1, temp2);
    
    // Check if moonset is valid
    int start = 0; // To skip the first day if invalid
    if (compare_datetime(&nm_moonset, &INVALID_DATETIME) == 0) {
        // Moonset doesn't exist, for extreme latitudes
        q_values[0] = -997.0;
        best_jds[0] = jd_new_moon;
        start++; // Skip day 0
    }
    else if (compare_datetime(&nm_moonset, &new_moon_dt) == 0) {
        // Moon is not visibile before the new moon
        q_values[0] = -999.0;
        best_jds[0] = gregorian_to_jd(nm_moonset, utc_offset);
        start++; // Skip day 0
    }

    for (int i = start; i < days; i++) {
        // Set the day parameters
        double test_jd_new_moon = jd_new_moon + i;
        datetime test_ymd_new_moon = add_days(new_moon_dt, i);
        double test_deltaT_new_moon = delta_t_approx(test_ymd_new_moon.year, test_ymd_new_moon.month);

        // Sunset and Moonset calculations
        double deltaPsi[3]; double true_obliquity[3]; // To store the nutations and avoid unnecessary computations
        
        datetime testing = find_proper_suntime(test_jd_new_moon, utc_offset, lat, lon, elev, temp, press, 
            STANDARD_SET_AND_RISE, 's');
        datetime test_nm_sunset = find_proper_suntime_w_nutation(test_jd_new_moon, utc_offset, lat, lon, elev, temp, press, 
            STANDARD_SET_AND_RISE, 's', deltaPsi, true_obliquity);

        datetime test_nm_moonset;
        if (i == 0)
            test_nm_moonset = nm_moonset;
        else
            test_nm_moonset = find_proper_moontime(test_jd_new_moon, utc_offset, lat, lon, elev, temp, press, 's', deltaPsi, true_obliquity);

        // For extreme latitudes where the moonset or sunset don't exist
        // TODO
        
        // If moonset is before sunset, continue
        if (compare_datetime(&test_nm_sunset, &test_nm_moonset) == 1) {
            q_values[i] = -998.0;
            best_jds[i] = gregorian_to_jd(test_nm_moonset, utc_offset);
            continue;
        }

        // Find the best time which is four ninths the moonset-sunset lag after sunset 
        double test_nm_sunset_jd = gregorian_to_jd(test_nm_sunset, utc_offset);
        double test_nm_moonset_jd = gregorian_to_jd(test_nm_moonset, utc_offset);

        double lag_days = fabs(test_nm_sunset_jd - test_nm_moonset_jd);
        double best_time_jd = test_nm_sunset_jd + 4.0 / 9.0 * lag_days;
        best_jds[i] = best_time_jd;
        
        // Compute sun and moon position and parameters
        double best_time_jde = best_time_jd + test_deltaT_new_moon / SECONDS_IN_DAY;
        SunResult nm_sun_params;
        compute_sun_result(best_time_jde, test_deltaT_new_moon, lat, lon, elev, temp, press, &nm_sun_params);

        MoonResult nm_moon_params;
        compute_moon_result(best_time_jde, test_deltaT_new_moon, lat, lon, elev, temp, press, 
            nm_sun_params.nutation_longitude, nm_sun_params.true_obliquity, &nm_moon_params);
        
        // Compute visibility
        switch (criterion) {
            case 0:
                // Odeh uses topocentric apparent (i.e. corrected for standard refraction) horizontal coordinates
                q_values[i] = yallop_odeh(nm_sun_params.true_azimuth, nm_sun_params.apparent_altitude, 
                                          nm_moon_params.true_azimuth, nm_moon_params.apparent_altitude, nm_moon_params.eh_parallax, 
                                          criterion);
                break;
            case 1:
                // Yallop uses geocentric horizontal coordinates
                double sun_geo_alt, sun_geo_az, moon_geo_alt, moon_geo_az;
                geocentric_horizontal_coordinates(lat, nm_sun_params.apparent_declination, nm_sun_params.local_hour_angle, &sun_geo_alt, &sun_geo_az);
                geocentric_horizontal_coordinates(lat, nm_moon_params.declination, nm_moon_params.local_hour_angle, &moon_geo_alt, &moon_geo_az);
                q_values[i] = yallop_odeh(sun_geo_az, sun_geo_alt, moon_geo_az, moon_geo_alt, nm_moon_params.eh_parallax, criterion);
                break;
            case 2:
                // TODO
                q_values[i] = shaukat();
                break;
            default:
                break;
        }
    }

    for (int i = 0; i < days; i++)
        jd_to_gregorian(best_jds[i], utc_offset, &best_dts[i]);
}
    
/* Python wrapper */

PyObject* py_compute_visibilities(PyObject* self, PyObject* args) {
    ENSURE_PYDATETIME();
    PyObject *input_datetime;
    double utc_offset, lat, lon, elev, temp, press;
    int days, criterion;

    if (!PyArg_ParseTuple(args, "O!ddddddii", PyDateTimeAPI->DateTimeType, &input_datetime, &utc_offset, &lat, &lon,
                                                                        &elev, &temp, &press, &days, &criterion))
        return NULL;

    datetime date;
    fill_in_datetime_values(&date, input_datetime);

    double* q_values = malloc(sizeof(double) * days);
    datetime* best_dts = malloc(sizeof(datetime) * days);
    compute_visibilities(date, utc_offset, lat, lon, elev, temp, press, q_values, best_dts, days, criterion);

    PyObject* result = PyTuple_New(days);
    if (!result) return NULL;

    for (int i = 0; i < days; i++) {
        PyObject* py_q = PyFloat_FromDouble(q_values[i]);
        PyObject* py_dt = datetime_to_pydatetime(best_dts[i]);
        if (!py_q || !py_dt) {
            Py_XDECREF(py_q);
            Py_XDECREF(py_dt);
            Py_DECREF(result);
            return NULL;
        }

        PyObject* pair = PyTuple_New(2);
        if (!pair) {
            Py_DECREF(py_q);
            Py_DECREF(py_dt);
            Py_DECREF(result);
            return NULL;
        }

        PyTuple_SET_ITEM(pair, 0, py_q);  // Transfers ownership
        PyTuple_SET_ITEM(pair, 1, py_dt); // Transfers ownership

        PyTuple_SET_ITEM(result, i, pair);  // Transfers ownership
    }

    return result;
}
