from __future__ import annotations

import math
from datetime import datetime, timedelta, timezone

import pytest

from islamic_times import time_equations as te
from islamic_times.it_dataclasses import Angle


def test_fraction_of_day_for_midday() -> None:
    dt = datetime(2025, 6, 1, 12, 0, 0, tzinfo=timezone.utc)
    assert math.isclose(te.fraction_of_day(dt), 0.5, rel_tol=0.0, abs_tol=1e-12)


def test_gmst_returns_angle_with_deprecation_warning() -> None:
    with pytest.warns(DeprecationWarning):
        gmst = te.greenwich_mean_sidereal_time(2451545.0)
    assert isinstance(gmst, Angle)
    assert 0.0 <= gmst.decimal < 360.0


def test_gregorian_jd_roundtrip_close_enough() -> None:
    dt = datetime(2025, 6, 1, 12, 30, 45, tzinfo=timezone.utc)
    with pytest.warns(DeprecationWarning):
        jd = te.gregorian_to_jd(dt, zone=0)
    with pytest.warns(DeprecationWarning):
        roundtrip = te.jd_to_gregorian(jd, adjust_for_tz_diff=0)

    delta_seconds = abs((roundtrip - dt.replace(tzinfo=None)).total_seconds())
    assert delta_seconds < 2


def test_jd_to_gregorian_minus_one_returns_datetime_min() -> None:
    with pytest.warns(DeprecationWarning):
        result = te.jd_to_gregorian(-1)
    assert result == datetime.min


def test_gregorian_to_hijri_output_shape() -> None:
    hijri = te.gregorian_to_hijri(2025, 6, 1)
    assert len(hijri) == 3
    year, month, day = hijri
    assert year > 0
    assert 1 <= month <= 12
    assert 1 <= day <= 30


def test_get_islamic_month_and_day_valid_and_invalid() -> None:
    assert te.get_islamic_month(9) == "Ramaḍān"
    assert te.get_islamic_month(13) == "Invalid month number"
    assert te.get_islamic_day("Friday") == "al-Jumuʿah"
    assert te.get_islamic_day("Funday") == "Invalid day"


@pytest.mark.parametrize("year,month", [(-600, 1), (0, 6), (1500, 3), (1670, 1), (2025, 6), (2200, 1)])
def test_delta_t_approx_branches_return_float(year: int, month: int) -> None:
    with pytest.warns(DeprecationWarning):
        delta_t = te.delta_t_approx(year, month)
    assert isinstance(delta_t, float)


def test_find_utc_offset_returns_timezone_and_offset() -> None:
    tz_name, offset_hours = te.find_utc_offset(43.65107, -79.347015, datetime(2025, 1, 15))
    assert isinstance(tz_name, str)
    assert tz_name != ""
    assert isinstance(offset_hours, float)
    assert -12.0 <= offset_hours <= 14.0


def test_time_midpoint_happy_path() -> None:
    dt1 = datetime(2025, 6, 1, 18, 0, 0, tzinfo=timezone.utc)
    dt2 = datetime(2025, 6, 2, 0, 0, 0, tzinfo=timezone.utc)
    midpoint = te.time_midpoint(dt1, dt2)
    assert midpoint == dt1 + timedelta(hours=3)


def test_time_midpoint_invalid_types_raise_typeerror() -> None:
    with pytest.raises(TypeError):
        te.time_midpoint("2025-01-01", datetime(2025, 1, 1))  # type: ignore[arg-type]


def test_format_utc_offset_formats_hours_and_minutes() -> None:
    assert te.format_utc_offset(5.5) == "UTC+05:30"
    assert te.format_utc_offset(-4.0) == "UTC-04:00"
