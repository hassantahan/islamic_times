from __future__ import annotations

import math
from datetime import datetime, timedelta, timezone

import pytest
import pytz

from islamic_times import time_equations as te
from islamic_times._legacy_py_impl import time_equations as legacy_te
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


def test_delta_t_approx_uses_modern_usno_values_for_2026() -> None:
    with pytest.warns(DeprecationWarning):
        delta_t = te.delta_t_approx(2026, 3)
    assert delta_t == pytest.approx(69.0939, abs=0.02)


def test_resolve_time_scales_matches_toronto_reference_timestamp() -> None:
    dt_utc = datetime(2026, 3, 19, 23, 58, 0, tzinfo=timezone.utc)
    with pytest.warns(DeprecationWarning):
        jd = te.gregorian_to_jd(dt_utc, zone=0)
    tt_minus_utc, ut1_minus_utc, delta_t = te.resolve_time_scales(jd)

    assert tt_minus_utc == pytest.approx(69.184, abs=1e-6)
    assert delta_t == pytest.approx(69.0928431609, abs=0.02)
    assert ut1_minus_utc == pytest.approx(0.0911568391, abs=0.02)


def _huber_standard_error_seconds(decimal_year: float) -> float:
    n = decimal_year - 2005.0
    if n <= 0.0:
        return 0.0
    return 365.25 * n * math.sqrt((n * 0.058 / 3.0) * (1.0 + n / 2500.0)) / 1000.0


def test_future_delta_t_splice_stays_continuous_after_usno_table_end() -> None:
    before_dt = datetime(2033, 10, 1, 0, 0, 0, tzinfo=timezone.utc)
    after_dt = datetime(2033, 10, 2, 0, 0, 0, tzinfo=timezone.utc)

    with pytest.warns(DeprecationWarning):
        before_jd = te.gregorian_to_jd(before_dt, zone=0)
    with pytest.warns(DeprecationWarning):
        after_jd = te.gregorian_to_jd(after_dt, zone=0)

    _, before_ut1_utc, before_delta_t = te.resolve_time_scales(before_jd)
    _, after_ut1_utc, after_delta_t = te.resolve_time_scales(after_jd)

    assert abs(after_delta_t - before_delta_t) < 1.0
    assert abs(after_ut1_utc - before_ut1_utc) < 1.0


@pytest.mark.parametrize("year,month", [(2050, 1), (2090, 1)])
def test_future_delta_t_smoothing_stays_within_huber_envelope(year: int, month: int) -> None:
    with pytest.warns(DeprecationWarning):
        smoothed_delta_t = te.delta_t_approx(year, month)
    polynomial_delta_t = legacy_te.delta_t_approx(year, month)

    decimal_year = year + (month - 0.5) / 12.0
    huber_sigma = _huber_standard_error_seconds(decimal_year)

    assert smoothed_delta_t < polynomial_delta_t
    assert abs(smoothed_delta_t - polynomial_delta_t) <= huber_sigma


def test_find_utc_offset_returns_timezone_and_offset() -> None:
    tz_name, offset_hours = te.find_utc_offset(43.65107, -79.347015, datetime(2025, 1, 15))
    assert isinstance(tz_name, str)
    assert tz_name != ""
    assert isinstance(offset_hours, float)
    assert -12.0 <= offset_hours <= 14.0


def test_find_utc_offset_raises_clear_error_when_timezone_unresolved(monkeypatch: pytest.MonkeyPatch) -> None:
    class DummyFinder:
        def certain_timezone_at(self, lat: float, lng: float) -> None:
            return None

        def timezone_at(self, lat: float, lng: float) -> None:
            return None

    monkeypatch.setattr(te, "_TZ_NAME_CACHE", {})
    monkeypatch.setattr(te, "_UTC_OFFSET_CACHE", {})
    monkeypatch.setattr(te, "_get_timezone_finder", lambda: DummyFinder())

    with pytest.raises(ValueError, match="Unable to determine timezone"):
        te.find_utc_offset(0.0, -160.0, datetime(2025, 1, 15))


def test_find_utc_offset_caches_timezone_resolution(monkeypatch: pytest.MonkeyPatch) -> None:
    class DummyFinder:
        def __init__(self) -> None:
            self.calls = 0

        def certain_timezone_at(self, lat: float, lng: float) -> str:
            self.calls += 1
            return "UTC"

    finder = DummyFinder()

    monkeypatch.setattr(te, "_TZ_NAME_CACHE", {})
    monkeypatch.setattr(te, "_UTC_OFFSET_CACHE", {})
    monkeypatch.setattr(te, "_get_timezone_finder", lambda: finder)

    first = te.find_utc_offset(43.65107, -79.347015, datetime(2025, 1, 15))
    second = te.find_utc_offset(43.65107, -79.347015, datetime(2025, 1, 15))

    assert first == second
    assert finder.calls == 1


def test_find_timezone_name_returns_iana_identifier() -> None:
    tz_name = te.find_timezone_name(43.65107, -79.347015, datetime(2025, 6, 1))
    assert isinstance(tz_name, str)
    assert "/" in tz_name


def test_localize_or_convert_datetime_rejects_nonexistent_local_time() -> None:
    tz = pytz.timezone("America/Toronto")
    with pytest.raises(ValueError, match="does not exist"):
        te.localize_or_convert_datetime(datetime(2025, 3, 9, 2, 30, 0), tz)


def test_localize_or_convert_datetime_rejects_ambiguous_local_time() -> None:
    tz = pytz.timezone("America/Toronto")
    with pytest.raises(ValueError, match="ambiguous"):
        te.localize_or_convert_datetime(datetime(2025, 11, 2, 1, 30, 0), tz)


def test_localize_or_convert_datetime_converts_aware_datetimes() -> None:
    tz = pytz.timezone("America/Toronto")
    converted = te.localize_or_convert_datetime(datetime(2025, 1, 15, 12, 0, 0, tzinfo=timezone.utc), tz)

    assert converted.hour == 7
    assert converted.utcoffset() is not None
    assert converted.utcoffset().total_seconds() / 3600.0 == pytest.approx(-5.0)


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
