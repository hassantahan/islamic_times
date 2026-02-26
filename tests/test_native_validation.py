from __future__ import annotations

import importlib
from datetime import datetime, timezone

import numpy as np
import pytest

import islamic_times.astro_core as fast_astro


def test_compute_visibilities_rejects_invalid_days() -> None:
    with pytest.raises(ValueError):
        fast_astro.compute_visibilities(
            datetime(2025, 6, 1, 12, 0, tzinfo=timezone.utc),
            0.0,
            43.651070,
            -79.347015,
            10.0,
            15.0,
            101.325,
            0,
            1,
        )


def test_compute_visibilities_rejects_invalid_criterion() -> None:
    with pytest.raises(ValueError):
        fast_astro.compute_visibilities(
            datetime(2025, 6, 1, 12, 0, tzinfo=timezone.utc),
            0.0,
            43.651070,
            -79.347015,
            10.0,
            15.0,
            101.325,
            1,
            99,
        )

def test_compute_visibilities_rejects_unimplemented_criterion_two() -> None:
    with pytest.raises(ValueError, match="Criterion 2 \\(Shaukat\\) is not implemented"):
        fast_astro.compute_visibilities(
            datetime(2025, 6, 1, 12, 0, tzinfo=timezone.utc),
            0.0,
            43.651070,
            -79.347015,
            10.0,
            15.0,
            101.325,
            1,
            2,
        )


def test_compute_visibilities_batch_requires_datetime_input() -> None:
    lat = np.array([43.651070], dtype=np.float64)
    lon = np.array([-79.347015], dtype=np.float64)
    with pytest.raises(TypeError):
        fast_astro.compute_visibilities_batch(
            lat,
            lon,
            "2025-06-01T12:00:00Z",
            1,
            1,
            0.0,
            10.0,
            15.0,
            101.325,
            "r",
        )


def test_compute_visibilities_batch_rejects_invalid_type_flag() -> None:
    lat = np.array([43.651070], dtype=np.float64)
    lon = np.array([-79.347015], dtype=np.float64)
    with pytest.raises(ValueError):
        fast_astro.compute_visibilities_batch(
            lat,
            lon,
            datetime(2025, 6, 1, 12, 0, tzinfo=timezone.utc),
            1,
            1,
            0.0,
            10.0,
            15.0,
            101.325,
            "x",
        )


def test_compute_visibilities_batch_rejects_invalid_days() -> None:
    lat = np.array([43.651070], dtype=np.float64)
    lon = np.array([-79.347015], dtype=np.float64)
    with pytest.raises(ValueError):
        fast_astro.compute_visibilities_batch(
            lat,
            lon,
            datetime(2025, 6, 1, 12, 0, tzinfo=timezone.utc),
            0,
            1,
            0.0,
            10.0,
            15.0,
            101.325,
            "r",
        )

def test_compute_visibilities_batch_rejects_unimplemented_criterion_two() -> None:
    lat = np.array([43.651070], dtype=np.float64)
    lon = np.array([-79.347015], dtype=np.float64)
    with pytest.raises(ValueError, match="Criterion 2 \\(Shaukat\\) is not implemented"):
        fast_astro.compute_visibilities_batch(
            lat,
            lon,
            datetime(2025, 6, 1, 12, 0, tzinfo=timezone.utc),
            1,
            2,
            0.0,
            10.0,
            15.0,
            101.325,
            "r",
        )


def test_compute_visibilities_batch_codes_returns_uint8_codes() -> None:
    lat = np.array([43.651070, 40.7128], dtype=np.float64)
    lon = np.array([-79.347015, -74.0060], dtype=np.float64)
    result = fast_astro.compute_visibilities_batch_codes(
        lat,
        lon,
        datetime(2025, 6, 1, 12, 0, tzinfo=timezone.utc),
        2,
        1,
        0.0,
        10.0,
        15.0,
        101.325,
    )
    assert result.dtype == np.uint8
    assert result.shape == (4,)


def test_jd_to_gregorian_rejects_nonpositive_jd() -> None:
    with pytest.raises(ValueError):
        fast_astro.jd_to_gregorian(-1.0, 0.0)


def test_moon_wrappers_accept_integer_nutation_sequences() -> None:
    dt = datetime(2025, 6, 1, 12, 0, tzinfo=timezone.utc)
    jd = fast_astro.gregorian_to_jd(dt, 0.0)

    transit = fast_astro.find_moon_transit(
        jd,
        73.0,
        43.651070,
        -79.347015,
        10.0,
        15.0,
        101.325,
        0.0,
        [0, 0, 0],
        [23, 23, 23],
    )
    moontime = fast_astro.find_proper_moontime(
        jd,
        73.0,
        43.651070,
        -79.347015,
        10.0,
        15.0,
        101.325,
        0.0,
        [0, 0, 0],
        [23, 23, 23],
        "s",
    )

    assert isinstance(transit, datetime)
    assert isinstance(moontime, datetime)


def test_moon_wrappers_accept_none_for_auto_solar_terms() -> None:
    dt = datetime(2025, 6, 1, 12, 0, tzinfo=timezone.utc)
    jd = fast_astro.gregorian_to_jd(dt, 0.0)

    transit = fast_astro.find_moon_transit(
        jd,
        73.0,
        43.651070,
        -79.347015,
        10.0,
        15.0,
        101.325,
        0.0,
        None,
        None,
    )
    moontime = fast_astro.find_proper_moontime(
        jd,
        73.0,
        43.651070,
        -79.347015,
        10.0,
        15.0,
        101.325,
        0.0,
        None,
        None,
        "s",
    )

    assert isinstance(transit, datetime)
    assert isinstance(moontime, datetime)


def test_moon_wrappers_reject_mixed_none_and_sequence_nutation_inputs() -> None:
    dt = datetime(2025, 6, 1, 12, 0, tzinfo=timezone.utc)
    jd = fast_astro.gregorian_to_jd(dt, 0.0)

    with pytest.raises(TypeError, match="must both be None or both be sequences"):
        fast_astro.find_moon_transit(
            jd,
            73.0,
            43.651070,
            -79.347015,
            10.0,
            15.0,
            101.325,
            0.0,
            None,
            [23.0, 23.0, 23.0],
        )


def test_find_proper_suntime_rejects_invalid_event_code() -> None:
    dt = datetime(2025, 6, 1, 12, 0, tzinfo=timezone.utc)
    jd = fast_astro.gregorian_to_jd(dt, 0.0)
    with pytest.raises(ValueError):
        fast_astro.find_proper_suntime(
            jd,
            43.651070,
            -79.347015,
            10.0,
            15.0,
            101.325,
            0.0,
            0.8333333333,
            "x",
        )


def test_find_proper_moontime_rejects_invalid_event_code() -> None:
    dt = datetime(2025, 6, 1, 12, 0, tzinfo=timezone.utc)
    jd = fast_astro.gregorian_to_jd(dt, 0.0)
    with pytest.raises(ValueError):
        fast_astro.find_proper_moontime(
            jd,
            73.0,
            43.651070,
            -79.347015,
            10.0,
            15.0,
            101.325,
            0.0,
            [0.0, 0.0, 0.0],
            [23.0, 23.0, 23.0],
            "x",
        )


def test_astro_core_reload_keeps_type_backed_wrappers_operational() -> None:
    module = fast_astro
    observer_dt = datetime(2025, 6, 1, 12, 0, tzinfo=timezone.utc)

    for _ in range(3):
        module = importlib.reload(module)

        jd = module.gregorian_to_jd(observer_dt, 0.0)
        delta_t = module.delta_t_approx(observer_dt.year, observer_dt.month)

        sun = module.compute_sun(
            jd + delta_t / 86400.0,
            delta_t,
            43.651070,
            -79.347015,
            10.0,
            15.0,
            101.325,
        )
        moon = module.compute_moon(
            jd + delta_t / 86400.0,
            delta_t,
            43.651070,
            -79.347015,
            10.0,
            15.0,
            101.325,
            0.0,
            23.4,
        )
        vis = module.compute_visibilities(
            observer_dt,
            0.0,
            43.651070,
            -79.347015,
            10.0,
            15.0,
            101.325,
            1,
            1,
        )

        assert sun is not None
        assert moon is not None
        assert len(vis.q_values) == 1
