from __future__ import annotations

import importlib
from datetime import datetime, timezone

import islamic_times.astro_core as fast_astro
import numpy as np
import pytest

_SPECIAL_CLASSIFICATIONS = {
    "Moonset before the new moon.",
    "Moonset before sunset.",
    "Moonset & Sunset don't exist.",
    "Sunset doesn't exist.",
    "Moonset doesn't exist.",
}
_SHAUKAT_CLASSIFICATIONS = {
    "A: Easily visible",
    "B: Visible under perfect conditions.",
    "C: Optical aid needed to find the moon.",
    "D: Visible with optical aid only.",
    "F: Not visible.",
}


def _expected_shaukat_classification(q: float) -> str:
    if q == -999.0:
        return "Moonset before the new moon."
    if q == -998.0:
        return "Moonset before sunset."
    if q == -997.0:
        return "Moonset & Sunset don't exist."
    if q == -996.0:
        return "Sunset doesn't exist."
    if q == -995.0:
        return "Moonset doesn't exist."

    if q > 0.27:
        return "A: Easily visible"
    if q > -0.024:
        return "B: Visible under perfect conditions."
    if q > -0.212:
        return "C: Optical aid needed to find the moon."
    if q > -0.48:
        return "D: Visible with optical aid only."
    return "F: Not visible."


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

def test_compute_visibilities_supports_criterion_two_shaukat() -> None:
    result = fast_astro.compute_visibilities(
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
    assert result.criterion == "Shaukat"
    assert len(result.classifications) == 1
    assert result.classifications[0] in (_SPECIAL_CLASSIFICATIONS | _SHAUKAT_CLASSIFICATIONS)


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

def test_compute_visibilities_batch_supports_criterion_two_shaukat() -> None:
    lat = np.array([43.651070], dtype=np.float64)
    lon = np.array([-79.347015], dtype=np.float64)
    labels = fast_astro.compute_visibilities_batch(
        lat,
        lon,
        datetime(2025, 6, 1, 12, 0, tzinfo=timezone.utc),
        1,
        2,
        0.0,
        10.0,
        15.0,
        101.325,
        "c",
    )
    assert labels.shape == (1,)
    assert str(labels[0]) in (_SPECIAL_CLASSIFICATIONS | _SHAUKAT_CLASSIFICATIONS)


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


def test_compute_visibilities_batch_codes_support_criterion_two() -> None:
    lat = np.array([43.651070, 40.7128], dtype=np.float64)
    lon = np.array([-79.347015, -74.0060], dtype=np.float64)
    result = fast_astro.compute_visibilities_batch_codes(
        lat,
        lon,
        datetime(2025, 6, 1, 12, 0, tzinfo=timezone.utc),
        2,
        2,
        0.0,
        10.0,
        15.0,
        101.325,
    )
    assert result.dtype == np.uint8
    assert result.shape == (4,)
    assert int(result.max()) <= 9


def test_compute_visibilities_batch_shaukat_labels_match_q_thresholds() -> None:
    lat = np.array([43.651070, 69.6492], dtype=np.float64)
    lon = np.array([-79.347015, 18.9553], dtype=np.float64)
    date = datetime(2025, 6, 1, 12, 0, tzinfo=timezone.utc)
    q_values = fast_astro.compute_visibilities_batch(
        lat, lon, date, 3, 2, 0.0, 10.0, 15.0, 101.325, "r"
    )
    labels = fast_astro.compute_visibilities_batch(
        lat, lon, date, 3, 2, 0.0, 10.0, 15.0, 101.325, "c"
    )

    assert q_values.shape == labels.shape
    for q, label in zip(q_values.tolist(), labels.tolist()):
        assert str(label) == _expected_shaukat_classification(float(q))


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
