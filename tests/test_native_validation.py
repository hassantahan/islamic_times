from __future__ import annotations

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


def test_jd_to_gregorian_rejects_nonpositive_jd() -> None:
    with pytest.raises(ValueError):
        fast_astro.jd_to_gregorian(-1.0, 0.0)
