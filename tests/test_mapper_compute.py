from __future__ import annotations

from datetime import datetime, timezone

import numpy as np

from islamic_times.mapper.compute import ComputeProfile, compute_visibility_volume, create_grid
from islamic_times.mapper.config import ComputeConfig


def test_compute_visibility_volume_raw_shape() -> None:
    _, _, lon_centers, lat_centers = create_grid(12, -20, 20, 0, 40)
    cfg = ComputeConfig(days_to_generate=2, criterion=1, max_workers=1)
    volume = compute_visibility_volume(
        lon_centers, lat_centers, datetime(2025, 6, 1, 12, 0, tzinfo=timezone.utc), cfg, "raw"
    )
    assert volume.shape == (len(lat_centers), len(lon_centers), 2)
    assert np.issubdtype(volume.dtype, np.floating)


def test_compute_visibility_volume_category_shape_and_dtype() -> None:
    _, _, lon_centers, lat_centers = create_grid(10, -15, 15, 10, 40)
    cfg = ComputeConfig(days_to_generate=1, criterion=1, max_workers=1)
    volume = compute_visibility_volume(
        lon_centers,
        lat_centers,
        datetime(2025, 6, 1, 12, 0, tzinfo=timezone.utc),
        cfg,
        "category",
    )
    assert volume.shape == (len(lat_centers), len(lon_centers), 1)
    assert volume.dtype == np.uint8
    assert int(volume.max()) <= 10


def test_compute_visibility_volume_category_shaukat_code_range() -> None:
    _, _, lon_centers, lat_centers = create_grid(10, -15, 15, 10, 40)
    cfg = ComputeConfig(days_to_generate=1, criterion=2, max_workers=1)
    volume = compute_visibility_volume(
        lon_centers,
        lat_centers,
        datetime(2025, 6, 1, 12, 0, tzinfo=timezone.utc),
        cfg,
        "category",
    )
    assert volume.shape == (len(lat_centers), len(lon_centers), 1)
    assert volume.dtype == np.uint8
    assert int(volume.max()) <= 9


def test_compute_visibility_volume_returns_profile_when_requested() -> None:
    _, _, lon_centers, lat_centers = create_grid(8, -10, 10, 10, 30)
    cfg = ComputeConfig(days_to_generate=2, criterion=1, max_workers=1)
    volume, profile = compute_visibility_volume(
        lon_centers,
        lat_centers,
        datetime(2025, 6, 1, 12, 0, tzinfo=timezone.utc),
        cfg,
        "raw",
        return_profile=True,
    )

    assert volume.shape == (len(lat_centers), len(lon_centers), 2)
    assert isinstance(profile, ComputeProfile)
    assert profile.mode == "raw"
    assert profile.criterion == 1
    assert profile.days == 2
    assert profile.worker_count == 1
    assert profile.chunk_count == 1
    assert profile.location_count == len(lat_centers) * len(lon_centers)
    assert profile.total_cells == profile.location_count * 2
    assert profile.total_compute_elapsed_s >= 0.0
    serialized = profile.to_dict()
    assert serialized["chunk_count"] == 1
