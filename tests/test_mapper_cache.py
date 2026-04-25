from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np

import islamic_times.mapper.pipeline as mapper_pipeline
from islamic_times.mapper.config import ComputeConfig, MapperConfig, RenderConfig


def _sample_cache(values: np.ndarray, *, mode: str = "category") -> mapper_pipeline.VisibilityCache:
    cfg = MapperConfig(
        date=datetime(2025, 6, 1, 12, 0, tzinfo=timezone.utc),
        total_months=1,
        map_region="WORLD",
        resolution=8,
        compute=ComputeConfig(days_to_generate=int(values.shape[2]), criterion=1, max_workers=1),
        render=RenderConfig(map_mode=mode),
    )
    return mapper_pipeline.VisibilityCache(
        config=cfg,
        bbox=(-10.0, 10.0, -10.0, 10.0),
        lon_edges=np.array([-10.0, 0.0, 10.0]),
        lat_edges=np.array([-10.0, 0.0, 10.0]),
        lon_centers=np.array([-5.0, 5.0]),
        lat_centers=np.array([-5.0, 5.0]),
        months=[
            mapper_pipeline.CachedMonthVisibility(
                month_index=0,
                new_moon_date_utc=datetime(2025, 6, 2, 0, 0, tzinfo=timezone.utc),
                islamic_year=1446,
                islamic_month_name="Dhū al-Ḥijjah",
                visibilities=values,
                compute_profile=None,
            )
        ],
    )


def test_visibility_cache_round_trip(tmp_path: Path) -> None:
    cfg = MapperConfig(
        date=datetime(2025, 6, 1, 12, 0, tzinfo=timezone.utc),
        total_months=1,
        map_region="WORLD",
        resolution=10,
        compute=ComputeConfig(days_to_generate=1, criterion=1, max_workers=1),
        render=RenderConfig(map_mode="category"),
    )
    cache = mapper_pipeline.build_visibility_cache(cfg)
    cache_path = mapper_pipeline.save_visibility_cache(cache, tmp_path / "cache" / "world_cache")
    loaded = mapper_pipeline.load_visibility_cache(cache_path)

    assert loaded.config.map_region == cfg.map_region
    assert loaded.config.resolution == cfg.resolution
    assert len(loaded.months) == 1
    assert loaded.months[0].month_index == 0
    assert loaded.months[0].compute_profile is not None
    assert np.array_equal(loaded.months[0].visibilities, cache.months[0].visibilities)
    day = loaded.day_slice(month_index=0, day_index=0)
    assert day.shape == (len(loaded.lat_centers), len(loaded.lon_centers))


def test_visibility_cache_tabular_round_trip(tmp_path: Path) -> None:
    cfg = MapperConfig(
        date=datetime(2025, 6, 1, 12, 0, tzinfo=timezone.utc),
        total_months=1,
        map_region="WORLD",
        resolution=10,
        compute=ComputeConfig(days_to_generate=1, criterion=1, max_workers=1),
        render=RenderConfig(map_mode="category"),
    )
    cache = mapper_pipeline.build_visibility_cache(cfg)
    cache_dir = mapper_pipeline.save_visibility_cache_tabular(cache, tmp_path / "cache_tabular")
    loaded = mapper_pipeline.load_visibility_cache_tabular(cache_dir)

    assert loaded.config.map_region == cfg.map_region
    assert len(loaded.months) == 1
    assert np.array_equal(loaded.months[0].visibilities, cache.months[0].visibilities)


def test_render_maps_from_cache_uses_cached_tensors(monkeypatch: Any, tmp_path: Path) -> None:
    cfg = MapperConfig(
        date=datetime(2025, 6, 1, 12, 0, tzinfo=timezone.utc),
        total_months=1,
        map_region="WORLD",
        resolution=8,
        compute=ComputeConfig(days_to_generate=1, criterion=1, max_workers=1),
        render=RenderConfig(map_mode="category"),
    )
    visibilities = np.zeros((2, 2, 1), dtype=np.uint8)
    cache = mapper_pipeline.VisibilityCache(
        config=cfg,
        bbox=(-10.0, 10.0, -10.0, 10.0),
        lon_edges=np.array([-10.0, 0.0, 10.0]),
        lat_edges=np.array([-10.0, 0.0, 10.0]),
        lon_centers=np.array([-5.0, 5.0]),
        lat_centers=np.array([-5.0, 5.0]),
        months=[
            mapper_pipeline.CachedMonthVisibility(
                month_index=0,
                new_moon_date_utc=datetime(2025, 6, 2, 0, 0, tzinfo=timezone.utc),
                islamic_year=1446,
                islamic_month_name="Dhū al-Ḥijjah",
                visibilities=visibilities,
                compute_profile={
                    "location_count": 4,
                    "total_cells": 4,
                    "chunk_count": 1,
                    "worker_count": 1,
                    "total_compute_elapsed_s": 0.01,
                },
            )
        ],
    )

    rendered: dict[str, Any] = {}

    monkeypatch.setattr(mapper_pipeline, "require_mapper_dependencies", lambda: None)

    class DummyGeoDataCache:
        def __init__(self, states_path: str, places_path: str) -> None:
            self.states_path = states_path
            self.places_path = places_path

        def get_clipped(self, region: str, cities: list[str], bbox: tuple[float, float, float, float]) -> tuple[Any, Any]:
            rendered["region"] = region
            rendered["bbox"] = bbox
            return SimpleNamespace(), SimpleNamespace()

    def fake_render_visibility_map(**kwargs: Any) -> Path:
        rendered["shape"] = kwargs["visibilities"].shape
        out_dir = Path(kwargs["out_dir"])
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / "cached-render.jpg"
        out_path.write_text("ok", encoding="utf-8")
        return out_path

    monkeypatch.setattr(mapper_pipeline, "GeoDataCache", DummyGeoDataCache)
    monkeypatch.setattr(mapper_pipeline, "render_visibility_map", fake_render_visibility_map)

    result = mapper_pipeline.render_maps_from_cache(
        cache,
        master_path=str(tmp_path / "maps_cached"),
        save_logs=False,
    )

    assert len(result.month_results) == 1
    assert result.aggregate_timings_s["compute"] == 0.0
    assert rendered["shape"] == (2, 2, 1)
    assert rendered["region"] == "WORLD"
    assert Path(result.month_results[0].output_path).name == "cached-render.jpg"


def test_cache_month_accessor_returns_first_month(tmp_path: Path) -> None:
    cache = _sample_cache(np.array([[[0], [1]], [[2], [3]]], dtype=np.uint8))
    assert cache.month(0).month_index == 0
