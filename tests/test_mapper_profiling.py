from __future__ import annotations

from datetime import datetime, timezone

from islamic_times.mapper.config import ComputeConfig, MapperConfig, RenderConfig
from islamic_times.mapper.pipeline import MapRunResult, MonthRunResult, _build_perf_payload


def test_build_perf_payload_includes_compute_profiles() -> None:
    cfg = MapperConfig(
        date=datetime(2025, 6, 1, 12, 0, tzinfo=timezone.utc),
        map_region="WORLD",
        resolution=120,
        compute=ComputeConfig(days_to_generate=2, criterion=1, max_workers=2),
        render=RenderConfig(map_mode="category"),
    )
    month_result = MonthRunResult(
        month_index=0,
        new_moon_date_utc=datetime(2025, 6, 2, 0, 0, tzinfo=timezone.utc),
        output_path="maps/world/test.jpg",
        timings_s={"compute": 1.25, "render": 0.75},
        compute_profile={
            "mode": "category",
            "criterion": 1,
            "days": 2,
            "worker_count": 2,
            "chunk_count": 2,
            "chunk_rows": [60, 60],
            "chunk_elapsed_s": [0.5, 0.6],
            "backend": "batch_codes",
            "used_multiprocessing": True,
            "location_count": 14400,
            "total_cells": 28800,
            "total_compute_elapsed_s": 1.1,
        },
    )
    run_result = MapRunResult(
        config=cfg,
        month_results=[month_result],
        total_time_s=2.5,
        aggregate_timings_s={"compute": 1.25, "render": 0.75, "other": 0.5},
        aggregate_compute_profile={
            "months": 1,
            "location_count": 14400,
            "total_cells": 28800,
            "chunk_count": 2,
            "worker_count_max": 2,
            "compute_elapsed_sum_s": 1.1,
        },
    )

    payload = _build_perf_payload(run_result)

    assert payload["total_time_s"] == 2.5
    assert payload["aggregate_timings_s"]["compute"] == 1.25
    assert payload["aggregate_compute_profile"]["chunk_count"] == 2
    assert payload["month_results"][0]["timings_s"]["render"] == 0.75
    assert payload["month_results"][0]["compute_profile"]["backend"] == "batch_codes"
    assert payload["month_results"][0]["compute_profile"]["total_cells"] == 28800

