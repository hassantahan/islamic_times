"""End-to-end mapper orchestration."""

from __future__ import annotations

import gc
import json
import os
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timedelta
from multiprocessing import Pool
from pathlib import Path
from time import perf_counter
from typing import Any, Dict, List

import islamic_times.astro_core as fast_astro
from islamic_times.time_equations import get_islamic_month, gregorian_to_hijri

from ._deps import _MAPPER_IMPORT_ERROR, require_mapper_dependencies
from .compute import compute_visibility_volume, create_grid
from .config import MapperConfig
from .geodata import GeoDataCache
from .profiling import StageTimer
from .regions import REGION_CITIES, REGION_COORDINATES
from .render import render_visibility_map

AVERAGE_LUNAR_MONTH_DAYS: float = 29.53059


@dataclass(slots=True)
class MonthRunResult:
    month_index: int
    new_moon_date_utc: datetime
    output_path: str
    timings_s: Dict[str, float]
    compute_profile: Dict[str, Any] | None = None


@dataclass(slots=True)
class MapRunResult:
    config: MapperConfig
    month_results: List[MonthRunResult]
    total_time_s: float
    aggregate_timings_s: Dict[str, float]
    aggregate_compute_profile: Dict[str, Any] | None = None


def _print_ts(message: str) -> None:
    line = f"[{datetime.now().strftime('%X %d-%m-%Y')}] {message}"
    encoding = getattr(sys.stdout, "encoding", None) or "utf-8"
    try:
        print(line)
    except UnicodeEncodeError:
        safe = line.encode(encoding, errors="replace").decode(encoding, errors="replace")
        print(safe)


class Tee:
    """Simple stdout tee for mapper logs."""

    def __init__(self, filename: str, log_dir: str = "mapper_logs", mode: str = "w+", encoding: str = "utf-8"):
        os.makedirs(log_dir, exist_ok=True)
        self.file = open(f"{log_dir}/{filename}", mode, encoding=encoding)
        self.stdout = sys.stdout

    def write(self, message: str) -> None:
        self.stdout.write(message)
        self.file.write(message)

    def flush(self) -> None:
        self.stdout.flush()
        self.file.flush()


def _resolve_islamic_month_bucket(new_moon_date: datetime) -> tuple[int, str]:
    islamic_year, islamic_month, islamic_day = gregorian_to_hijri(new_moon_date.year, new_moon_date.month, new_moon_date.day)
    if islamic_day > 6:
        islamic_month += 1
        if islamic_month > 12:
            islamic_month = 1
            islamic_year += 1
    return islamic_year, get_islamic_month(islamic_month)


def _next_new_moon_date(date: datetime, month_index: int) -> datetime:
    return fast_astro.next_phases_of_moon_utc(date + timedelta(days=month_index * AVERAGE_LUNAR_MONTH_DAYS))[0]


def _build_perf_payload(result: MapRunResult) -> Dict[str, Any]:
    """Build machine-readable perf-report payload from run result."""
    return {
        "config": {
            "date": result.config.date.isoformat(),
            "total_months": result.config.total_months,
            "map_region": result.config.map_region,
            "resolution": result.config.resolution,
            "master_path": result.config.master_path,
            "compute": asdict(result.config.compute),
            "render": asdict(result.config.render),
        },
        "total_time_s": result.total_time_s,
        "aggregate_timings_s": result.aggregate_timings_s,
        "aggregate_compute_profile": result.aggregate_compute_profile,
        "month_results": [
            {
                "month_index": month.month_index,
                "new_moon_date_utc": month.new_moon_date_utc.isoformat(),
                "output_path": month.output_path,
                "timings_s": month.timings_s,
                "compute_profile": month.compute_profile,
            }
            for month in result.month_results
        ],
        "mapper_import_error": str(_MAPPER_IMPORT_ERROR) if _MAPPER_IMPORT_ERROR else None,
    }


def generate_maps(config: MapperConfig, perf_report_path: str | None = None) -> MapRunResult:
    """Generate one or more monthly visibility maps for configured region/workload."""
    require_mapper_dependencies()
    if config.map_region == "WORLD_FULL":
        raise NotImplementedError("WORLD_FULL is not currently supported in mapper pipeline.")

    start_total = perf_counter()
    original_stdout = sys.stdout
    if config.save_logs:
        sys.stdout = Tee(f"mapper_{datetime.now().strftime('%Y-%m-%d_%H%M%S')}.log")

    try:
        _print_ts(f"Mapper start: region={config.map_region}, months={config.total_months}, resolution={config.resolution}")
        cities = REGION_CITIES[config.map_region]
        bbox = REGION_COORDINATES[config.map_region]
        lon_edges, lat_edges, lon_centers, lat_centers = create_grid(config.resolution, *bbox)

        geodata_cache = GeoDataCache(config.states_path, config.places_path)
        states_clip, places_clip = geodata_cache.get_clipped(config.map_region, cities, bbox)

        pool: Pool | None = None
        worker_count = min(max(1, config.compute.resolved_workers), len(lat_centers))
        if worker_count > 1:
            pool = Pool(worker_count)

        month_results: list[MonthRunResult] = []
        aggregate = {"compute": 0.0, "render": 0.0, "other": 0.0}
        aggregate_compute_profile = {
            "months": 0,
            "location_count": 0,
            "total_cells": 0,
            "chunk_count": 0,
            "worker_count_max": 0,
            "compute_elapsed_sum_s": 0.0,
        }

        try:
            for month_idx in range(config.total_months):
                timer = StageTimer()
                new_moon_date = _next_new_moon_date(config.date, month_idx)
                islamic_year, islamic_month_name = _resolve_islamic_month_bucket(new_moon_date)
                out_dir = config.output_root / config.map_region.replace("_", " ").title() / str(islamic_year)

                _print_ts(f"===Generating map for {islamic_month_name}, {islamic_year}===")
                timer.start("compute")
                visibilities, compute_profile = compute_visibility_volume(
                    lon_centers,
                    lat_centers,
                    new_moon_date,
                    cfg=config.compute,
                    mode=config.render.map_mode,
                    pool=pool,
                    return_profile=True,
                )
                compute_dt = timer.stop("compute")
                _print_ts(f"Compute done in {compute_dt:.2f}s")

                timer.start("render")
                output_path = render_visibility_map(
                    lon_edges=lon_edges,
                    lat_edges=lat_edges,
                    lon_centers=lon_centers,
                    lat_centers=lat_centers,
                    visibilities=visibilities,
                    states_clip=states_clip,
                    places_clip=places_clip,
                    start_date=new_moon_date,
                    islamic_month_name=islamic_month_name,
                    islamic_year=islamic_year,
                    criterion=config.compute.criterion,
                    render_cfg=config.render,
                    out_dir=out_dir,
                )
                render_dt = timer.stop("render")
                _print_ts(f"Render done in {render_dt:.2f}s -> {output_path}")
                _print_ts(f"RSS after month clean-up: {__import__('psutil').Process(os.getpid()).memory_info().rss // (1024*1024)} MB")

                aggregate["compute"] += compute_dt
                aggregate["render"] += render_dt
                aggregate_compute_profile["months"] += 1
                aggregate_compute_profile["location_count"] += int(compute_profile.location_count)
                aggregate_compute_profile["total_cells"] += int(compute_profile.total_cells)
                aggregate_compute_profile["chunk_count"] += int(compute_profile.chunk_count)
                aggregate_compute_profile["worker_count_max"] = max(
                    int(aggregate_compute_profile["worker_count_max"]), int(compute_profile.worker_count)
                )
                aggregate_compute_profile["compute_elapsed_sum_s"] += float(compute_profile.total_compute_elapsed_s)
                month_results.append(
                    MonthRunResult(
                        month_index=month_idx,
                        new_moon_date_utc=new_moon_date,
                        output_path=str(output_path),
                        timings_s=dict(timer.elapsed),
                        compute_profile=compute_profile.to_dict(),
                    )
                )
                gc.collect()
        finally:
            if pool is not None:
                pool.close()
                pool.join()

        total_time = perf_counter() - start_total
        aggregate["other"] = max(0.0, total_time - aggregate["compute"] - aggregate["render"])
        _print_ts(f"~~~ --- === Total time taken: {total_time:.2f}s === --- ~~~")

        result = MapRunResult(
            config=config,
            month_results=month_results,
            total_time_s=total_time,
            aggregate_timings_s=aggregate,
            aggregate_compute_profile=aggregate_compute_profile,
        )

        if perf_report_path is not None:
            payload = _build_perf_payload(result)
            Path(perf_report_path).parent.mkdir(parents=True, exist_ok=True)
            Path(perf_report_path).write_text(json.dumps(payload, indent=2), encoding="utf-8")

        return result
    finally:
        if config.save_logs:
            assert isinstance(sys.stdout, Tee)
            sys.stdout.file.close()
            sys.stdout = original_stdout


def main(
    date: datetime | None = None,
    master_path: str = "maps/",
    total_months: int = 1,
    map_region: str = "WORLD",
    map_mode: str = "category",
    resolution: int = 300,
    days_to_generate: int = 3,
    criterion: int = 1,
    save_logs: bool = False,
    max_workers: int | None = None,
    perf_report_path: str | None = None,
) -> MapRunResult:
    """Backward-compatible convenience entrypoint for mapper workflows."""
    from .config import ComputeConfig, RenderConfig

    date = date or datetime.now()
    cfg = MapperConfig(
        date=date,
        total_months=total_months,
        map_region=map_region,
        resolution=resolution,
        master_path=master_path,
        save_logs=save_logs,
        compute=ComputeConfig(days_to_generate=days_to_generate, criterion=criterion, max_workers=max_workers),
        render=RenderConfig(map_mode=map_mode),
    )
    return generate_maps(cfg, perf_report_path=perf_report_path)
