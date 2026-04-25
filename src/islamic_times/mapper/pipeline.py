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
import numpy as np
from islamic_times.time_equations import get_islamic_month, gregorian_to_hijri

from ._deps import _MAPPER_IMPORT_ERROR, require_mapper_dependencies
from .compute import compute_visibility_volume, create_grid
from .config import ComputeConfig, MapperConfig, RenderConfig
from .geodata import GeoDataCache
from .profiling import StageTimer
from .regions import REGION_CITIES, REGION_COORDINATES
from .render import render_visibility_map

AVERAGE_LUNAR_MONTH_DAYS: float = 29.53059
CACHE_TABULAR_SCHEMA_VERSION: str = "1.0"


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


@dataclass(slots=True)
class CachedMonthVisibility:
    month_index: int
    new_moon_date_utc: datetime
    islamic_year: int
    islamic_month_name: str
    visibilities: np.ndarray
    compute_profile: Dict[str, Any] | None = None

    def day_slice(self, day_index: int = 0) -> np.ndarray:
        days = int(self.visibilities.shape[2])
        if day_index < 0 or day_index >= days:
            raise IndexError(f"'day_index' must be in [0, {days - 1}].")
        return self.visibilities[:, :, day_index]


@dataclass(slots=True)
class VisibilityCache:
    config: MapperConfig
    bbox: tuple[float, float, float, float]
    lon_edges: np.ndarray
    lat_edges: np.ndarray
    lon_centers: np.ndarray
    lat_centers: np.ndarray
    months: List[CachedMonthVisibility]

    def month(self, month_index: int = 0) -> CachedMonthVisibility:
        if month_index < 0 or month_index >= len(self.months):
            raise IndexError(f"'month_index' must be in [0, {len(self.months) - 1}].")
        return self.months[month_index]

    def day_slice(self, month_index: int = 0, day_index: int = 0) -> np.ndarray:
        return self.month(month_index).day_slice(day_index)


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


def _config_payload(config: MapperConfig) -> Dict[str, Any]:
    return {
        "date": config.date.isoformat(),
        "total_months": config.total_months,
        "map_region": config.map_region,
        "resolution": config.resolution,
        "master_path": config.master_path,
        "save_logs": config.save_logs,
        "states_path": config.states_path,
        "places_path": config.places_path,
        "compute": asdict(config.compute),
        "render": asdict(config.render),
    }


def _config_from_payload(payload: Dict[str, Any]) -> MapperConfig:
    return MapperConfig(
        date=datetime.fromisoformat(str(payload["date"])),
        total_months=int(payload.get("total_months", 1)),
        map_region=str(payload.get("map_region", "WORLD")),
        resolution=int(payload.get("resolution", 300)),
        master_path=str(payload.get("master_path", "maps/")),
        save_logs=bool(payload.get("save_logs", False)),
        states_path=str(payload.get("states_path", "map_shp_files/combined_polygons.shp")),
        places_path=str(payload.get("places_path", "map_shp_files/combined_points.shp")),
        compute=ComputeConfig(**dict(payload.get("compute", {}))),
        render=RenderConfig(**dict(payload.get("render", {}))),
    )


def write_perf_report(result: MapRunResult, perf_report_path: str | Path) -> None:
    payload = _build_perf_payload(result)
    out_path = Path(perf_report_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def build_visibility_cache(config: MapperConfig) -> VisibilityCache:
    """Compute and retain visibility tensors for later reuse without recompute."""
    if config.map_region == "WORLD_FULL":
        raise NotImplementedError("WORLD_FULL is not currently supported in mapper pipeline.")

    bbox = REGION_COORDINATES[config.map_region]
    lon_edges, lat_edges, lon_centers, lat_centers = create_grid(config.resolution, *bbox)

    pool: Pool | None = None
    worker_count = min(max(1, config.compute.resolved_workers), len(lat_centers))
    if worker_count > 1:
        pool = Pool(worker_count)

    months: list[CachedMonthVisibility] = []
    try:
        for month_idx in range(config.total_months):
            new_moon_date = _next_new_moon_date(config.date, month_idx)
            islamic_year, islamic_month_name = _resolve_islamic_month_bucket(new_moon_date)
            visibilities, compute_profile = compute_visibility_volume(
                lon_centers,
                lat_centers,
                new_moon_date,
                cfg=config.compute,
                mode=config.render.map_mode,
                pool=pool,
                return_profile=True,
            )
            months.append(
                CachedMonthVisibility(
                    month_index=month_idx,
                    new_moon_date_utc=new_moon_date,
                    islamic_year=islamic_year,
                    islamic_month_name=islamic_month_name,
                    visibilities=visibilities,
                    compute_profile=compute_profile.to_dict(),
                )
            )
    finally:
        if pool is not None:
            pool.close()
            pool.join()

    return VisibilityCache(
        config=config,
        bbox=bbox,
        lon_edges=lon_edges,
        lat_edges=lat_edges,
        lon_centers=lon_centers,
        lat_centers=lat_centers,
        months=months,
    )


def save_visibility_cache(cache: VisibilityCache, cache_path: str | Path) -> Path:
    """Persist precomputed visibility tensors and metadata as tabular cache."""
    return save_visibility_cache_tabular(cache, cache_path)


def save_visibility_cache_tabular(cache: VisibilityCache, cache_path: str | Path) -> Path:
    """Persist cache as a directory with JSON manifest and month-level CSV tables."""
    out_dir = Path(cache_path)
    out_dir.mkdir(parents=True, exist_ok=True)

    lat_grid, lon_grid = np.meshgrid(cache.lat_centers, cache.lon_centers, indexing="ij")
    lat_flat = lat_grid.ravel()
    lon_flat = lon_grid.ravel()

    months_payload: list[Dict[str, Any]] = []
    is_category = cache.config.render.map_mode == "category"
    for month in cache.months:
        csv_name = f"month_{month.month_index:03d}.csv"
        csv_path = out_dir / csv_name
        days = int(month.visibilities.shape[2])
        values_flat = month.visibilities.reshape(-1, days)
        matrix = np.column_stack([lat_flat, lon_flat, values_flat])
        header = ",".join(["lat", "lon"] + [f"day_{idx}" for idx in range(days)])
        day_fmt = "%d" if is_category else "%.10f"
        fmt = ["%.8f", "%.8f"] + [day_fmt] * days
        np.savetxt(csv_path, matrix, delimiter=",", header=header, comments="", fmt=fmt)
        months_payload.append(
            {
                "month_index": month.month_index,
                "new_moon_date_utc": month.new_moon_date_utc.isoformat(),
                "islamic_year": month.islamic_year,
                "islamic_month_name": month.islamic_month_name,
                "compute_profile": month.compute_profile,
                "days": days,
                "csv_file": csv_name,
            }
        )

    manifest = {
        "schema_version": CACHE_TABULAR_SCHEMA_VERSION,
        "cache_format": "tabular",
        "config": _config_payload(cache.config),
        "bbox": list(cache.bbox),
        "lon_edges": cache.lon_edges.tolist(),
        "lat_edges": cache.lat_edges.tolist(),
        "lon_centers": cache.lon_centers.tolist(),
        "lat_centers": cache.lat_centers.tolist(),
        "months": months_payload,
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return out_dir


def load_visibility_cache(cache_path: str | Path) -> VisibilityCache:
    """Load tabular visibility cache created by ``save_visibility_cache_tabular``."""
    return load_visibility_cache_tabular(cache_path)


def load_visibility_cache_tabular(cache_path: str | Path) -> VisibilityCache:
    """Load tabular CSV export data from a directory or manifest JSON path."""
    input_path = Path(cache_path)
    manifest_path = input_path / "manifest.json" if input_path.is_dir() else input_path
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != CACHE_TABULAR_SCHEMA_VERSION:
        raise ValueError(
            f"Unsupported tabular cache schema version: {manifest.get('schema_version')!r}. "
            f"Expected {CACHE_TABULAR_SCHEMA_VERSION!r}."
        )
    if manifest.get("cache_format") != "tabular":
        raise ValueError(f"Unsupported tabular cache format marker: {manifest.get('cache_format')!r}.")

    base_dir = manifest_path.parent
    config = _config_from_payload(dict(manifest["config"]))
    lon_edges = np.array(manifest["lon_edges"], dtype=np.float64)
    lat_edges = np.array(manifest["lat_edges"], dtype=np.float64)
    lon_centers = np.array(manifest["lon_centers"], dtype=np.float64)
    lat_centers = np.array(manifest["lat_centers"], dtype=np.float64)
    lat_count = len(lat_centers)
    lon_count = len(lon_centers)
    total_cells = lat_count * lon_count

    months: list[CachedMonthVisibility] = []
    for month_payload in manifest["months"]:
        csv_file = base_dir / str(month_payload["csv_file"])
        rows = np.loadtxt(csv_file, delimiter=",", skiprows=1)
        if rows.ndim == 1:
            rows = rows.reshape(1, -1)
        if rows.shape[0] != total_cells:
            raise ValueError(
                f"Cache table cell count mismatch for {csv_file}: got {rows.shape[0]}, expected {total_cells}."
            )

        values_flat = rows[:, 2:]
        if values_flat.ndim == 1:
            values_flat = values_flat.reshape(-1, 1)
        days = int(month_payload["days"])
        visibilities = values_flat.reshape(lat_count, lon_count, days)
        if config.render.map_mode == "category":
            visibilities = np.rint(visibilities).astype(np.uint8)
        else:
            visibilities = visibilities.astype(np.float64)

        months.append(
            CachedMonthVisibility(
                month_index=int(month_payload["month_index"]),
                new_moon_date_utc=datetime.fromisoformat(str(month_payload["new_moon_date_utc"])),
                islamic_year=int(month_payload["islamic_year"]),
                islamic_month_name=str(month_payload["islamic_month_name"]),
                visibilities=visibilities,
                compute_profile=dict(month_payload["compute_profile"])
                if month_payload.get("compute_profile") is not None
                else None,
            )
        )

    return VisibilityCache(
        config=config,
        bbox=tuple(float(v) for v in manifest["bbox"]),
        lon_edges=lon_edges,
        lat_edges=lat_edges,
        lon_centers=lon_centers,
        lat_centers=lat_centers,
        months=months,
    )


def render_maps_from_cache(
    cache: VisibilityCache,
    master_path: str | None = None,
    save_logs: bool | None = None,
    perf_report_path: str | None = None,
) -> MapRunResult:
    """Render maps from precomputed cache tensors without rerunning compute."""
    require_mapper_dependencies()
    output_root = Path(master_path) if master_path is not None else cache.config.output_root
    log_enabled = cache.config.save_logs if save_logs is None else bool(save_logs)
    result_config = MapperConfig(
        date=cache.config.date,
        total_months=len(cache.months),
        map_region=cache.config.map_region,
        resolution=cache.config.resolution,
        master_path=str(output_root),
        save_logs=log_enabled,
        states_path=cache.config.states_path,
        places_path=cache.config.places_path,
        compute=cache.config.compute,
        render=cache.config.render,
    )

    start_total = perf_counter()
    original_stdout = sys.stdout
    if result_config.save_logs:
        sys.stdout = Tee(f"mapper_render_cache_{datetime.now().strftime('%Y-%m-%d_%H%M%S')}.log")

    try:
        cities = REGION_CITIES[result_config.map_region]
        geodata_cache = GeoDataCache(result_config.states_path, result_config.places_path)
        states_clip, places_clip = geodata_cache.get_clipped(result_config.map_region, cities, cache.bbox)
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

        for month in cache.months:
            timer = StageTimer()
            out_dir = result_config.output_root / result_config.map_region.replace("_", " ").title() / str(month.islamic_year)
            timer.start("render")
            output_path = render_visibility_map(
                lon_edges=cache.lon_edges,
                lat_edges=cache.lat_edges,
                lon_centers=cache.lon_centers,
                lat_centers=cache.lat_centers,
                visibilities=month.visibilities,
                states_clip=states_clip,
                places_clip=places_clip,
                start_date=month.new_moon_date_utc,
                islamic_month_name=month.islamic_month_name,
                islamic_year=month.islamic_year,
                criterion=result_config.compute.criterion,
                render_cfg=result_config.render,
                out_dir=out_dir,
            )
            render_dt = timer.stop("render")
            aggregate["render"] += render_dt

            if month.compute_profile is not None:
                aggregate_compute_profile["months"] += 1
                aggregate_compute_profile["location_count"] += int(month.compute_profile.get("location_count", 0))
                aggregate_compute_profile["total_cells"] += int(month.compute_profile.get("total_cells", 0))
                aggregate_compute_profile["chunk_count"] += int(month.compute_profile.get("chunk_count", 0))
                aggregate_compute_profile["worker_count_max"] = max(
                    int(aggregate_compute_profile["worker_count_max"]),
                    int(month.compute_profile.get("worker_count", 0)),
                )
                aggregate_compute_profile["compute_elapsed_sum_s"] += float(
                    month.compute_profile.get("total_compute_elapsed_s", 0.0)
                )

            month_results.append(
                MonthRunResult(
                    month_index=month.month_index,
                    new_moon_date_utc=month.new_moon_date_utc,
                    output_path=str(output_path),
                    timings_s={"compute": 0.0, "render": render_dt},
                    compute_profile=month.compute_profile,
                )
            )

        total_time = perf_counter() - start_total
        aggregate["other"] = max(0.0, total_time - aggregate["compute"] - aggregate["render"])
        result = MapRunResult(
            config=result_config,
            month_results=month_results,
            total_time_s=total_time,
            aggregate_timings_s=aggregate,
            aggregate_compute_profile=aggregate_compute_profile,
        )

        if perf_report_path is not None:
            write_perf_report(result, perf_report_path)
        return result
    finally:
        if result_config.save_logs:
            assert isinstance(sys.stdout, Tee)
            sys.stdout.file.close()
            sys.stdout = original_stdout


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
            write_perf_report(result, perf_report_path)

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
