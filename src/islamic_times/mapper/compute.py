"""Numerical visibility-grid computation utilities for mapper."""

from __future__ import annotations

from dataclasses import dataclass
from multiprocessing import Pool, cpu_count
from time import perf_counter
from typing import Any

import numpy as np

import islamic_times.astro_core as fast_astro

from .config import ComputeConfig
from .palette import category_code_map


@dataclass(frozen=True, slots=True)
class ComputeProfile:
    """Machine-readable compute profiling summary for one volume run."""

    mode: str
    criterion: int
    days: int
    worker_count: int
    chunk_count: int
    chunk_rows: tuple[int, ...]
    chunk_elapsed_s: tuple[float, ...]
    backend: str
    used_multiprocessing: bool
    location_count: int
    total_cells: int
    total_compute_elapsed_s: float

    def to_dict(self) -> dict[str, Any]:
        """Serialize profile for perf-report payloads."""
        return {
            "mode": self.mode,
            "criterion": self.criterion,
            "days": self.days,
            "worker_count": self.worker_count,
            "chunk_count": self.chunk_count,
            "chunk_rows": list(self.chunk_rows),
            "chunk_elapsed_s": list(self.chunk_elapsed_s),
            "backend": self.backend,
            "used_multiprocessing": self.used_multiprocessing,
            "location_count": self.location_count,
            "total_cells": self.total_cells,
            "total_compute_elapsed_s": self.total_compute_elapsed_s,
        }


def create_grid(
    resolution: int, minx: float = -179.0, maxx: float = 180.0, miny: float = -61.0, maxy: float = 61.0
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Create edge and center arrays for regular lon/lat grid."""
    lon_edges = np.linspace(minx, maxx, resolution + 1)
    lat_edges = np.linspace(miny, maxy, resolution + 1)
    lon_centers = 0.5 * (lon_edges[:-1] + lon_edges[1:])
    lat_centers = 0.5 * (lat_edges[:-1] + lat_edges[1:])
    return lon_edges, lat_edges, lon_centers, lat_centers


def _map_category_labels_to_codes(labels: np.ndarray, criterion: int) -> np.ndarray:
    code_map = category_code_map(criterion)
    mapped = np.zeros(labels.shape, dtype=np.uint8)
    for label, idx in code_map.items():
        mask = labels == label
        if mask.any():
            mapped[mask] = idx
    return mapped


def _compute_category_points(
    lats: np.ndarray,
    lons: np.ndarray,
    new_moon_date: Any,
    cfg: ComputeConfig,
) -> np.ndarray:
    point_count = int(lats.size)
    if point_count == 0:
        return np.empty((0, cfg.days_to_generate), dtype=np.uint8)
    lats_flat = np.ascontiguousarray(lats, dtype=np.float64)
    lons_flat = np.ascontiguousarray(lons, dtype=np.float64)

    if hasattr(fast_astro, "compute_visibilities_batch_codes"):
        result_flat = fast_astro.compute_visibilities_batch_codes(
            lats_flat,
            lons_flat,
            new_moon_date,
            cfg.days_to_generate,
            cfg.criterion,
            cfg.utc_offset,
            cfg.elevation_m,
            cfg.temperature_c,
            cfg.pressure_kpa,
        )
        return np.ascontiguousarray(result_flat, dtype=np.uint8).reshape(point_count, cfg.days_to_generate)

    labels_flat = fast_astro.compute_visibilities_batch(
        lats_flat,
        lons_flat,
        new_moon_date,
        cfg.days_to_generate,
        cfg.criterion,
        cfg.utc_offset,
        cfg.elevation_m,
        cfg.temperature_c,
        cfg.pressure_kpa,
        "c",
    ).reshape(point_count, cfg.days_to_generate)
    return _map_category_labels_to_codes(labels_flat, cfg.criterion)


def _compute_category_block(
    lat_values: np.ndarray,
    lon_values: np.ndarray,
    new_moon_date: Any,
    cfg: ComputeConfig,
) -> np.ndarray:
    lat_grid, lon_grid = np.meshgrid(lat_values, lon_values, indexing="ij")
    block = _compute_category_points(lat_grid.ravel(), lon_grid.ravel(), new_moon_date, cfg)
    return block.reshape(len(lat_values), len(lon_values), cfg.days_to_generate)


def _probe_lines(start: int, end: int) -> tuple[int, ...]:
    size = end - start
    indices = {start, end - 1, start + (size // 2)}
    if size >= 4:
        indices.add(start + (size // 4))
        indices.add(start + ((3 * size) // 4))
    return tuple(sorted(index for index in indices if start <= index < end))


def _build_probe_indices(y0: int, y1: int, x0: int, x1: int) -> tuple[np.ndarray, np.ndarray]:
    rows = _probe_lines(y0, y1)
    cols = _probe_lines(x0, x1)
    probe_points: set[tuple[int, int]] = set()
    for row in rows:
        for col in range(x0, x1):
            probe_points.add((row, col))
    for col in cols:
        for row in range(y0, y1):
            probe_points.add((row, col))
    if not probe_points:
        return np.empty(0, dtype=np.int32), np.empty(0, dtype=np.int32)
    probe_rows = np.fromiter((pair[0] for pair in probe_points), dtype=np.int32, count=len(probe_points))
    probe_cols = np.fromiter((pair[1] for pair in probe_points), dtype=np.int32, count=len(probe_points))
    return probe_rows, probe_cols


def _is_uniform_probe(probe_codes: np.ndarray) -> bool:
    return probe_codes.shape[0] > 0 and bool(np.all(probe_codes == probe_codes[0]))


def _compute_chunk_category_adaptive(
    lat_chunk: np.ndarray,
    lon_centers: np.ndarray,
    new_moon_date: Any,
    cfg: ComputeConfig,
) -> np.ndarray:
    ny = len(lat_chunk)
    nx = len(lon_centers)
    output = np.empty((ny, nx, cfg.days_to_generate), dtype=np.uint8)
    stack: list[tuple[int, int, int, int, int]] = [(0, ny, 0, nx, 0)]

    while stack:
        y0, y1, x0, x1, depth = stack.pop()
        height = y1 - y0
        width = x1 - x0
        if height <= 0 or width <= 0:
            continue
        cells = height * width
        if (
            cells <= cfg.adaptive_min_block_cells
            or depth >= cfg.adaptive_max_depth
            or height < 2
            or width < 2
        ):
            output[y0:y1, x0:x1, :] = _compute_category_block(
                lat_chunk[y0:y1], lon_centers[x0:x1], new_moon_date, cfg
            )
            continue

        probe_rows, probe_cols = _build_probe_indices(y0, y1, x0, x1)
        probe_codes = _compute_category_points(lat_chunk[probe_rows], lon_centers[probe_cols], new_moon_date, cfg)
        if _is_uniform_probe(probe_codes):
            output[y0:y1, x0:x1, :] = probe_codes[0]
            continue

        y_mid = y0 + (height // 2)
        x_mid = x0 + (width // 2)
        if y_mid <= y0 or y_mid >= y1 or x_mid <= x0 or x_mid >= x1:
            output[y0:y1, x0:x1, :] = _compute_category_block(
                lat_chunk[y0:y1], lon_centers[x0:x1], new_moon_date, cfg
            )
            continue

        next_depth = depth + 1
        stack.append((y0, y_mid, x0, x_mid, next_depth))
        stack.append((y0, y_mid, x_mid, x1, next_depth))
        stack.append((y_mid, y1, x0, x_mid, next_depth))
        stack.append((y_mid, y1, x_mid, x1, next_depth))

    return output


def _compute_chunk_raw(
    lat_chunk: np.ndarray,
    lon_centers: np.ndarray,
    new_moon_date: Any,
    cfg: ComputeConfig,
) -> np.ndarray:
    lat_grid, lon_grid = np.meshgrid(lat_chunk, lon_centers, indexing="ij")
    lats_flat = np.ascontiguousarray(lat_grid.ravel(), dtype=np.float64)
    lons_flat = np.ascontiguousarray(lon_grid.ravel(), dtype=np.float64)
    ny, nx = lat_grid.shape
    return fast_astro.compute_visibilities_batch(
        lats_flat,
        lons_flat,
        new_moon_date,
        cfg.days_to_generate,
        cfg.criterion,
        cfg.utc_offset,
        cfg.elevation_m,
        cfg.temperature_c,
        cfg.pressure_kpa,
        "r",
    ).reshape(ny, nx, cfg.days_to_generate)


def _compute_chunk(
    lat_chunk: np.ndarray,
    lon_centers: np.ndarray,
    new_moon_date: Any,
    cfg: ComputeConfig,
    mode: str,
) -> np.ndarray:
    if mode == "category":
        if cfg.adaptive_category:
            return _compute_chunk_category_adaptive(lat_chunk, lon_centers, new_moon_date, cfg)
        return _compute_category_block(lat_chunk, lon_centers, new_moon_date, cfg)
    return _compute_chunk_raw(lat_chunk, lon_centers, new_moon_date, cfg)


def _compute_chunk_profiled(
    lat_chunk: np.ndarray,
    lon_centers: np.ndarray,
    new_moon_date: Any,
    cfg: ComputeConfig,
    mode: str,
) -> tuple[np.ndarray, float]:
    """Compute one chunk and return elapsed wall-clock seconds."""
    t0 = perf_counter()
    output = _compute_chunk(lat_chunk, lon_centers, new_moon_date, cfg, mode)
    return output, perf_counter() - t0


def _split_chunk_ranges(
    total_rows: int,
    workers: int,
    chunk_multiplier: int,
    min_chunk_rows: int,
) -> list[tuple[int, int]]:
    if total_rows < 1:
        return []
    target_chunks = min(total_rows, max(1, workers * max(1, chunk_multiplier)))
    max_chunks_for_min_rows = max(1, total_rows // max(1, min_chunk_rows))
    target_chunks = min(target_chunks, max_chunks_for_min_rows)
    target_chunks = max(1, target_chunks)

    base, remainder = divmod(total_rows, target_chunks)
    ranges: list[tuple[int, int]] = []
    start = 0
    for idx in range(target_chunks):
        rows = base + (1 if idx < remainder else 0)
        end = start + rows
        if end > start:
            ranges.append((start, end))
        start = end
    return ranges


def _compute_chunk_profiled_task(
    task: tuple[int, int, np.ndarray, np.ndarray, Any, ComputeConfig, str],
) -> tuple[int, int, np.ndarray, float]:
    chunk_idx, row_start, lat_chunk, lon_centers, new_moon_date, cfg, mode = task
    output, elapsed = _compute_chunk_profiled(lat_chunk, lon_centers, new_moon_date, cfg, mode)
    return chunk_idx, row_start, output, float(elapsed)


def _assemble_volume_from_results(
    chunk_results: list[tuple[int, int, np.ndarray, float]],
    total_rows: int,
    lon_count: int,
    days: int,
    chunk_count: int,
) -> tuple[np.ndarray, tuple[float, ...]]:
    if not chunk_results:
        raise ValueError("No chunk results were produced.")

    volume = np.empty((total_rows, lon_count, days), dtype=chunk_results[0][2].dtype)
    chunk_elapsed: list[float] = [0.0] * chunk_count
    for chunk_idx, row_start, output, elapsed in chunk_results:
        row_end = row_start + output.shape[0]
        volume[row_start:row_end, :, :] = output
        chunk_elapsed[chunk_idx] = float(elapsed)
    return volume, tuple(chunk_elapsed)


def compute_visibility_volume(
    lon_centers: np.ndarray,
    lat_centers: np.ndarray,
    new_moon_date: Any,
    cfg: ComputeConfig,
    mode: str,
    pool: Pool | None = None,
    return_profile: bool = False,
) -> np.ndarray | tuple[np.ndarray, ComputeProfile]:
    """Compute (lat, lon, day) visibility tensor.

    When ``return_profile=True``, returns ``(volume, ComputeProfile)``.
    """
    total_rows = len(lat_centers)
    workers = min(max(1, cfg.max_workers or cpu_count()), total_rows)
    chunk_ranges = _split_chunk_ranges(
        total_rows=total_rows,
        workers=workers,
        chunk_multiplier=cfg.chunk_multiplier,
        min_chunk_rows=cfg.min_chunk_rows,
    )
    chunk_rows = tuple(row_end - row_start for row_start, row_end in chunk_ranges)
    chunks = [lat_centers[row_start:row_end] for row_start, row_end in chunk_ranges]
    effective_workers = min(workers, len(chunk_ranges))

    backend = "batch_raw"
    if mode == "category":
        backend = "batch_codes" if hasattr(fast_astro, "compute_visibilities_batch_codes") else "batch_labels"

    if len(chunks) == 1:
        chunk_output, chunk_elapsed = _compute_chunk_profiled(chunks[0], lon_centers, new_moon_date, cfg, mode)
        if not return_profile:
            return chunk_output

        location_count = len(lat_centers) * len(lon_centers)
        profile = ComputeProfile(
            mode=mode,
            criterion=cfg.criterion,
            days=cfg.days_to_generate,
            worker_count=1,
            chunk_count=1,
            chunk_rows=(chunk_rows[0],),
            chunk_elapsed_s=(chunk_elapsed,),
            backend=backend,
            used_multiprocessing=False,
            location_count=location_count,
            total_cells=location_count * cfg.days_to_generate,
            total_compute_elapsed_s=chunk_elapsed,
        )
        return chunk_output, profile

    tasks = [
        (chunk_idx, row_start, lat_centers[row_start:row_end], lon_centers, new_moon_date, cfg, mode)
        for chunk_idx, (row_start, row_end) in enumerate(chunk_ranges)
    ]
    if pool is not None:
        chunk_results = list(pool.imap_unordered(_compute_chunk_profiled_task, tasks))
    elif effective_workers > 1:
        with Pool(effective_workers) as local_pool:
            chunk_results = list(local_pool.imap_unordered(_compute_chunk_profiled_task, tasks))
    else:
        chunk_results = [_compute_chunk_profiled_task(task) for task in tasks]

    volume, chunk_elapsed_s = _assemble_volume_from_results(
        chunk_results=chunk_results,
        total_rows=total_rows,
        lon_count=len(lon_centers),
        days=cfg.days_to_generate,
        chunk_count=len(chunk_ranges),
    )
    if not return_profile:
        return volume

    location_count = len(lat_centers) * len(lon_centers)
    used_multiprocessing = pool is not None or effective_workers > 1
    profile_workers = effective_workers if used_multiprocessing else 1
    profile = ComputeProfile(
        mode=mode,
        criterion=cfg.criterion,
        days=cfg.days_to_generate,
        worker_count=profile_workers,
        chunk_count=len(chunk_ranges),
        chunk_rows=chunk_rows,
        chunk_elapsed_s=chunk_elapsed_s,
        backend=backend,
        used_multiprocessing=used_multiprocessing,
        location_count=location_count,
        total_cells=location_count * cfg.days_to_generate,
        total_compute_elapsed_s=sum(chunk_elapsed_s),
    )
    return volume, profile
