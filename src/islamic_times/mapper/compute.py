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


def _split_chunks(lat_centers: np.ndarray, workers: int) -> list[np.ndarray]:
    return [chunk for chunk in np.array_split(lat_centers, workers) if chunk.size]


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
    workers = min(max(1, cfg.max_workers or cpu_count()), len(lat_centers))
    chunks = _split_chunks(lat_centers, workers)

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
            chunk_rows=(len(chunks[0]),),
            chunk_elapsed_s=(chunk_elapsed,),
            backend=backend,
            used_multiprocessing=False,
            location_count=location_count,
            total_cells=location_count * cfg.days_to_generate,
            total_compute_elapsed_s=chunk_elapsed,
        )
        return chunk_output, profile

    args = [(chunk, lon_centers, new_moon_date, cfg, mode) for chunk in chunks]
    if pool is not None:
        outputs_profiled = pool.starmap(_compute_chunk_profiled, args)
    else:
        with Pool(workers) as local_pool:
            outputs_profiled = local_pool.starmap(_compute_chunk_profiled, args)

    outputs = [out for out, _ in outputs_profiled]
    chunk_elapsed_s = tuple(float(elapsed) for _, elapsed in outputs_profiled)
    volume = np.concatenate(outputs, axis=0)
    if not return_profile:
        return volume

    location_count = len(lat_centers) * len(lon_centers)
    profile = ComputeProfile(
        mode=mode,
        criterion=cfg.criterion,
        days=cfg.days_to_generate,
        worker_count=workers,
        chunk_count=len(chunks),
        chunk_rows=tuple(len(chunk) for chunk in chunks),
        chunk_elapsed_s=chunk_elapsed_s,
        backend=backend,
        used_multiprocessing=True,
        location_count=location_count,
        total_cells=location_count * cfg.days_to_generate,
        total_compute_elapsed_s=sum(chunk_elapsed_s),
    )
    return volume, profile
