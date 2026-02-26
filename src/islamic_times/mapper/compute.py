"""Numerical visibility-grid computation utilities for mapper."""

from __future__ import annotations

from multiprocessing import Pool, cpu_count
from typing import Any

import numpy as np

import islamic_times.astro_core as fast_astro

from .config import ComputeConfig
from .palette import category_code_map


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


def _compute_chunk(
    lat_chunk: np.ndarray,
    lon_centers: np.ndarray,
    new_moon_date: Any,
    cfg: ComputeConfig,
    mode: str,
) -> np.ndarray:
    lat_grid, lon_grid = np.meshgrid(lat_chunk, lon_centers, indexing="ij")
    lats_flat = np.ascontiguousarray(lat_grid.ravel(), dtype=np.float64)
    lons_flat = np.ascontiguousarray(lon_grid.ravel(), dtype=np.float64)
    ny, nx = lat_grid.shape

    if mode == "category" and hasattr(fast_astro, "compute_visibilities_batch_codes"):
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
        return result_flat.reshape(ny, nx, cfg.days_to_generate)

    raw_or_labels = fast_astro.compute_visibilities_batch(
        lats_flat,
        lons_flat,
        new_moon_date,
        cfg.days_to_generate,
        cfg.criterion,
        cfg.utc_offset,
        cfg.elevation_m,
        cfg.temperature_c,
        cfg.pressure_kpa,
        "c" if mode == "category" else "r",
    ).reshape(ny, nx, cfg.days_to_generate)

    if mode == "category":
        return _map_category_labels_to_codes(raw_or_labels, cfg.criterion)
    return raw_or_labels


def _split_chunks(lat_centers: np.ndarray, workers: int) -> list[np.ndarray]:
    return [chunk for chunk in np.array_split(lat_centers, workers) if chunk.size]


def compute_visibility_volume(
    lon_centers: np.ndarray,
    lat_centers: np.ndarray,
    new_moon_date: Any,
    cfg: ComputeConfig,
    mode: str,
    pool: Pool | None = None,
) -> np.ndarray:
    """Compute (lat, lon, day) visibility tensor."""
    workers = min(max(1, cfg.max_workers or cpu_count()), len(lat_centers))
    chunks = _split_chunks(lat_centers, workers)
    if len(chunks) == 1:
        return _compute_chunk(chunks[0], lon_centers, new_moon_date, cfg, mode)

    args = [(chunk, lon_centers, new_moon_date, cfg, mode) for chunk in chunks]
    if pool is not None:
        outputs = pool.starmap(_compute_chunk, args)
    else:
        with Pool(workers) as local_pool:
            outputs = local_pool.starmap(_compute_chunk, args)
    return np.concatenate(outputs, axis=0)
