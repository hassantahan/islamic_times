"""Typed mapper configuration with explicit validation."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from multiprocessing import cpu_count
from pathlib import Path

from .regions import REGION_COORDINATES


@dataclass(slots=True, frozen=True)
class ComputeConfig:
    """Numerical compute parameters for visibility grid generation."""

    days_to_generate: int = 3
    criterion: int = 1
    utc_offset: float = 0.0
    elevation_m: float = 0.0
    temperature_c: float = 20.0
    pressure_kpa: float = 101.325
    max_workers: int | None = None

    def __post_init__(self) -> None:
        if self.days_to_generate < 1:
            raise ValueError("'days_to_generate' must be >= 1.")
        if self.criterion not in (0, 1, 2):
            raise ValueError("'criterion' must be 0 (Odeh), 1 (Yallop), or 2 (Shaukat).")
        if self.max_workers is not None and self.max_workers < 1:
            raise ValueError("'max_workers' must be >= 1 when provided.")

    @property
    def resolved_workers(self) -> int:
        return self.max_workers or cpu_count()


@dataclass(slots=True, frozen=True)
class RenderConfig:
    """Rendering parameters for map output."""

    map_mode: str = "category"
    dpi: int = 300
    annotate_cities: bool = True
    image_format: str = "jpg"
    jpeg_quality: int | None = None

    def __post_init__(self) -> None:
        mode = self.map_mode.lower()
        if mode not in ("raw", "category"):
            raise ValueError("'map_mode' must be either 'raw' or 'category'.")
        object.__setattr__(self, "map_mode", mode)
        if self.dpi < 72:
            raise ValueError("'dpi' must be >= 72.")
        if self.jpeg_quality is not None and not (1 <= self.jpeg_quality <= 100):
            raise ValueError("'jpeg_quality' must be in [1, 100] when provided.")


@dataclass(slots=True, frozen=True)
class MapperConfig:
    """Top-level mapper execution configuration."""

    date: datetime
    total_months: int = 1
    map_region: str = "WORLD"
    resolution: int = 300
    master_path: str = "maps/"
    save_logs: bool = False
    states_path: str = "map_shp_files/combined_polygons.shp"
    places_path: str = "map_shp_files/combined_points.shp"
    compute: ComputeConfig = ComputeConfig()
    render: RenderConfig = RenderConfig()

    def __post_init__(self) -> None:
        if self.total_months < 1:
            raise ValueError("'total_months' must be >= 1.")
        if self.resolution < 8:
            raise ValueError("'resolution' must be >= 8.")
        region = self.map_region.upper()
        if region not in REGION_COORDINATES:
            raise ValueError(f"'map_region' must be one of: {sorted(REGION_COORDINATES)}.")
        object.__setattr__(self, "map_region", region)

    @property
    def output_root(self) -> Path:
        return Path(self.master_path)
