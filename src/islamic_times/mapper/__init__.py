"""Public mapper API."""

from __future__ import annotations

from ._deps import _MAPPER_IMPORT_ERROR, require_mapper_dependencies
from .config import ComputeConfig, MapperConfig, RenderConfig
from .pipeline import (
    CachedMonthVisibility,
    MapRunResult,
    MonthRunResult,
    VisibilityCache,
    build_visibility_cache,
    generate_maps,
    load_visibility_cache_tabular,
    load_visibility_cache,
    main,
    render_maps_from_cache,
    save_visibility_cache,
    save_visibility_cache_tabular,
    write_perf_report,
)

__all__ = [
    "_MAPPER_IMPORT_ERROR",
    "require_mapper_dependencies",
    "ComputeConfig",
    "RenderConfig",
    "MapperConfig",
    "MonthRunResult",
    "MapRunResult",
    "CachedMonthVisibility",
    "VisibilityCache",
    "build_visibility_cache",
    "save_visibility_cache",
    "save_visibility_cache_tabular",
    "load_visibility_cache",
    "load_visibility_cache_tabular",
    "render_maps_from_cache",
    "write_perf_report",
    "generate_maps",
    "main",
]

