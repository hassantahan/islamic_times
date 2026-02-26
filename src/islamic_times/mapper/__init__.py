"""Public mapper API."""

from __future__ import annotations

from ._deps import _MAPPER_IMPORT_ERROR, require_mapper_dependencies
from .config import ComputeConfig, MapperConfig, RenderConfig
from .pipeline import MapRunResult, MonthRunResult, generate_maps, main

__all__ = [
    "_MAPPER_IMPORT_ERROR",
    "require_mapper_dependencies",
    "ComputeConfig",
    "RenderConfig",
    "MapperConfig",
    "MonthRunResult",
    "MapRunResult",
    "generate_maps",
    "main",
]

