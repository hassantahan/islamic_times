"""Optional dependency probing and guardrails for mapper workflows."""

from __future__ import annotations

_MAPPER_IMPORT_ERROR: Exception | None = None


def _probe_mapping_dependencies() -> None:
    import geopandas  # noqa: F401
    import matplotlib  # noqa: F401
    import psutil  # noqa: F401
    import shapely  # noqa: F401


try:
    _probe_mapping_dependencies()
except ImportError as exc:  # pragma: no cover - depends on local environment
    _MAPPER_IMPORT_ERROR = exc


def require_mapper_dependencies() -> None:
    """Fail fast with actionable guidance when mapper-only deps are missing."""
    if _MAPPER_IMPORT_ERROR is not None:
        raise ImportError(
            "Mapping dependencies are not installed. Install with: "
            "pip install \"islamic_times[map]\" (or pip install -e \".[map]\" for local development)."
        ) from _MAPPER_IMPORT_ERROR

