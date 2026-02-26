"""Geodata loading, clipping, and reuse cache for mapper runs."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Sequence, Tuple

from ._deps import require_mapper_dependencies


def _make_valid(geometry: Any) -> Any:
    """Best-effort geometry repair across Shapely 1/2."""
    try:
        from shapely.validation import make_valid as shapely_make_valid  # type: ignore

        return shapely_make_valid(geometry)
    except Exception:
        try:
            return geometry.buffer(0)
        except Exception:
            return geometry


def _clean_geometries(gdf: Any) -> Any:
    gdf = gdf.copy()
    gdf["geometry"] = gdf.geometry.apply(lambda g: _make_valid(g) if not g.is_valid else g)
    gdf = gdf.explode(index_parts=False, ignore_index=True)
    return gdf[~gdf.geometry.is_empty]


def load_shapefiles(states_path: str, places_path: str, cities: Sequence[str]) -> tuple[Any, Any]:
    """Load state and places layers and keep one place feature per city name."""
    require_mapper_dependencies()
    import geopandas as gpd

    states_gdf = gpd.read_file(states_path)
    places_gdf = gpd.read_file(places_path)
    places_gdf = places_gdf[places_gdf["NAME"].isin(cities)]
    places_gdf = places_gdf.loc[places_gdf.groupby("NAME")["POP_MAX"].idxmax()]
    return states_gdf, places_gdf


def clip_map(states_gdf: Any, places_gdf: Any, bbox: tuple[float, float, float, float]) -> tuple[Any, Any]:
    """Clip state/place geodata to a rectangular bbox."""
    require_mapper_dependencies()
    from shapely.geometry import box

    minx, maxx, miny, maxy = bbox
    bbox_geom = box(minx, miny, maxx, maxy)

    if places_gdf.crs != states_gdf.crs:
        places_gdf = places_gdf.to_crs(states_gdf.crs)

    idx = list(states_gdf.sindex.query(bbox_geom))
    states_sub = states_gdf.iloc[idx]
    states_sub = _clean_geometries(states_sub)

    states_clip = states_sub.clip(bbox_geom)
    places_clip = places_gdf.clip(bbox_geom)
    return states_clip, places_clip


@dataclass(slots=True)
class GeoDataCache:
    """Process-local cache for loaded and clipped geodata."""

    states_path: str
    places_path: str
    _states: Any | None = None
    _places: Any | None = None
    _clips: Dict[Tuple[str, Tuple[str, ...], tuple[float, float, float, float]], tuple[Any, Any]] = field(
        default_factory=dict
    )

    def get_clipped(self, region: str, cities: Sequence[str], bbox: tuple[float, float, float, float]) -> tuple[Any, Any]:
        key = (region, tuple(cities), bbox)
        cached = self._clips.get(key)
        if cached is not None:
            return cached

        if self._states is None or self._places is None:
            self._states, self._places = load_shapefiles(self.states_path, self.places_path, cities)

        clipped = clip_map(self._states, self._places, bbox)
        self._clips[key] = clipped
        return clipped

