from __future__ import annotations

from datetime import datetime, timezone

import pytest

from islamic_times.mapper.config import ComputeConfig, MapperConfig, RenderConfig


def test_render_config_normalizes_mode_case() -> None:
    cfg = RenderConfig(map_mode="RAW")
    assert cfg.map_mode == "raw"


def test_mapper_config_rejects_unknown_region() -> None:
    with pytest.raises(ValueError, match="map_region"):
        MapperConfig(date=datetime(2025, 6, 1, tzinfo=timezone.utc), map_region="MARS")


def test_compute_config_resolved_workers_respects_override() -> None:
    cfg = ComputeConfig(max_workers=3)
    assert cfg.resolved_workers == 3

