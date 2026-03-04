from __future__ import annotations

from datetime import datetime, timezone

import pytest

from islamic_times.mapper.config import ComputeConfig, MapperConfig, RenderConfig
from islamic_times.mapper.palette import category_labels


def test_render_config_normalizes_mode_case() -> None:
    cfg = RenderConfig(map_mode="RAW")
    assert cfg.map_mode == "raw"


def test_mapper_config_rejects_unknown_region() -> None:
    with pytest.raises(ValueError, match="map_region"):
        MapperConfig(date=datetime(2025, 6, 1, tzinfo=timezone.utc), map_region="MARS")


def test_compute_config_resolved_workers_respects_override() -> None:
    cfg = ComputeConfig(max_workers=3)
    assert cfg.resolved_workers == 3


def test_compute_config_accepts_shaukat_criterion() -> None:
    cfg = ComputeConfig(criterion=2)
    assert cfg.criterion == 2


def test_compute_config_rejects_unknown_criterion() -> None:
    with pytest.raises(ValueError, match="criterion"):
        ComputeConfig(criterion=9)


def test_compute_config_rejects_nonpositive_adaptive_min_block_cells() -> None:
    with pytest.raises(ValueError, match="adaptive_min_block_cells"):
        ComputeConfig(adaptive_min_block_cells=0)


def test_compute_config_rejects_nonpositive_adaptive_max_depth() -> None:
    with pytest.raises(ValueError, match="adaptive_max_depth"):
        ComputeConfig(adaptive_max_depth=0)


def test_shaukat_palette_contains_expected_terminal_labels() -> None:
    labels = category_labels(2)
    assert labels[-5:] == [
        "F: Not visible.",
        "D: Visible with optical aid only.",
        "C: Optical aid needed to find the moon.",
        "B: Visible under perfect conditions.",
        "A: Easily visible",
    ]
