from __future__ import annotations

from datetime import datetime
from types import SimpleNamespace
from typing import Any

import pytest

import islamic_times.mapper.cli as mapper_cli


def test_list_methods_prints_primary_keys(capsys: pytest.CaptureFixture[str]) -> None:
    exit_code = mapper_cli.main(["--list-methods"])
    captured = capsys.readouterr()

    assert exit_code == 0
    assert "Available prayer methods" in captured.out
    assert "JAFARI" in captured.out
    assert "FRANCE" in captured.out
    assert "RUSSIA" in captured.out
    assert "SINGAPORE" in captured.out


def test_list_criteria_prints_supported_codes(capsys: pytest.CaptureFixture[str]) -> None:
    exit_code = mapper_cli.main(["--list-criteria"])
    captured = capsys.readouterr()

    assert exit_code == 0
    assert "Available visibility criteria" in captured.out
    assert "0: Odeh (2006)" in captured.out
    assert "1: Yallop (1997)" in captured.out
    assert "2: Shaukat (n.d.)" in captured.out


def test_explain_criterion_prints_detailed_text(capsys: pytest.CaptureFixture[str]) -> None:
    exit_code = mapper_cli.main(["--explain-criterion", "2"])
    captured = capsys.readouterr()

    assert exit_code == 0
    assert "Criterion 2: Shaukat (n.d.)" in captured.out
    assert "Yallop q-values with Shaukat-specific classification thresholds." in captured.out


def test_invalid_criterion_emits_actionable_hint(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit):
        mapper_cli.main(["generate", "--criterion", "9"])
    captured = capsys.readouterr()

    assert "Hint: supported values are 0 (Odeh), 1 (Yallop), and 2 (Shaukat)." in captured.err
    assert "--list-criteria" in captured.err


def test_invalid_map_mode_emits_actionable_hint(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit):
        mapper_cli.main(["generate", "--map_mode", "heat"])
    captured = capsys.readouterr()

    assert "Hint: --map_mode must be one of: raw, category." in captured.err


def test_generate_path_still_dispatches_to_pipeline(monkeypatch: pytest.MonkeyPatch) -> None:
    called: dict[str, Any] = {}

    def fake_generate_maps(cfg: Any, perf_report_path: str | None = None) -> SimpleNamespace:
        called["cfg"] = cfg
        called["perf_report_path"] = perf_report_path
        return SimpleNamespace()

    monkeypatch.setattr(mapper_cli, "generate_maps", fake_generate_maps)
    exit_code = mapper_cli.main(
        [
            "generate",
            "--date",
            "2025-06-01T00:00:00",
            "--criterion",
            "2",
            "--map_mode",
            "category",
            "--perf_report",
            "mapper_logs/test.json",
        ]
    )

    assert exit_code == 0
    assert called["cfg"].compute.criterion == 2
    assert called["cfg"].render.map_mode == "category"
    assert called["cfg"].date == datetime(2025, 6, 1, 0, 0, 0)
    assert called["perf_report_path"] == "mapper_logs/test.json"


def test_cache_command_builds_and_saves_cache(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Any,
    capsys: pytest.CaptureFixture[str],
) -> None:
    called: dict[str, Any] = {}

    def fake_build_visibility_cache(cfg: Any) -> Any:
        called["cfg"] = cfg
        return SimpleNamespace()

    def fake_save_visibility_cache_tabular(cache: Any, cache_path: str) -> Any:
        called["cache"] = cache
        called["cache_path"] = cache_path
        return tmp_path / "cache_dir"

    monkeypatch.setattr(mapper_cli, "build_visibility_cache", fake_build_visibility_cache)
    monkeypatch.setattr(mapper_cli, "save_visibility_cache_tabular", fake_save_visibility_cache_tabular)

    exit_code = mapper_cli.main(
        [
            "cache",
            "--date",
            "2025-06-01T00:00:00",
            "--criterion",
            "1",
            "--cache_path",
            "cache/test_cache",
        ]
    )

    assert exit_code == 0
    assert called["cfg"].compute.criterion == 1
    assert called["cache_path"] == "cache/test_cache"
    assert f"Saved cache to {tmp_path / 'cache_dir'}" in capsys.readouterr().out


def test_export_csv_command_builds_and_saves_csv(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Any,
    capsys: pytest.CaptureFixture[str],
) -> None:
    called: dict[str, Any] = {}

    def fake_build_visibility_cache(cfg: Any) -> Any:
        called["cfg"] = cfg
        return SimpleNamespace()

    def fake_save_visibility_cache_tabular(cache: Any, csv_path: str) -> Any:
        called["cache"] = cache
        called["csv_path"] = csv_path
        return tmp_path / "csv_export"

    monkeypatch.setattr(mapper_cli, "build_visibility_cache", fake_build_visibility_cache)
    monkeypatch.setattr(mapper_cli, "save_visibility_cache_tabular", fake_save_visibility_cache_tabular)

    exit_code = mapper_cli.main(
        [
            "export-csv",
            "--date",
            "2025-06-01T00:00:00",
            "--criterion",
            "1",
            "--csv_path",
            "exports/test_csv",
        ]
    )

    assert exit_code == 0
    assert called["cfg"].compute.criterion == 1
    assert called["csv_path"] == "exports/test_csv"
    assert f"Saved CSV export to {tmp_path / 'csv_export'}" in capsys.readouterr().out


def test_render_cache_command_loads_and_renders(monkeypatch: pytest.MonkeyPatch) -> None:
    called: dict[str, Any] = {}
    fake_cache = SimpleNamespace()

    def fake_load_visibility_cache_tabular(cache_path: str) -> Any:
        called["load_path"] = cache_path
        return fake_cache

    def fake_render_maps_from_cache(
        cache: Any,
        master_path: str | None = None,
        save_logs: bool | None = None,
        perf_report_path: str | None = None,
    ) -> Any:
        called["cache"] = cache
        called["master_path"] = master_path
        called["save_logs"] = save_logs
        called["perf_report_path"] = perf_report_path
        return SimpleNamespace()

    monkeypatch.setattr(mapper_cli, "load_visibility_cache_tabular", fake_load_visibility_cache_tabular)
    monkeypatch.setattr(mapper_cli, "render_maps_from_cache", fake_render_maps_from_cache)

    exit_code = mapper_cli.main(
        [
            "render-cache",
            "--cache_path",
            "cache/test",
            "--master_path",
            "maps_cached/",
            "--save_logs",
            "--perf_report",
            "mapper_logs/cache_render.json",
        ]
    )

    assert exit_code == 0
    assert called["load_path"] == "cache/test"
    assert called["cache"] is fake_cache
    assert called["master_path"] == "maps_cached/"
    assert called["save_logs"] is True
    assert called["perf_report_path"] == "mapper_logs/cache_render.json"
