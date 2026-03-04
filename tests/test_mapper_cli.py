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

